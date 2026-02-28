
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include </opt/local/dune_210/dune-alugrid-2.10.dev20221009/dune/alugrid/3d/alugrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/common/scsgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscvec.h"
#include "petscis.h"
#include "slepc.h"
#include "slepceps.h"
#include "slepcst.h"

#include "shapefunctions.hh"
#include "qudsim_timer.hh"
#include "qudsim_physics.hh"

using namespace Dune;
using namespace QuDSim;


struct CoordTuple {
    double c[3];
    bool operator<(const CoordTuple& o) const {
        for (int d = 0; d < 3; ++d) {
            if (c[d] < o.c[d] - 1e-14) return true;
            if (c[d] > o.c[d] + 1e-14) return false;
        }
        return false;
    }
    bool operator==(const CoordTuple& o) const {
        for (int d = 0; d < 3; ++d)
            if (std::abs(c[d] - o.c[d]) > 1e-14) return false;
        return true;
    }
};

template<class GV>
void buildGlobalIndexMap(const GV& gv, std::vector<PetscInt>& l2g, PetscInt& globalN,
                         int mpiRank, int mpiSize)
{
    static const int dim = GV::dimension;
    const auto& set = gv.indexSet();
    const int N = gv.size(dim);

    // Step 1: Collect local vertex coordinates
    std::vector<double> myCoords(N * dim);
    for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
        int localIdx = set.index(*it);
        auto pos = it->geometry().center();
        for (int d = 0; d < dim; ++d)
            myCoords[localIdx * dim + d] = pos[d];
    }

    // Step 2: Gather all coordinates from all ranks
    int myCount = N * dim;
    std::vector<int> allCounts(mpiSize);
    MPI_Allgather(&myCount, 1, MPI_INT, allCounts.data(), 1, MPI_INT, PETSC_COMM_WORLD);

    std::vector<int> displs(mpiSize, 0);
    for (int p = 1; p < mpiSize; ++p)
        displs[p] = displs[p-1] + allCounts[p-1];
    int totalDoubles = displs[mpiSize-1] + allCounts[mpiSize-1];

    std::vector<double> allCoords(totalDoubles);
    MPI_Allgatherv(myCoords.data(), myCount, MPI_DOUBLE,
                   allCoords.data(), allCounts.data(), displs.data(), MPI_DOUBLE,
                   PETSC_COMM_WORLD);

    // Step 3: Build sorted unique vertex list
    int totalVerts = totalDoubles / dim;
    std::vector<CoordTuple> allVerts(totalVerts);
    for (int v = 0; v < totalVerts; ++v) {
        allVerts[v].c[0] = allVerts[v].c[1] = allVerts[v].c[2] = 0.0;
        for (int d = 0; d < dim; ++d)
            allVerts[v].c[d] = allCoords[v * dim + d];
    }
    std::sort(allVerts.begin(), allVerts.end());
    allVerts.erase(std::unique(allVerts.begin(), allVerts.end()), allVerts.end());
    globalN = (PetscInt)allVerts.size();

    // Step 4: Coord -> global index map
    std::map<CoordTuple, PetscInt> coordToGlobal;
    for (PetscInt i = 0; i < globalN; ++i)
        coordToGlobal[allVerts[i]] = i;

    // Step 5: Build local -> global
    l2g.resize(N);
    for (int i = 0; i < N; ++i) {
        CoordTuple ct;
        ct.c[0] = ct.c[1] = ct.c[2] = 0.0;
        for (int d = 0; d < dim; ++d)
            ct.c[d] = myCoords[i * dim + d];
        auto findIt = coordToGlobal.find(ct);
        if (findIt != coordToGlobal.end()) {
            l2g[i] = findIt->second;
        } else {
            bool found = false;
            for (PetscInt g = 0; g < globalN; ++g) {
                bool match = true;
                for (int d = 0; d < dim; ++d)
                    if (std::abs(ct.c[d] - allVerts[g].c[d]) > 1e-12) { match = false; break; }
                if (match) { l2g[i] = g; found = true; break; }
            }
            if (!found) {
                std::cerr << "ERROR: vertex " << i << " on rank " << mpiRank
                          << " not found in global list!" << std::endl;
                MPI_Abort(PETSC_COMM_WORLD, 1);
            }
        }
    }
    if (mpiRank == 0)
        std::cout << "Global index map: " << globalN << " unique vertices" << std::endl;
}


template<class GV>
void gatherVecToLocal(Vec petscVec, const GV& gv, const std::vector<PetscInt>& l2g,
                      std::vector<double>& localData)
{
    static const int dim = GV::dimension;
    const int N = gv.size(dim);
    localData.resize(N, 0.0);

    Vec u_all = nullptr;
    VecScatter scatterAll = nullptr;
    VecScatterCreateToAll(petscVec, &scatterAll, &u_all);
    VecScatterBegin(scatterAll, petscVec, u_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatterAll, petscVec, u_all, INSERT_VALUES, SCATTER_FORWARD);

    const PetscScalar* arr;
    VecGetArrayRead(u_all, &arr);
    for (int i = 0; i < N; ++i)
        localData[i] = PetscRealPart(arr[l2g[i]]);
    VecRestoreArrayRead(u_all, &arr);

    VecScatterDestroy(&scatterAll);
    VecDestroy(&u_all);
}

template<class GV>
class ParallelPoissonFEM
{
public:
    static const int dim = GV::dimension;
    typedef typename GV::ctype ctype;
    using Vector = Dune::BlockVector<Dune::FieldVector<ctype, 1>>;

private:
    const GV& gv_;
    const std::vector<PetscInt>& l2g_;
    PetscInt globalN_;

public:
    ParallelPoissonFEM(const GV& gv, const std::vector<PetscInt>& l2g, PetscInt globalN)
        : gv_(gv), l2g_(l2g), globalN_(globalN) {}

    void solve(Vector& charge, Vector& potential, double Vss,
               Vector& deltaRho, Vector& deltaPotential, int it0,
               const Vector& epsi, const Vector& regionid)
    {
        const auto& set = gv_.indexSet();
        const int N = gv_.size(dim);
        P1ShapeFunctionSet<ctype,ctype,dim> basis =
            P1ShapeFunctionSet<ctype,ctype,dim>::instance();
        Dune::SingleCodimSingleGeomTypeMapper<GV, 0> elementmap(gv_);

        PetscMPIInt rank, poissonMpiSize;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &poissonMpiSize);

        // Create ISLocalToGlobalMapping
        ISLocalToGlobalMapping lgmap;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, N, l2g_.data(),
                                     PETSC_COPY_VALUES, &lgmap);

        // Create parallel vectors
        Vec b, x;
        VecCreate(PETSC_COMM_WORLD, &b);
        VecSetSizes(b, PETSC_DECIDE, globalN_);
        VecSetFromOptions(b);
        VecSetLocalToGlobalMapping(b, lgmap);
        VecDuplicate(b, &x);
        VecSet(b, 0.0);

        // Create parallel matrix
        Mat A;
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, globalN_, globalN_);
        MatSetFromOptions(A);
        MatSetLocalToGlobalMapping(A, lgmap, lgmap);

        // Preallocation
        int maxNnz = 0;
        std::vector<std::set<int>> adj(N);
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto ref = Dune::referenceElement(it->geometry());
            int nv = ref.size(dim);
            for (int i = 0; i < nv; ++i) {
                int ii = set.subIndex(*it, i, dim);
                for (int j = 0; j < nv; ++j)
                    adj[ii].insert(set.subIndex(*it, j, dim));
            }
        }
        for (int i = 0; i < N; ++i)
            if ((int)adj[i].size() > maxNnz) maxNnz = (int)adj[i].size();
        MatMPIAIJSetPreallocation(A, maxNnz, NULL, maxNnz, NULL);
        MatSeqAIJSetPreallocation(A, maxNnz, NULL);
        MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        // Element assembly
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto geo = it->geometry();
            auto ref = Dune::referenceElement(geo);
            int nv = ref.size(dim);
            int elemIdx = elementmap.index(*it);
            double eps_elem = (double)epsi[elemIdx];

            std::vector<PetscInt> localRows(nv);
            for (int i = 0; i < nv; ++i)
                localRows[i] = (PetscInt)set.subIndex(*it, i, dim);

            std::vector<PetscScalar> Ke(nv*nv, 0.0);
            std::vector<PetscScalar> be(nv, 0.0);

            // Stiffness
            const auto& rule = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 1);
            for (auto r = rule.begin(); r != rule.end(); ++r) {
                auto JinvT = geo.jacobianInverseTransposed(r->position());
                ctype w = r->weight();
                ctype detJ = geo.integrationElement(r->position());
                for (int i = 0; i < nv; ++i) {
                    Dune::FieldVector<ctype,dim> gi;
                    JinvT.mv(basis[i].evaluateGradient(r->position()), gi);
                    for (int j = 0; j < nv; ++j) {
                        Dune::FieldVector<ctype,dim> gj;
                        JinvT.mv(basis[j].evaluateGradient(r->position()), gj);
                        Ke[i*nv+j] += (PetscScalar)(eps_elem * (gi * gj) * w * detJ);
                    }
                }
            }

            // RHS + NR mass term
            const auto& rule2 = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 2);
            for (auto r = rule2.begin(); r != rule2.end(); ++r) {
                ctype w = r->weight();
                ctype detJ = geo.integrationElement(r->position());
                double rho_qp = 0.0, dRho_qp = 0.0;
                for (int i = 0; i < nv; ++i) {
                    ctype phi = basis[i].evaluateFunction(r->position());
                    int idx = set.subIndex(*it, i, dim);
                    rho_qp  += phi * (double)charge[idx];
                    dRho_qp += phi * (double)deltaRho[idx];
                }
                for (int i = 0; i < nv; ++i) {
                    ctype phi_i = basis[i].evaluateFunction(r->position());
                    be[i] -= (PetscScalar)(phi_i * rho_qp * w * detJ);
                    for (int j = 0; j < nv; ++j) {
                        ctype phi_j = basis[j].evaluateFunction(r->position());
                        Ke[i*nv+j] -= (PetscScalar)(dRho_qp * phi_i * phi_j * w * detJ);
                    }
                }
            }

            MatSetValuesLocal(A, nv, localRows.data(), nv, localRows.data(),
                              Ke.data(), ADD_VALUES);
            VecSetValuesLocal(b, nv, localRows.data(), be.data(), ADD_VALUES);
        }

        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b);
        VecAssemblyEnd(b);

       
        {
            // Step 1: Each rank finds its local boundary nodes
            std::vector<PetscInt> myBcRows;
            std::vector<PetscScalar> myBcVals;
            for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
                auto ref = Dune::referenceElement(it->geometry());
                for (auto is = gv_.ibegin(*it); is != gv_.iend(*it); ++is) {
                    if (!is->boundary()) continue;
                    int nv_f = ref.size(is->indexInInside(), 1, dim);
                    for (int k = 0; k < nv_f; ++k) {
                        int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
                        int li = set.subIndex(*it, lv, dim);
                        PetscInt gi = l2g_[li];
                        auto pos = it->geometry().corner(lv);
                        double xcoord = pos[0];

                        PetscScalar bcval = 0.0;
                        if (xcoord <= -6.9) {
                            // Gate contact
                            bcval = (it0 == 0) ? (PetscScalar)(Vss - (double)potential[li])
                                               : (PetscScalar)0.0;
                        } else if (xcoord >= 49.9) {
                            // Bulk silicon contact
                            bcval = (it0 == 0) ? (PetscScalar)(0.0 - (double)potential[li])
                                               : (PetscScalar)0.0;
                        } else {
                            // Lateral boundaries (y=0 or y=50)
                            bcval = (it0 == 0) ? (PetscScalar)(0.0 - (double)potential[li])
                                               : (PetscScalar)0.0;
                        }
                        myBcRows.push_back(gi);
                        myBcVals.push_back(bcval);
                    }
                }
            }

            // Step 2: Gather ALL BC rows to ALL ranks (MatZeroRows is collective!)
            int myBcCount = (int)myBcRows.size();
            std::vector<int> allBcCounts(poissonMpiSize);
            MPI_Allgather(&myBcCount, 1, MPI_INT, allBcCounts.data(), 1, MPI_INT,
                          PETSC_COMM_WORLD);

            std::vector<int> bcDispls(poissonMpiSize, 0);
            for (int p = 1; p < poissonMpiSize; ++p)
                bcDispls[p] = bcDispls[p-1] + allBcCounts[p-1];
            int totalBcCount = bcDispls[poissonMpiSize-1] + allBcCounts[poissonMpiSize-1];

            std::vector<PetscInt> allBcRows(totalBcCount);
            {
                std::vector<int> myBcRows_int(myBcRows.begin(), myBcRows.end());
                std::vector<int> allBcRows_int(totalBcCount);
                MPI_Allgatherv(myBcRows_int.data(), myBcCount, MPI_INT,
                               allBcRows_int.data(), allBcCounts.data(), bcDispls.data(),
                               MPI_INT, PETSC_COMM_WORLD);
                for (int i = 0; i < totalBcCount; ++i)
                    allBcRows[i] = (PetscInt)allBcRows_int[i];
            }
            std::sort(allBcRows.begin(), allBcRows.end());
            allBcRows.erase(std::unique(allBcRows.begin(), allBcRows.end()), allBcRows.end());

            // Step 3: Set RHS values for BC rows
            for (size_t k = 0; k < myBcRows.size(); ++k)
                VecSetValue(b, myBcRows[k], myBcVals[k], INSERT_VALUES);
            VecAssemblyBegin(b); VecAssemblyEnd(b);

            // Step 4: Zero rows in matrix — ALL ranks pass the SAME complete list
            // Equivalent to serial: A[i]=0; A[i][i]=1; b[i]=bcval
            MatZeroRows(A, (PetscInt)allBcRows.size(), allBcRows.data(),
                        (PetscScalar)1.0, NULL, NULL);

            if (rank == 0)
                std::cout << "  [Poisson] Dirichlet BCs on ALL boundaries: "
                          << allBcRows.size() << " DOFs\n";
        }

        // Solve
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, A, A);
        KSPSetType(ksp, KSPGMRES);
        KSPGMRESSetRestart(ksp, 100);
        PC pc; KSPGetPC(ksp, &pc);
        PCSetType(pc, PCASM);
        PCASMSetOverlap(pc, 1);
        KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b, x);

        PetscInt its; KSPConvergedReason reason;
        KSPGetIterationNumber(ksp, &its);
        KSPGetConvergedReason(ksp, &reason);
        if (rank == 0)
            std::cout << "  Poisson KSP: " << its << " iters (reason=" << reason << ")" << std::endl;

        // Gather solution back to local
        std::vector<double> localSol;
        gatherVecToLocal(x, gv_, l2g_, localSol);
        for (int i = 0; i < N; ++i)
            deltaPotential[i] = localSol[i];

        ISLocalToGlobalMappingDestroy(&lgmap);
        KSPDestroy(&ksp);
        MatDestroy(&A);
        VecDestroy(&b);
        VecDestroy(&x);
    }
};


template<class GV>
class ParallelSchrodingerFEM
{
public:
    static const int dim = GV::dimension;
    typedef typename GV::ctype ctype;
    using Vector = Dune::BlockVector<Dune::FieldVector<ctype, 1>>;

private:
    const GV& gv_;
    const std::vector<PetscInt>& l2g_;
    PetscInt globalN_;
    int carrier_type_;
    double Ef_;
    double beta_;
    int num_requested_;
    int nconv_;
    std::vector<double> eigenvalues_;
    std::vector<std::vector<double>> wavefunctions_global_;

public:
    ParallelSchrodingerFEM(const GV& gv, const std::vector<PetscInt>& l2g,
                           PetscInt globalN, int carrier_type, double Ef,
                           int num_eigenvalues = 10, double T = 300.0)
        : gv_(gv), l2g_(l2g), globalN_(globalN),
          carrier_type_(carrier_type), Ef_(Ef),
          num_requested_(num_eigenvalues), nconv_(0)
    {
        beta_ = Constants::Hartree_to_eV / (Constants::kB_eV * T);
    }

    void solve(const Vector& potential, const Vector& emass,
               Vector& carrier_density)
    {
        const auto& set = gv_.indexSet();
        const int N = gv_.size(dim);
        P1ShapeFunctionSet<ctype,ctype,dim> basis =
            P1ShapeFunctionSet<ctype,ctype,dim>::instance();
        Dune::SingleCodimSingleGeomTypeMapper<GV, 0> elementmap(gv_);

        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        // Create ISLocalToGlobalMapping
        ISLocalToGlobalMapping lgmap;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, N, l2g_.data(),
                                     PETSC_COPY_VALUES, &lgmap);

        // Create parallel H and M matrices
        Mat H, M;

        // Preallocation
        int maxNnz = 0;
        std::vector<std::set<int>> adj(N);
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto ref = Dune::referenceElement(it->geometry());
            int nv = ref.size(dim);
            for (int i = 0; i < nv; ++i) {
                int ii = set.subIndex(*it, i, dim);
                for (int j = 0; j < nv; ++j)
                    adj[ii].insert(set.subIndex(*it, j, dim));
            }
        }
        for (int i = 0; i < N; ++i)
            if ((int)adj[i].size() > maxNnz) maxNnz = (int)adj[i].size();

        MatCreate(PETSC_COMM_WORLD, &H);
        MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, globalN_, globalN_);
        MatSetFromOptions(H);
        MatSetLocalToGlobalMapping(H, lgmap, lgmap);
        MatMPIAIJSetPreallocation(H, maxNnz, NULL, maxNnz, NULL);
        MatSeqAIJSetPreallocation(H, maxNnz, NULL);
        MatSetOption(H, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        MatCreate(PETSC_COMM_WORLD, &M);
        MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, globalN_, globalN_);
        MatSetFromOptions(M);
        MatSetLocalToGlobalMapping(M, lgmap, lgmap);
        MatMPIAIJSetPreallocation(M, maxNnz, NULL, maxNnz, NULL);
        MatSeqAIJSetPreallocation(M, maxNnz, NULL);
        MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        // Element assembly
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto geo = it->geometry();
            auto ref = Dune::referenceElement(geo);
            int nv = ref.size(dim);
            int elemIdx = elementmap.index(*it);
            double meff = (double)emass[elemIdx];
            if (meff <= 0.0) meff = 1.0;

            std::vector<PetscInt> localRows(nv);
            for (int i = 0; i < nv; ++i)
                localRows[i] = (PetscInt)set.subIndex(*it, i, dim);

            std::vector<PetscScalar> He(nv*nv, 0.0), Me(nv*nv, 0.0);

            // Kinetic energy
            const auto& rule = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 1);
            for (auto r = rule.begin(); r != rule.end(); ++r) {
                auto JinvT = geo.jacobianInverseTransposed(r->position());
                ctype w = r->weight();
                ctype detJ = geo.integrationElement(r->position());
                for (int i = 0; i < nv; ++i) {
                    Dune::FieldVector<ctype,dim> gi;
                    JinvT.mv(basis[i].evaluateGradient(r->position()), gi);
                    for (int j = 0; j < nv; ++j) {
                        Dune::FieldVector<ctype,dim> gj;
                        JinvT.mv(basis[j].evaluateGradient(r->position()), gj);
                        He[i*nv+j] += (PetscScalar)((1.0/(2.0*meff)) * (gi*gj) * w * detJ);
                    }
                }
            }

            // Potential + mass
            const auto& rule2 = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 2);
            for (auto r = rule2.begin(); r != rule2.end(); ++r) {
                ctype w = r->weight();
                ctype detJ = geo.integrationElement(r->position());
                double V_qp = 0.0;
                for (int i = 0; i < nv; ++i) {
                    ctype phi = basis[i].evaluateFunction(r->position());
                    V_qp += phi * (double)potential[set.subIndex(*it, i, dim)];
                }
                for (int i = 0; i < nv; ++i) {
                    ctype phi_i = basis[i].evaluateFunction(r->position());
                    for (int j = 0; j < nv; ++j) {
                        ctype phi_j = basis[j].evaluateFunction(r->position());
                        He[i*nv+j] += (PetscScalar)(V_qp * phi_i * phi_j * w * detJ);
                        Me[i*nv+j] += (PetscScalar)(phi_i * phi_j * w * detJ);
                    }
                }
            }

            MatSetValuesLocal(H, nv, localRows.data(), nv, localRows.data(),
                              He.data(), ADD_VALUES);
            MatSetValuesLocal(M, nv, localRows.data(), nv, localRows.data(),
                              Me.data(), ADD_VALUES);
        }

        MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

        // Mark as Hermitian (required for complex PETSc with EPS_GHEP)
        MatSetOption(H, MAT_HERMITIAN, PETSC_TRUE);
        MatSetOption(M, MAT_HERMITIAN, PETSC_TRUE);
        MatSetOption(M, MAT_SPD, PETSC_TRUE);

        // BC: psi=0 at domain boundaries
        // Gather ALL BC rows to ALL ranks (MatZeroRowsColumns is collective!)
        std::vector<PetscInt> myBcRows;
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto ref = Dune::referenceElement(it->geometry());
            for (auto is = gv_.ibegin(*it); is != gv_.iend(*it); ++is) {
                if (!is->boundary()) continue;
                int nv_f = ref.size(is->indexInInside(), 1, dim);
                for (int k = 0; k < nv_f; ++k) {
                    int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
                    int li = set.subIndex(*it, lv, dim);
                    myBcRows.push_back(l2g_[li]);
                }
            }
        }

        // Allgather BC rows so every rank has the complete list
        PetscMPIInt schMpiSize;
        MPI_Comm_size(PETSC_COMM_WORLD, &schMpiSize);
        int myBcCount = (int)myBcRows.size();
        std::vector<int> allBcCounts(schMpiSize);
        MPI_Allgather(&myBcCount, 1, MPI_INT, allBcCounts.data(), 1, MPI_INT, PETSC_COMM_WORLD);

        std::vector<int> bcDispls(schMpiSize, 0);
        for (int p = 1; p < schMpiSize; ++p)
            bcDispls[p] = bcDispls[p-1] + allBcCounts[p-1];
        int totalBcCount = bcDispls[schMpiSize-1] + allBcCounts[schMpiSize-1];

        std::vector<int> allBcRows_int(totalBcCount);
        {
            std::vector<int> myBcRows_int(myBcRows.begin(), myBcRows.end());
            MPI_Allgatherv(myBcRows_int.data(), myBcCount, MPI_INT,
                           allBcRows_int.data(), allBcCounts.data(), bcDispls.data(),
                           MPI_INT, PETSC_COMM_WORLD);
        }
        std::vector<PetscInt> bcRows(totalBcCount);
        for (int i = 0; i < totalBcCount; ++i)
            bcRows[i] = (PetscInt)allBcRows_int[i];
        std::sort(bcRows.begin(), bcRows.end());
        bcRows.erase(std::unique(bcRows.begin(), bcRows.end()), bcRows.end());

        // ALL ranks call MatZeroRowsColumns with the SAME complete list
        {
            Vec diagH, diagM;
            MatCreateVecs(H, NULL, &diagH);
            MatCreateVecs(M, NULL, &diagM);
            VecSet(diagH, 0.0);
            VecSet(diagM, 0.0);
            PetscScalar bigVal = 1.0e6;
            MatZeroRowsColumns(H, (PetscInt)bcRows.size(), bcRows.data(), bigVal, diagH, diagH);
            MatZeroRowsColumns(M, (PetscInt)bcRows.size(), bcRows.data(), (PetscScalar)1.0, diagM, diagM);
            VecDestroy(&diagH);
            VecDestroy(&diagM);
        }

        // SLEPc eigenvalue solve (parallel)
        EPS eps;
        EPSCreate(PETSC_COMM_WORLD, &eps);
        EPSSetOperators(eps, H, M);
        EPSSetProblemType(eps, EPS_GHEP);
        EPSSetType(eps, EPSKRYLOVSCHUR);
        EPSSetWhichEigenpairs(eps, EPS_TARGET_REAL);
        EPSSetTarget(eps, (PetscScalar)0.0);
        EPSSetDimensions(eps, num_requested_, PETSC_DEFAULT, PETSC_DEFAULT);
        EPSSetTolerances(eps, 1e-8, 1000);

//EPSSetTrueResidual(eps, PETSC_TRUE);
//EPSSetBalance(eps, EPS_BALANCE_NONE);

        // Use shift-and-invert with iterative solver (LU not available for mpiaij)
        ST st;
        EPSGetST(eps, &st);
        STSetType(st, STSINVERT);
        KSP st_ksp;
        STGetKSP(st, &st_ksp);
       // KSPSetType(st_ksp, KSPGMRES);
KSPSetType(st_ksp, KSPCG);
       // KSPGMRESSetRestart(st_ksp, 100);
        PC st_pc;
        KSPGetPC(st_ksp, &st_pc);
       // PCSetType(st_pc, PCASM);
PCSetType(st_pc, PCGAMG);
PCGAMGSetType(st_pc, PCGAMGAGG);
//PCGAMGSetThreshold(st_pc, 0.02);
PCGAMGSetNSmooths(st_pc, 1);
        //PCASMSetOverlap(st_pc, 1);
        KSPSetTolerances(st_ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
KSPSetNormType(st_ksp, KSP_NORM_UNPRECONDITIONED);
        EPSSetFromOptions(eps);
        EPSSolve(eps);

        EPSGetConverged(eps, &nconv_);
        if (rank == 0)
            std::cout << "  Schrodinger (" << (carrier_type_==0?"elec":"hole")
                      << "): " << nconv_ << " eigenvalues" << std::endl;

        // Extract eigenvalues and wavefunctions (gather to all ranks)
        eigenvalues_.resize(nconv_);
        wavefunctions_global_.resize(nconv_);
        Vec xr, xi;
        MatCreateVecs(H, NULL, &xr);
        VecDuplicate(xr, &xi);

        for (int i = 0; i < nconv_; ++i) {
            PetscScalar kr, ki;
            EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
            eigenvalues_[i] = PetscRealPart(kr);

            // Gather full eigenvector to all ranks
            Vec xr_all = nullptr;
            VecScatter scatterAll = nullptr;
            VecScatterCreateToAll(xr, &scatterAll, &xr_all);
            VecScatterBegin(scatterAll, xr, xr_all, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterAll, xr, xr_all, INSERT_VALUES, SCATTER_FORWARD);

            const PetscScalar* arr;
            VecGetArrayRead(xr_all, &arr);
            wavefunctions_global_[i].resize(globalN_);
            for (PetscInt g = 0; g < globalN_; ++g)
                wavefunctions_global_[i][g] = PetscRealPart(arr[g]);
            VecRestoreArrayRead(xr_all, &arr);

            VecScatterDestroy(&scatterAll);
            VecDestroy(&xr_all);

            if (rank == 0 && i < 5)
                std::cout << "    E[" << i << "] = " << eigenvalues_[i]
                          << " Ha = " << Ha_to_eV(eigenvalues_[i]) << " eV" << std::endl;
        }

        // Compute local carrier density from global wavefunctions
        computeCarrierDensity(carrier_density);

        VecDestroy(&xr); VecDestroy(&xi);
        EPSDestroy(&eps);
        ISLocalToGlobalMappingDestroy(&lgmap);
        MatDestroy(&H); MatDestroy(&M);
    }

    void computeCarrierDensity(Vector& carrier_density)
    {
        const int N = gv_.size(dim);
        carrier_density = 0.0;
        for (int i = 0; i < nconv_; ++i) {
            double Ei = eigenvalues_[i];
            double occ;
            if (carrier_type_ == 0) {
                double arg = beta_ * (Ei - Ef_);
                occ = (arg > 500.0) ? 0.0 : (arg < -500.0) ? 1.0 : 1.0/(1.0+std::exp(arg));
            } else {
                double arg = beta_ * (Ei - Ef_);
                occ = (arg > 500.0) ? 1.0 : (arg < -500.0) ? 0.0 : std::exp(arg)/(1.0+std::exp(arg));
            }
            for (int j = 0; j < N; ++j) {
                double psi = wavefunctions_global_[i][l2g_[j]];
                carrier_density[j] += 2.0 * occ * psi * psi;
            }
        }
        double maxn = cm3_to_au3(1e25);
        for (int j = 0; j < N; ++j)
            if ((double)carrier_density[j] > maxn) carrier_density[j] = maxn;
    }

    int getnconv() const { return nconv_; }
    double geteigen(int i) const { return (i < nconv_) ? eigenvalues_[i] : 0.0; }
    double getwave_local(int eigIdx, int localIdx) const {
        return wavefunctions_global_[eigIdx][l2g_[localIdx]];
    }
    const std::vector<double>& getEigenvalues() const { return eigenvalues_; }
};



int main(int argc, char** argv)
{
    try {
        
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
        SlepcInitialize(&argc, &argv, NULL, NULL);

        PetscMPIInt rank, mpiSize;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);

        QuDSimTimer timer(PETSC_COMM_WORLD);

        if (rank == 0) {
            std::cout << "========================================================\n"
                      << "  QuDSim: PARALLEL Self-Consistent Schrodinger-Poisson\n"
                      << "  MOSCAP Structure - 2D FEM (" << mpiSize << " MPI ranks)\n"
                      << "========================================================\n";
        }

        // ============================================================
        // Simulation Parameters (rank 0 reads, then broadcast)
        // ============================================================
        SimulationParams params;
        params.Na = 1e18;
        params.ni = 1.5e10;
        params.T  = 300.0;
        params.Vfb = -0.88;
        params.tolerance = 1e-6;
        params.mixing_alpha = 0.3;
        params.num_eigenvalues = 15;

        int ACCINV = 2;
        double Vs = 1.0;
        int IMAX = 50;

        if (rank == 0) {
            std::cout << "\n  MOSCAP Simulation Setup\n"
                      << "  1 = Accumulation (negative Vg)\n"
                      << "  2 = Inversion (positive Vg)\n"
                      << "  Select mode: ";
            std::cin >> ACCINV;
            std::cout << "  Enter gate voltage (V): ";
            std::cin >> Vs;
            std::cout << "  Enter max iterations: ";
            std::cin >> IMAX;
        }

        // Broadcast inputs to all ranks
        MPI_Bcast(&ACCINV, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(&Vs, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
        MPI_Bcast(&IMAX, 1, MPI_INT, 0, PETSC_COMM_WORLD);

        params.Vg = Vs;
        params.max_iterations = IMAX;
        params.initialize();

        double Vss = params.getVss();
        double Ef  = params.Ef;
        double beta = params.beta;
        double Eg_Ha = eV_to_Ha(Materials::Silicon.Eg);

        if (rank == 0) {
            std::cout << "\n  Parameters (Hartree a.u.):\n"
                      << "    Vss  = " << Vss << " Ha (" << Ha_to_eV(Vss) << " eV)\n"
                      << "    Ef   = " << Ef  << " Ha (" << Ha_to_eV(Ef)  << " eV)\n"
                      << "    Eg   = " << Eg_Ha << " Ha\n"
                      << "    beta = " << beta << "\n"
                      << "    Na   = " << params.Na << " cm^-3\n\n";
        }

        
        timer.start("Mesh Loading");

        constexpr int dim = 2;
        using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming>;
        using GV   = Grid::LeafGridView;
        using ctype = Grid::ctype;
        using ElementMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;
        using VertexMapper  = Dune::SingleCodimSingleGeomTypeMapper<GV, dim>;
        using Vector = Dune::BlockVector<Dune::FieldVector<ctype, 1>>;

        // Read mesh and loadBalance for MPI distribution
        const std::string gridfile = "/home/athira/dune-projects/dune-sc1/src/moscap_40k.msh";

        auto gridptr = Dune::GmshReader<Grid>::read(gridfile);
        gridptr->loadBalance();

        Grid& grid = *gridptr;
        const GV& gv = grid.leafGridView();

        const int N = gv.size(dim);
        const auto& set = gv.indexSet();
        ElementMapper elementmap(gv);
        VertexMapper vertexmap(gv);

        std::cout << "Rank " << rank << ": " << gv.size(0)
                  << " elements, " << N << " vertices (local)" << std::endl;

        // Build coordinate-based global index map
        std::vector<PetscInt> l2g;
        PetscInt globalN;
        buildGlobalIndexMap(gv, l2g, globalN, rank, mpiSize);

        timer.stop("Mesh Loading");

        int nElem = gv.size(0);
        Vector charge(N);          charge = 0.0;
        Vector potential(N);       potential = 0.0;
        Vector potential_elec(N);  potential_elec = 0.0;
        Vector potential_hole(N);  potential_hole = 0.0;
        Vector deltapotential(N);  deltapotential = 0.0;
        Vector deltaRho(N);        deltaRho = 0.0;
        Vector Nen(N);             Nen = 0.0;
        Vector Nep(N);             Nep = 0.0;

        Vector epsi(nElem);    Vector emass(nElem);
        Vector hmass(nElem);   Vector regionid_vec(nElem);

    
        timer.start("Material Setup");

        double Na_au = cm3_to_au3(params.Na);

        for (auto it = gv.begin<0>(); it != gv.end<0>(); ++it) {
            int ei = elementmap.index(*it);
            auto center = it->geometry().center();
            double xc = center[0];
            int rid;
            if (xc < -2.0)       rid = 300;  // Metal gate
            else if (xc < 0.0)   rid = 400;  // SiO2
            else if (xc < 5.0)   rid = 200;  // Si near interface
            else if (xc < 20.0)  rid = 201;  // Si mid
            else                  rid = 202;  // Si bulk

            regionid_vec[ei] = rid;
            Material mat = Materials::getMaterialByRegion(rid);
            epsi[ei]  = Constants::eps0 * mat.eps_r;
            emass[ei] = mat.me_eff;
            hmass[ei] = mat.mh_eff;
        }

        timer.stop("Material Setup");
        if (rank == 0) std::cout << "Materials assigned.\n";

        ParallelSchrodingerFEM<GV> schElec(gv, l2g, globalN, 0, -Ef,
                                            params.num_eigenvalues, params.T);
        ParallelSchrodingerFEM<GV> schHole(gv, l2g, globalN, 1, Ef,
                                            params.num_eigenvalues, params.T);
        ParallelPoissonFEM<GV> poisson(gv, l2g, globalN);

        std::ofstream res_plot;
        if (rank == 0) res_plot.open("Results.dat");

    
        timer.start("Self-Consistent Loop");

        int imax = params.max_iterations;
        double tolerance = params.tolerance;
        double alpha = params.mixing_alpha;
        int it0 = 0;

        if (rank == 0) {
            std::cout << "\n=============================================\n"
                      << "  Self-Consistent Iteration\n"
                      << "  Vg = " << Vs << " V, max_iter = " << imax << "\n"
                      << "=============================================\n";
        }

        for (int iter = 0; iter < imax; ++iter) {
            if (rank == 0)
                std::cout << "\n--- Iteration " << iter << " ---\n";

            // Step 1: Build Schrodinger potentials
            timer.start("Potential Setup");

            for (int p = 0; p < N; ++p) {
                potential_elec[p] = Eg_Ha / 2.0 - potential[p];
                potential_hole[p] = Eg_Ha / 2.0 + potential[p];
            }

            for (auto it = gv.begin<0>(); it != gv.end<0>(); ++it) {
                int ei = elementmap.index(*it);
                int rid = (int)regionid_vec[ei];
                if (rid == 200 || rid == 201 || rid == 202) continue;

                Material mat = Materials::getMaterialByRegion(rid);
                auto ref = Dune::referenceElement(it->geometry());
                int nv = ref.size(dim);
                for (int pt = 0; pt < nv; ++pt) {
                    int idx = set.subIndex(*it, pt, dim);
                    if (rid == 400) {
                        double CBO = eV_to_Ha(mat.VBO_electron);
                        double VBO = eV_to_Ha(mat.VBO_hole);
                        potential_elec[idx] = Eg_Ha/2.0 + CBO - potential[idx];
                        potential_hole[idx] = Eg_Ha/2.0 + VBO + potential[idx];
                    } else if (rid == 300) {
                        potential_elec[idx] = 100.0;
                        potential_hole[idx] = 100.0;
                    }
                }
            }
            timer.stop("Potential Setup");

            // Step 2: Solve Schrodinger for electrons
            timer.start("Schrodinger Electrons");
            if (rank == 0) std::cout << "Solving Schrodinger (electrons)...\n";
            schElec.solve(potential_elec, emass, Nen);
            timer.stop("Schrodinger Electrons");

            // Step 3: Solve Schrodinger for holes
            timer.start("Schrodinger Holes");
            if (rank == 0) std::cout << "Solving Schrodinger (holes)...\n";
            schHole.solve(potential_hole, hmass, Nep);
            timer.stop("Schrodinger Holes");

            // Step 4: Compute charge density (only in Silicon)
            timer.start("Charge Density");
            for (auto it = gv.begin<0>(); it != gv.end<0>(); ++it) {
                int ei = elementmap.index(*it);
                int rid = (int)regionid_vec[ei];
                if (rid != 200 && rid != 201 && rid != 202) continue;

                auto ref = Dune::referenceElement(it->geometry());
                int nv = ref.size(dim);
                for (int pt = 0; pt < nv; ++pt) {
                    int idx = set.subIndex(*it, pt, dim);
                    double Np = (double)Nep[idx];
                    double Nn = (double)Nen[idx];
                    double maxn = cm3_to_au3(1e25);
                    if (Np > maxn) Np = maxn;
                    if (Nn > maxn) Nn = maxn;
                    charge[idx] = Constants::q0 * (Np - Nn - Na_au);
                    deltaRho[idx] = Constants::q0 * (-beta*Np - beta*Nn);
                }
            }
            timer.stop("Charge Density");

            // Step 5: VTK output (collective - all ranks must call)
            {
                std::ostringstream fn;
                fn << "Results_iter_" << iter;
                Dune::VTKWriter<GV> vw(gv);
                std::vector<double> pot_vec(N), chg_vec(N), nen_vec(N), nep_vec(N);
                for (int i = 0; i < N; ++i) {
                    pot_vec[i] = (double)potential[i];
                    chg_vec[i] = (double)charge[i];
                    nen_vec[i] = (double)Nen[i];
                    nep_vec[i] = (double)Nep[i];
                }
                vw.addVertexData(pot_vec, "potential");
                vw.addVertexData(chg_vec, "charge");
                vw.addVertexData(nen_vec, "Nen");
                vw.addVertexData(nep_vec, "Nep");
                vw.write(fn.str(), Dune::VTK::appendedraw);
            }

            // Step 6: Solve Poisson
            timer.start("Poisson Solve");
            if (rank == 0) std::cout << "Solving Poisson...\n";
            poisson.solve(charge, potential, Vss, deltaRho, deltapotential, it0,
                         epsi, regionid_vec);
            timer.stop("Poisson Solve");

            // Step 7: Update potential with mixing
            double maxDelta_local = 0.0;
            for (int p = 0; p < N; ++p) {
                double dp = alpha * (double)deltapotential[p];
                potential[p] += dp;
                if (std::fabs(dp) > maxDelta_local) maxDelta_local = std::fabs(dp);
            }
            it0 = 1;

            // Global max across all ranks
            double maxDelta_global = 0.0;
            MPI_Allreduce(&maxDelta_local, &maxDelta_global, 1, MPI_DOUBLE, MPI_MAX,
                         PETSC_COMM_WORLD);

            if (rank == 0)
                std::cout << "  Max |delta_phi| = " << maxDelta_global
                          << " Ha = " << Ha_to_eV(maxDelta_global) << " eV\n";

            // Step 8: Convergence check (global)
            if (maxDelta_global < tolerance) {
                if (rank == 0)
                    std::cout << "\n*** CONVERGENCE at iteration " << iter << " ***\n\n";
                break;
            }
            if (iter == imax-1 && rank == 0)
                std::cout << "\nWARNING: Max iterations reached!\n";

            if (rank == 0)
                res_plot << iter << " " << Vss << " " << maxDelta_global << "\n";
        }

        timer.stop("Self-Consistent Loop");

      
        timer.start("Post-Processing");

        // Final VTK (collective)
        {
            std::vector<double> pot_vec(N), chg_vec(N), nen_vec(N), nep_vec(N);
            std::vector<double> pelec_vec(N), phole_vec(N);
            for (int i = 0; i < N; ++i) {
                pot_vec[i] = (double)potential[i];
                chg_vec[i] = (double)charge[i];
                nen_vec[i] = (double)Nen[i];
                nep_vec[i] = (double)Nep[i];
                pelec_vec[i] = (double)potential_elec[i];
                phole_vec[i] = (double)potential_hole[i];
            }
            Dune::VTKWriter<GV> vtkw(gv);
            vtkw.addVertexData(pot_vec, "potential");
            vtkw.addVertexData(chg_vec, "charge");
            vtkw.addVertexData(nen_vec, "Nen");
            vtkw.addVertexData(nep_vec, "Nep");
            vtkw.addVertexData(pelec_vec, "potential_elec");
            vtkw.addVertexData(phole_vec, "potential_hole");
            vtkw.write("Results_final", Dune::VTK::appendedraw);
        }

        // Eigenvalue output (rank 0 only)
        if (rank == 0) {
            std::ofstream mplot("moscap_eigenvalues.dat");
            mplot << "# Electron eigenvalues\n";
            for (int i = 0; i < schElec.getnconv(); ++i)
                mplot << i << " " << schElec.geteigen(i) << " "
                      << Ha_to_eV(schElec.geteigen(i)) << " eV\n";
            mplot << "\n# Hole eigenvalues\n";
            for (int i = 0; i < schHole.getnconv(); ++i)
                mplot << i << " " << schHole.geteigen(i) << " "
                      << Ha_to_eV(schHole.geteigen(i)) << " eV\n";
            mplot.close();
        }

        // Wavefunction output (rank 0 only)
        if (rank == 0) {
            std::ofstream wfp("wavefunctions.dat");
            wfp << "# x y psi_e0 psi_e1 psi_e2 potential\n";
            for (auto vt = gv.begin<dim>(); vt != gv.end<dim>(); ++vt) {
                auto pos = vt->geometry().corner(0);
                int vi = vertexmap.index(*vt);
                if (std::fabs(pos[1] - 25.0) < 3.0) {
                    wfp << pos[0] << " " << pos[1];
                    for (int i = 0; i < std::min(3, schElec.getnconv()); ++i)
                        wfp << " " << schElec.getwave_local(i, vi);
                    wfp << " " << (double)potential[vi] << "\n";
                }
            }
            wfp.close();
        }

        timer.stop("Post-Processing");
        if (rank == 0) res_plot.close();

        timer.report();
        timer.reportCSV("timings.csv");

        SlepcFinalize();

        if (rank == 0)
            std::cout << "\nSimulation complete.\n"
                      << "Output: Results_final.pvtu, moscap_eigenvalues.dat\n";

        return 0;
    }
    catch (Dune::Exception& e) { std::cerr << "DUNE: " << e.what() << "\n"; }
    catch (std::exception& e)  { std::cerr << "Error: " << e.what() << "\n"; }
    catch (...)                 { std::cerr << "Unknown exception!\n"; }
    return 1;
}
