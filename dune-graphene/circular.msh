/**
 * QuDSim - PARALLEL Graphene Quantum Dot Solver
 * 2D FEM: Schrodinger with PEP (quadratic polynomial eigenvalue)
 * Computes: eigenvalues, wavefunctions, LDOS, tunneling current, dI/dV
 * MPI parallel: DUNE ALUGrid loadBalance + PETSc ISLocalToGlobalMapping
 * Uses: DUNE (grid/FEM), PETSc (linear algebra), SLEPc (PEP eigenvalues)
 */
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
#include "slepcpep.h"
#include "slepcst.h"

#include "shapefunctions.hh"
#include "qudsim_timer.hh"

using namespace Dune;

// ============================================================================
//  Coordinate-based Global Index Map (same pattern as MOSCAP solver)
// ============================================================================
struct CoordTuple {
    double c[3];
    bool operator<(const CoordTuple& o) const {
        for (int d = 0; d < 3; ++d) {
            if (c[d] < o.c[d] - 1e-12) return true;
            if (c[d] > o.c[d] + 1e-12) return false;
        }
        return false;
    }
    bool operator==(const CoordTuple& o) const {
        for (int d = 0; d < 3; ++d)
            if (std::abs(c[d] - o.c[d]) > 1e-12) return false;
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

    std::vector<double> myCoords(N * dim);
    for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
        int localIdx = set.index(*it);
        auto pos = it->geometry().center();
        for (int d = 0; d < dim; ++d)
            myCoords[localIdx * dim + d] = pos[d];
    }

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

    std::map<CoordTuple, PetscInt> coordToGlobal;
    for (PetscInt i = 0; i < globalN; ++i)
        coordToGlobal[allVerts[i]] = i;

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
                    if (std::abs(ct.c[d] - allVerts[g].c[d]) > 1e-10) { match = false; break; }
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


// ============================================================================
//  Create a parallel PETSc matrix with ISLocalToGlobalMapping
// ============================================================================
Mat createParallelMat(PetscInt globalN, ISLocalToGlobalMapping lgmap, int maxNnz)
{
    Mat M;
    MatCreate(PETSC_COMM_WORLD, &M);
    MatSetSizes(M, PETSC_DECIDE, PETSC_DECIDE, globalN, globalN);
    MatSetFromOptions(M);
    MatSetLocalToGlobalMapping(M, lgmap, lgmap);
    MatMPIAIJSetPreallocation(M, maxNnz, NULL, maxNnz, NULL);
    MatSeqAIJSetPreallocation(M, maxNnz, NULL);
    MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    return M;
}


// ============================================================================
//  Gather distributed PETSc Vec to local std::vector via l2g
// ============================================================================
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


// ============================================================================
//                              MAIN
// ============================================================================
int main(int argc, char** argv)
{
    try {
        // ============================================================
        // Initialize MPI + PETSc + SLEPc
        // ============================================================
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
        SlepcInitialize(&argc, &argv, NULL, NULL);

        PetscMPIInt rank, mpiSize;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &mpiSize);

        QuDSimTimer timer(PETSC_COMM_WORLD);

        if (rank == 0) {
            std::cout << "========================================================\n"
                      << "  QuDSim: PARALLEL Graphene Quantum Dot Solver\n"
                      << "  2D FEM + PEP Eigenvalue (" << mpiSize << " MPI ranks)\n"
                      << "========================================================\n";
        }

        // ============================================================
        // Grid Setup with loadBalance
        // ============================================================
        timer.start("Mesh Loading");

        constexpr int dim = 2;
        using Grid = Dune::ALUGrid<dim, 2, Dune::simplex, Dune::conforming>;
        using GV   = Grid::LeafGridView;
        using ctype = Grid::ctype;
        using ElementMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;
        using VertexMapper  = Dune::SingleCodimSingleGeomTypeMapper<GV, dim>;
        using Vector = Dune::BlockVector<Dune::FieldVector<ctype, 1>>;

        const std::string gridfile = "/home/athira/dune-projects/dune-graphene1/src/gqd_refined.msh";
       // const std::string gridfile = "/home/local/dune_210/dune-graphene/src/gqd_hr2.msh"; 

        auto grid = Dune::GmshReader<Grid>::read(gridfile);
        grid->loadBalance();

        const GV& gv = grid->leafGridView();
        const int N = gv.size(dim);
        const auto& iset = gv.indexSet();
        ElementMapper elementmap(gv);
        VertexMapper vertexmap(gv);

        std::cout << "Rank " << rank << ": " << gv.size(0)
                  << " elements, " << N << " vertices (local)" << std::endl;

        // Build global index map
        std::vector<PetscInt> l2g;
        PetscInt globalN;
        buildGlobalIndexMap(gv, l2g, globalN, rank, mpiSize);

        timer.stop("Mesh Loading");

        // ============================================================
        // Compute max adjacency for preallocation
        // ============================================================
        int maxNnz = 0;
        {
            std::vector<std::set<int>> adj(N);
            for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
                auto ref = Dune::referenceElement(it->geometry());
                int nv = ref.size(dim);
                for (int i = 0; i < nv; ++i) {
                    int ii = iset.subIndex(*it, i, dim);
                    for (int j = 0; j < nv; ++j)
                        adj[ii].insert(iset.subIndex(*it, j, dim));
                }
            }
            for (int i = 0; i < N; ++i)
                if ((int)adj[i].size() > maxNnz) maxNnz = (int)adj[i].size();
        }

        // ============================================================
        // Create ISLocalToGlobalMapping
        // ============================================================
        ISLocalToGlobalMapping lgmap;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, N, l2g.data(),
                                     PETSC_COPY_VALUES, &lgmap);

        // ============================================================
        // Create parallel matrices
        // ============================================================
        timer.start("Assembly");

        Mat A_p  = createParallelMat(globalN, lgmap, maxNnz);
        Mat B_p  = createParallelMat(globalN, lgmap, maxNnz);
        Mat C1_p = createParallelMat(globalN, lgmap, maxNnz);
        Mat C2_p = createParallelMat(globalN, lgmap, maxNnz);
        Mat H_p;  // will be built from MatDuplicate after A_p is assembled

        // ============================================================
        // FEM Assembly: stiffness A and mass B using MatSetValuesLocal
        // ============================================================
        P1ShapeFunctionSet<ctype,ctype,dim> basis = P1ShapeFunctionSet<ctype,ctype,dim>::instance();

        const double R_dot = 54.0;       // Quantum dot radius (a.u.)
        const double V_barrier = 0.016;  // Barrier potential outside dot (Hartree)

        for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
            auto geo = it->geometry();
            auto ref = Dune::referenceElement(geo);
            int nv = ref.size(dim);

            std::vector<PetscInt> localRows(nv);
            for (int i = 0; i < nv; ++i)
                localRows[i] = (PetscInt)iset.subIndex(*it, i, dim);

            std::vector<PetscScalar> Ae(nv*nv, 0.0), Be(nv*nv, 0.0);

            const auto& rule = Dune::QuadratureRules<ctype,dim>::rule(it->type(), 2);
            for (auto r = rule.begin(); r != rule.end(); ++r) {
                auto JinvT = geo.jacobianInverseTransposed(r->position());
                ctype w = r->weight();
                ctype detJ = geo.integrationElement(r->position());

                for (int i = 0; i < nv; ++i) {
                    Dune::FieldVector<ctype,dim> gi;
                    JinvT.mv(basis[i].evaluateGradient(r->position()), gi);
                    ctype phii = basis[i].evaluateFunction(r->position());

                    // Radial check at corner i for potential assignment
                    double x_i = geo.corner(i)[0];
                    double y_i = geo.corner(i)[1];
                    double rr = std::sqrt(x_i*x_i + y_i*y_i);
                    double fval = (rr < R_dot) ? 0.0 : V_barrier;

                    for (int j = 0; j < nv; ++j) {
                        Dune::FieldVector<ctype,dim> gj;
                        JinvT.mv(basis[j].evaluateGradient(r->position()), gj);
                        ctype phij = basis[j].evaluateFunction(r->position());

                        // A = 0.5*grad*grad + V*phi_i*phi_j
                        Ae[i*nv+j] += (PetscScalar)((gi*gj)*0.5*w*detJ
                                     + phii*fval*phij*w*detJ);
                        // B = phi_i * phi_j (mass)
                        Be[i*nv+j] += (PetscScalar)(phii*phij*w*detJ);
                    }
                }
            }

            MatSetValuesLocal(A_p, nv, localRows.data(), nv, localRows.data(),
                              Ae.data(), ADD_VALUES);
            MatSetValuesLocal(B_p, nv, localRows.data(), nv, localRows.data(),
                              Be.data(), ADD_VALUES);
        }

        MatAssemblyBegin(A_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(B_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(B_p, MAT_FINAL_ASSEMBLY);

        if (rank == 0) std::cout << "A and B assembled.\n";

        // ============================================================
        // Boundary Conditions
        //   At domain boundary: A[i]=0, A[i][i]=1, B[i]=0
        //   C1[i][i] = 0.25/94.6,  C2[i][i] = i*(-0.5)
        // ============================================================
        {
            std::vector<PetscInt> bcRows;
            for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
                auto ref = Dune::referenceElement(it->geometry());
                for (auto is = gv.ibegin(*it); is != gv.iend(*it); ++is) {
                    if (!is->boundary()) continue;
                    int nv_f = ref.size(is->indexInInside(), 1, dim);
                    for (int k = 0; k < nv_f; ++k) {
                        int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
                        int li = iset.subIndex(*it, lv, dim);
                        bcRows.push_back(l2g[li]);
                    }
                }
            }
            std::sort(bcRows.begin(), bcRows.end());
            bcRows.erase(std::unique(bcRows.begin(), bcRows.end()), bcRows.end());

            if (rank == 0) std::cout << "Boundary DOFs: " << bcRows.size() << std::endl;

            // A: zero rows+cols, diagonal=1 (maintains symmetry for complex PETSc)
            if (!bcRows.empty()) {
                Vec zeroVec;
                MatCreateVecs(A_p, NULL, &zeroVec);
                VecSet(zeroVec, 0.0);
                MatZeroRowsColumns(A_p, (PetscInt)bcRows.size(), bcRows.data(),
                                   (PetscScalar)1.0, zeroVec, zeroVec);
                VecDestroy(&zeroVec);

                // B: zero rows+cols, diagonal=0
                Vec zeroVec2;
                MatCreateVecs(B_p, NULL, &zeroVec2);
                VecSet(zeroVec2, 0.0);
                MatZeroRowsColumns(B_p, (PetscInt)bcRows.size(), bcRows.data(),
                                   (PetscScalar)0.0, zeroVec2, zeroVec2);
                VecDestroy(&zeroVec2);
            }

            // C1 and C2: boundary absorbing terms
            for (size_t k = 0; k < bcRows.size(); ++k) {
                PetscInt gi = bcRows[k];
                PetscScalar c1_val = (PetscScalar)(0.25 / 94.6);
                PetscScalar c2_val = PETSC_i * (PetscScalar)(-0.5);
                MatSetValue(C1_p, gi, gi, c1_val, INSERT_VALUES);
                MatSetValue(C2_p, gi, gi, c2_val, INSERT_VALUES);
            }
        }

        MatAssemblyBegin(C1_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(C1_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(C2_p, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(C2_p, MAT_FINAL_ASSEMBLY);

        // H = A + C1 - fr*B  (build from duplicate of A_p which is already assembled)
        const double fr = -0.03976;
        MatDuplicate(A_p, MAT_COPY_VALUES, &H_p);
        MatAXPY(H_p, 1.0, C1_p, DIFFERENT_NONZERO_PATTERN);
        MatAXPY(H_p, -fr, B_p, DIFFERENT_NONZERO_PATTERN);

        // D[0] = -B for PEP
        Mat negB;
        MatDuplicate(B_p, MAT_COPY_VALUES, &negB);
        MatScale(negB, -1.0);

        timer.stop("Assembly");
        if (rank == 0) std::cout << "Assembly complete. H = A + C1 - fr*B\n";

        // ============================================================
        // PEP Solve:  D[0] + D[1]*lambda + D[2]*lambda^2 = 0
        //   D[0] = -B,  D[1] = C2,  D[2] = H
        // ============================================================
        timer.start("PEP Solve");

        Mat D[3];
        D[0] = negB;
        D[1] = C2_p;
        D[2] = H_p;

        PEP pep;
        PEPCreate(PETSC_COMM_WORLD, &pep);
        PEPSetOperators(pep, 3, D);

        // ST/KSP for parallel mpiaij
        ST st;
        PEPGetST(pep, &st);
        KSP ksp;
        STGetKSP(st, &ksp);
        KSPSetType(ksp, KSPGMRES);
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCBJACOBI);

        int nev_requested = 7;
        PEPSetDimensions(pep, nev_requested, PETSC_DECIDE, PETSC_DECIDE);
        PEPSetFromOptions(pep);

        if (rank == 0) std::cout << "Solving PEP eigenvalue problem...\n";
        PEPSolve(pep);

        PetscInt its;
        PEPGetIterationNumber(pep, &its);
        if (rank == 0) std::cout << "PEP iterations: " << its << std::endl;

        PetscInt nconv;
        PEPGetConverged(pep, &nconv);
        if (rank == 0) std::cout << "Converged eigenpairs: " << nconv << std::endl;

        timer.stop("PEP Solve");

        // ============================================================
        // Extract eigenvalues and wavefunctions
        // ============================================================
        timer.start("Post-Processing");

        if (nconv > 0) {
            Vec xr, xi;
            MatCreateVecs(B_p, &xr, &xi);

            std::vector<double> eigenvalues_re(nconv), eigenvalues_im(nconv);
            std::vector<std::vector<double>> wavefunctions(nconv);
            std::vector<Vector> wf_vectors(nconv);

            for (int i = 0; i < nconv; ++i) {
                PetscScalar kr, ki_val;
                PEPGetEigenpair(pep, i, &kr, &ki_val, xr, xi);

                double re = PetscRealPart(kr);
                double im = PetscImaginaryPart(kr);
                eigenvalues_re[i] = re;
                eigenvalues_im[i] = im;

                PetscReal error;
                PEPComputeError(pep, i, PEP_ERROR_BACKWARD, &error);

                if (rank == 0) {
                    if (im != 0.0)
                        std::cout << "  E[" << i << "] = " << re*27.21
                                  << " + " << im*27.21 << "i eV"
                                  << "  (err=" << error << ")\n";
                    else
                        std::cout << "  E[" << i << "] = " << re*27.21 << " eV"
                                  << "  (err=" << error << ")\n";
                }

                // Gather wavefunction to all ranks
                gatherVecToLocal(xr, gv, l2g, wavefunctions[i]);
                wf_vectors[i].resize(N);
                for (int j = 0; j < N; ++j)
                    wf_vectors[i][j] = wavefunctions[i][j];
            }

            // VTK output (collective — all ranks must call)
            {
                Dune::VTKWriter<GV> vtkwriter(gv, Dune::VTK::conforming);
                for (int i = 0; i < nconv; ++i) {
                    std::string name = "wavefunction" + std::to_string(i);
                    vtkwriter.addVertexData(wf_vectors[i], name);
                }
                vtkwriter.write("graphene_gqd", Dune::VTK::appendedraw);
            }
            if (rank == 0) std::cout << "VTK output: graphene_gqd.pvtu\n";

            // ============================================================
            // LDOS Calculation (rank 0 output)
            // ============================================================
            if (rank == 0) {
                const double he = 27.2;
                const double eta = 0.001 / he;
                std::vector<double> r_points = {0.0, 37.0, 75.0, 100.0};

                std::ofstream ldos_file("ldos_output.dat");
                ldos_file << "# r_target  LDOS\n";

                for (size_t ri = 0; ri < r_points.size(); ++ri) {
                    double r_target = r_points[ri];
                    double ldos_total = 0.0;

                    for (int n = 0; n < nconv; ++n) {
                        double E_n = eigenvalues_re[n];

                        // Sum |psi|^2 near r_target
                        double sum_psi2 = 0.0;
                        for (auto vt = gv.template begin<dim>(); vt != gv.template end<dim>(); ++vt) {
                            auto pos = vt->geometry().corner(0);
                            double r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
                            if (std::fabs(r - r_target) < 2.0) {
                                int idx = vertexmap.index(*vt);
                                sum_psi2 += wavefunctions[n][idx] * wavefunctions[n][idx];
                            }
                        }

                        // Energy scan +/-50 meV around eigenvalue
                        double E_min = E_n - 0.05;
                        double E_max = E_n + 0.05;
                        double dE = 0.001;
                        for (double E_probe = E_min; E_probe <= E_max; E_probe += dE) {
                            double deltaE = E_probe - E_n;
                            double lorentz = (1.0/M_PI) * (eta / (deltaE*deltaE + eta*eta));
                            ldos_total += sum_psi2 * lorentz * dE;
                        }
                    }

                    double ldos_final = ldos_total * he * he;
                    std::cout << "LDOS at r=" << r_target << " : " << ldos_final << std::endl;
                    ldos_file << r_target << "  " << ldos_final << "\n";
                }
                ldos_file.close();
            }

            // ============================================================
            // Tunneling Current & dI/dV (rank 0)
            // ============================================================
            if (rank == 0) {
                const double he = 27.2;
                const double hbar = 1.0;
                const double eta = 0.001 / he;
                const double S = 0.0;

                // LDOS at center (r=0) for tunneling
                double ldos_center = 0.0;
                for (int n = 0; n < nconv; ++n) {
                    double E_n = eigenvalues_re[n];
                    double sum_psi2 = 0.0;
                    for (auto vt = gv.template begin<dim>(); vt != gv.template end<dim>(); ++vt) {
                        auto pos = vt->geometry().corner(0);
                        double r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
                        if (r < 2.0) {
                            int idx = vertexmap.index(*vt);
                            sum_psi2 += wavefunctions[n][idx] * wavefunctions[n][idx];
                        }
                    }
                    double deltaE = 0.067/he - E_n;
                    double lorentz = (1.0/M_PI) * (eta / (deltaE*deltaE + eta*eta));
                    ldos_center += sum_psi2 * std::exp(-2.0 * S / hbar) * lorentz;
                }
                ldos_center *= he * he;

                // Bias sweep
                std::ofstream iv_file("tunneling_current.dat");
                iv_file << "# V_bias(meV)  I(a.u.)  dI/dV\n";

                double V_start = -700.0, V_end = 0.0, V_step = 50.0;
                int bias_size = (int)((V_end - V_start) / V_step) + 1;
                double dE_bias = V_step / he;
                double conv_factor = 0.067 / 27.2;

                std::vector<double> current(bias_size, 0.0);
                for (int i = 0; i < bias_size; ++i) {
                    for (int j = 0; j <= i; ++j)
                        current[i] += ldos_center * dE_bias;

                    double V_bias = V_start + i * V_step;
                    double dIdV = current[i] / conv_factor;

                    std::cout << "V=" << V_bias << " meV, I=" << current[i]
                              << ", dI/dV=" << dIdV << std::endl;
                    iv_file << V_bias << "  " << current[i] << "  " << dIdV << "\n";
                }
                iv_file.close();
                std::cout << "Tunneling current: tunneling_current.dat\n";
            }

            // Eigenvalues file (rank 0)
            if (rank == 0) {
                std::ofstream ef("eigenvalues.dat");
                ef << "# index  E_re(Ha)  E_im(Ha)  E_re(eV)  E_im(eV)\n";
                for (int i = 0; i < nconv; ++i)
                    ef << i << "  " << eigenvalues_re[i] << "  " << eigenvalues_im[i]
                       << "  " << eigenvalues_re[i]*27.21 << "  " << eigenvalues_im[i]*27.21 << "\n";
                ef.close();
            }

            VecDestroy(&xr);
            VecDestroy(&xi);
        }

        timer.stop("Post-Processing");

        // Cleanup
        PEPDestroy(&pep);
        MatDestroy(&negB);
        ISLocalToGlobalMappingDestroy(&lgmap);
        MatDestroy(&A_p);
        MatDestroy(&B_p);
        MatDestroy(&C1_p);
        MatDestroy(&C2_p);
        MatDestroy(&H_p);

        timer.report();
        timer.reportCSV("timings.csv");
        timer.appendScalingLine("scaling_results.csv");

        SlepcFinalize();

        if (rank == 0)
            std::cout << "\nSimulation complete.\n"
                      << "Output: graphene_gqd.pvtu, eigenvalues.dat, "
                      << "ldos_output.dat, tunneling_current.dat\n";

        return 0;
    }
    catch (Dune::Exception& e) { std::cerr << "DUNE: " << e.what() << "\n"; }
    catch (std::exception& e)  { std::cerr << "Error: " << e.what() << "\n"; }
    catch (...)                 { std::cerr << "Unknown exception!\n"; }
    return 1;
}
