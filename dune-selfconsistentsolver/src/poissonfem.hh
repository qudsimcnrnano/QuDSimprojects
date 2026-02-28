/**
 * QuDSim - Poisson FEM Solver (Newton-Raphson) for MOSCAP
 * Solves: -div(eps * grad(phi)) = rho(phi)
 * NR linearization: (K - dRho*M)*delta_phi = -(K*phi + rho)
 * Uses DUNE P1 FEM + PETSc KSP
 * Authors: QuDSim Team
 */
#ifndef QUDSIM_POISSONFEM_HH
#define QUDSIM_POISSONFEM_HH

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include </opt/local/dune_210/dune-grid-2.10.dev20221009/dune/grid/common/scsgmapper.hh>

#include "shapefunctions.hh"
#include "qudsim_physics.hh"

#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

namespace QuDSim {

template<class GV, class Vector>
class PoissonFEM
{
public:
    static const int dim = GV::dimension;
    typedef typename GV::ctype ctype;

private:
    using IndexSet      = typename GV::IndexSet;
    using ElementMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;
    using VertexMapper  = Dune::SingleCodimSingleGeomTypeMapper<GV, dim>;

    const GV&     gv_;
    Vector&       charge_;
    Vector&       potential_;
    double        Vss_;
    Vector&       deltaRho_;
    Vector&       deltaPotential_;
    int           it0_;
    const Vector& epsi_;
    const Vector& regionid_;
    const Vector& boundaryid_;

    Mat A_;
    Vec b_, x_;

public:
    PoissonFEM(const GV& gv, Vector& charge, Vector& potential,
               double Vss, Vector& deltaRho, Vector& deltaPotential,
               int it0, const Vector& epsi, const Vector& regionid,
               const Vector& boundaryid)
        : gv_(gv), charge_(charge), potential_(potential),
          Vss_(Vss), deltaRho_(deltaRho), deltaPotential_(deltaPotential),
          it0_(it0), epsi_(epsi), regionid_(regionid), boundaryid_(boundaryid),
          A_(nullptr), b_(nullptr), x_(nullptr)
    {}

    ~PoissonFEM() {
        if (A_) MatDestroy(&A_);
        if (b_) VecDestroy(&b_);
        if (x_) VecDestroy(&x_);
    }

    void apply()
    {
        const IndexSet& set = gv_.indexSet();
        const int N = gv_.size(dim);
        P1ShapeFunctionSet<ctype,ctype,dim> basis =
            P1ShapeFunctionSet<ctype,ctype,dim>::instance();
        ElementMapper elementmap(gv_);

        // Sparsity
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

        // PETSc setup
        if (A_) MatDestroy(&A_);
        if (b_) VecDestroy(&b_);
        if (x_) VecDestroy(&x_);

        std::vector<PetscInt> nnz(N);
        for (int i = 0; i < N; ++i) nnz[i] = (PetscInt)adj[i].size();
        MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 0, nnz.data(), &A_);
        MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        VecCreateSeq(PETSC_COMM_SELF, N, &b_);
        VecDuplicate(b_, &x_);
        VecSet(b_, 0.0);

        // Element assembly
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto geo = it->geometry();
            auto ref = Dune::referenceElement(geo);
            int nv = ref.size(dim);
            int elemIdx = elementmap.index(*it);
            double eps_elem = epsi_[elemIdx];

            std::vector<PetscInt> rows(nv);
            for (int i = 0; i < nv; ++i)
                rows[i] = (PetscInt)set.subIndex(*it, i, dim);

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
                        Ke[i*nv+j] += eps_elem * (gi * gj) * w * detJ;
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
                    rho_qp  += phi * charge_[idx];
                    dRho_qp += phi * deltaRho_[idx];
                }
                for (int i = 0; i < nv; ++i) {
                    ctype phi_i = basis[i].evaluateFunction(r->position());
                    be[i] -= phi_i * rho_qp * w * detJ;
                    for (int j = 0; j < nv; ++j) {
                        ctype phi_j = basis[j].evaluateFunction(r->position());
                        Ke[i*nv+j] -= dRho_qp * phi_i * phi_j * w * detJ;
                    }
                }
            }

            MatSetValues(A_, nv, rows.data(), nv, rows.data(), Ke.data(), ADD_VALUES);
            VecSetValues(b_, nv, rows.data(), be.data(), ADD_VALUES);
        }

        MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b_);
        VecAssemblyEnd(b_);

        // Boundary conditions
        std::vector<PetscInt> bcRows;
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto ref = Dune::referenceElement(it->geometry());
            for (auto is = gv_.ibegin(*it); is != gv_.iend(*it); ++is) {
                if (!is->boundary()) continue;
                int nv_f = ref.size(is->indexInInside(), 1, dim);
                for (int k = 0; k < nv_f; ++k) {
                    int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
                    int idx = set.subIndex(*it, lv, dim);
                    auto pos = it->geometry().corner(lv);
                    double x = pos[0];
                    if (x <= -6.9) {
                        bcRows.push_back(idx);
                        if (it0_ == 0)
                            VecSetValue(b_, idx, (PetscScalar)(Vss_ - (double)potential_[idx]), INSERT_VALUES);
                        else
                            VecSetValue(b_, idx, (PetscScalar)0.0, INSERT_VALUES);
                    } else if (x >= 49.9) {
                        bcRows.push_back(idx);
                        if (it0_ == 0)
                            VecSetValue(b_, idx, (PetscScalar)(0.0 - (double)potential_[idx]), INSERT_VALUES);
                        else
                            VecSetValue(b_, idx, (PetscScalar)0.0, INSERT_VALUES);
                    }
                }
            }
        }
        std::sort(bcRows.begin(), bcRows.end());
        bcRows.erase(std::unique(bcRows.begin(), bcRows.end()), bcRows.end());

        VecAssemblyBegin(b_); VecAssemblyEnd(b_);
        if (!bcRows.empty())
            MatZeroRows(A_, (PetscInt)bcRows.size(), bcRows.data(), (PetscScalar)1.0, NULL, NULL);

        // Solve
        KSP ksp;
        KSPCreate(PETSC_COMM_SELF, &ksp);
        KSPSetOperators(ksp, A_, A_);
        KSPSetType(ksp, KSPGMRES);
        PC pc; KSPGetPC(ksp, &pc); PCSetType(pc, PCILU);
        KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, b_, x_);

        PetscInt its; KSPConvergedReason reason;
        KSPGetIterationNumber(ksp, &its);
        KSPGetConvergedReason(ksp, &reason);
        std::cout << "  Poisson KSP: " << its << " iters (reason=" << reason << ")" << std::endl;

        const PetscScalar* arr;
        VecGetArrayRead(x_, &arr);
        for (int i = 0; i < N; ++i) deltaPotential_[i] = PetscRealPart(arr[i]);
        VecRestoreArrayRead(x_, &arr);
        KSPDestroy(&ksp);
    }
};

} // namespace QuDSim
#endif
