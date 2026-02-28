/**
 * QuDSim - Schrodinger FEM Solver using SLEPc
 * Solves: -1/(2*m_eff) * laplacian(psi) + V(x)*psi = E*psi
 * Generalized eigenvalue: (H + V*M)*psi = E*M*psi
 * Uses DUNE P1 FEM + SLEPc EPS
 * Authors: QuDSim Team
 */
#ifndef QUDSIM_SCHRODINGERFEM_HH
#define QUDSIM_SCHRODINGERFEM_HH

#include <vector>
#include <set>
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
#include "petscmat.h"
#include "petscvec.h"
#include "slepceps.h"

namespace QuDSim {

template<class GV, class Vector>
class SchrodingerFEM
{
public:
    static const int dim = GV::dimension;
    typedef typename GV::ctype ctype;

private:
    using IndexSet      = typename GV::IndexSet;
    using ElementMapper = Dune::SingleCodimSingleGeomTypeMapper<GV, 0>;

    const GV&     gv_;
    Vector&       carrier_density_;
    const Vector& potential_;
    const Vector& wf_bc_;
    int           carrier_type_;  // 0=electron, 1=hole
    double        Ef_;
    const Vector& emass_;
    const Vector& regionid_;
    int           nconv_;
    int           num_requested_;
    std::vector<double>  eigenvalues_;
    std::vector<Vector>  wavefunctions_;
    double beta_;

public:
    SchrodingerFEM(const GV& gv, Vector& carrier_density,
                   const Vector& potential, const Vector& wf_bc,
                   int carrier_type, double Ef,
                   const Vector& emass, const Vector& regionid,
                   int num_eigenvalues = 10, double T = 300.0)
        : gv_(gv), carrier_density_(carrier_density), potential_(potential),
          wf_bc_(wf_bc), carrier_type_(carrier_type), Ef_(Ef),
          emass_(emass), regionid_(regionid),
          nconv_(0), num_requested_(num_eigenvalues)
    {
        beta_ = Constants::Hartree_to_eV / (Constants::kB_eV * T);
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

        // Create H and M matrices
        Mat H, M;
        std::vector<PetscInt> nnz(N);
        for (int i = 0; i < N; ++i) nnz[i] = (PetscInt)adj[i].size();
        MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 0, nnz.data(), &H);
        MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 0, nnz.data(), &M);
        MatSetOption(H, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        // Element assembly
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto geo = it->geometry();
            auto ref = Dune::referenceElement(geo);
            int nv = ref.size(dim);
            int elemIdx = elementmap.index(*it);
            double meff = emass_[elemIdx];
            if (meff <= 0.0) meff = 1.0;

            std::vector<PetscInt> rows(nv);
            for (int i = 0; i < nv; ++i)
                rows[i] = (PetscInt)set.subIndex(*it, i, dim);

            std::vector<PetscScalar> He(nv*nv, 0.0), Me(nv*nv, 0.0);

            // Kinetic energy (stiffness)
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
                        He[i*nv+j] += (1.0/(2.0*meff)) * (gi*gj) * w * detJ;
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
                    V_qp += phi * potential_[set.subIndex(*it, i, dim)];
                }
                for (int i = 0; i < nv; ++i) {
                    ctype phi_i = basis[i].evaluateFunction(r->position());
                    for (int j = 0; j < nv; ++j) {
                        ctype phi_j = basis[j].evaluateFunction(r->position());
                        He[i*nv+j] += V_qp * phi_i * phi_j * w * detJ;
                        Me[i*nv+j] += phi_i * phi_j * w * detJ;
                    }
                }
            }

            MatSetValues(H, nv, rows.data(), nv, rows.data(), He.data(), ADD_VALUES);
            MatSetValues(M, nv, rows.data(), nv, rows.data(), Me.data(), ADD_VALUES);
        }

        MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

        // BC: psi=0 at all domain boundaries
        std::vector<PetscInt> bcRows;
        for (auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it) {
            auto ref = Dune::referenceElement(it->geometry());
            for (auto is = gv_.ibegin(*it); is != gv_.iend(*it); ++is) {
                if (!is->boundary()) continue;
                int nv_f = ref.size(is->indexInInside(), 1, dim);
                for (int k = 0; k < nv_f; ++k) {
                    int lv = ref.subEntity(is->indexInInside(), 1, k, dim);
                    bcRows.push_back(set.subIndex(*it, lv, dim));
                }
            }
        }
        std::sort(bcRows.begin(), bcRows.end());
        bcRows.erase(std::unique(bcRows.begin(), bcRows.end()), bcRows.end());

        if (!bcRows.empty()) {
            PetscScalar bigVal = 1.0e6;
            MatZeroRows(H, (PetscInt)bcRows.size(), bcRows.data(), bigVal, NULL, NULL);
            MatZeroRows(M, (PetscInt)bcRows.size(), bcRows.data(), (PetscScalar)1.0, NULL, NULL);
        }

        // SLEPc eigenvalue solve
        EPS eps;
        EPSCreate(PETSC_COMM_SELF, &eps);
        EPSSetOperators(eps, H, M);
        EPSSetProblemType(eps, EPS_GHEP);
        EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
        EPSSetDimensions(eps, num_requested_, PETSC_DEFAULT, PETSC_DEFAULT);
        EPSSetType(eps, EPSKRYLOVSCHUR);
        EPSSetTolerances(eps, 1e-8, 1000);
        EPSSetFromOptions(eps);
        EPSSolve(eps);

        EPSGetConverged(eps, &nconv_);
        std::cout << "  Schrodinger (" << (carrier_type_==0?"elec":"hole")
                  << "): " << nconv_ << " eigenvalues" << std::endl;

        // Extract results
        eigenvalues_.resize(nconv_);
        wavefunctions_.resize(nconv_);
        Vec xr, xi;
        MatCreateVecs(H, NULL, &xr);
        VecDuplicate(xr, &xi);

        for (int i = 0; i < nconv_; ++i) {
            PetscScalar kr, ki;
            EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
            eigenvalues_[i] = PetscRealPart(kr);
            wavefunctions_[i].resize(gv_.size(dim));
            wavefunctions_[i] = 0.0;
            const PetscScalar* arr;
            VecGetArrayRead(xr, &arr);
            for (int j = 0; j < (int)gv_.size(dim); ++j)
                wavefunctions_[i][j] = PetscRealPart(arr[j]);
            VecRestoreArrayRead(xr, &arr);
            if (i < 5)
                std::cout << "    E[" << i << "] = " << eigenvalues_[i]
                          << " Ha = " << Ha_to_eV(eigenvalues_[i]) << " eV" << std::endl;
        }

        // Carrier density
        computeCarrierDensity();

        VecDestroy(&xr); VecDestroy(&xi);
        EPSDestroy(&eps);
        MatDestroy(&H); MatDestroy(&M);
    }

    void computeCarrierDensity()
    {
        const int N = gv_.size(dim);
        carrier_density_ = 0.0;
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
                double psi = wavefunctions_[i][j];
                carrier_density_[j] += 2.0 * occ * psi * psi;
            }
        }
        double maxn = cm3_to_au3(1e25);
        for (int j = 0; j < N; ++j)
            if (carrier_density_[j] > maxn) carrier_density_[j] = maxn;
    }

    int getnconv() const { return nconv_; }
    double geteigen(int i) const { return (i<nconv_) ? eigenvalues_[i] : 0.0; }
    const Vector& getwave(int i) const { return wavefunctions_[i]; }
    const std::vector<double>& getEigenvalues() const { return eigenvalues_; }
};

} // namespace QuDSim
#endif
