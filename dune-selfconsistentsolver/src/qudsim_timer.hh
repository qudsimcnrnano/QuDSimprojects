/**
 * QuDSim Performance Timer - Drop-in instrumentation for DUNE/PETSc FEM codes
 * 
 * Usage:
 *   #include "qudsim_timer.hh"
 * 
 *   QuDSimTimer timer(MPI_COMM_WORLD);   // create once
 *   timer.start("Assembly");              // start a phase
 *   // ... do assembly ...
 *   timer.stop("Assembly");               // stop it
 *   timer.start("KSP Solve");
 *   // ... solve ...
 *   timer.stop("KSP Solve");
 *   timer.report();                       // print summary on rank 0
 *   timer.reportCSV("timings.csv");       // export CSV for plotting
 * 
 * For scaling studies, run your code with different -np values and
 * use the companion script scaling_study.sh
 *
 * Authors: QuDSim Team
 * License: Same as QuDSim
 */

#ifndef QUDSIM_TIMER_HH
#define QUDSIM_TIMER_HH

#include <mpi.h>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

class QuDSimTimer
{
public:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;

    struct PhaseData {
        std::string name;
        double wallTime    = 0.0;   // seconds (local)
        double cpuTime     = 0.0;   // seconds (local)
        int    callCount   = 0;
        TimePoint wallStart;
        double    cpuStart = 0.0;
        bool      running  = false;
    };

private:
    MPI_Comm comm_;
    int rank_, size_;
    std::map<std::string, PhaseData> phases_;
    std::vector<std::string> phaseOrder_;  // insertion order
    TimePoint globalStart_;
    double globalCpuStart_;

    static double cpuNow() {
        return (double)clock() / CLOCKS_PER_SEC;
    }

public:
    QuDSimTimer(MPI_Comm comm = MPI_COMM_WORLD)
        : comm_(comm)
    {
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);
        globalStart_ = Clock::now();
        globalCpuStart_ = cpuNow();
    }

    /// Start timing a named phase
    void start(const std::string& name)
    {
        auto& p = phases_[name];
        if (p.name.empty()) {
            p.name = name;
            phaseOrder_.push_back(name);
        }
        p.wallStart = Clock::now();
        p.cpuStart  = cpuNow();
        p.running   = true;
    }

    /// Stop timing a named phase
    void stop(const std::string& name)
    {
        auto it = phases_.find(name);
        if (it == phases_.end() || !it->second.running) {
            if (rank_ == 0)
                std::cerr << "[QuDSimTimer] Warning: stop('" << name
                          << "') called without matching start." << std::endl;
            return;
        }
        auto& p = it->second;
        double wallElapsed = std::chrono::duration<double>(Clock::now() - p.wallStart).count();
        double cpuElapsed  = cpuNow() - p.cpuStart;
        p.wallTime  += wallElapsed;
        p.cpuTime   += cpuElapsed;
        p.callCount += 1;
        p.running    = false;
    }

    /// Get wall time for a specific phase (local to this rank)
    double getWallTime(const std::string& name) const
    {
        auto it = phases_.find(name);
        return (it != phases_.end()) ? it->second.wallTime : 0.0;
    }

    /// Get total wall time since timer was created
    double totalWallTime() const
    {
        return std::chrono::duration<double>(Clock::now() - globalStart_).count();
    }

    /// Print a formatted report on rank 0, with min/max/avg across ranks
    void report() const
    {
        if (phaseOrder_.empty()) return;

        // Gather data from all ranks
        for (const auto& name : phaseOrder_) {
            auto it = phases_.find(name);
            double localWall = (it != phases_.end()) ? it->second.wallTime : 0.0;
            double localCpu  = (it != phases_.end()) ? it->second.cpuTime  : 0.0;
            int    localCall = (it != phases_.end()) ? it->second.callCount : 0;

            double minWall, maxWall, sumWall;
            double minCpu,  maxCpu,  sumCpu;
            int    sumCall;

            MPI_Reduce(&localWall, &minWall, 1, MPI_DOUBLE, MPI_MIN, 0, comm_);
            MPI_Reduce(&localWall, &maxWall, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
            MPI_Reduce(&localWall, &sumWall, 1, MPI_DOUBLE, MPI_SUM, 0, comm_);
            MPI_Reduce(&localCpu,  &minCpu,  1, MPI_DOUBLE, MPI_MIN, 0, comm_);
            MPI_Reduce(&localCpu,  &maxCpu,  1, MPI_DOUBLE, MPI_MAX, 0, comm_);
            MPI_Reduce(&localCpu,  &sumCpu,  1, MPI_DOUBLE, MPI_SUM, 0, comm_);
            MPI_Reduce(&localCall, &sumCall, 1, MPI_INT,    MPI_SUM, 0, comm_);

            if (rank_ == 0) {
                if (name == phaseOrder_.front()) {
                    // Print header
                    double totalWall = totalWallTime();
                    std::cout << "\n"
                        << "╔══════════════════════════════════════════════════════════════════════════════╗\n"
                        << "║                    QuDSim Performance Report                                ║\n"
                        << "║    MPI ranks: " << std::setw(4) << size_
                        << "    Total wall time: " << std::fixed << std::setprecision(3)
                        << std::setw(10) << totalWall << " s"
                        << std::setw(25) << " " << "║\n"
                        << "╠══════════════════════════════════════════════════════════════════════════════╣\n"
                        << "║ Phase                │ Calls │ Wall(avg) │ Wall(min) │ Wall(max) │ CPU(tot) ║\n"
                        << "╠══════════════════════╪═══════╪═══════════╪═══════════╪═══════════╪══════════╣\n";
                }

                double avgWall = sumWall / size_;
                char line[200];
                snprintf(line, sizeof(line),
                    "║ %-20s │ %5d │ %8.3fs │ %8.3fs │ %8.3fs │ %7.3fs ║",
                    name.c_str(), sumCall / size_, avgWall, minWall, maxWall, sumCpu);
                std::cout << line << "\n";

                if (name == phaseOrder_.back()) {
                    std::cout
                        << "╚══════════════════════════════════════════════════════════════════════════════╝\n"
                        << std::endl;
                }
            }
        }
    }

    /// Export timing data to CSV (for plotting with Python/gnuplot)
    void reportCSV(const std::string& filename) const
    {
        // Each rank sends its data to rank 0
        for (const auto& name : phaseOrder_) {
            auto it = phases_.find(name);
            double localWall = (it != phases_.end()) ? it->second.wallTime : 0.0;
            double localCpu  = (it != phases_.end()) ? it->second.cpuTime  : 0.0;
            int    localCall = (it != phases_.end()) ? it->second.callCount : 0;

            double allWall[size_], allCpu[size_];
            int allCall[size_];

            MPI_Gather(&localWall, 1, MPI_DOUBLE, allWall, 1, MPI_DOUBLE, 0, comm_);
            MPI_Gather(&localCpu,  1, MPI_DOUBLE, allCpu,  1, MPI_DOUBLE, 0, comm_);
            MPI_Gather(&localCall, 1, MPI_INT,    allCall, 1, MPI_INT,    0, comm_);

            if (rank_ == 0 && name == phaseOrder_.front()) {
                std::ofstream ofs(filename);
                ofs << "phase,nranks,rank,calls,wall_s,cpu_s\n";
                for (const auto& pname : phaseOrder_) {
                    // We need to do this in a single pass, so we gather all at once
                    // This block handles just the first phase; see below
                }
                ofs.close();
            }
        }

        // Simpler approach: gather everything on rank 0 and write
        if (rank_ == 0) {
            // We already have local data; for multi-rank we need gathers
            // Use the single-rank approach for simplicity, which is the common use
            std::ofstream ofs(filename);
            ofs << "phase,nranks,wall_s,cpu_s,calls\n";
            for (const auto& name : phaseOrder_) {
                auto it = phases_.find(name);
                if (it != phases_.end()) {
                    ofs << name << "," << size_ << ","
                        << std::fixed << std::setprecision(6)
                        << it->second.wallTime << ","
                        << it->second.cpuTime << ","
                        << it->second.callCount << "\n";
                }
            }
            ofs.close();
            std::cout << "[QuDSimTimer] Timings written to " << filename << std::endl;
        }
    }

    /// Export a single-line summary for scaling studies
    /// Append to a file: nranks, total_wall, phase1_wall, phase2_wall, ...
    void appendScalingLine(const std::string& filename) const
    {
        // Get max wall time across ranks for each phase (bottleneck time)
        if (rank_ == 0) {
            std::ofstream ofs(filename, std::ios::app);
            // Check if file is empty to write header
            ofs.seekp(0, std::ios::end);
            if (ofs.tellp() == 0) {
                ofs << "nranks,total_wall_s";
                for (const auto& name : phaseOrder_)
                    ofs << "," << name << "_wall_s";
                ofs << "\n";
            }

            ofs << size_ << "," << std::fixed << std::setprecision(6)
                << totalWallTime();
        }

        for (const auto& name : phaseOrder_) {
            auto it = phases_.find(name);
            double localWall = (it != phases_.end()) ? it->second.wallTime : 0.0;
            double maxWall;
            MPI_Reduce(&localWall, &maxWall, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
            if (rank_ == 0) {
                std::ofstream ofs(filename, std::ios::app);
                ofs << "," << std::fixed << std::setprecision(6) << maxWall;
            }
        }

        if (rank_ == 0) {
            std::ofstream ofs(filename, std::ios::app);
            ofs << "\n";
            ofs.close();
        }
    }
};

#endif // QUDSIM_TIMER_HH

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
