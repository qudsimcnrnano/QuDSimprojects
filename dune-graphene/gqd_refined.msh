#ifndef QUDSIM_TIMER_HH
#define QUDSIM_TIMER_HH

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <mpi.h>

class QuDSimTimer {
    struct TimerEntry {
        double wall_start = 0;
        double wall_total = 0;
        double cpu_total  = 0;
        double wall_min   = 1e30;
        double wall_max   = 0;
        int    calls      = 0;
    };
    std::map<std::string, TimerEntry> timers_;
    std::vector<std::string> order_;
    MPI_Comm comm_;
    double global_start_;
    int rank_, size_;

public:
    QuDSimTimer(MPI_Comm comm = MPI_COMM_WORLD) : comm_(comm) {
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);
        global_start_ = MPI_Wtime();
    }

    void start(const std::string& name) {
        auto& t = timers_[name];
        if (t.calls == 0) order_.push_back(name);
        t.wall_start = MPI_Wtime();
    }

    void stop(const std::string& name) {
        double elapsed = MPI_Wtime() - timers_[name].wall_start;
        auto& t = timers_[name];
        t.wall_total += elapsed;
        t.cpu_total  += elapsed;
        t.wall_min = std::min(t.wall_min, elapsed);
        t.wall_max = std::max(t.wall_max, elapsed);
        t.calls++;
    }

    void report() {
        double total_wall = MPI_Wtime() - global_start_;
        if (rank_ != 0) return;
        std::cout << "\n";
        std::cout << std::string(78, '=') << "\n";
        std::cout << "  QuDSim Performance Report  |  MPI ranks: " << size_
                  << "  |  Total wall time: " << std::fixed << std::setprecision(3)
                  << total_wall << " s\n";
        std::cout << std::string(78, '=') << "\n";
        std::cout << std::left << std::setw(24) << " Phase"
                  << std::setw(8) << "Calls"
                  << std::setw(12) << "Wall(avg)"
                  << std::setw(12) << "Wall(min)"
                  << std::setw(12) << "Wall(max)" << "\n";
        std::cout << std::string(78, '-') << "\n";
        for (auto& name : order_) {
            auto& t = timers_[name];
            double avg = t.wall_total / std::max(t.calls, 1);
            std::cout << std::left << std::setw(24) << (" " + name)
                      << std::setw(8) << t.calls
                      << std::setw(12) << (std::to_string(avg).substr(0,7) + "s")
                      << std::setw(12) << (std::to_string(t.wall_min).substr(0,7) + "s")
                      << std::setw(12) << (std::to_string(t.wall_max).substr(0,7) + "s")
                      << "\n";
        }
        std::cout << std::string(78, '=') << "\n\n";
    }

    void reportCSV(const std::string& filename) {
        if (rank_ != 0) return;
        std::ofstream f(filename);
        f << "phase,calls,wall_avg,wall_min,wall_max\n";
        for (auto& name : order_) {
            auto& t = timers_[name];
            double avg = t.wall_total / std::max(t.calls, 1);
            f << name << "," << t.calls << "," << avg << "," << t.wall_min << "," << t.wall_max << "\n";
        }
        f.close();
        if (rank_ == 0) std::cout << "[QuDSimTimer] Timings written to " << filename << "\n";
    }

    void appendScalingLine(const std::string& filename) {
        double total_wall = MPI_Wtime() - global_start_;
        if (rank_ != 0) return;
        bool exists = std::ifstream(filename).good();
        std::ofstream f(filename, std::ios::app);
        if (!exists) f << "nranks,total_wall_s\n";
        f << size_ << "," << total_wall << "\n";
        f.close();
    }
};

#endif // QUDSIM_TIMER_HH
