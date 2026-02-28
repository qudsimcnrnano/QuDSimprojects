// =============================================================================
// QDSimTimer.hh - MPI-Aware Simulation Timer for GAA Nanowire FET
// Profiles: Grid setup, Assembly, BC application, Eigenvalue solve,
//           Poisson solve, Communication, Total runtime
// Compatible with scaling.sh and plot_scaling.py for performance analysis
// =============================================================================

#ifndef QDSIMTIMER_HH
#define QDSIMTIMER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <mpi.h>

class QDSimTimer
{
public:

  // Timer entry: stores accumulated time and call count
  struct TimerEntry {
    double total_time;
    double start_time;
    int call_count;
    bool running;

    TimerEntry() : total_time(0.0), start_time(0.0), call_count(0), running(false) {}
  };

private:
  std::map<std::string, TimerEntry> timers;
  double wall_start;
  int rank;
  int nprocs;

  // Order in which timers were first created (for consistent output)
  std::vector<std::string> timer_order;

public:

  QDSimTimer()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    wall_start = MPI_Wtime();
  }

  // Start a named timer
  void start(const std::string& name)
  {
    if (timers.find(name) == timers.end()) {
      timers[name] = TimerEntry();
      timer_order.push_back(name);
    }
    timers[name].start_time = MPI_Wtime();
    timers[name].running = true;
  }

  // Stop a named timer
  void stop(const std::string& name)
  {
    if (timers.find(name) != timers.end() && timers[name].running) {
      double elapsed = MPI_Wtime() - timers[name].start_time;
      timers[name].total_time += elapsed;
      timers[name].call_count++;
      timers[name].running = false;
    }
  }

  // Get elapsed time for a timer
  double getTime(const std::string& name) const
  {
    auto it = timers.find(name);
    if (it != timers.end())
      return it->second.total_time;
    return 0.0;
  }

  // Get total wall time since construction
  double getTotalWallTime() const
  {
    return MPI_Wtime() - wall_start;
  }

  // Print summary table (rank 0 only, with min/max/avg across ranks)
  void printSummary() const
  {
    double total_wall = getTotalWallTime();

    if (rank == 0) {
      std::cout << "\n";
      std::cout << "================================================================\n";
      std::cout << "  TIMING SUMMARY (" << nprocs << " MPI processes)\n";
      std::cout << "================================================================\n";
      std::cout << std::left << std::setw(28) << "  Phase"
                << std::right << std::setw(10) << "Min(s)"
                << std::setw(10) << "Max(s)"
                << std::setw(10) << "Avg(s)"
                << std::setw(8) << "Calls"
                << std::setw(8) << "%Total"
                << "\n";
      std::cout << "  " << std::string(72, '-') << "\n";
    }

    for (const auto& name : timer_order) {
      auto it = timers.find(name);
      if (it == timers.end()) continue;

      double local_time = it->second.total_time;
      int local_calls = it->second.call_count;

      double min_time, max_time, sum_time;
      int total_calls;

      MPI_Reduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&local_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&local_calls, &total_calls, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      if (rank == 0) {
        double avg_time = sum_time / nprocs;
        double pct = (max_time / total_wall) * 100.0;
        std::cout << "  " << std::left << std::setw(26) << name
                  << std::right << std::fixed << std::setprecision(4)
                  << std::setw(10) << min_time
                  << std::setw(10) << max_time
                  << std::setw(10) << avg_time
                  << std::setw(8) << (total_calls / nprocs)
                  << std::setw(7) << std::setprecision(1) << pct << "%"
                  << "\n";
      }
    }

    // Total wall time
    double min_wall, max_wall, sum_wall;
    MPI_Reduce(&total_wall, &min_wall, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_wall, &max_wall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_wall, &sum_wall, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      std::cout << "  " << std::string(72, '-') << "\n";
      std::cout << "  " << std::left << std::setw(26) << "TOTAL WALL TIME"
                << std::right << std::fixed << std::setprecision(4)
                << std::setw(10) << min_wall
                << std::setw(10) << max_wall
                << std::setw(10) << (sum_wall / nprocs)
                << std::setw(8) << ""
                << std::setw(7) << "100%"
                << "\n";
      std::cout << "================================================================\n\n";
    }
  }

  // Write timing data to file for scaling analysis
  // Format: nprocs, total_time, grid_time, assembly_time, solve_time, comm_time
  void writeScalingData(const std::string& filename = "scaling_data.csv") const
  {
    double total_wall = getTotalWallTime();

    // Collect max times across ranks (bottleneck = slowest rank)
    double max_wall;
    MPI_Reduce(&total_wall, &max_wall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      // Check if file exists to write header
      bool write_header = true;
      {
        std::ifstream check(filename);
        if (check.good()) {
          std::string line;
          std::getline(check, line);
          if (line.find("nprocs") != std::string::npos)
            write_header = false;
        }
      }

      std::ofstream out(filename, std::ios::app);
      if (write_header) {
        out << "nprocs";
        for (const auto& name : timer_order)
          out << "," << name;
        out << ",total_wall_time\n";
      }

      out << nprocs;
      for (const auto& name : timer_order) {
        auto it = timers.find(name);
        double local_time = (it != timers.end()) ? it->second.total_time : 0.0;
        double max_time;
        MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        out << "," << std::fixed << std::setprecision(6) << max_time;
      }
      out << "," << std::fixed << std::setprecision(6) << max_wall << "\n";
      out.close();
    } else {
      // Non-root ranks still participate in MPI_Reduce
      for (const auto& name : timer_order) {
        auto it = timers.find(name);
        double local_time = (it != timers.end()) ? it->second.total_time : 0.0;
        double max_time;
        MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      }
    }
  }
};

#endif // QDSIMTIMER_HH
