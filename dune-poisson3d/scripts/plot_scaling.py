#!/usr/bin/env python3
"""
QuDSim Scaling Study Plotter
Generates publication-quality plots for IEEE TCAD submission.

Usage:
    python3 plot_scaling.py                           # default: scaling_results.csv
    python3 plot_scaling.py --input my_results.csv
    python3 plot_scaling.py --phase-csv timings.csv   # per-phase breakdown

Output:
    scaling_speedup.pdf      - Speedup vs. #cores
    scaling_efficiency.pdf   - Parallel efficiency vs. #cores
    scaling_walltime.pdf     - Wall time vs. #cores
    phase_breakdown.pdf      - Per-phase time breakdown (if --phase-csv given)
"""

import argparse
import numpy as np
import sys

try:
    import matplotlib
    matplotlib.use('Agg')  # non-interactive backend for HPC
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
except ImportError:
    print("ERROR: matplotlib not found. Install with: pip install matplotlib")
    sys.exit(1)

# IEEE-compatible figure style
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (3.5, 2.8),  # IEEE single-column width
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'lines.linewidth': 1.5,
    'lines.markersize': 6,
    'axes.grid': True,
    'grid.alpha': 0.3,
})


def load_scaling_data(filename):
    """Load scaling_results.csv: nranks,run,wall_time_s"""
    data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=None, encoding='utf-8')
    
    results = {}  # nranks -> list of wall times
    for row in data:
        np_val = int(row[0])
        wall = float(row[2])
        if np_val not in results:
            results[np_val] = []
        results[np_val].append(wall)
    
    # Take minimum wall time for each core count
    nranks = sorted(results.keys())
    best_times = [min(results[n]) for n in nranks]
    avg_times = [np.mean(results[n]) for n in nranks]
    std_times = [np.std(results[n]) for n in nranks]
    
    return np.array(nranks), np.array(best_times), np.array(avg_times), np.array(std_times)


def load_phase_data(filename):
    """Load per-phase timing CSV: phase,nranks,wall_s,cpu_s,calls"""
    phases = {}
    with open(filename) as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            phase = parts[0]
            nranks = int(parts[1])
            wall = float(parts[2])
            cpu = float(parts[3])
            if phase not in phases:
                phases[phase] = {'nranks': [], 'wall': [], 'cpu': []}
            phases[phase]['nranks'].append(nranks)
            phases[phase]['wall'].append(wall)
            phases[phase]['cpu'].append(cpu)
    return phases


def plot_speedup(nranks, best_times, output='scaling_speedup.pdf'):
    """Speedup vs. number of cores"""
    T1 = best_times[0]
    speedup = T1 / best_times
    ideal = nranks / nranks[0]
    
    fig, ax = plt.subplots()
    ax.plot(nranks, ideal, 'k--', label='Ideal', alpha=0.5)
    ax.plot(nranks, speedup, 'bo-', label='QuDSim', markerfacecolor='white')
    ax.set_xlabel('Number of MPI Ranks')
    ax.set_ylabel('Speedup')
    ax.set_title('Strong Scaling: Speedup')
    ax.legend()
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig(output)
    plt.close()
    print(f"  Saved: {output}")


def plot_efficiency(nranks, best_times, output='scaling_efficiency.pdf'):
    """Parallel efficiency vs. number of cores"""
    T1 = best_times[0]
    efficiency = (T1 / best_times) / (nranks / nranks[0]) * 100
    
    fig, ax = plt.subplots()
    ax.plot(nranks, efficiency, 'rs-', markerfacecolor='white')
    ax.axhline(y=100, color='k', linestyle='--', alpha=0.3, label='Ideal (100%)')
    ax.axhline(y=80, color='gray', linestyle=':', alpha=0.3, label='80% threshold')
    ax.set_xlabel('Number of MPI Ranks')
    ax.set_ylabel('Parallel Efficiency (%)')
    ax.set_title('Strong Scaling: Efficiency')
    ax.set_ylim(0, 110)
    ax.legend()
    ax.set_xscale('log', base=2)
    plt.savefig(output)
    plt.close()
    print(f"  Saved: {output}")


def plot_walltime(nranks, best_times, avg_times, std_times, output='scaling_walltime.pdf'):
    """Wall time vs. number of cores"""
    fig, ax = plt.subplots()
    ax.errorbar(nranks, avg_times, yerr=std_times, fmt='go-', capsize=3,
                markerfacecolor='white', label='Avg ± std')
    ax.plot(nranks, best_times, 'b^--', markerfacecolor='white', label='Best')
    
    # Ideal scaling line
    T1 = best_times[0]
    ideal = T1 * nranks[0] / nranks
    ax.plot(nranks, ideal, 'k:', alpha=0.4, label='Ideal')
    
    ax.set_xlabel('Number of MPI Ranks')
    ax.set_ylabel('Wall Time (s)')
    ax.set_title('Strong Scaling: Wall Time')
    ax.legend()
    ax.set_xscale('log', base=2)
    ax.set_yscale('log')
    plt.savefig(output)
    plt.close()
    print(f"  Saved: {output}")


def plot_phase_breakdown(phases, output='phase_breakdown.pdf'):
    """Stacked bar chart of per-phase times"""
    phase_names = list(phases.keys())
    
    fig, ax = plt.subplots(figsize=(4.5, 3.0))
    colors = ['#2196F3', '#FF9800', '#4CAF50', '#F44336', '#9C27B0', '#795548']
    
    x = np.arange(len(phase_names))
    walls = [phases[p]['wall'][0] for p in phase_names]
    cpus = [phases[p]['cpu'][0] for p in phase_names]
    
    width = 0.35
    ax.bar(x - width/2, walls, width, label='Wall Time', color=colors[0], edgecolor='white')
    ax.bar(x + width/2, cpus,  width, label='CPU Time',  color=colors[1], edgecolor='white')
    
    ax.set_xlabel('Solver Phase')
    ax.set_ylabel('Time (s)')
    ax.set_title('Time Breakdown by Phase')
    ax.set_xticks(x)
    ax.set_xticklabels(phase_names, rotation=30, ha='right')
    ax.legend()
    plt.savefig(output)
    plt.close()
    print(f"  Saved: {output}")


def plot_grid_convergence(grid_sizes, eigenvalues, times, output='grid_convergence.pdf'):
    """
    Grid refinement convergence study.
    Call this function from your own script with measured data, e.g.:
    
        grid_sizes = [1000, 5000, 10000, 50000, 100000]
        eigenvalues = [0.2345, 0.2351, 0.2353, 0.2354, 0.2354]
        times = [0.5, 2.1, 8.3, 45.2, 210.5]
        plot_grid_convergence(grid_sizes, eigenvalues, times)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.8))
    
    # Eigenvalue convergence
    ax1.semilogx(grid_sizes, eigenvalues, 'bo-', markerfacecolor='white')
    ax1.set_xlabel('Number of Grid Points')
    ax1.set_ylabel('Eigenvalue (eV)')
    ax1.set_title('Eigenvalue Convergence')
    
    # Computation time
    ax2.loglog(grid_sizes, times, 'rs-', markerfacecolor='white')
    ax2.set_xlabel('Number of Grid Points')
    ax2.set_ylabel('Wall Time (s)')
    ax2.set_title('Computational Cost')
    
    plt.tight_layout()
    plt.savefig(output)
    plt.close()
    print(f"  Saved: {output}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='QuDSim Scaling Study Plotter')
    parser.add_argument('--input', default='scaling_results.csv', help='Scaling CSV file')
    parser.add_argument('--phase-csv', default=None, help='Per-phase timing CSV')
    args = parser.parse_args()
    
    print("QuDSim Scaling Study Plotter")
    print("============================")
    
    try:
        nranks, best, avg, std = load_scaling_data(args.input)
        print(f"\nLoaded {len(nranks)} data points from {args.input}")
        print(f"Core counts: {nranks.tolist()}")
        print(f"Best times:  {best.tolist()}\n")
        
        plot_speedup(nranks, best)
        plot_efficiency(nranks, best)
        plot_walltime(nranks, best, avg, std)
        
    except FileNotFoundError:
        print(f"\nWARNING: {args.input} not found. Run scaling_study.sh first.")
        print("Generating example plots with synthetic data...\n")
        
        # Synthetic data for testing
        nranks = np.array([1, 2, 4, 8, 16, 32, 64, 128])
        T1 = 240.0
        # Realistic Amdahl's law with 5% serial fraction
        serial_frac = 0.05
        best = T1 * (serial_frac + (1 - serial_frac) / nranks)
        avg = best * 1.05  # 5% overhead
        std = best * 0.02
        
        plot_speedup(nranks, best)
        plot_efficiency(nranks, best)
        plot_walltime(nranks, best, avg, std)
        
        # Example grid convergence
        grid_sizes = np.array([500, 1000, 5000, 10000, 50000, 100000])
        evals = np.array([0.230, 0.234, 0.2352, 0.2354, 0.23545, 0.23546])
        times = np.array([0.3, 0.8, 5.2, 18.5, 120.3, 580.1])
        plot_grid_convergence(grid_sizes, evals, times)
    
    if args.phase_csv:
        try:
            phases = load_phase_data(args.phase_csv)
            plot_phase_breakdown(phases)
        except FileNotFoundError:
            print(f"WARNING: {args.phase_csv} not found.")
    
    print("\nDone! Import these PDFs into your LaTeX manuscript.")
