#!/usr/bin/env python3
"""
Plot scaling results from QuDSim scaling study.
Reads: scaling_results.csv OR timings.csv files
Produces: scaling_plot.png
"""

import csv
import sys
import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available, generating text report only")


def read_scaling_csv(filename="scaling_results.csv"):
    """Read scaling_results.csv: nranks, run, wall_time_s"""
    data = {}
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            np_val = int(row['nranks'])
            t = float(row['wall_time_s'])
            if np_val not in data:
                data[np_val] = []
            data[np_val].append(t)
    # Take minimum time per rank count
    result = {}
    for np_val, times in sorted(data.items()):
        result[np_val] = min(times)
    return result


def read_timings_csv(directory="."):
    """Read multiple timings.csv files if scaling_results.csv doesn't exist."""
    # Look for timings files with rank info
    data = {}
    for f in os.listdir(directory):
        if f.startswith("timings") and f.endswith(".csv"):
            with open(os.path.join(directory, f)) as fh:
                reader = csv.DictReader(fh)
                total = 0
                for row in reader:
                    total += float(row.get('wall_avg', 0))
            # Try to extract np from filename
            # Default: just use as single entry
            data[1] = total
    return data


def plot_scaling(data, output="scaling_plot.png"):
    """Generate scaling plot."""
    ranks = sorted(data.keys())
    times = [data[r] for r in ranks]
    baseline = times[0]
    baseline_np = ranks[0]

    speedups = [baseline / t for t in times]
    ideal = [r / baseline_np for r in ranks]
    efficiencies = [100.0 * baseline / (t * r / baseline_np) for r, t in zip(ranks, times)]

    # Text report
    print("\nScaling Summary:")
    print(f"{'Ranks':>6} {'Time(s)':>10} {'Speedup':>10} {'Efficiency':>12}")
    print("-" * 40)
    for r, t, s, e in zip(ranks, times, speedups, efficiencies):
        print(f"{r:>6} {t:>10.3f} {s:>10.2f}x {e:>11.1f}%")

    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Wall time
    axes[0].plot(ranks, times, 'bo-', linewidth=2, markersize=8)
    axes[0].set_xlabel('MPI Ranks')
    axes[0].set_ylabel('Wall Time (s)')
    axes[0].set_title('Strong Scaling: Wall Time')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_xticks(ranks)

    # Speedup
    axes[1].plot(ranks, speedups, 'ro-', linewidth=2, markersize=8, label='Actual')
    axes[1].plot(ranks, ideal, 'k--', linewidth=1.5, label='Ideal')
    axes[1].set_xlabel('MPI Ranks')
    axes[1].set_ylabel('Speedup')
    axes[1].set_title('Strong Scaling: Speedup')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xticks(ranks)

    # Efficiency
    axes[2].bar(ranks, efficiencies, color='green', alpha=0.7, width=[r*0.3 for r in ranks])
    axes[2].axhline(y=100, color='k', linestyle='--', linewidth=1)
    axes[2].set_xlabel('MPI Ranks')
    axes[2].set_ylabel('Efficiency (%)')
    axes[2].set_title('Parallel Efficiency')
    axes[2].grid(True, alpha=0.3, axis='y')
    axes[2].set_xticks(ranks)
    axes[2].set_ylim(0, 120)

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {output}")


def main():
    csv_file = "scaling_results.csv"
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]

    if os.path.exists(csv_file):
        data = read_scaling_csv(csv_file)
    else:
        print(f"File {csv_file} not found.")
        print("Run the scaling study first: ./scaling_study.sh")
        sys.exit(1)

    if not data:
        print("No data found!")
        sys.exit(1)

    plot_scaling(data)


if __name__ == '__main__':
    main()
