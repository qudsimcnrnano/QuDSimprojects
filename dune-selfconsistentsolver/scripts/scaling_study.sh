#!/bin/bash
###############################################################################
# QuDSim Scaling Study Script
#
# Runs the QuDSim executable with varying numbers of MPI ranks
# and collects timing data for strong/weak scaling analysis.
#
# Usage:
#   ./scaling_study.sh <executable> [mesh_file] [extra_args...]
#
# Example:
#   ./scaling_study.sh ./qudsim_poisson mesh.msh
#   ./scaling_study.sh ./qudsim_gaa device.geo -ksp_type gmres
#
# Output:
#   scaling_results.csv          - timing data
#   scaling_report.txt           - human-readable summary
#   scaling_plot.py              - Python script to generate plots
###############################################################################

set -e

EXE="${1:?Usage: $0 <executable> [mesh_file] [extra_args...]}"
shift
EXTRA_ARGS="$@"

# Configuration - modify these for your system
CORE_LIST="1 2 4 8 16 32"    # Number of MPI ranks to test
REPEATS=3                      # Runs per configuration (take minimum)
RESULTS_FILE="scaling_results.csv"
REPORT_FILE="scaling_report.txt"
HOSTFILE=""                    # Set to your MPI hostfile if needed, e.g. "hostfile.txt"

# MPI launcher detection
if command -v mpirun &> /dev/null; then
    MPI_LAUNCHER="mpirun"
elif command -v mpiexec &> /dev/null; then
    MPI_LAUNCHER="mpiexec"
elif command -v srun &> /dev/null; then
    MPI_LAUNCHER="srun"  # SLURM
else
    echo "ERROR: No MPI launcher found (mpirun/mpiexec/srun)"
    exit 1
fi

echo "=============================================="
echo " QuDSim Strong Scaling Study"
echo "=============================================="
echo " Executable:  $EXE"
echo " Extra args:  $EXTRA_ARGS"
echo " Core counts: $CORE_LIST"
echo " Repeats:     $REPEATS per config"
echo " MPI:         $MPI_LAUNCHER"
echo "=============================================="
echo ""

# Initialize CSV
echo "nranks,run,wall_time_s" > "$RESULTS_FILE"

# Initialize report
cat > "$REPORT_FILE" << EOF
QuDSim Scaling Study Report
============================
Date: $(date)
Executable: $EXE
Args: $EXTRA_ARGS
Repeats per config: $REPEATS

EOF

declare -A BEST_TIMES

for NP in $CORE_LIST; do
    echo "--- Testing with $NP MPI rank(s) ---"
    BEST=999999.0

    for RUN in $(seq 1 $REPEATS); do
        echo -n "  Run $RUN/$REPEATS ... "

        # Build MPI command
        MPI_CMD="$MPI_LAUNCHER"
        if [ "$MPI_LAUNCHER" = "srun" ]; then
            MPI_CMD="$MPI_CMD -n $NP"
        else
            MPI_CMD="$MPI_CMD -np $NP"
            if [ -n "$HOSTFILE" ]; then
                MPI_CMD="$MPI_CMD --hostfile $HOSTFILE"
            fi
            # Allow oversubscription for testing on local machine
            MPI_CMD="$MPI_CMD --oversubscribe"
        fi

        # Run and capture wall time
        START_TIME=$(date +%s%N)
        $MPI_CMD $EXE $EXTRA_ARGS > /dev/null 2>&1
        END_TIME=$(date +%s%N)

        # Calculate elapsed time in seconds (with nanosecond precision)
        ELAPSED=$(echo "scale=6; ($END_TIME - $START_TIME) / 1000000000" | bc)

        echo "${ELAPSED}s"
        echo "$NP,$RUN,$ELAPSED" >> "$RESULTS_FILE"

        # Track best (minimum) time
        BETTER=$(echo "$ELAPSED < $BEST" | bc -l)
        if [ "$BETTER" -eq 1 ]; then
            BEST=$ELAPSED
        fi
    done

    BEST_TIMES[$NP]=$BEST
    echo "  Best: ${BEST}s"
    echo ""
done

# Generate report
echo "Results Summary" >> "$REPORT_FILE"
echo "===============" >> "$REPORT_FILE"
printf "%-10s %-15s %-12s %-12s\n" "Ranks" "Best Wall(s)" "Speedup" "Efficiency" >> "$REPORT_FILE"
printf "%-10s %-15s %-12s %-12s\n" "-----" "-----------" "-------" "----------" >> "$REPORT_FILE"

# Get baseline (1 core)
BASELINE=${BEST_TIMES[1]:-${BEST_TIMES[2]}}
BASELINE_NP=1
if [ -z "${BEST_TIMES[1]}" ]; then
    BASELINE_NP=2
fi

for NP in $CORE_LIST; do
    T=${BEST_TIMES[$NP]}
    if [ -n "$T" ] && [ -n "$BASELINE" ]; then
        SPEEDUP=$(echo "scale=2; $BASELINE / $T" | bc -l)
        EFFICIENCY=$(echo "scale=1; 100 * $BASELINE / ($T * $NP / $BASELINE_NP)" | bc -l)
        printf "%-10s %-15s %-12s %-12s\n" "$NP" "${T}" "${SPEEDUP}x" "${EFFICIENCY}%" >> "$REPORT_FILE"
    fi
done

echo "" >> "$REPORT_FILE"
echo "Baseline: $BASELINE_NP rank(s) = ${BASELINE}s" >> "$REPORT_FILE"

# Print report to screen
echo ""
cat "$REPORT_FILE"

echo ""
echo "Data saved to: $RESULTS_FILE"
echo "Report saved to: $REPORT_FILE"
echo ""
echo "To generate plots, run: python3 plot_scaling.py"
