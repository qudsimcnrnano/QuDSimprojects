#!/bin/bash
###############################################################################
# QuDSim Graphene QD Scaling Study
# Runs with varying MPI ranks and collects timing data
#
# Usage:
#   ./scaling_study.sh
#
# Prerequisites:
#   - dune-graphene1 must be compiled
#   - Adjust EXE and MPIEXEC paths below
###############################################################################

set -e

# ---- Configuration ----
EXE="./dune-graphene1/build-cmake/src/dune-graphene1"
MPIEXEC="/opt/local/petsc-3.15.2/arch-linux-cxx-debug/bin/mpiexec"
CORE_LIST="1 2 4 8"
REPEATS=2
RESULTS_FILE="scaling_results.csv"
REPORT_FILE="scaling_report.txt"

echo "=============================================="
echo " QuDSim Graphene QD Scaling Study"
echo "=============================================="
echo " Executable:  $EXE"
echo " Core counts: $CORE_LIST"
echo " Repeats:     $REPEATS per config"
echo "=============================================="
echo ""

# Initialize CSV
echo "nranks,run,wall_time_s" > "$RESULTS_FILE"

# Initialize report
cat > "$REPORT_FILE" << EOF
QuDSim Graphene QD Scaling Study Report
========================================
Date: $(date)
Executable: $EXE
Repeats per config: $REPEATS

EOF

declare -A BEST_TIMES

for NP in $CORE_LIST; do
    echo "--- Testing with $NP MPI rank(s) ---"
    BEST=999999.0

    for RUN in $(seq 1 $REPEATS); do
        echo -n "  Run $RUN/$REPEATS ... "

        START_TIME=$(date +%s%N)
        $MPIEXEC -np $NP $EXE > /dev/null 2>&1 || true
        END_TIME=$(date +%s%N)

        ELAPSED=$(echo "scale=6; ($END_TIME - $START_TIME) / 1000000000" | bc)
        echo "${ELAPSED}s"
        echo "$NP,$RUN,$ELAPSED" >> "$RESULTS_FILE"

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

echo ""
cat "$REPORT_FILE"
echo ""
echo "Data: $RESULTS_FILE"
echo "Report: $REPORT_FILE"
echo "Plot:  python3 plot_scaling.py"
