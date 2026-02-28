#!/bin/bash
# =============================================================================
# scaling.sh - MPI Scaling Test for GAA Nanowire FET Simulator
# Runs the simulation with different numbers of MPI processes
# and collects timing data for strong scaling analysis.
#
# Usage:
#   ./scaling.sh                    # Default: test 1,2,4,8 cores
#   ./scaling.sh 1 2 4 8 16 32     # Custom core counts
#   ./scaling.sh --weak 1000 2000 4000  # Weak scaling with different mesh sizes
# =============================================================================

EXEC="./nanowire_parallel"
SLEPC_OPTS="-eps_type krylovschur -eps_nev 3 -eps_smallest_real -eps_tol 1e-10"
SCALING_FILE="scaling_data.csv"
LOG_DIR="scaling_logs"
REPEATS=3  # Number of repetitions for averaging

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =================== Parse Arguments ========================================
WEAK_SCALING=0
CORE_COUNTS=()

if [ "$1" == "--weak" ]; then
    WEAK_SCALING=1
    shift
    MESH_SIZES=("$@")
    if [ ${#MESH_SIZES[@]} -eq 0 ]; then
        echo "Usage: $0 --weak mesh_size1 mesh_size2 ..."
        exit 1
    fi
elif [ $# -gt 0 ]; then
    CORE_COUNTS=("$@")
else
    # Default: test with 1, 2, 4, 8 cores
    CORE_COUNTS=(1 2 4 8)
fi

# =================== Pre-flight Checks ======================================
if [ ! -f "$EXEC" ]; then
    echo -e "${RED}ERROR: Executable '$EXEC' not found.${NC}"
    echo "Please build first: make"
    exit 1
fi

# Detect available cores
MAX_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
echo -e "${BLUE}======================================================"
echo "  MPI Scaling Test - GAA Nanowire FET Simulator"
echo "  Available cores: $MAX_CORES"
echo "  Repetitions per config: $REPEATS"
echo -e "======================================================${NC}"

# Create log directory
mkdir -p "$LOG_DIR"

# Remove old scaling data
rm -f "$SCALING_FILE"

# =================== Strong Scaling =========================================
if [ $WEAK_SCALING -eq 0 ]; then

    echo ""
    echo -e "${YELLOW}Strong Scaling Test${NC}"
    echo "Core counts: ${CORE_COUNTS[*]}"
    echo ""

    for NP in "${CORE_COUNTS[@]}"; do

        if [ "$NP" -gt "$MAX_CORES" ]; then
            echo -e "${YELLOW}  WARNING: $NP cores requested but only $MAX_CORES available."
            echo -e "  Using --oversubscribe flag.${NC}"
            OVERSUBSCRIBE="--oversubscribe"
        else
            OVERSUBSCRIBE=""
        fi

        echo -e "${GREEN}--- Testing with $NP MPI process(es) ---${NC}"

        BEST_TIME=999999
        for REP in $(seq 1 $REPEATS); do
            echo -n "  Run $REP/$REPEATS ... "

            LOGFILE="$LOG_DIR/run_np${NP}_rep${REP}.log"

            # Run simulation
            START_TIME=$(date +%s%N)
            mpirun $OVERSUBSCRIBE -np $NP $EXEC $SLEPC_OPTS > "$LOGFILE" 2>&1
            EXIT_CODE=$?
            END_TIME=$(date +%s%N)

            ELAPSED=$(echo "scale=4; ($END_TIME - $START_TIME) / 1000000000" | bc)

            if [ $EXIT_CODE -eq 0 ]; then
                echo -e "${GREEN}OK${NC} (${ELAPSED}s)"

                # Track best time
                BETTER=$(echo "$ELAPSED < $BEST_TIME" | bc)
                if [ "$BETTER" -eq 1 ]; then
                    BEST_TIME=$ELAPSED
                fi
            else
                echo -e "${RED}FAILED${NC} (exit code $EXIT_CODE)"
                echo "  See: $LOGFILE"
            fi
        done

        echo "  Best time with $NP core(s): ${BEST_TIME}s"
        echo ""
    done
fi

# =================== Results Summary =========================================
echo ""
echo -e "${BLUE}======================================================"
echo "  Scaling Test Complete"
echo "======================================================${NC}"
echo ""

if [ -f "$SCALING_FILE" ]; then
    echo "Scaling data saved to: $SCALING_FILE"
    echo ""
    echo "Contents:"
    cat "$SCALING_FILE"
    echo ""
    echo "Generate plots with:"
    echo "  python3 plot_scaling.py"
else
    echo -e "${YELLOW}No scaling data file generated.${NC}"
    echo "Check logs in $LOG_DIR/ for errors."
fi

echo ""
echo "Log files in: $LOG_DIR/"
