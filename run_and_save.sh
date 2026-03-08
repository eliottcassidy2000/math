#!/bin/bash
# Usage: ./run_and_save.sh SCRIPT.py [timeout_seconds]
# Runs a script from 04-computation/ and saves output to 05-knowledge/results/
#
# Examples:
#   ./run_and_save.sh reversal_proof_attempt.py
#   ./run_and_save.sh M_n7_structure.py 300

SCRIPT="$1"
TIMEOUT="${2:-120}"

if [ -z "$SCRIPT" ]; then
    echo "Usage: $0 SCRIPT.py [timeout_seconds]"
    exit 1
fi

BASE=$(basename "$SCRIPT" .py)
SCRIPT_PATH="04-computation/$SCRIPT"
OUT_PATH="05-knowledge/results/${BASE}.out"

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: $SCRIPT_PATH not found"
    exit 1
fi

echo "Running $SCRIPT_PATH (timeout: ${TIMEOUT}s)..."
echo "Output: $OUT_PATH"
echo "---"

timeout "$TIMEOUT" python3 "$SCRIPT_PATH" 2>&1 | tee "$OUT_PATH"
EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 124 ]; then
    echo "" >> "$OUT_PATH"
    echo "*** TIMED OUT after ${TIMEOUT}s ***" >> "$OUT_PATH"
    echo "*** Script timed out after ${TIMEOUT}s ***"
elif [ $EXIT_CODE -ne 0 ]; then
    echo "" >> "$OUT_PATH"
    echo "*** EXIT CODE: $EXIT_CODE ***" >> "$OUT_PATH"
fi

echo "---"
echo "Saved to $OUT_PATH"
