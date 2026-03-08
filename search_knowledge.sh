#!/bin/bash
# Search the knowledge web for a term
# Usage: ./search_knowledge.sh TERM [area]
#
# Areas: variables, hypotheses, results, all (default)
#
# Examples:
#   ./search_knowledge.sh "transfer matrix"
#   ./search_knowledge.sh t_3 variables
#   ./search_knowledge.sh REFUTED hypotheses

TERM="$1"
AREA="${2:-all}"

if [ -z "$TERM" ]; then
    echo "Usage: $0 TERM [variables|hypotheses|results|all]"
    exit 1
fi

echo "=== Searching for '$TERM' ==="
echo ""

case "$AREA" in
    variables)
        echo "--- Variables ---"
        grep -rn "$TERM" 05-knowledge/variables/ 2>/dev/null
        ;;
    hypotheses)
        echo "--- Hypotheses ---"
        grep -rn "$TERM" 05-knowledge/hypotheses/ 2>/dev/null
        ;;
    results)
        echo "--- Results ---"
        grep -rn "$TERM" 05-knowledge/results/ 2>/dev/null | head -50
        ;;
    all)
        echo "--- Variables ---"
        grep -rn "$TERM" 05-knowledge/variables/ 2>/dev/null
        echo ""
        echo "--- Hypotheses ---"
        grep -rn "$TERM" 05-knowledge/hypotheses/ 2>/dev/null
        echo ""
        echo "--- Results (first 30 matches) ---"
        grep -rn "$TERM" 05-knowledge/results/ 2>/dev/null | head -30
        echo ""
        echo "--- Navigation ---"
        grep -n "$TERM" 00-navigation/*.md 2>/dev/null
        echo ""
        echo "--- Theorems ---"
        grep -rn "$TERM" 01-canon/theorems/ 2>/dev/null
        echo ""
        echo "--- Computation scripts (filenames) ---"
        ls 04-computation/ | grep -i "$(echo $TERM | tr ' ' '_')" 2>/dev/null
        ;;
esac
