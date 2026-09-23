#!/bin/bash
# AMM 12592 procgen lane 2026-09-23: reproduce every number in
#   05-knowledge/results/amm12592_procgen_20260923_uniform_frontier.md
# Output: 05-knowledge/results/amm12592_procgen_20260923.out
# Memory: every step stays below ~400 MB; steps run sequentially.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO" || exit 1
OUT="05-knowledge/results/amm12592_procgen_20260923.out"
DRAFT="scratch/procgen_amm/long_draft"
UA="Mozilla/5.0 (research; math-repo)"
BASE="https://raw.githubusercontent.com/long-mathematics/deterministic-von-neumann-fair-extractor/main"

fetch_draft() {
  mkdir -p "$DRAFT/scripts/certificates"
  for f in README.md deterministic_von_neumann_fair_extractor.tex scripts/verify_glazer_critical.py \
           scripts/supplementary_checks.py scripts/certificates/finite_blocks.json \
           scripts/verification_output.txt scripts/supplementary_verification_output.txt \
           scripts/finite_certificate_check.txt; do
    [ -s "$DRAFT/$f" ] || curl -sS -A "$UA" -o "$DRAFT/$f" "$BASE/$f" || { echo "fetch failed: $f"; return 1; }
  done
}

{
  echo "# AMM 12592 procgen 2026-09-23 -- run $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# python: $(python3 --version 2>&1); numpy $(python3 -c 'import numpy;print(numpy.__version__)'); scipy $(python3 -c 'import scipy;print(scipy.__version__)'); mpmath $(python3 -c 'import mpmath;print(mpmath.__version__)'); ortools $(python3 -c 'import ortools;print(ortools.__version__)')"
  for s in 04-computation/experiments/amm12592_procgen_20260923_*.py 04-computation/experiments/amm12592_procgen_20260923_run.sh; do
    echo "# sha256 $(shasum -a 256 "$s" | cut -d' ' -f1)  $s"
  done
  echo
  echo "################ Part 0. Long's draft: fetch (if absent) and his own exact verifier"
  fetch_draft
  echo "# draft commit inspected: 6c5b43349cf96973d434d4cf43d6f92568b6c0ef (2026-09-22T23:16:43Z)"
  echo "# sha256 verifier: $(shasum -a 256 $DRAFT/scripts/verify_glazer_critical.py | cut -d' ' -f1)"
  (cd "$DRAFT" && python3 scripts/verify_glazer_critical.py && \
     python3 scripts/verify_glazer_critical.py --check-certificate scripts/certificates/finite_blocks.json && \
     python3 scripts/supplementary_checks.py && \
     python3 scripts/verify_glazer_critical.py --extra 128 256 512 1024 --export /tmp/amm12592_long_reexport.json && \
     echo "re-exported certificate identical to distributed one: $(cmp -s /tmp/amm12592_long_reexport.json scripts/certificates/finite_blocks.json && echo yes || echo NO)")
  echo
  echo "################ Part 1. Independent audit of Long's draft"
  python3 04-computation/experiments/amm12592_procgen_20260923_long_audit.py --nmax 512 --nmax-lac 128
  echo
  echo "################ Part 2. Polya-capacity lower bound (obstruction side)"
  python3 04-computation/experiments/amm12592_procgen_20260923_polya_capacity.py
  echo
  echo "################ Part 3. Super-blocks (construction side)"
  python3 04-computation/experiments/amm12592_procgen_20260923_superblocks.py --tlimit 120 --fold-nmax 1024
  echo
  echo "# sha256 subordination certificate: $(shasum -a 256 05-knowledge/results/amm12592_procgen_20260923_subordination_certificate.json | cut -d' ' -f1)"
  echo "# end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} 2>&1 | tee "$OUT"
