#!/bin/bash
# AMM 12592 procgen lane 2026-09-23, HYP-9128 / HYP-9129 follow-up. This script reproduces every number in
#   05-knowledge/results/amm12592_procgen_20260923_hyp9128_proof.md
# Output: 05-knowledge/results/amm12592_procgen_20260923_hyp9128.out
# The steps run one at a time. Each step runs under a resident-memory watchdog (limit 880 MB, polled every 0.5 s),
# and the observed peak is printed. The total is about 15 minutes.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO" || exit 1
OUT="05-knowledge/results/amm12592_procgen_20260923_hyp9128.out"
E="04-computation/experiments"
CERT="05-knowledge/results/amm12592_procgen_20260923_hyp9128_gamma377_certificate.json"
LONGCERT="scratch/procgen_amm/long_draft/scripts/certificates/finite_blocks.json"
UA="Mozilla/5.0 (research; math-repo)"
BASE="https://raw.githubusercontent.com/long-mathematics/deterministic-von-neumann-fair-extractor/main"
LIMIT_MB=880

step() {   # step "title" cmd args...
  local title="$1"; shift
  local tmp; tmp="$(mktemp)"
  echo "################ $title"
  echo "# \$ $*"
  local t0; t0=$(date +%s)
  "$@" > "$tmp" 2>&1 &
  local pid=$! peak=0 rss
  while kill -0 "$pid" 2>/dev/null; do
    rss=$(ps -o rss= -p "$pid" 2>/dev/null | awk '{print int($1/1024)}')
    [ "${rss:-0}" -gt "$peak" ] && peak=${rss:-0}
    if [ "${rss:-0}" -gt "$LIMIT_MB" ]; then echo "WATCHDOG: killing step at ${rss} MB" >> "$tmp"; kill -9 "$pid"; fi
    sleep 0.5
  done
  wait "$pid"; local rc=$?
  cat "$tmp"; rm -f "$tmp"
  echo "# step exit=$rc, wall $(( $(date +%s) - t0 )) s, peak RSS (0.5 s polling) ${peak} MB"
  echo
}

{
  echo "# AMM 12592 procgen 2026-09-23, HYP-9128/9129 follow-up -- run $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# python: $(python3 --version 2>&1); numpy $(python3 -c 'import numpy;print(numpy.__version__)'); scipy $(python3 -c 'import scipy;print(scipy.__version__)'); mpmath $(python3 -c 'import mpmath;print(mpmath.__version__)')"
  for s in "$E"/amm12592_procgen_20260923_hyp9128_*.py "$E"/amm12592_procgen_20260923_hyp9128_run.sh "$CERT"; do
    echo "# sha256 $(shasum -a 256 "$s" | cut -d' ' -f1)  $s"
  done
  if [ ! -s "$LONGCERT" ]; then
    mkdir -p "$(dirname "$LONGCERT")"
    curl -sS -A "$UA" -o "$LONGCERT" "$BASE/scripts/certificates/finite_blocks.json" || echo "# fetch of Long's certificate failed"
  fi
  echo "# sha256 $(shasum -a 256 "$LONGCERT" | cut -d' ' -f1)  $LONGCERT (Long's distributed certificate, used by the majorant audit)"
  echo

  echo "======================= PROOF OF HYP-9128 (C* <= 159/100 < C_*) ======================="
  step "Part 1. Analytic certificates, N >= N_A = 4096 (interval arithmetic, uniform in v = i/N)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_contours.py" --NA 4096
  step "Part 2. Finite exact certificates 16 <= N <= 2048 (margins) and independent Lemma-R rounding (N <= 256)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_finite.py" --nmin 16 --nmax 2048 --round-max 256
  step "Part 3. Exact checks of the fold, level-0 and pairing identities used in the proof" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_lemmas.py"

  echo "======================= THM-4467 extension (gamma = 377/1000) ======================="
  step "Part 4. Stored subordination certificate at gamma = 377/1000 (verify-only)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_gamma377.py" --verify-only "$CERT"

  echo "======================= Structured majorant (THM-4467 triangle step) ======================="
  step "Part 5. Box sharpness, parity skeleton, defect tables, exact kernel saturation" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_majorant.py"

  echo "======================= HYP-9129 numerics (NUMERICAL) ======================="
  step "Part 6a. Fold-sufficient thresholds of the golden-zero families (measure level)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" family --B 4 8 --steps 10
  step "Part 6b. H2: realizable 2-power cyclotomic states (alternating LP)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" lp --dict cyclo --B 4 --lo 1.54 --hi 1.60 --steps 9
  step "Part 6c. H2: realizable integer-factor states, degree <= 4 (alternating LP)" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" lp --dict factor --B 4 --lo 1.55 --hi 1.58 --steps 8 --ntheta 121
  step "Part 6d. H2: real zero-measure relaxation (not realizable), B = 4" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" lp --dict atoms --B 4 --lo 1.45 --hi 1.60 --steps 10 --ntheta 121 --full-rows
  step "Part 6e. H2: real zero-measure relaxation (not realizable), B = 8" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" lp --dict atoms --B 8 --lo 1.48 --hi 1.62 --steps 10 --ntheta 121 --full-rows
  step "Part 6f. H4: necessary-condition LP (single-block evaluation + mass + Mahler), B = 2, 4, 8" \
       python3 "$E/amm12592_procgen_20260923_hyp9128_h2.py" h4 --B 2 4 8
  echo "# end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} 2>&1 | tee "$OUT"
