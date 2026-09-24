#!/usr/bin/env bash
# Regenerates the implication-atlas output (session collatz-procgen-20260922, lane atlas, 2026-09-24).
#   bash 04-computation/experiments/procgen_atlas_20260924_run.sh > 05-knowledge/results/procgen_atlas_20260924.out
# One process at a time; peak memory well under 100 MB; about 1 minute on the mac-mini.
# Also (re)writes 05-knowledge/results/procgen_atlas_20260924_implications.dot.
set -euo pipefail
cd "$(dirname "$0")/../.."
echo "procgen_atlas_20260924: implication atlas output (deterministic; timing on stderr)"
for s in graph cycles mahler_bridge collatz_permutation martingale primegame; do
  echo
  echo "################################################################################"
  echo "# procgen_atlas_20260924_${s}.py"
  echo "################################################################################"
  start=$(date +%s)
  python3 "04-computation/experiments/procgen_atlas_20260924_${s}.py"
  echo "[${s}: $(( $(date +%s) - start )) s]" >&2
done
