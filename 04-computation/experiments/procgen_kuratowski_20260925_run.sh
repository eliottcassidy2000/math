#!/usr/bin/env bash
# procgen_kuratowski_20260925: Kohl's group as a Tait-coloured graph (Kempe chains, quotients, boundary,
# switches), the sign fold and Althofer's 3n+-1 union graph (planarity 52, K5 68, projective 76, linkless 92,
# Petersen 104), the triple dictionary, and the excluded-minor calculus for Collatz-type maps.
# One process at a time; peak memory < 500 MB; about 6-10 minutes (CP-SAT infeasibility proofs dominate).
# Usage (repository root):
#   bash 04-computation/experiments/procgen_kuratowski_20260925_run.sh > 05-knowledge/results/procgen_kuratowski_20260925.out
set -euo pipefail
cd "$(dirname "$0")/../.."
for part in kohl_tait union_minors triples excluded_minors; do
  echo "#### procgen_kuratowski_20260925_${part}.py"
  python3 "04-computation/experiments/procgen_kuratowski_20260925_${part}.py"
  echo
done
