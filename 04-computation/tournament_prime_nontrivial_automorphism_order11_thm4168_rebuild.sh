#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C

REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
CXX_BIN=${CXX_BIN:-clang++}
LABELG_BIN=${LABELG_BIN:-/opt/homebrew/bin/labelg}
SOURCE=${SOURCE:-$REPO_ROOT/04-computation/tournament_prime_nontrivial_automorphism_order11_thm4168_census.cpp}
OUTPUT_D6=${OUTPUT_D6:-$REPO_ROOT/05-knowledge/results/tournament_prime_nontrivial_automorphism_order11_thm4168.d6}
RUN_DIR=${RUN_DIR:-$(mktemp -d /tmp/thm4168-prime-symmetry.XXXXXX)}
BINARY="$RUN_DIR/prime-symmetry-census"

"$CXX_BIN" -std=c++17 -O3 -Wall -Wextra -Werror "$SOURCE" -o "$BINARY"

for residue in $(seq 0 15); do
  (
    "$BINARY" --type 3x2 --res "$residue" --mod 16 --emit-digraph6 \
      | tee >(wc -l > "$RUN_DIR/3x2.raw.$residue.count") \
      | "$LABELG_BIN" -q \
      | sort -u > "$RUN_DIR/3x2.canon.$residue.d6"
  ) &
done

for type in 3x3 5x2 11; do
  (
    "$BINARY" --type "$type" --emit-digraph6 \
      | tee >(wc -l > "$RUN_DIR/$type.raw.count") \
      | "$LABELG_BIN" -q \
      | sort -u > "$RUN_DIR/$type.canon.d6"
  ) &
done
wait

sort -u "$RUN_DIR"/3x2.canon.*.d6 > "$RUN_DIR/3x2.union.d6"
sort -u \
  "$RUN_DIR/3x2.union.d6" \
  "$RUN_DIR/3x3.canon.d6" \
  "$RUN_DIR/5x2.canon.d6" \
  "$RUN_DIR/11.canon.d6" \
  > "$RUN_DIR/all.union.d6"

echo "run_dir=<temporary>"
echo "raw_3x2_by_residue"
for residue in $(seq 0 15); do
  printf '%s ' "$residue"
  tr -d ' ' < "$RUN_DIR/3x2.raw.$residue.count"
done
echo "raw_other_types"
for type in 3x3 5x2 11; do
  printf '%s ' "$type"
  tr -d ' ' < "$RUN_DIR/$type.raw.count"
done
echo "canonical_counts"
for type in 3x2.union 3x3.canon 5x2.canon 11.canon all.union; do
  count=$(wc -l < "$RUN_DIR/$type.d6" | tr -d ' ')
  printf '%s %s\n' "$type" "$count"
done
printf 'census_sha256 '
sha256sum "$SOURCE" | cut -d' ' -f1
printf 'canonical_d6_sha256 '
sha256sum "$RUN_DIR/all.union.d6" | cut -d' ' -f1

if [[ -e "$OUTPUT_D6" ]]; then
  cmp "$RUN_DIR/all.union.d6" "$OUTPUT_D6"
else
  cp "$RUN_DIR/all.union.d6" "$OUTPUT_D6"
fi
cmp "$RUN_DIR/all.union.d6" "$OUTPUT_D6"
echo "frozen_d6_byte_match=1"
