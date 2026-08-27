#!/bin/sh
set -eu

script_dir=$(CDPATH= cd "$(dirname "$0")" && pwd)
repo=${1:-$(dirname "$script_dir")}
cxx=${CXX:-clang++}
opt_level=${OPT_LEVEL:--O3}
source_file="$repo/04-computation/tournament_plus_min_obstruction_thm4224.cpp"
generator="$repo/04-computation/tournament_two_reversal_words_thm4239.py"
expected_source=4c0b1ac8f80cf41d9091420c852e1bd7595fed94105ff577535d4939fc88afb6
expected_generator=450504c9cb05720f389cf02758c32bbd80fc01eb10c5c483bbe2960063405238
run_dir=$(mktemp -d /tmp/two-reversal-order18.XXXXXX)
trap 'rm -rf -- "$run_dir"' EXIT HUP INT TERM

source_hash=$(sha256sum "$source_file" | awk '{print $1}')
generator_hash=$(sha256sum "$generator" | awk '{print $1}')
if [ "$source_hash" != "$expected_source" ]; then
  echo "ERROR: inherited C++ source hash mismatch" >&2
  exit 1
fi
if [ "$generator_hash" != "$expected_generator" ]; then
  echo "ERROR: word-generator hash mismatch" >&2
  exit 1
fi

generated="$run_dir/two_reversal_order18.cpp"
binary="$run_dir/two_reversal_order18"
LC_ALL=C sed \
  -e 's/n > 12/n > 20/' \
  -e 's/census scope is 3\.\.12/census scope is 3..20/' \
  -e 's/strong_classes=/strong_presentations=/' \
  "$source_file" > "$generated"
"$cxx" -x c++ -std=c++20 "$opt_level" -DNDEBUG "$generated" -o "$binary"

printf 'dependency_source_sha256=%s\n' "$source_hash"
printf 'generator_sha256=%s\n' "$generator_hash"
printf 'compile_transform=order_cap_20,presentation_label; optimization=%s\n' "$opt_level"
for order in 15 16 17 18; do
  python3 -B "$generator" "$order" | "$binary" --census
done
