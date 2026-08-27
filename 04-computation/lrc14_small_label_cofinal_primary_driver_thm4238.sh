#!/bin/sh
set -eu

repo=${1:-/Users/e/Documents/GitHub/math}
cxx=${CXX:-clang++}
libomp_prefix=${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}
jobs=${JOBS:-3}
opt_level=${OPT_LEVEL:--O3}
ledger=${LEDGER_PATH:-/tmp/lrc14-small-label-primary-ledger.out}
run_dir=$(mktemp -d /tmp/lrc14-small-label-primary.XXXXXX)
trap 'rm -rf -- "$run_dir"' EXIT HUP INT TERM
binary="$run_dir/runner"

"$cxx" "$opt_level" -std=c++20 -Wno-unknown-pragmas \
  -I"$libomp_prefix/include" \
  "$repo/04-computation/lrc14_fixed_one_outsider_cofinal_tail_primary_thm4231.cpp" \
  -o "$binary"

printf '%s\n' 2 3 4 5 6 7 9 11 12 13 14 17 18 19 21 22 23 24 25 26 \
  27 28 29 31 32 33 34 35 36 37 38 39 41 43 44 45 46 47 48 49 |
xargs -P "$jobs" -n 1 sh -c '
  binary="$1"
  run_dir="$2"
  q="$3"
  "$binary" "$q" > "$run_dir/q${q}.out"
' sh "$binary" "$run_dir"

awk '
  /^Q / {q=$2; D=$4; cells=$6; mc=$8; cc=$10; bodies=$12}
  /^LIMIT_SURPLUS/ {np=$3; mint=$5; minw=$7}
  /^TAIL_MAX/ {K=$3; cnt=$5; w=$7; mass=$9; comp=$11; surp=$13}
  /^FINGERPRINT/ {fp=$2}
  /^VERDICT/ {print q,D,cells,mc,cc,bodies,np,mint,minw,K,cnt,w,mass,comp,surp,fp,$2}
' "$run_dir"/q*.out | sort -n > "$ledger"
