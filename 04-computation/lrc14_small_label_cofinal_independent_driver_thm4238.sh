#!/bin/sh
set -eu

repo=${1:-/Users/e/Documents/GitHub/math}
cxx=${CXX:-clang++}
libomp_prefix=${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}
jobs=${JOBS:-3}
opt_level=${OPT_LEVEL:--O3}
ledger=${LEDGER_PATH:-/tmp/lrc14-small-label-independent-ledger.out}
template="$repo/04-computation/lrc14_fixed_one_outsider_cofinal_tail_independent_audit_thm4231.cpp"
run_dir=$(mktemp -d /tmp/lrc14-small-label-independent.XXXXXX)
trap 'rm -rf -- "$run_dir"' EXIT HUP INT TERM

printf '%s\n' \
  '2 563' '3 587' '4 589' '5 528' '6 614' '7 557' '9 547' \
  '11 555' '12 586' '13 557' '14 538' '17 581' '18 560' '19 572' \
  '21 545' '22 579' '23 565' '24 589' '25 598' '26 561' '27 552' \
  '28 550' '29 546' '31 571' '32 569' '33 562' '34 551' '35 541' \
  '36 546' '37 552' '38 569' '39 535' '41 540' '43 549' '44 568' \
  '45 546' '46 572' '47 577' '48 580' '49 538' |
xargs -P "$jobs" -n 2 sh -c '
  repo="$1"
  cxx="$2"
  libomp_prefix="$3"
  template="$4"
  opt_level="$5"
  run_dir="$6"
  q="$7"
  k="$8"
  km1=$((k-1))
  d=$((14*q))
  src="$run_dir/q${q}.cpp"
  bin="$run_dir/q${q}"
  out="$run_dir/q${q}.out"
  LC_ALL=C sed \
    -e "s/q=1/q=${q}/g" \
    -e "s/q1/q${q}/g" \
    -e "s/Q1/Q${q}/g" \
    -e "s/542/${km1}/g" \
    -e "s/543/${k}/g" \
    -e "s/append_speed(1, 30)/append_speed(${q}, 30)/" \
    -e "s/push_back(1)/push_back(${q})/g" \
    -e "s/i64{14});/i64{${d}});/" \
    "$template" > "$src"
  "$cxx" "$opt_level" -std=c++20 -Xpreprocessor -fopenmp \
    -I"$libomp_prefix/include" -L"$libomp_prefix/lib" \
    -Wl,-rpath,"$libomp_prefix/lib" -lomp "$src" -o "$bin"
  OMP_NUM_THREADS=4 "$bin" > "$out"
' sh "$repo" "$cxx" "$libomp_prefix" "$template" "$opt_level" "$run_dir"

awk '
  /^LRC14_FIXED_Q/ {q=$1; sub("LRC14_FIXED_Q","",q); sub("_RAY.*","",q)}
  /^GEOMETRY/ {D=$3; walls=$5; cells=$7; keys=$9}
  /^CENSUS/ {b=$3; pos=$5; kmin=$7; minw=$9; kmax=$11; cnt=$13; maxw=$15}
  /^SECOND_THRESHOLD/ {k2=$2; c2=$4}
  /^EXTREMAL_FIXED_GEOMETRY/ {mass=$3; comp=$5; surp=$7}
  /^COMMUTATIVE_DIGEST/ {dx=$3; ds=$5}
  /^VERDICT/ {print q,D,walls,cells,keys,b,pos,kmin,minw,kmax,cnt,maxw,k2,c2,mass,comp,surp,dx,ds,$3}
' "$run_dir"/q*.out | sort -n > "$ledger"
