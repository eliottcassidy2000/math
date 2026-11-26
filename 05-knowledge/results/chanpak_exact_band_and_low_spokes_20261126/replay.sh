#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT
cxx="${CXX:-g++}"
"$cxx" -std=c++17 -O3 -Wall -Wextra band.cpp -o "$scratch/band"
"$cxx" -std=c++17 -O3 -Wall -Wextra verify.cpp -o "$scratch/verify"
awk '{print $1,$2,$3,$4}' used_graphs.tsv >"$scratch/profiles"
lo=34359738369
hi=34379738367
"$scratch/band" "$lo" "$hi" 2048 "$scratch/witnesses" "$scratch/profiles"
"$scratch/verify" "$lo" "$hi" "$scratch/witnesses" used_graphs.tsv
echo '9e06c1ae9aca7bb758136756071ec1a28d79e1329ed230366ae53b4bbc84469c  witnesses' | (cd "$scratch"; sha256sum -c -)
