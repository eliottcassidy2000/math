#!/bin/sh
cd "$(dirname "$0")"
for d in $(seq 25 35); do
  for e2 in $(seq 1 $((d-9))); do echo "$d $e2"; done
done | xargs -P 8 -n 2 sh -c './lrc14_k11_band_mu certify "$0" "$1" "$1"'
echo "ALLDONE"
