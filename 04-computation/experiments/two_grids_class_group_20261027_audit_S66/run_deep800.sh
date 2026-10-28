#!/bin/sh
cd "$(dirname "$0")"
for n in 6 54; do
  ./a3_search2.exe deep $n 800 5 | grep -E "^deep n = .*processed|hit:" >> a3_deep800.out
done
echo DONE800 >> a3_deep800.out
