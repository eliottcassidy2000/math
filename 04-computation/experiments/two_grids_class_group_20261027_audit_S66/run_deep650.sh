#!/bin/sh
# After the n = 6 least-box runs: B = 650 for the other persistent misses
cd "$(dirname "$0")"
until grep -q DONEMIN a3_six_minbox.out 2>/dev/null; do sleep 5; done
for n in 78 60 124 140 75; do
  ./a3_search2.exe deep $n 650 5 | grep -E "^deep n = .*processed|hit:" >> a3_deep650.out
done
echo DONE650 >> a3_deep650.out
