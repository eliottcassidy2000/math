#!/bin/sh
# After the B = 800 batch: confirm the least box for 6 (expect 0 hits at B = 540, the representation at B = 541)
cd "$(dirname "$0")"
until grep -q DONE800 a3_deep800.out 2>/dev/null; do sleep 5; done
for B in 540 541; do
  ./a3_search2.exe deep 6 $B 5 | grep -E "^deep n = .*processed|hit:" >> a3_six_minbox.out
done
echo DONEMIN >> a3_six_minbox.out
