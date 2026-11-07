#!/bin/sh
# level-floor minimisers: order, level, value
for spec in "3 1 3" "4 2 5" "5 3 9" "5 4 11" "5 5 15" "6 6 23" "6 7 33" "6 8 41" "7 9 51" "7 10 69" "8 11 97" "8 12 107" "8 13 141" "7 14 171" "9 15 213" "9 16 227" "9 17 321"; do
  set -- $spec
  gentourng -q $1 | ./hspec3 $1 witness $3 $2
done
