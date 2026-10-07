#!/bin/bash
cd "$(dirname "$0")"
until [ -f queue_done.txt ]; do sleep 5; done
mkdir -p A
for R in 1 2 3 4 6 8; do ./qsurv trans $R 200000 1000000 $((100+R)) > A/trans_R$R.txt; done
for R in 12 16 24 32 48 64; do ./qsurv trans $R 100000 1000000 $((100+R)) > A/trans_R$R.txt; done
./qsurv mers 1 200000 1000000 201 > A/mers_D1.txt
./qsurv mers 61 100000 1000000 261 > A/mers_D61.txt
./qsurv mers 241 40000 1000000 441 > A/mers_D241.txt
./qsurv mers 961 10000 100000 961 > A/mers_D961.txt
echo DONE > A/done.txt
