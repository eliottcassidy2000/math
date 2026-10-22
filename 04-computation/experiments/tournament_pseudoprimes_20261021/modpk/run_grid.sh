#!/bin/bash
# grid for THM-4748 tests
cd "$(dirname "$0")"
{
for k in 6 7 8 9 10; do ./pk_period 3 $k; done
for k in 4 5 6; do ./pk_period 5 $k; done
for k in 4 5; do ./pk_period 7 $k; done
./pk_period 11 4; ./pk_period 13 4
for p in 17 19 23; do ./pk_period $p 3; done
for p in 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97; do ./pk_period $p 2; done
} > pk_period_grid.out 2>&1
echo GRID_DONE >> pk_period_grid.out
