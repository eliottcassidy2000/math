#!/bin/sh
# Regenerates every census and receipt from nauty's gentourng (nauty 2.9.3, /opt/homebrew).  ~5 minutes on 8 cores.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); W=${W:-/tmp/a568_work}; mkdir -p "$W"; cd "$W"
cc -O2 -I/opt/homebrew/include/nauty "$HERE/tourn_enum.c" -L/opt/homebrew/lib -lnauty -o tourn_enum
cc -O2 "$HERE/c3class.c" -o c3class; cc -O2 "$HERE/c3score.c" -o c3score
for n in 3 4 5 6 7 8 9; do gentourng $n 2>/dev/null | ./tourn_enum $n --fas --hpaths --c3canon --table t$n.tsv --primes p$n.txt > s$n.txt; done
for n in 3 4 5 6 7; do gentourng $n 2>/dev/null | ./c3class $n > cc$n.txt; gentourng $n 2>/dev/null | ./c3score $n | tail -1; done
for r in 0 1 2 3 4 5 6 7; do (gentourng 10 $r/8 2>/dev/null | ./tourn_enum 10 --c3canon --noreg --compact --table c10_$r.tsv > s10_$r.txt) & done; wait
cat s3.txt s4.txt s5.txt s6.txt s7.txt s8.txt s9.txt
DATA="$W" python3 "$HERE/species.py" 9 p9.txt; DATA="$W" python3 "$HERE/species10.py"; DATA="$W" python3 "$HERE/analyze.py"; python3 "$HERE/conn_check.py"
