#!/bin/sh
# Regenerates the packet data (nauty 2.9.3 under /opt/homebrew; python3 with numpy for the nuclear norms).
set -e
cd "$(dirname "$0")"
cc -O2 -I/opt/homebrew/include/nauty tourn_types.c -L/opt/homebrew/lib -lnauty -o tourn_types
for k in 3 4 5 6 7 8 9; do /opt/homebrew/bin/gentourng -q $k | ./tourn_types $k --types t$k.txt; done
for i in 0 1 2 3 4 5 6 7; do /opt/homebrew/bin/gentourng -q 10 $i/8 | ./tourn_types 10 --types prime_types_n10_part$i.txt > census_n10_part$i.txt & done; wait
python3 merge_types.py 'prime_types_n10_part*.txt' > t10.txt
python3 species11.py 11 > species11.out
python3 species11b.py > species11b.out
python3 nuclear.py > nuclear.out
for n in 4 5 6 7 8; do python3 nucmax.py $n; done > nucmax_n4_n8.out
python3 nucmax.py 9 > nucmax_n9.out; python3 nucmax9sv.py > nucmax_n9_maximizers.out
python3 parse_bound.py > parse_bound.out; python3 misc_checks.py > misc_checks.out
/opt/homebrew/bin/gentourng -q 7 | python3 blackcore.py 7 gentourng 0 all > bc7.out 2>/dev/null
/opt/homebrew/bin/gentourng -q 9 | python3 blackcore.py 9 gentourng 60 first > bc9.out 2>/dev/null
python3 blackcore.py 11 paley 40 > bc11p.out 2>/dev/null; python3 blackcore.py 11 randreg 40 12 > bc11r.out 2>/dev/null
