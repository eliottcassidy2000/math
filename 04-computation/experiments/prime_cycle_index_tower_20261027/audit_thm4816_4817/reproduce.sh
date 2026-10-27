#!/bin/sh
# Independent audit of THM-4816 / THM-4817 (prime cycle-index tower; order-12..14 ledgers).
# Run from this directory.  Timings on the shared 8-core Mac mini, one core each:
#   orders 3-11 < 1 s, order 12 33 s, order 13 ~17 min (~270 MB), order-14 p = 5, 7 partial census ~1 s.
set -e
NAUTY=/opt/homebrew/opt/nauty
cc -O3 -Wall -I$NAUTY/include audit_census.c   $NAUTY/lib/libnauty.a -o audit_census
cc -O3 -Wall -I$NAUTY/include audit_census_p.c $NAUTY/lib/libnauty.a -o audit_census_p
cc -O3 -Wno-unused-function -I$NAUTY/include audit_reclass.c  $NAUTY/lib/libnauty.a -o audit_reclass
cc -O3 -Wno-unused-function -I$NAUTY/include prime_fraction.c $NAUTY/lib/libnauty.a -o prime_fraction
# 1. non-rigid prime census (auditor's own code), orders 3-13
for n in 3 4 5 6 7 8 9 10 11 12; do ./audit_census $n my_nr$n.txt my_cls$n.txt my_stats$n.txt; done
./audit_census 13 my_nr13.txt - my_stats13.txt
# 2. class-by-class comparison with the producer's class lists (re-canonicalised by the auditor's code)
for n in 3 5 6 7 8 9 10 11 12; do ./audit_reclass $n < ../out/nrclasses$n.txt > re_cls$n.txt; done
python3 - <<'EOF' > compare_classes.out
for n in [3,5,6,7,8,9,10,11,12]:
    mine = {l.split()[0]: tuple(l.split()[1:3]) for l in open(f'my_cls{n}.txt')}
    re = {l.split()[0]: tuple(l.split()[1:3]) for l in open(f're_cls{n}.txt')}
    print(f"n={n}: mine={len(mine)} producer={len(re)} same_class_set={set(mine)==set(re)} aut/sc agree={all(mine[k]==re[k] for k in mine if k in re)}")
EOF
# 3. Burnside A(n), SC(n) by constraint propagation; exact ledger and all comparisons
python3 audit_burnside.py 16 > my_burnside.out
python3 audit_main.py mine > audit_main.out
python3 narrative_checks.py > narrative_checks.out
# 4. circulants, primality cross-check, MISTAKE-683 checks, order-14 partial census, prime fractions
for p in 11 13; do python3 circulants.py $p burnside; python3 circulants.py $p arcs | ./audit_reclass $p 2>/dev/null | sort | uniq -c; done > circulants.out
python3 primality_crosscheck.py > primality_crosscheck.out
python3 check_s4s4.py > check_s4s4.out; python3 check_mistake683_examples.py > check_mistake683_examples.out
./audit_census_p 14 my_nr14_p7.txt - my_stats14_p7.txt 7; ./audit_census_p 14 my_nr14_p5.txt - my_stats14_p5.txt 5
for n in 10 11 12 13; do ./prime_fraction $n; done > prime_fraction.out
