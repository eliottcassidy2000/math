#!/bin/sh
# Regeneration recipe (macOS, nauty 2.9.3 in /opt/homebrew, python3 with numpy and scipy).  All timings on an 8-core Apple machine.
set -e
NAUTY=/opt/homebrew
# 1. order-11 prime cycle-index types, symmetry first (candidates: 1,316,960 lines, ~75 MB; 1 min to generate, 7 s to census)
for t in '3^2 1^5' '3^3 1^2' '5^2 1' '11' '2^5 1' '6 2^2 1' '10 1'; do python3 gen_sym11.py "$t" "cand_$(echo $t | tr ' ^' '_p').txt"; done
cat cand_3p2_1p5.txt cand_3p3_1p2.txt cand_5p2_1.txt cand_11.txt cand_2p5_1.txt cand_6_2p2_1.txt cand_10_1.txt > cand_all.txt
cc -O2 -I$NAUTY/include/nauty tourn_types_dedup.c -L$NAUTY/lib -lnauty -o tourn_types_dedup
./tourn_types_dedup 11 --types t11_sym.txt < cand_all.txt > s11_sym.txt
python3 anti_types_check.py 12 > anti_types_check.out            # the 2 mod 4 law for anti-automorphism cycles
# 2. the order-12 ledger: data/t3.txt..t10.txt are the order-10 packet's prime_types_n3..n10.txt; t11.txt = t11_sym.txt + the rigid non-SC line
mkdir -p data; (cat t11_sym.txt; echo "TYPE count=796605430 k=11 aut=1 AUT: 1^11,:1 SC=0 ANTI:") > data/t11.txt
DATA=data python3 species11_fast.py 12 > species11_fast_12.out    # A, S (non-prime side), c_s, c, P_lab, l   (seconds)
DATA=data python3 species12b.py 12 > species12b.out               # SC decomposition at 12, s_lab(12), P_lab(12)
# 3. the dipole census on isomorphism classes, orders 4-8 (mg8.txt is 1.8 M lines, regenerated in seconds)
cc -O2 -I$NAUTY/include/nauty movegraph.c -L$NAUTY/lib -lnauty -o movegraph
for n in 4 5 6 7 8; do $NAUTY/bin/gentourng -q $n | ./movegraph $n > mg$n.txt; done
for n in 5 6 7 8; do python3 analyze_mg.py $n > an$n.txt; done
for n in 6 7 8; do python3 singleton_failures.py $n > singleton_failures_n$n.out; done
python3 lp_dipole.py 5 > lp5.txt; python3 lp_dipole.py 6 > lp6.txt   # n=6: ~3 min (80,640 LPs)
cc -O2 lbfs7.c -o lbfs7; ./lbfs7 paley > lbfs7_paley.out; ./lbfs7 tt > lbfs7_tt.out   # 2^21 labelled states, 15 s each
# 4. all 1,223 regular 11-tournaments: an all-black word of length 5 to a singleton plus a strong core (~8 min)
$NAUTY/bin/gentourng -q -d5D5 11 | python3 blackcore.py 11 gentourng 60 first > bc11all.out 2> bc11all.err
