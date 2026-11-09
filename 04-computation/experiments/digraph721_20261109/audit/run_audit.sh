#!/bin/bash
# AUDIT D1 (THM-5055, HYP-9625; session mac-mini-2026-11-09-digraph721): regenerate every audit output.
# Run from this directory.  No third-party code (only the repo's own code and the auditor's code).
# Each step was run individually during the audit; total wall time about 25 minutes (n = 9 enumerations ~5 min each).
set -e
cd "$(dirname "$0")"

## 1. Session pipeline (PROOF.md section 6), compared byte for byte with the committed files
gcc -O2 -o dom_even_audit ../../paleyeven_20261109/dom_even.c
( time ./dom_even_audit 511 16 > dom_K511_regen.out 2> dom_K511_regen.err ) > census_time.txt 2>&1
shasum -a 256 dom_K511_regen.out >> census_time.txt              # = sha256.txt: 6307dfa1...
diff dom_K511_regen.err ../dom_K511.err
python3 ../mkbase.py dom_K511_regen.out 9 > base9f_regen.txt && cmp base9f_regen.txt ../base9f.txt
awk '{print $1,$2,$3,$4,$5}' base9f_regen.txt > base9_regen.txt && cmp base9_regen.txt ../base9.txt
gcc -O2 -o evenclos_a ../evenclos.c && gcc -O2 -o closure4_a ../closure4.c && gcc -O2 -o wcheck_a ../wcheck.c
./evenclos_a 26 base9f_regen.txt 9 12 > evenclos_full_regen.out && cmp evenclos_full_regen.out ../evenclos_full.out
./evenclos_a 26 base9f_regen.txt 9 12 noR > evenclos_noR_regen.out && cmp evenclos_noR_regen.out ../evenclos_noR.out
./closure4_a 26 base9_regen.txt 9 > closure4_regen.out && cmp closure4_regen.out ../closure4.out
./evenclos_a 24 base9f_regen.txt 9 12 full 81 > gamma_dump_regen.out && cmp gamma_dump_regen.out ../gamma_dump.out
python3 ../seedcases.py gamma_dump_regen.out                    # same counts as ../seedcases.out
( time python3 ../witness.py dom_K511_regen.out 9 22 witnesses_regen.txt 12 > witness_regen.out ) > witness_time.txt 2>&1
shasum -a 256 witnesses_regen.txt >> witness_time.txt && cmp witnesses_regen.txt ../witnesses.txt
./wcheck_a < witnesses_regen.txt                                 # ok=6147 bad=0
./wcheck_a < negctl.txt || true                                  # negative controls: K_2, P_4, wrong value -> 3 BAD
( cd .. && python3 lemmacheck.py && python3 chaincheck.py )      # 706 / 217 tests, 0 failures

## 2. Independent verifier (d by inclusion-exclusion; parity from a refinement-based stabiliser chain)
gcc -O2 -Wall -o aw_verify aw_verify.c
./aw_verify selftest 30000                                       # vs brute force over all permutations: 0 failures
python3 structured_tests.py && ( ./aw_verify bf < structured_tests.txt > structured_tests.out || true )   # 200 AGREE
./aw_verify < ../witnesses.txt > aw_verify_witnesses.out         # ok=6147 bad=0
python3 aw_coverage.py aw_verify_witnesses.out > aw_coverage.out # ALL COVERED
./aw_verify < negctl.txt || true

## 3. Exhaustive enumeration of all graphs on <= 9 vertices (degree-sorted labelled representatives)
gcc -O3 -Wall -Wno-unused-function -o aw_enum aw_enum.c
./aw_enum 8 > aw_enum8.out 2> aw_enum8.err
( time ./aw_enum 9 9 > aw_enum9.out 2> aw_enum9.err ) > aw_enum9.time 2>&1
diff <(cat aw_enum8.out aw_enum9.out | sort) <(sort ../base9f.txt) > aw_enum9.diff   # empty: base9f.txt reproduced
./aw_enum 8 1 rep > aw_enum8_rep.out 2>/dev/null                 # representatives for the closure

## 4. Auditor's closure from the auditor's base on <= 8 vertices
python3 aw_closure.py aw_enum8_rep.out 8 22 12 aw_closure8_22_w.txt 12 > aw_closure8_22.out
./aw_verify < aw_closure8_22_w.txt > aw_closure8_22_verify.out
python3 aw_closure.py aw_enum8_rep.out 8 26 12 aw_closure8_26_w.txt 4 > aw_closure8_26.out
./aw_verify < aw_closure8_26_w.txt > aw_closure8_26_verify.out

## 5. Graphs built exactly as the proof prescribes (section 5, Lemma 8 cases at J = 13)
python3 aw_construct.py ../witnesses.txt aw_construct_graphs.txt 20 25 > aw_construct.out
./aw_verify < aw_construct_graphs.txt > aw_construct_verify.out

## 6. Arithmetic of the case analyses (J, n <= 400)
python3 aw_arith.py 400 > aw_arith.out

## 7. Prop 3.3 over its candidate class, odd n <= 15
python3 aw_prop33.py 15 > aw_prop33_graphs.txt && ./aw_verify < aw_prop33_graphs.txt > aw_prop33_verify.out
python3 aw_prop33_summary.py aw_prop33_verify.out > aw_prop33.out

## 8. Small-d census (Prop 3.2, Remark 7.1, HYP-9625 evidence)
gcc -O3 -Wall -o aw_dump aw_dump.c
( time ./aw_dump 9 1 3 5 7 11 13 23 25 37 47 > aw_dump_smalld.out 2> aw_dump_smalld.err ) > aw_dump.time 2>&1
python3 aw_smalld.py aw_dump_smalld.out > aw_smalld.out

## 9. Remark 7.1: four-move closure from the auditor's base on <= 9 vertices
cat aw_enum8.out aw_enum9.out > aw_base9_tuples.txt
python3 aw_four.py aw_base9_tuples.txt 9 26 > aw_four.out

## 10. Remark 7.2: top-band values realised by unions
python3 aw_union.py ../witnesses.txt > /dev/null && ./aw_verify < aw_union_topband.txt > aw_union_topband.out

## 11. HYP-9625
python3 aw_hyp.py > aw_hyp_facts.out && ./aw_verify < aw_hyp_graphs.txt > aw_hyp_verify.out || true   # parity column: ODD
gcc -O3 -Wall -Wno-unused-function -o aw_oddvals aw_oddvals.c
( time ./aw_oddvals 9 > aw_oddvals.out ) > aw_oddvals.time 2>&1

## cleanup of binaries and large regenerable files
rm -f dom_even_audit evenclos_a closure4_a wcheck_a aw_verify aw_enum aw_dump aw_oddvals dom_K511_regen.out \
      witnesses_regen.txt base9f_regen.txt base9_regen.txt evenclos_full_regen.out evenclos_noR_regen.out \
      closure4_regen.out gamma_dump_regen.out aw_closure8_2?_w.txt aw_construct_graphs.txt aw_prop33_graphs.txt \
      aw_base9_tuples.txt
