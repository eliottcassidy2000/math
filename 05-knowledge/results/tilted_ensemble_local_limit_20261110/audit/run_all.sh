#!/bin/sh
# AUDIT: rerun every independent check of audit/AUDIT.md.  Scratch files and binaries go to $W (default
# /tmp/claude-audit-ensemble); text outputs are rewritten next to this script.  Needs cc, python3 and nauty's gentourng.
# Total time on 4 cores: about 10 minutes (the order-10 passes dominate).  Usage: sh run_all.sh [quick]
#   quick: skip every order-10 pass (then d_window_analyze.py is run to order 9 only).
set -e
A=$(cd "$(dirname "$0")" && pwd); W=${W:-/tmp/claude-audit-ensemble}; G=${GENTOURNG:-gentourng}; QUICK=$1
mkdir -p "$W"; cd "$A"
for p in a_labjoint b_thmC b_classdisc b_rankonly d_window d_witness d_single_realiser e_oneapex; do cc -O2 -o "$W/$p" "$p.c"; done

# (a) Theorems A and B
python3 -B a_thmAB.py > a_thmAB.out
for n in 3 4 5 6 7; do "$W/a_labjoint" $n > "$W/labjoint_$n.txt"; done
python3 -B a_labjoint_analyze.py "$W" > a_labjoint_analyze.out
python3 -B a_prior_art_links.py > a_prior_art_links.out

# proofs (insertion argument, Lemma C4, Propositions 5.1 / 6.1, alternating matrices)
python3 -B p_proofchecks.py > p_proofchecks.out

# (b) Theorem C
(for m in 2 3 4 5 6 7; do $G -q $m | "$W/b_thmC" $m direct; done) > b_thmC_classes_2_7.out
(for m in 2 3 4 5 6; do "$W/b_thmC" $m labelled; done) > b_thmC_labelled_2_6.out
$G -q 8 | "$W/b_thmC" 8 direct > b_thmC_classes_8.out
$G -q 9 | "$W/b_thmC" 9 fast > b_thmC_classes_9.out
(for n in 3 4 5 6 7 8 9; do $G -q $n | "$W/b_classdisc" $n; done) > b_classdisc_3_9.out
python3 -B b_digit1_link.py > b_digit1_link.out

# (c), (d) histograms and deletion pairs, orders 2..9
(for m in 2 3 4 5 6 7 8 9; do $G -q $m | "$W/d_window" $m "$W/win_$m"; done) > d_window_runs.out
python3 -B c_hist_analyze.py 9 "$W/win_9.hist" > c_hist_9.out
for m in 2 3 4 5 6 7 8 9; do $G -q $m | "$W/d_single_realiser" $m > "$W/single_$m.txt"; done
python3 -c "
import random
rng = random.Random(20261110); out = []
for _ in range(5):
    h = 237 + 2 * rng.randrange((400 - 237) // 2 + 1); n = 3 * h + 2 * rng.randrange(h)
    out.append((h, n))
open('d_targets_random.txt', 'w').write(''.join('%d %d\n' % t for t in out))"
$G -q 9 | "$W/d_witness" 9 d_targets_random.txt > d_witnesses_random.txt
$G -q 9 | "$W/d_witness" 9 d_targets_extra.txt > "$W/wit_extra_9.txt" 2> "$W/wit_extra_9.err" || true

if [ "$QUICK" = quick ]; then
  python3 -B d_window_analyze.py "$W" 9 > d_window_analyze.out
  cp "$W/wit_extra_9.txt" d_witnesses_extra.txt
else
  for r in 0 1 2 3; do ( $G -q 10 $r/4 | "$W/d_window" 10 "$W/win_10_$r" > "$W/win_10_$r.log" ) & done; wait
  cat "$W"/win_10_?.log >> d_window_runs.out
  python3 -B d_window_analyze.py "$W" 10 > d_window_analyze.out
  python3 -B c_hist_analyze.py 10 "$W"/win_10_?.hist > c_hist_10.out
  python3 -B d_spotchecks.py "$W" > d_spotchecks.out
  for r in 0 1 2 3; do ( $G -q 10 $r/4 | "$W/d_single_realiser" 10 > "$W/single_10_$r.txt" ) & done; wait
  sed -n 's/NOT FOUND at order 9: //p' "$W/wit_extra_9.err" > "$W/targets_extra_10.txt"
  $G -q 10 | "$W/d_witness" 10 "$W/targets_extra_10.txt" > "$W/wit_extra_10.txt"
  cat "$W/wit_extra_9.txt" "$W/wit_extra_10.txt" > d_witnesses_extra.txt
  for r in 0 1 2 3; do ( $G -q 10 $r/4 | "$W/b_classdisc" 10 > "$W/classdisc_10_$r.out" ) & done; wait
  cat "$W"/classdisc_10_?.out > b_classdisc_10_parts.out
  for r in 0 1 2 3; do ( $G -q 10 $r/4 | "$W/b_rankonly" 10 > "$W/rank10_$r.out" ) & done; wait
  cat "$W"/rank10_?.out > b_rankonly_10_parts.out
fi
python3 -B d_single_analyze.py "$W" > d_single_realiser.out
(python3 -B d_witness_verify.py < d_witnesses_random.txt; python3 -B d_witness_verify.py < d_witnesses_extra.txt) > d_witness_verify.out

# (e) section 4.2
"$W/e_oneapex" 8 > e_oneapex.out
echo "done; see the *.out files in $A"
