#!/usr/bin/env bash
# collatz_procgen_20260923_sweep_run.sh -- reproduces collatz_procgen_20260923_sweep.out
# (hypothesis sweep HYP-9120..9127, 9132; note collatz_procgen_20260923_hypothesis_sweep.md)
#
#   bash 04-computation/experiments/collatz_procgen_20260923_sweep_run.sh > 05-knowledge/results/collatz_procgen_20260923_sweep.out
#
# REUSE=1 : do not recompute the record objects; print their logs and re-verify the stored objects
#           (objects in $WORKDIR/objects, produced by an earlier full run).
# DEEP=1  : also run the deeper certificate searches for the open generation-1 points (about 1 h).
# QUICK=1 : records only up to q=15601 (pos) / n=25781 (neg), E_S DFS to m=55, generation-1 sweep to i=120.
# Peak memory < 150 MB per process; checkpoint files (up to ~1 GB) under $WORKDIR/ck, deleted per object.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/04-computation/experiments"
W="${WORKDIR:-$ROOT/scratch/procgen_sweep}"
mkdir -p "$W/ck" "$W/objects"
export TIGHTDP="$W/tightdp_sweep" CKDIR="$W/ck" OBJDIR="$W/objects" FB=2.0 FC=3.0 RT=64
cc -O2 -o "$W/tightdp_sweep" "$EXP/collatz_procgen_20260923_sweep_tightdp.c" -lm
cc -O2 -o "$W/es_forward_sweep" "$EXP/collatz_procgen_20260923_sweep_es_forward.c" -lm
cc -O2 -o "$W/rn_cover_sweep" "$EXP/collatz_procgen_20260923_sweep_rn_cover.c" -lm
QMAXP=79335; QMAXN=50508; ESM=64; GI=300
if [ "${QUICK:-0}" = "1" ]; then QMAXP=15601; QMAXN=25781; ESM=55; GI=120; fi

echo "### S1. HYP-9122 record objects (sign +): B_q, C_q for every record q <= $QMAXP (tight DP, reconstruction, verification)"
if [ "${REUSE:-0}" = "1" ]; then cat "$OBJDIR/pos_log.jsonl"
else rm -f "$OBJDIR/pos_log.jsonl"; python3 "$EXP/collatz_procgen_20260923_sweep_records.py" pos $QMAXP --recon > /dev/null; cat "$OBJDIR/pos_log.jsonl"; fi
echo "### S2. HYP-9125 record objects (sign -): B_n (n>=11), C_n for every record n <= $QMAXN"
if [ "${REUSE:-0}" = "1" ]; then cat "$OBJDIR/neg_log.jsonl"
else rm -f "$OBJDIR/neg_log.jsonl"; python3 "$EXP/collatz_procgen_20260923_sweep_records.py" neg $QMAXN --recon > /dev/null; cat "$OBJDIR/neg_log.jsonl"; fi
echo "### S3. independent re-verification of every stored object; tightness lemmas T1-T3 (records <= 6000); hub hypothesis"
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" pos reverify
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" neg reverify
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" pos tightcheck
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" neg tightcheck
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" pos hubcheck
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" neg hubcheck
echo "### S4. explicit spliced loops (Theorem 4.2 / 2.2) at sample lengths, each verified exactly"
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" pos 2,3,40,305,6001,6290,12345,15600,20000,47466,47467,79333,79334,120000,190535
python3 "$EXP/collatz_procgen_20260923_sweep_splice.py" neg 10,18,4001,25779,25780,40000,50506,50507,100000,176249
echo "### S5. HYP-9121: E control (J=1,S={0}) and E_S (6 mod 8) thread counts and projections"
"$W/es_forward_sweep" 26 50 1 0 | grep -E "^m=(30|31|32|35|40|43|46|50) |^P m=(16|20|24|25|28|30|32|35|36|40|44|45) "
"$W/es_forward_sweep" 26 $ESM 3 6 | grep -E "^m="
ES_DUMP="$W/es_bad_m55.txt" "$W/es_forward_sweep" 26 55 3 6 | grep -E "^P m=(20|25|30|35|40|45|50) "
"$W/es_forward_sweep" 26 59 3 6 | grep -E "^P m=(20|25|30|35|40|45|50) "
"$W/es_forward_sweep" 26 50 3 6 | grep -E "^P m=(20|25|30|35|40|45) "
echo "### S6. HYP-9121 refutation: Lemma F, Lemma G, loops through -1 in E_S, Theorem R, cross-checks"
ES_FILTER=0 python3 "$EXP/collatz_procgen_20260923_sweep_es_refute.py" "$W/es_bad_m55.txt"
echo "### S7. HYP-9126 generation-1 points above 3/2 (i <= $GI)"
python3 "$EXP/collatz_procgen_20260923_sweep_gen1.py" $GI 3000
if [ "${DEEP:-0}" = "1" ]; then
  echo "### S7b. deeper searches for the open generation-1 points"
  python3 "$EXP/collatz_procgen_20260923_sweep_gen1.py" deep-back 159 40000
  python3 "$EXP/collatz_procgen_20260923_sweep_gen1.py" deep-fwd 45,64,148,213,232,278,297 1500 60000
fi
echo "### S8. chain-sum coverage R_s (note section 3.5)"
"$W/rn_cover_sweep" 16 3.0
