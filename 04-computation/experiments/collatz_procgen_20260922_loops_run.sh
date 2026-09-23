#!/usr/bin/env bash
# Reproduce 05-knowledge/results/collatz_procgen_20260922_loops.out
#   bash 04-computation/experiments/collatz_procgen_20260922_loops_run.sh > 05-knowledge/results/collatz_procgen_20260922_loops.out
# About 15 minutes on one core; peak memory about 700 MB (section 1: 360 MB; the s=2300 base loop: 650 MB).
set -euo pipefail
cd "$(dirname "$0")"
P=collatz_procgen_20260922
T=$(mktemp -d)
for f in loops_through_one loops_dp loops_height loops_recon1 loops_hubs loops_seesaw; do
  clang -O3 -o $T/$f ${P}_$f.c -lm
done
echo "### 0. baseline (lane one program, values <= 2e7), extended to s <= 80"
$T/loops_through_one 80 | tail -40
echo
echo "### 1. lexicographic DP (min K, then min height) for loops through 1, s <= 6000, values <= 3e7"
$T/loops_height 6000 30000000 > $T/height.txt
python3 - "$T/height.txt" << 'EOF'
import re, sys, math, statistics
rows = []
for line in open(sys.argv[1]):
    m = re.match(r's=\s*(\d+) K0=\s*(\d+) minK=\s*(\d+) c=([\d.]+) H=(\d+)', line)
    if m: rows.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4)), int(m.group(5))))
bad = [r[0] for r in rows if r[0] >= 2 and r[2] != r[1]]
print(f"rows {len(rows)}; s>=2 with minK != K0(s): {bad if bad else 'none'};  s=1: minK={rows[0][2]} (trivial loop)")
print("largest minimal heights H(s) (all s <= 6000), with eps = log2(c/1.5) and H*eps*3ln2/s (proved lower bound ~ 1):")
for s, K0, K, c, H in sorted(rows, key=lambda r: -r[4])[:12]:
    eps = math.log2(c / 1.5)
    print(f"  s={s:5d} c={c:.6f} eps={eps:.6f} H={H:9d}  H*eps*3ln2/s={H*eps*3*math.log(2)/s:.2f}")
for lo in range(0, 6000, 1000):
    hs = [r[4] for r in rows if lo < r[0] <= lo + 1000]
    print(f"  s in ({lo},{lo+1000}]: median H = {statistics.median(hs):.0f}, max H = {max(hs)}")
with open(sys.argv[1] + '.map', 'w') as f:
    for r in rows: f.write(f"{r[0]} {r[4]}\n")
EOF
echo
echo "### 2. minimal-height loops for s <= 80 (checkpointed reconstruction at value cap H(s)), independent check, table"
: > $T/minh80.txt
for s in $(seq 1 80); do H=$(awk -v s=$s '$1==s {print $2}' $T/height.txt.map); $T/loops_recon1 $s $H 8 >> $T/minh80.txt; done
python3 ${P}_loops_verify.py $T/minh80.txt --quiet
echo
echo "### 3. hub cycles: minimal (q,p) through climb points; records (smallest eps = p - q log2 3 so far)"
for h in 4 13 40 121; do
  $T/loops_hubs $h 120 3000000 > $T/hub_$h.txt
  python3 - $h $T/hub_$h.txt << 'EOF'
import re, sys
best = 9; recs = []
for line in open(sys.argv[2]):
    m = re.search(r'q=\s*(\d+) p=\s*(\d+) ratio=([\d.]+) eps=p-q\*log2\(3\)=([-\d.]+)', line)
    q, p, r, e = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
    if e < best: best = e; recs.append((q, p, round(r, 6)))
print(f"  hub {sys.argv[1]}: record cycles (q,p,ratio) for q <= 120: {recs}")
EOF
done
python3 - << 'EOF' > $T/records.txt
best = 10.0; import math
for q in range(1, 2600):
    p = (3 ** q).bit_length(); e = p - q * math.log2(3)
    if e < best: best = e; print(q, p, e)
EOF
echo "  record denominators q <= 2600 (upper best approximations of log2 3): $(awk '{printf "%s ", $1}' $T/records.txt)"
: > $T/hubcyc_1093.txt
for q in 17 29 41 94 147 200 253 306; do
  B=$(python3 -c "import math;print(max(4,int(math.sqrt($q))))")
  $T/loops_recon1 $q 4000000 $B 1093 | sed "s/^s=/q=/" >> $T/hubcyc_1093.txt
done
python3 - $T/hubcyc_1093.txt << 'EOF'
import re, sys
for line in open(sys.argv[1]):
    m = re.search(r'q=(\d+).*minK=(\d+).*k:(\S+)', line)
    q, p = int(m.group(1)), int(m.group(2)); ks = [int(t) for t in m.group(3).split(',')]
    x = 1093; mx = x
    for k in ks:
        t = (1 << k) * x - 1; assert t % 3 == 0; x = t // 3; assert x >= 1 and x % 3; mx = max(mx, x)
    assert x == 1093 and sum(ks) == p and len(ks) == q and p == (3 ** q).bit_length()
    print(f"  cycle through 1093 verified: (q,p)=({q},{p}), p = ceil(q log2 3), ratio 2^p/3^q = {2**p/3**q:.6f}, max value {mx}")
EOF
echo
echo "### 4. record base loops (s = q-1 for records q <= 2301, and s = 2) and the splicing family, 2 <= s <= 2500"
echo "s=2 k:2,2" > $T/base.txt
for q in $(awk '$1>=5 && $1<=2301 {print $1}' $T/records.txt); do
  s=$((q-1)); H=$(awk -v s=$s '$1==s {print $2}' $T/height.txt.map)
  B=$(python3 -c "import math;print(max(8,int(math.sqrt($s))))")
  $T/loops_recon1 $s $H $B >> $T/base.txt
done
python3 ${P}_loops_verify.py $T/base.txt
python3 ${P}_loops_family.py 2500 $T/base.txt $T/family2500r.txt $T/hubcyc_1093.txt
python3 ${P}_loops_verify.py $T/family2500r.txt --quiet
echo "loop table, s <= 80 (minimal-height loop; construction of the family loop):"
python3 ${P}_loops_table.py $T/minh80.txt $T/family2500r.txt 80
echo "family constructions for s = 81..160:"
awk -F' ' '{split($1,a,"="); if (a[2]>80 && a[2]<=160) printf "%s:%s ", a[2], $2} END {print ""}' $T/family2500r.txt | fold -w 200
echo
echo "### 5. hostile points 1/2 and dyadic threads: Bad_inf proofs, transfer lemmas on integers, prices"
$T/loops_seesaw 12 $T/half12.txt > $T/seesaw12.txt
python3 ${P}_loops_escapes.py $T/family2500r.txt $T/half12.txt
echo
echo "### 6. see-saw: exact certificates at level r = 12 (classes mod 3^13)"
cat $T/seesaw12.txt
rm -rf $T
