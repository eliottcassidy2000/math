#!/usr/bin/env bash
# Reproduce 05-knowledge/results/collatz_procgen_20260922_mirror.out (Q1 mirror lane: loops through -1,
# the -1 escape, hostile chains, the level-24 see-saw, the endgame reduction, the forward/backward duality).
#   bash 04-computation/experiments/collatz_procgen_20260922_mirror_run.sh > 05-knowledge/results/collatz_procgen_20260922_mirror.out
# FULL (default): 38 min on one core (measured); peak memory 0.55 GB (the loop DP, section 1); section 9
#   (divide-and-conquer reconstruction of B_1054 and C_1054, family to K <= 4000) uses 0.23-0.36 GB.
# QUICK=1: K <= 1200 with cap 2e6, chains to 3e7, no section 9 (a few minutes, < 0.1 GB); the duality
#   section still needs the dimension lane's dumps.
# DIM: directory with the dimension lane's dumps fwd30/{fwd_bad_m61,fwd_bad_m63,hostile_z13_M61_J25}.txt and
#   bwd17/bwd_bad_r41.txt (default: scratch/procgen_dim at the repository root; regenerate them with
#   OUTDIR=<dir> FULL=1 bash 04-computation/experiments/collatz_procgen_20260922_dim_run.sh).
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
P=collatz_procgen_20260922
T=$(mktemp -d)
DIM=${DIM:-$ROOT/scratch/procgen_dim}
for f in mirror_loops_dp mirror_recon1 mirror_seesaw mirror_chain dim_alpha; do
  clang -O3 -o $T/$f ${P}_$f.c -lm
done
if [ "${QUICK:-0}" = 1 ]; then KMAX=1200; HCAP=2000000; CHAIN_N=30000000; else KMAX=4000; HCAP=30000000; CHAIN_N=4294967296; fi

echo "### 0. loop equation through -1, exhaustive for K <= 12, a <= 12 (Prop 1.1, Prop 1.2, exceptions)"
python3 ${P}_mirror_loops_verify.py
echo
echo "### 1. layered DP (min a, then min height) for loops through -1, K <= $KMAX, multiplication points <= $HCAP, pruning slack 2"
$T/mirror_loops_dp $KMAX $HCAP 2 > $T/dp.txt
python3 - "$T/dp.txt" << 'EOF'
import re, sys, math, statistics
rows = {}
for line in open(sys.argv[1]):
    m = re.match(r'K=\s*(\d+) a0=\s*(\d+) amin=\s*(\d+) ratio=([\d.]+) eta=([\d.]+) H=(\d+)', line)
    if m: rows[int(m.group(1))] = tuple(float(g) if '.' in g else int(g) for g in m.groups()[1:])
K1 = max(rows)
print(f"rows {len(rows)} (K = 1..{K1}); K with amin != a0(K): {[K for K in rows if rows[K][1] != rows[K][0]]}")
print("  (K=1: the basic loop MH, a=1 < a0; K=5..9: a=K > a0, proved optimal in section 0)")
print(f"  max height H(K) = {max(r[4] for r in rows.values())}; the cap is not binding (all H <= 2/3 of the state cap)")
ln3 = math.log(3)
print("hard cases eta(K+1) < 0.006: K+1, K, a0, ratio 3^a0/2^K, eta, H(K), proved lower bound (Prop 3.1), H*eta*3ln3/a0:")
def lb(a, eta):
    H = 1
    while True:
        rhs = (a - 2 - math.log(2 * H, 3)) / (3 * eta * ln3)
        if H > rhs: return H
        H = max(H + 1, int(rhs))
for K in sorted(rows):
    a0, am, r, eta, H = rows[K]
    if K >= 10 and eta < 0.006:
        print(f"  {K+1:5d} {K:5d} {a0:5d} {r:.6f} {eta:.7f} {H:10d} {lb(a0, eta):9d}  {H*eta*3*ln3/a0:.3f}")
for lo in range(0, K1, 1000):
    hs = [rows[K][4] for K in rows if lo < K <= lo + 1000 and K >= 10]
    if hs: print(f"  K in ({lo},{lo+1000}]: median H = {statistics.median(hs):.0f}, max H = {max(hs)}")
print("full table K <= 100 (K a0 amin ratio eta H):")
for K in range(1, 101):
    a0, am, r, eta, H = rows[K]; print(f"  {K:3d} {a0:3d} {am:3d} {r:.6f} {eta:.5f} {H}")
with open(sys.argv[1] + '.H', 'w') as f:
    for K in sorted(rows): f.write(f"{K} {rows[K][4]}\n")
EOF
echo
echo "### 2. minimal-height loops for K <= 80 (checkpointed reconstruction at cap H(K)), independent verification"
: > $T/minh.txt
for K in $(seq 2 80); do
  H=$(awk -v k=$K '$1==k {print $2}' $T/dp.txt.H)
  if [ $K -ge 5 ] && [ $K -le 9 ]; then SL=-1; else SL=2; fi
  $T/mirror_recon1 1 $K $H $SL 0 recon >> $T/minh.txt
done
grep -c 'RECON_CHECK.*OK' $T/minh.txt | sed 's/^/  reconstructions OK: /'
python3 - "$T/minh.txt" << 'EOF'
import re, sys, importlib.util
spec = importlib.util.spec_from_file_location("v", "collatz_procgen_20260922_mirror_loops_verify.py")
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
ok = n = 0
for line in open(sys.argv[1]):
    m = re.match(r'RECON h=-1 P=(\d+) q=(\d+) word:(.*)', line)
    if not m: continue
    K = int(m.group(1)); w = v.parse_word(m.group(3)); n += 1
    g, res = v.verify_word(w, K)
    exc = 5 <= K <= 9 and res['legal'] and res['closed'] and res['ident'] and res['cost'] and res['a'] == K
    ok += (g or exc)
    print(f"  K={K:2d} a={res['a']:2d} H={res.get('H')}: {m.group(3).strip()}")
print(f"  independently verified: {ok}/{n} (a = a0(K) for K >= 10 and K = 2,3,4; a = K for K = 5..9)")
EOF
echo
echo "### 3. hub cycles through the climb points -(3^j+1)/2; eta-records (lower best approximations of log2 3)"
for h in 1 5 14 41 122 365 1094; do $T/mirror_recon1 $h 100 300000 3 1 scan > $T/hub_$h.txt; done
for h in 3281 9842; do $T/mirror_recon1 $h 600 2000000 2 1 scan > $T/hub_$h.txt; done
python3 - $T << 'EOF'
import re, sys, math
th = math.log(2) / math.log(3)
for h in (1, 5, 14, 41, 122, 365, 1094, 3281, 9842):
    best = 9; out = []
    for line in open(f"{sys.argv[1]}/hub_{h}.txt"):
        m = re.search(r'p=\s*(\d+) q=\s*(\d+) ratio=([\d.]+)', line)
        if not m: continue
        p, q = int(m.group(1)), int(m.group(2)); d = q - p * th
        if d < best - 1e-12: best = d; out.append(f"({p},{q},{float(m.group(3)):.6f})")
    print(f"  hub -{h}: successive minima of eta over cycles (p halvings, q = mult.): {' '.join(out)}")
recs = []; best = None; a = 0; p3 = 1; p2 = 1
for n in range(1, 60001):
    p2 *= 2
    while p3 <= p2: p3 *= 3; a += 1          # a = ceil(n log3 2), incrementally
    if best is None or p3 * best[2] < best[1] * p2: recs.append((n, a)); best = (a, p3, p2)
print("  eta-records n <= 60000 (n/a with a = ceil(n log3 2); lower best approximations of log2 3): "
      + " ".join(f"{n}/{a}" for n, a in recs))
EOF
echo
echo "### 4. record base loops B_n (K = n-1) and record cycles C_n; the splicing family (Theorem 2.2), K <= 1100"
( $T/mirror_recon1 1 10 43 2 0 recon; $T/mirror_recon1 1 18 547 2 0 recon; $T/mirror_recon1 1 83 23947 2 0 recon
  $T/mirror_recon1 1 568 458899 2 0 recon 24
  $T/mirror_recon1 1 1 10 -1 1 recon; $T/mirror_recon1 5 3 100 -1 1 recon; $T/mirror_recon1 14 11 100000 2 1 recon
  $T/mirror_recon1 41 19 100000 2 1 recon; $T/mirror_recon1 365 84 300000 2 1 recon
  $T/mirror_recon1 1094 569 2000000 2 1 recon ) > $T/recon_all.txt
grep RECON_CHECK $T/recon_all.txt | sed 's/^/  /'
grep 'word:' $T/recon_all.txt | cut -c1-240 | sed 's/^/  /'
python3 ${P}_mirror_loops_family.py $T/recon_all.txt 1100
echo
echo "### 5. the escape near -1: exact prices, exit lemma on integers, landing classes, canonical chains"
python3 ${P}_mirror_escape.py 40 | tee $T/escape.txt
args=""; for m in $(seq 2 40); do args="$args $(( (1<<m) - 1 )) 1"; done
$T/dim_alpha 22 41 $args > $T/dim_alpha.txt
python3 - "$T/dim_alpha.txt" "$T/escape.txt" << 'EOF'
import re, sys
esc = {int(m.group(1)): int(m.group(2)) for m in (re.search(r'm=\s*(\d+) a0\(m\)=\s*\d+ A\*=\s*(\d+)', l) for l in open(sys.argv[2])) if m}
bad = n = 0
for l in open(sys.argv[1]):
    if not l.startswith('h='): continue
    num = int(l.split('=')[1].split('/')[0]); vals = list(map(int, l.split(':')[1].split())); m = (num + 1).bit_length() - 1
    n += 1; bad += (vals[m] != esc.get(m))
print(f"  cross-check with the dimension lane's dim_alpha (f_(m+1)(2^m - 1), independent code): {n} values, mismatches {bad}")
EOF
echo
echo "### 6. level-24 forward certificates (mirror of the level-12 see-saw): classes mod 2^24"
ML=$(grep 'seesaw argument list' $T/escape.txt | sed 's/.*list: //')
$T/mirror_seesaw 24 "$ML"
echo
echo "### 7. shortest descents from the -1 thread (mirror of half_chain.c): all n = 7 mod 8, n <= $CHAIN_N"
$T/mirror_chain $CHAIN_N 1000 120 22
echo
echo "### 8. backward prover, forward/backward numerator duality, generation-1 clocks"
python3 ${P}_mirror_bwd_prover.py 1/2 1 43/32 59/64 145/128 209/256 371/256 499/512 1/4 5/4 7/8
python3 ${P}_mirror_duality.py $DIM
echo
if [ "${QUICK:-0}" != 1 ]; then
echo "### 9. the record base loop B_1054 (height 1.9e7) and cycle C_1054 by divide-and-conquer reconstruction; family to K <= 4000"
clang -O3 -o $T/mirror_mitm ${P}_mirror_mitm.c -lm
$T/mirror_mitm 1 1053 19012579 2 0 2> $T/mitm_b.log > $T/B1054.txt
grep RECON_CHECK $T/B1054.txt | sed 's/^/  /'
A1=$(python3 -c "import re,sys; w=open('$T/B1054.txt').read().split('word:')[1].split()[0]; print(int(re.match(r'M(\d+)',w).group(1)))")
echo "  B_1054 begins with a climb of $A1 multiplications; scanning climb hubs for a (1054,665) cycle"
HUB=""
for j in $(seq $((A1 - 1)) -1 6); do
  h=$(( (3**j + 1) / 2 ))
  $T/mirror_recon1 $h 1054 30000000 2 1 scan > $T/hub1054_$h.txt
  if grep -qE 'p=1054 q= 665 ' $T/hub1054_$h.txt; then HUB=$h; echo "  hub -$h: (1054,665) realized"; break; else echo "  hub -$h: no (1054,665) cycle within the cap"; fi
done
if [ -n "$HUB" ]; then
  $T/mirror_mitm $HUB 1054 30000000 2 1 2> $T/mitm_c.log > $T/C1054.txt
  grep RECON_CHECK $T/C1054.txt | sed 's/^/  /'
  cat $T/recon_all.txt $T/B1054.txt $T/C1054.txt > $T/recon_all2.txt
  python3 ${P}_mirror_loops_family.py $T/recon_all2.txt 4000
fi
fi
echo "done"
