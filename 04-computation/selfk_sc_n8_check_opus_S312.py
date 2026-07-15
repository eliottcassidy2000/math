# opus-2026-07-15-S312 -- HYP-6895: the 2*selfK = SC law at n=8.
# Law (S311): #{non-grid-sym tilings t : kappa(t) iso t} = SC(n) for n >= 5.
# Prediction at n=8: black quasi-fixed tilings = SC(8) = 176 (A000570).
# Also computes the Aut-WEIGHTED count Sum_{qf t} |Aut(t)| (the pair-counting
# quantity M/n!) to decide which version is the true law.
# Reuses the S308 subH-certified classifier (6880 classes exact).
import sys, time
from collections import defaultdict
sys.path.insert(0, '04-computation')

t0 = time.time()
n = 8
m = n*(n-1)//2 - (n-1)  # 21 tiles
code = open('04-computation/coarse_concordance_n8_bothsides_opus_S308.py',
            encoding='utf-8', errors='surrogateescape').read()
core = code.split('# discordances by level distance')[0]
ns = {'__file__': '04-computation/coarse_concordance_n8_bothsides_opus_S308.py'}
exec(core, ns)
cls_of = ns['cls_of']
NC = max(cls_of) + 1
print(f"classifier reused: {NC} classes ({time.time()-t0:.0f}s)", flush=True)
assert NC == 6880, "classifier not certified!"

# tiles and maps -----------------------------------------------------------
tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
tidx = {t: i for i, t in enumerate(tiles)}
assert len(tiles) == m
sig = [tidx[(n+1-y, n+1-x)] for (x, y) in tiles]     # grid reflection sigma
FULL = (1 << m) - 1

def sig_t(t):
    s = 0
    for i in range(m):
        if (t >> i) & 1: s |= 1 << sig[i]
    return s

# quasi-fixed scan ---------------------------------------------------------
qf_total = qf_black = qf_blue = 0
qf_by_class = defaultdict(int)          # class -> #qf tilings
black_by_class = defaultdict(int)
for t in range(1 << m):
    if cls_of[t] == cls_of[t ^ FULL]:
        qf_total += 1
        c = cls_of[t]
        qf_by_class[c] += 1
        if sig_t(t) == t: qf_blue += 1
        else:
            qf_black += 1
            black_by_class[c] += 1
print(f"scan done ({time.time()-t0:.0f}s)", flush=True)

# SC(8) via sigma-realizes-reversal: rev(cls(t)) = cls(sig_t(t)) ------------
rcls = {}
seen = set()
for t in range(1 << m):
    c = cls_of[t]
    if c not in seen:
        seen.add(c)
        rcls[c] = cls_of[sig_t(t)]
SC = sum(1 for c, rc in rcls.items() if rc == c)
print(f"SC(8) via sigma = {SC} (A000570 says 176)", flush=True)

# Aut weights for carrier classes: |Aut| = H / fiber ------------------------
fiber = defaultdict(int)
rep = {}
for t in range(1 << m):
    c = cls_of[t]
    fiber[c] += 1
    if c not in rep: rep[c] = t

def ham_paths(t):
    # arcs: base path k -> k-1; tile (x,y) bit1 = x -> y else y -> x
    adj = [[False]*n for _ in range(n)]
    for k in range(2, n+1): adj[k-1][k-2] = True
    for i, (x, y) in enumerate(tiles):
        if (t >> i) & 1: adj[x-1][y-1] = True
        else: adj[y-1][x-1] = True
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            d = dp[S][v]
            if not d: continue
            for w in range(n):
                if not (S >> w) & 1 and adj[v][w]:
                    dp[S | (1 << w)][w] += d
    return sum(dp[(1 << n) - 1])

weighted = 0
carrier_rows = []
for c, k in sorted(qf_by_class.items()):
    H = ham_paths(rep[c])
    aut = H // fiber[c]
    assert aut * fiber[c] == H
    weighted += k * aut
    carrier_rows.append((c, k, black_by_class.get(c, 0), H, aut,
                         'SC' if rcls[c] == c else 'NS'))

print(f"\n=== n=8 RESULTS ===")
print(f"quasi-fixed total  = {qf_total}")
print(f"  black (non-gs)   = {qf_black}   [LAW predicts SC = {SC}]")
print(f"  blue  (grid-sym) = {qf_blue}   (selfB lines = {qf_blue//2})")
print(f"LAW 2*selfK = SC:  {qf_black} == {SC}  ->  {qf_black == SC}")
print(f"total = SC + 2*selfB: {qf_total} == {SC + qf_blue}  ->  {qf_total == SC + qf_blue}")
print(f"Aut-WEIGHTED qf count = {weighted}  (vs unweighted {qf_total})")
print(f"\ncarrier classes: {len(qf_by_class)}; "
      f"SC carriers: {sum(1 for r in carrier_rows if r[5]=='SC')}")
mult = defaultdict(int)
for c, k in qf_by_class.items(): mult[k] += 1
print(f"qf-per-class histogram: {dict(sorted(mult.items()))}")
aut_hist = defaultdict(int)
for r in carrier_rows: aut_hist[r[4]] += 1
print(f"carrier Aut histogram: {dict(sorted(aut_hist.items()))}")
print(f"\nfirst 30 carriers (class, qf, black-qf, H, Aut, SC?):")
for r in carrier_rows[:30]: print("  ", r)
print(f"\ntotal {time.time()-t0:.0f}s")
