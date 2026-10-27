"""The transitive-regular slab at odd orders, from the class move graphs dumped by dipole_exactness --dump.

For a tournament with sorted scores s_(0) <= ... <= s_(n-1) and m = (n-1)/2:
    L = (1/2) sum |s_(i) - i|  (THM-4810: distance to the transitive class >= L)
    R = (1/2) sum |s_(i) - m|  (THM-4815 B8: distance to the regular class  = R exactly)
Termwise triangle inequality: L + R >= (1/2) sum |i - m| = (n^2-1)/8, with equality iff every s_(i) lies between i and m.
The slab G_n is the equality set.  Checked here, class by class:
  (1) every class on a geodesic between TT and the regular class lies in G_n (d(TT,T) + d(T,REG) = (n^2-1)/8 => T in G_n);
  (2) inside G_n, d(T, TT) = L(T) exactly for the classes reachable from TT by R-descending moves, and these are exactly
      the slab classes with d(T,TT) = L(T);
  (3) which slab classes fail (d > L), and why (sources/sinks, cut vertices), level by level.
Usage: python3 slab_census.py n  (n odd; reads ../dipole_exactness_20261027/graph_n{n}.txt)"""
import sys
from collections import deque, Counter
n = int(sys.argv[1]); m = (n - 1) // 2; c = (n * n - 1) // 8
cls = []
for line in open(f"../dipole_exactness_20261027/graph_n{n}.txt"):
    head, anyn, blk, arcs = line.rstrip("\n").split("|")
    f = head.split()
    s = sorted(int(ch) for ch in f[2])
    L = sum(abs(s[i] - i) for i in range(n)) // 2
    R = sum(abs(x - m) for x in s) // 2
    cls.append(dict(scores=f[2], comps=f[3], aut=int(f[4]), sc=int(f[5]), lam=int(f[6]), L=L, R=R,
                    nb=[int(x) for x in anyn.split()]))
C = len(cls)
def bfs(pred):
    d = [-1] * C; q = deque()
    for i in range(C):
        if pred(cls[i]): d[i] = 0; q.append(i)
    while q:
        u = q.popleft()
        for t in cls[u]['nb']:
            if d[t] < 0: d[t] = d[u] + 1; q.append(t)
    return d
dTT = bfs(lambda x: x['L'] == 0)
dREG = bfs(lambda x: x['R'] == 0)
assert all(dREG[i] == cls[i]['R'] for i in range(C)), "B8 failed"
slab = [i for i in range(C) if cls[i]['L'] + cls[i]['R'] == c]
assert all(cls[i]['L'] + cls[i]['R'] >= c for i in range(C)), "slab inequality failed"
geo = [i for i in range(C) if dTT[i] + dREG[i] == c]
assert set(geo) <= set(slab), "a geodesic class outside the slab"
# R-descending reachability from TT inside the slab
reach = set(i for i in range(C) if cls[i]['L'] == 0); frontier = list(reach)
while frontier:
    nxt = []
    for u in frontier:
        for t in cls[u]['nb']:
            if cls[t]['R'] == cls[u]['R'] - 1 and t not in reach: reach.add(t); nxt.append(t)
    frontier = nxt
exact_slab = set(i for i in slab if dTT[i] == cls[i]['L'])
assert reach == exact_slab == set(geo), (len(reach), len(exact_slab), len(geo))
reg = [i for i in range(C) if cls[i]['R'] == 0]
print(f"n={n}: classes {C}; slab classes {len(slab)}; on a TT-REG geodesic (= reachable from TT by R-descent = L-exact in the slab) {len(geo)}; regular classes {len(reg)}, at distance {c} from TT: {sum(1 for i in reg if dTT[i] == c)}")
print("  level L: slab classes / L-exact / failing (failing: with source or sink, strong with lambda = 1, other)")
for Lv in range(c + 1):
    lev = [i for i in slab if cls[i]['L'] == Lv]
    ex = [i for i in lev if i in exact_slab]
    bad = [i for i in lev if i not in exact_slab]
    srcsink = sum(1 for i in bad if '0' in cls[i]['scores'] or str(n - 1) in cls[i]['scores'])
    lam1 = sum(1 for i in bad if not ('0' in cls[i]['scores'] or str(n - 1) in cls[i]['scores']) and cls[i]['lam'] == 1)
    print(f"   L={Lv:2d} R={c-Lv:2d}: {len(lev):5d} / {len(ex):5d} / {len(bad):5d}  ({srcsink}, {lam1}, {len(bad)-srcsink-lam1})")
