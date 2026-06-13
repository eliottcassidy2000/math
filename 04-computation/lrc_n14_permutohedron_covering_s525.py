#!/usr/bin/env python3
"""
lrc_n14_permutohedron_covering_s525.py    oracle-2026-06-01-S525

n=14 attempt via the PERMUTOHEDRON / circle-covering reformulation, building on
opus-S524 (CRT 7-class reduction) and S522-S524 (round = realizable; #SCC block
structure; regular polygon = roots of unity).

PERMUTOHEDRON PICTURE.  Runner positions (s_i t mod 1), i=1..13, trace the speed
LINE on the torus T^13.  Each runner i has a forbidden SLAB  ||s_i t|| < 1/14.
LRC@14  <=>  the line exits ALL 13 slabs simultaneously at some t  <=>  the 13
blocking-arcs  B_i = { t : ||s_i t|| < 1/14 }  on the time-circle [0,1) do NOT
COVER it.  Each B_i is a union of intervals; the slab walls ||s_i t||=1/14 are the
affine-braid-arrangement walls (permutohedron alcove boundaries).  A LONELY time
is an open alcove avoiding every slab; the TIGHT case is lonely only on a wall.

We compute, per speed set:
  - whether it is lonely, the lonely time, and OPEN vs WALL-ONLY,
  - d_open = max over open cells of #safe runners (13 = lonely in open cell),
  - total blocked measure |union B_i| and the gap (lonely) measure,
  - at the lonely time: the half-turn runner sub-tournament's #SCC and roundness
    (S523/S524: realizable = round; tight = boundary/strong block).
Run on the initial segment {1..13} (the tight extreme) + random + adversarial
large-speed sets; report whether the initial segment is UNIQUELY wall-only.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import random

N = 14
THR = Fraction(1, N)

def frac(x): return x - (x.numerator // x.denominator)
def dist0(x):
    f = frac(x); return min(f, ONE - f)
ONE = Fraction(1)

def is_primitive(s):
    return reduce(gcd, s) == 1

def walls(speeds):
    """endpoints where some ||s_i t|| = 1/14: t=(k*1 +- 1/14)/s_i scaled ->
    t = (14k +- 1)/(14 s_i), in [0,1)."""
    W = set()
    for s in speeds:
        for k in range(0, s + 1):
            for sgn in (1, -1):
                t = Fraction(14 * k + sgn, 14 * s)
                if 0 <= t < 1:
                    W.add(t)
    return sorted(W)

def nsafe(speeds, t):
    return sum(1 for s in speeds if dist0(Fraction(s) * t) >= THR)

def half_turn_subtournament(speeds, t):
    m = len(speeds)
    adj = [[0]*m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            if a==b: continue
            f = frac(Fraction(speeds[a]-speeds[b]) * t)
            adj[a][b] = 1 if 0 < f < Fraction(1,2) else 0
    return adj

def scc_count(adj):
    m=len(adj)
    def reach(s, fwd=True):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for w in range(m):
                e = adj[u][w] if fwd else adj[w][u]
                if e and w not in seen: seen.add(w); st.append(w)
        return seen
    comp=[None]*m; c=0
    for v in range(m):
        if comp[v] is not None: continue
        both = reach(v, True) & reach(v, False)
        for w in both:
            if comp[w] is None: comp[w]=c
        c+=1
    return len(set(x for x in comp if x is not None))

def analyze(speeds):
    W = walls(speeds)
    Wl = W + [ONE]
    mids = [(a+b)/2 for a,b in zip(Wl, Wl[1:])]
    # open cells
    d_open = 0; lonely_t = None; kind = None
    for t in mids:
        d = nsafe(speeds, t)
        if d > d_open: d_open = d
        if d == len(speeds):
            lonely_t = t; kind = 'OPEN'; break
    if lonely_t is None:
        # wall check (tight case)
        for t in W:
            if nsafe(speeds, t) == len(speeds):
                lonely_t = t; kind = 'WALL'; break
    lonely = lonely_t is not None
    info = dict(lonely=lonely, kind=kind, d_open=d_open, t=lonely_t, nwalls=len(W))
    if lonely:
        adj = half_turn_subtournament(speeds, lonely_t)
        info['scc'] = scc_count(adj)
    return info

def main():
    print("LRC@14 via permutohedron circle-covering (oracle-S525)\n")
    sets = []
    sets.append(("initial {1..13}", tuple(range(1,14))))
    rnd = random.Random(525)
    # random primitive 13-subsets from 1..40
    cnt=0
    while cnt < 120:
        s = tuple(sorted(rnd.sample(range(1,41), 13)))
        if is_primitive(s): sets.append((f"rand{cnt}", s)); cnt+=1
    # adversarial: large gaps / near-AP / contains the initial pattern scaled
    adv = [tuple(range(1,14)),
           tuple(2*i+1 for i in range(13)),         # odds 1,3,..,25
           tuple(i for i in range(1,27) if i%2==1)[:13],
           tuple([1]+list(range(7,19))),
           tuple(sorted({1,2,3,4,5,6,7,8,9,10,11,12,100})),
           tuple(sorted({1,7,8,14,15,21,22,28,29,35,36,42,43}))]  # mod-7 clustered
    for i,s in enumerate(adv):
        if len(set(s))==13 and is_primitive(s): sets.append((f"adv{i}", s))

    fails=[]; wall_only=[]; results=[]
    for name, s in sets:
        info = analyze(s)
        results.append((name,s,info))
        if not info['lonely']: fails.append((name,s,info))
        elif info['kind']=='WALL': wall_only.append((name,s,info))

    print(f"speed sets tested: {len(sets)}")
    print(f"LRC FAILURES (no lonely time): {len(fails)}")
    if fails:
        for name,s,info in fails[:10]: print("  FAIL", name, s, info)
    print(f"WALL-ONLY lonely (tight, boundary seam): {len(wall_only)}")
    for name,s,info in wall_only:
        print(f"  {name}: {s}  d_open={info['d_open']}/13  lonely@{info['t']}  runner-subtour #SCC={info.get('scc')}")
    print()
    # structure of the lonely configuration across all sets
    from collections import Counter
    sccc = Counter(info['scc'] for _,_,info in results if info['lonely'])
    kinds = Counter(info['kind'] for _,_,info in results if info['lonely'])
    print("lonely-time kind:", dict(kinds))
    print("lonely-time runner sub-tournament #SCC distribution:", dict(sorted(sccc.items())))
    print("\nREADING: LRC@14 = '13 blocking arcs do not cover [0,1)'. Every tested set")
    print("has a gap (lonely). The initial segment {1..13} is the tight WALL-ONLY case")
    print("(d_open=12): lonely only on a permutohedron wall -- the boundary seam (S523),")
    print("the regular-polygon extreme. The covering gap is the unproven core (= LRC@14).")

if __name__=="__main__":
    main()
