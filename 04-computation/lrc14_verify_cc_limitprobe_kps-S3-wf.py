#!/usr/bin/env python3
"""
LIMIT PROBE: does best_margin approach 1 (the claimed limit infimum) for real covering S3
sets as V0 grows along a FIXED bad shape?  If margin -> 1 from above and covering primitive
sets exist arbitrarily far out on the shape, then C(S) has NO uniform margin and the
"realized floor 1.33" is NOT a theorem.  We hunt the worst shape and push V0 up.

Also: directly test the claim that for a FIXED shape the finite margin converges to the
limit FROM ABOVE (monotone), which the package flagged as VERIFIED-not-proved.
"""
from fractions import Fraction as F
from math import gcd
import random

H = F(1,14)
def nrm(x):
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1,2) else 1 - r
def safe_components(A, h=H):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe
def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)
def best_margin(S):
    S = sorted(set(S)); best = F(-1); arg = None
    for v in S:
        m = Wwidth([u for u in S if u != v]) * 7 * v
        if m > best: best = m; arg = v
    return best, arg
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2,15))
def is_primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1
def case_of(S):
    S=sorted(set(S)); Vmin=S[0]; Vmax=S[-1]; k=sum(1 for v in S if v>13)
    if k<=1: return "S1"
    if Vmax < 13*Vmin: return "S2"
    return "S3"

print("="*70)
print("PART 1: fixed-shape margin vs V0 (look for descent toward 1)")
print("="*70)
# A shape: small part covers 2..13, cluster = {V0+off}, offsets fixed, includes mult-of-14 anchor.
# We slide V0 over multiples of 14 (so offset 0 = mult of 14) so the cluster always covers q=14.
small = [7,8,9,10,11,12,13]      # covers 2..13
offset_sets = [
    [0,14,28,42,45,46],          # 6 cluster speeds (spread 46)
    [0,2,4,6,8,10],
    [0,1,2,30,31,45],
    [0,7,14,21,28,35],
    [0,15,16,30,31,45],
]
for offs in offset_sets:
    print(f"\noffsets={offs} small={small}")
    rows=[]
    for k in [1,2,5,10,50,200,1000,5000,25000,260260//14*1]:  # V0 = 14*k (mult of 14)
        V0 = 14*k
        S = sorted(set(small + [V0+o for o in offs]))
        if len(S)!=13:
            rows.append((V0,"sizebad")); continue
        prim = is_primitive(S); cov = is_covering(S); cs = case_of(S)
        bm,arg = best_margin(S)
        rows.append((V0, float(bm), arg, prim, cov, cs))
    prev=None
    for r in rows:
        if len(r)==2: print("  V0=",r[0], r[1]); continue
        V0,bm,arg,prim,cov,cs = r
        mono = "" if prev is None else ("DOWN" if bm<prev else ("UP" if bm>prim else "flat"))
        mono = "" if prev is None else ("v" if bm<prev else ("^" if bm>prev else "="))
        print(f"  V0={V0:8d} margin={bm:.6f} {mono} arg={arg} prim={prim} cov={cov} {cs}")
        prev=bm

print()
print("="*70)
print("PART 2: ADVERSARIAL minimize best_margin over covering primitive S3, all scales")
print("="*70)
# Broad randomized hunt focused on getting margin as LOW as possible. Accept only covering+primitive S3.
rng = random.Random(20260618)
small_pool = [
    [7,8,9,10,11,12,13],
    [5,7,8,9,10,11,12,13],
    [6,7,8,9,10,11,12,13],
    [4,7,8,9,10,11,12,13],
    [3,6,7,8,9,10,11,12,13],
    [9,10,11,12,13],
    [8,9,10,11,12,13],
    [11,12,13],
]
worst = F(10**9); worst_set=None; tested=0; cfails=0
for _ in range(400000):
    small = rng.choice(small_pool)
    csize = 13 - len(small)
    if csize < 2: continue
    scale = rng.choice([1,1,1,2,3,5,10,30,100,1000,10000])
    V0 = 14*rng.randint(1, 26026)*1  # multiple of 14 anchor for q=14 coverage
    V0 = V0 // 14 * 14
    spread = rng.randint(csize, rng.choice([14,25,45,80,150]))
    cand14 = [x for x in range(0,spread+1) if (V0+x)%14==0]
    if not cand14: continue
    offs = set([rng.choice(cand14)])
    while len(offs) < csize:
        offs.add(rng.randint(0, spread))
    S = sorted(set(small + [V0+o for o in offs]))
    if len(S)!=13: continue
    if not is_primitive(S) or not is_covering(S): continue
    if case_of(S)!="S3": continue
    tested += 1
    bm,arg = best_margin(S)
    if bm <= 1: cfails += 1
    if bm < worst:
        worst = bm; worst_set = (S[:], arg)
print(f"tested covering primitive S3: {tested}")
print(f"C(S) failures (best_margin<=1): {cfails}")
print(f"global min best_margin: {worst} = {float(worst):.6f}")
print(f"  at {worst_set}")
print(f"  compared to claimed realized floor ~1.336 and 1/14 threshold (margin>1)")
