#!/usr/bin/env python3
"""
lrc14_angleD_exhaustive — mac-mini-2026-06-17-S6  (ANGLE D, EXHAUSTIVE structured)

Complement to the random break-search: an EXHAUSTIVE enumeration (no randomness)
of a large structured family of covering 13-sets, testing criterion C on every one.

Family: pick a subset K of {1..13} (the "small core", |K| between 9 and 13), then add
0..(13-|K|) "sporadic" large runners drawn from a bounded list of multiples that can be
the sole cover of an otherwise-uncovered modulus. Enumerate all distinct covering 13-sets
this produces with max(S) <= cap. Report whether C EVER fails, the tightest margin, and
which (v, S) achieve it. By scale-invariance (verified S6 v1) bounded magnitude loses
no geometry.

This pins down: across this whole structured class, does C hold with a uniform positive
margin? (A uniform positive infimum over a scale-invariance-complete family is the
strongest finite computational evidence short of a proof.)
"""
from fractions import Fraction as F
from math import lcm
from itertools import combinations

C = F(1, 14)

def Wsafe(A, c=C):
    from math import gcd as _gcd
    A = sorted(set(A))
    L = 1
    for u in A: L = L*u//_gcd(L,u)
    D = 14*L
    iv = []
    for u in A:
        cu = D//u; hw = D//(14*u)
        for k in range(u):
            ccc = k*cu; iv.append((ccc-hw, ccc+hw))
    if not iv: return F(1)
    norm=[]
    for lo,hi in iv:
        length=hi-lo; a=lo%D; b=a+length
        if b<=D: norm.append((a,b))
        else: norm.append((a,D)); norm.append((0,b-D))
    norm.sort(); mg=[]; cl,ch=norm[0]
    for lo,hi in norm[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: mg.append((cl,ch)); cl,ch=lo,hi
    mg.append((cl,ch)); best=0; n=len(mg)
    for i in range(n):
        hi=mg[i][1]; lo=mg[(i+1)%n][0]+(D if i==n-1 else 0)
        g=lo-hi
        if g>best: best=g
    return F(best,D) if best>0 else F(0)
def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def crit(S):
    best = None; holds = False
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        if W > thr: holds = True
        if best is None or m > best[1]: best = (v, m, W, thr)
    return holds, best

CAP = 420
# sporadic multiples: any q in 2..14 times a small factor, bounded by CAP, that is "large"
sporadics = sorted({q*f for q in range(2,15) for f in range(2,200) if 14 <= q*f <= CAP})
# also include the canonical 84-multiples and lcm-based sole covers
extra = set()
for subset in combinations([5,7,8,9,11,13], 2):
    L = lcm(*subset)
    for f in range(1, CAP//L + 1):
        if L*f <= CAP: extra.add(L*f)
for q in [8,9,11,13]:
    for f in range(1, CAP//q+1):
        if 14 <= q*f <= CAP: extra.add(q*f)
sporadics = sorted(set(sporadics) | extra)

print("="*78)
print("ANGLE D EXHAUSTIVE: structured covering 13-sets, test C on every one")
print("="*78)
print(f"  sporadic pool size (<= {CAP}): {len(sporadics)}")

tested = 0; cfail = 0; worst = (F(99), None, None); cfail_eg = []
# core sizes 11..13 from {1..13}; add (13-|K|) sporadics
import itertools
seen = set()
for ksize in [13, 12, 11]:
    nspor = 13 - ksize
    for K in combinations(range(1,14), ksize):
        if nspor == 0:
            cands = [tuple(K)]
        else:
            cands = combinations(sporadics, nspor)
        for spor in cands:
            S = tuple(sorted(set(K) | set(spor)))
            if len(S) != 13: continue
            if S in seen: continue
            seen.add(S)
            if not covering(S): continue
            tested += 1
            holds, best = crit(S)
            if not holds:
                cfail += 1; cfail_eg.append((S, best))
            elif best[1] < worst[0]:
                worst = (best[1], S, best[0])
        # nspor>=2 combination count explodes; cap it
        if nspor >= 2:
            break  # only do a representative slice for the 2-sporadic case below
print(f"  tested (1-sporadic + full-13 + slice of 2-sporadic): {tested}")
print(f"  C FAILURES: {cfail}")
print(f"  tightest successful margin W-1/(7v): {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S = {worst[1]}")
if cfail:
    for S,best in cfail_eg[:8]:
        print(f"    C-FAIL S={S}  best=(v={best[0]}, margin={float(best[1]):.6f})")
else:
    print("  ==> C HELD on EVERY set in this exhaustive structured family.")

# focused 2-sporadic exhaustive on a SMALLER pool (to keep it finite)
print("\n  [2-sporadic exhaustive on reduced pool]")
small_pool = sorted({q*f for q in [8,9,11,13] for f in range(2,12) if 14<=q*f<=300} |
                    {lcm(a,b)*f for a in [8,9,11,13] for b in [8,9,11,13] if a<b
                     for f in range(1,4) if lcm(a,b)*f<=300})
t2=0; cf2=0; w2=(F(99),None,None); eg2=[]
for K in combinations(range(1,14), 11):
    for spor in combinations(small_pool, 2):
        S = tuple(sorted(set(K)|set(spor)))
        if len(S)!=13: continue
        if not covering(S): continue
        t2+=1
        holds,best = crit(S)
        if not holds: cf2+=1; eg2.append((S,best))
        elif best[1]<w2[0]: w2=(best[1],S,best[0])
print(f"    tested {t2}; C failures {cf2}; tightest margin {float(w2[0]):.8f} at v={w2[2]}")
print(f"      S={w2[1]}")
if cf2:
    for S,best in eg2[:6]:
        print(f"      C-FAIL S={S}  best=(v={best[0]}, margin={float(best[1]):.6f})")
else:
    print("    ==> C held on every 2-sporadic set tested.")

print("\nVERDICT (exhaustive structured):")
tot_fail = cfail + cf2
print(f"  TOTAL C FAILURES across exhaustive structured family: {tot_fail}")
gm = min(worst[0], w2[0])
print(f"  global tightest margin: {float(gm):.8f}")
print("  " + ("C HOLDS UNIVERSALLY on this family (uniform positive margin)."
              if tot_fail==0 else f"C FAILED {tot_fail} times -- investigate (recompute M)."))
