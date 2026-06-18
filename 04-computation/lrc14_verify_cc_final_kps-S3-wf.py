#!/usr/bin/env python3
"""
CONSOLIDATED adversarial verification (unbuffered) of cluster-collapse-threedistance, LRC14 S3.
Exact Fractions throughout. Prints flushed.

Sections:
 0. Tool self-consistency vs prompt-provided routines.
 1. CLUSTER-COLLAPSE LEMMA: exact re-derivation, formula, threshold.
 2. ARC-WIDTH IMPLICATION soundness spot-checks (C(S) => M(S)>=1/14).
 3. BROAD adversarial hunt (moderate V0 <= 4000) for C(S) failures / M<1/14. Hundreds of thousands.
 4. LARGE-V0 targeted batch (few sets, V0 up to ~360360) for the carry-phase limit regime.
 5. v=max determinism sub-claim re-check.
"""
import sys
from fractions import Fraction as F
from math import gcd
import random

def out(*a): print(*a); sys.stdout.flush()

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
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C
def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2,15))
def is_primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1
def best_margin(S):
    S = sorted(set(S)); best = F(-1); arg = None
    for v in S:
        m = Wwidth([u for u in S if u != v]) * 7 * v
        if m > best: best = m; arg = v
    return best, arg
def case_of(S):
    S=sorted(set(S)); Vmin=S[0]; Vmax=S[-1]; k=sum(1 for v in S if v>13)
    if k<=1: return "S1"
    if Vmax < 13*Vmin: return "S2"
    return "S3"

# ---------- SECTION 0 ----------
out("="*70); out("SECTION 0: tool self-consistency"); out("="*70)
S0=[1,2,3,4,5,6,7,8,9,10,11,12,13]
out("S={1..13} M=", Mval(S0), float(Mval(S0)), ">=1/14?", Mval(S0)>=H, "case", case_of(S0))
bm,arg=best_margin(S0); out("  best_margin=",bm,float(bm),"arg=",arg)

# ---------- SECTION 1 ----------
out(); out("="*70); out("SECTION 1: cluster-collapse lemma (exact)"); out("="*70)
def cc(V0,Vmax):
    lo=F(1,14*V0); hi=F(13,14*Vmax)
    return (hi-lo), (hi-lo)*7*Vmax, (F(13)-F(Vmax,V0))/2, (F(1,Vmax)-F(1,14*Vmax)==hi)
ok_all=True; thr_all=True
for V0 in [14,15,29,100,1000,360360]:
    for Vmax in [V0+1,V0+14,2*V0,10*V0,11*V0-1,11*V0,12*V0]:
        if 13*V0<=Vmax: continue
        w,m,f,hiok=cc(V0,Vmax)
        ok_all = ok_all and (m==f) and hiok
        thr_all = thr_all and ((m>1)==(Vmax<11*V0))
out("formula width*7*Vmax == (13 - Vmax/V0)/2 AND hi==1/Vmax-1/(14Vmax):", ok_all)
out("threshold margin>1 <=> Vmax<11*V0:", thr_all)
# exhaustive safety check: random L in [V0,Vmax], band inside safe_components(L)?
random.seed(1); ccfail=0
for _ in range(4000):
    V0=random.randint(14,3000); spread=random.randint(1,min(12*V0-1,50))
    Vmax=V0+spread; c=random.randint(2,8)
    L=set([V0,Vmax])
    while len(L)<c: L.add(V0+random.randint(0,spread))
    L=sorted(L)
    lo=F(1,14*V0); hi=F(13,14*Vmax)
    sc=safe_components(L)
    contained=any(a<=lo and hi<=b for a,b in sc)
    mid=(lo+hi)/2; allsafe=all(nrm(u*mid)>=H for u in L)
    if not (contained and allsafe): ccfail+=1
out("cluster-collapse safety failures in 4000 random L:", ccfail)

# ---------- SECTION 2 ----------
out(); out("="*70); out("SECTION 2: arc-width implication spot-check (C => M>=1/14)"); out("="*70)
# For sets where C holds via v, confirm M(S)>=1/14 directly. Pull a few covering sets.
def quick_cov(rng):
    small=[7,8,9,10,11,12,13]
    V0=14*rng.randint(2,80); spread=rng.randint(6,45)
    c14=[x for x in range(0,spread+1) if (V0+x)%14==0]
    if not c14: return None
    offs=set([rng.choice(c14)])
    while len(offs)<6: offs.add(rng.randint(0,spread))
    S=sorted(set(small+[V0+o for o in offs]))
    if len(S)!=13 or not is_primitive(S) or not is_covering(S): return None
    return S
rng=random.Random(5); chk=0; bad=0
while chk<300:
    S=quick_cov(rng)
    if S is None or case_of(S)!="S3": continue
    chk+=1
    bm,arg=best_margin(S); M=Mval(S)
    if bm>1 and M<H: bad+=1; out("  IMPLICATION VIOLATION:",S,"bm",float(bm),"M",M)
out(f"checked {chk} S3 sets with C holding; arc-width implication violations: {bad}")

# ---------- SECTION 3 ----------
out(); out("="*70); out("SECTION 3: BROAD hunt (V0<=4000) for C-failures / M<1/14"); out("="*70)
small_pool=[
 [7,8,9,10,11,12,13],[5,7,8,9,10,11,12,13],[6,7,8,9,10,11,12,13],
 [4,7,8,9,10,11,12,13],[3,6,7,8,9,10,11,12,13],[9,10,11,12,13],
 [8,9,10,11,12,13],[11,12,13],[1,7,8,9,10,11,12,13],[2,7,8,9,10,11,12,13],
]
rng=random.Random(20260618)
worst=F(10**9); worst_set=None; tested=0; cfail=0; mfail=0; cfx=[]; mfx=[]
for _ in range(250000):
    small=rng.choice(small_pool); csize=13-len(small)
    if csize<2: continue
    V0=14*rng.randint(1,280)  # V0 up to ~3920, anchored to mult of 14
    spread=rng.randint(csize, rng.choice([14,25,45,90]))
    c14=[x for x in range(0,spread+1) if (V0+x)%14==0]
    if not c14: continue
    offs=set([rng.choice(c14)])
    while len(offs)<csize: offs.add(rng.randint(0,spread))
    S=sorted(set(small+[V0+o for o in sorted(offs)]))
    if len(S)!=13 or not is_primitive(S) or not is_covering(S): continue
    if case_of(S)!="S3": continue
    tested+=1
    bm,arg=best_margin(S)
    if bm<worst: worst=bm; worst_set=(S[:],arg)
    if bm<=1:
        cfail+=1; M=Mval(S)
        if len(cfx)<10: cfx.append((S[:],bm,arg,M))
        if M<H: mfail+=1;
        if M<H and len(mfx)<10: mfx.append((S[:],M))
out(f"tested covering primitive S3: {tested}")
out(f"C(S) failures (best_margin<=1): {cfail}")
out(f"TRUE M<1/14 counterexamples: {mfail}")
out(f"min best_margin: {worst} = {float(worst):.6f} at {worst_set}")
for s,b,a,M in cfx[:10]: out("  C-FAIL",s,float(b),"arg",a,"M",float(M))
for s,M in mfx[:10]: out("  M<1/14",s,M,float(M))

# ---------- SECTION 4 ----------
out(); out("="*70); out("SECTION 4: LARGE-V0 targeted batch (carry-phase limit regime)"); out("="*70)
# A handful of hand-picked bad shapes pushed to large V0. Exact, slow, bounded count.
small=[7,8,9,10,11,12,13]
shapes=[[0,14,28,42,45,46],[0,1,2,30,31,45],[0,7,14,28,42,45],[0,2,16,30,44,45],[0,15,30,31,44,45]]
for offs in shapes:
    out(f" shape offsets={offs}")
    for V0 in [140, 1400, 14000, 140000, 364000]:
        V0=(V0//14)*14
        S=sorted(set(small+[V0+o for o in offs]))
        if len(S)!=13: out(f"   V0={V0} sizebad"); continue
        prim=is_primitive(S); cov=is_covering(S); cs=case_of(S)
        bm,arg=best_margin(S)
        out(f"   V0={V0:7d} margin={float(bm):.6f} arg={arg} prim={prim} cov={cov} {cs}")

# ---------- SECTION 5 ----------
out(); out("="*70); out("SECTION 5: v=max determinism sub-claim"); out("="*70)
out("claimed v=max counterexample margin 242015/245024 =", float(F(242015,245024)), "(<1)")
out("claimed M 123/989 =", float(F(123,989)), ">=1/14?", F(123,989)>=H)
rng=random.Random(7); vmax_fail=0; ex=[]
for _ in range(120000):
    small=rng.choice(small_pool); csize=13-len(small)
    if csize<2: continue
    V0=14*rng.randint(1,200); spread=rng.randint(csize, rng.choice([14,25,45]))
    c14=[x for x in range(0,spread+1) if (V0+x)%14==0]
    if not c14: continue
    offs=set([rng.choice(c14)])
    while len(offs)<csize: offs.add(rng.randint(0,spread))
    S=sorted(set(small+[V0+o for o in sorted(offs)]))
    if len(S)!=13 or not is_primitive(S) or not is_covering(S): continue
    if case_of(S)!="S3": continue
    vmax=max(S)
    mvmax=Wwidth([u for u in S if u!=vmax])*7*vmax
    if mvmax<=1:
        bm,arg=best_margin(S); vmax_fail+=1
        if len(ex)<6: ex.append((S[:],mvmax,bm,arg,Mval(S)))
out(f"S3 sets where v=max margin<=1 (determinism FALSE): {vmax_fail}")
for s,mv,bm,arg,M in ex:
    out(f"  S={s} vmax_margin={float(mv):.5f} best_margin={float(bm):.5f} arg={arg} M={M}={float(M):.6f}")
out(); out("DONE")
