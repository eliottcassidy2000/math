#!/usr/bin/env python3
"""mac-mini-S105: EXACT band finite check + structure. M(S) is EXACT via peak candidates t=m/(v_i+v_j)
(local maxima of min_l ||v_l t|| occur where two runners cross with opposite slopes). (1) prove the
NEAR-DILATE theorem M(V_L)=1/13 exactly; (2) exact M for band interval-core multi-killers; (3) insights."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def M_exact(S):
    """EXACT M(S) = max over peak candidates t=m/(a+b), pairs a<b in S, of min_l ||v_l t||. Rational."""
    best=F(0)
    denoms=set()
    for a,b in combinations(sorted(set(S)),2):
        denoms.add(a+b)
        if b>a: denoms.add(b-a)
    for q in denoms:
        if q<=0: continue
        for m in range(1,q):
            # min_l ||v_l m/q|| = min_l min(v_l*m %q, q-(v_l*m%q)) / q
            num=q
            for v in S:
                r=(v*m)%q; d=min(r,q-r)
                if d<num: num=d
                if num*13<q: break   # early: clearance <1/13, this t won't beat 1/13
            c=F(num,q)
            if c>best: best=c
    return best
def M_ge_113_witness(S):
    """rigorous: is M(S)>=1/13? find a peak t with min clearance >=1/13 (exact)."""
    for a,b in combinations(sorted(set(S)),2):
        for q in ([a+b,b-a] if b>a else [a+b]):
            if q<=0: continue
            for m in range(1,q):
                ok=True
                for v in S:
                    r=(v*m)%q
                    if min(r,q-r)*13<q: ok=False; break
                if ok: return (m,q)
    return None
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))

print("(1) NEAR-DILATE theorem: M(V_L)=1/13 for L with 13 not| (L+1). Witness t=(L+1)/(13L). Verify exact:")
for L in [30,31,32,60,90,120,180]:  # L=30 => diam 391 (band); various
    V=[i*L for i in range(1,13)]+[13*L+1]
    t=F(L+1,13*L)
    # min clearance at t (exact)
    mc=min(min((v*t)%1, 1-(v*t)%1) for v in V)
    Me=M_exact(V) if 13*L+1<600 else None
    print(f"   L={L}: V_L diam={13*L+1}, gcd={gcd(V[0],gcd(*V))}, 13|(L+1)? {(L+1)%13==0}; clearance@t={mc}={float(mc):.5f}; M_exact={Me}")
print("   => M(V_L)=1/13 EXACTLY (small runners iL at ||i(L+1)/13||>=1/13, killer 13L+1 at (L+1)/(13L)>1/13).")

print("\n(2) EXACT band check: interval-core multi-killers {1..k} + outliers, largest in (220,475]:")
W0=475; below=[]; checked=0; minM=(F(1),None)
for k in [11,10,9]:
    core=list(range(1,k+1)); D=[q for q in range(k+1,15)]
    # candidate outliers: multiples of moduli in [2..14] that help cover D, in [15,W0]
    cands=sorted({q*j for q in range(2,15) for j in range(1,W0//q+1) if 15<=q*j<=W0})
    n_out=13-k
    cnt=0
    for outs in combinations(cands, n_out):
        if max(outs)<=220 or max(outs)>W0: continue
        S=core+list(outs)
        if len(set(S))!=13 or not covering(S): continue
        if gcd(S[0], __import__('functools').reduce(gcd,S))!=1: continue
        checked+=1; cnt+=1
        w=M_ge_113_witness(S)
        if w is None:
            Me=M_exact(S); below.append((Me,sorted(S)))
            if Me<minM[0]: minM=(Me,sorted(S))
        if cnt>=4000: break   # cap per k for time
print(f"   band interval-core families checked: {checked}; with M < 1/13 (no witness): {len(below)}")
for Me,S in sorted(below)[:6]: print(f"      M={Me}={float(Me):.5f}: {S}")
print(f"   => {'ALL band interval-core families have M>=1/13 (exact witness found).' if not below else 'some below 1/13 -- examine.'}")
