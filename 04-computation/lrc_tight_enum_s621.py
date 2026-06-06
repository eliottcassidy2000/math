#!/usr/bin/env python3
"""
S621 — enumerate the EXACTLY-TIGHT Lonely Runner instances (the "collapse family"
of HYP-2190: primitive speed sets with free measure p_0 = 0, equivalently
loneliness gap gamma(S) = exactly 1/(n+1), equivalently the width-2/(n+1)
forbidden arcs COVER the circle).

Fast path: a non-tight set has lonely measure ~ (1-2/(n+1))^n ~ 0.13 (an open set),
so a cheap float sample at K times rules it out w.h.p.; survivors are verified
EXACTLY over Q with gap_exact.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, time

def norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= Fr(1,2) else 1-f

def not_tight_fast(speeds, K, df):
    """Return True if a sample time is lonelier than 1/(n+1) (=> p_0>0 => not tight)."""
    margin = df + 1e-9
    for k in range(1, K+1):
        t = (k + 0.5)/(2*K+1)           # irrational-ish offset, avoids exact grid resonance
        m = 1.0
        for v in speeds:
            f = (v*t) % 1.0
            d = f if f <= 0.5 else 1-f
            if d < m:
                m = d
                if m <= margin: break
        if m > margin:
            return True
    return False

def gap_exact(speeds):
    V=[abs(v) for v in speeds]; cands=set()
    for i in range(len(V)):
        vi=V[i]
        for k in range(0,2*vi+1):
            t=Fr(2*k+1,2*vi)
            if 0<t<=Fr(1,2): cands.add(t)
        for j in range(i):
            vj=V[j]
            for d in (vi+vj,abs(vi-vj)):
                if d==0: continue
                kk=1
                while Fr(kk,d)<=Fr(1,2):
                    cands.add(Fr(kk,d)); kk+=1
    best=Fr(0)
    for t in cands:
        m=min(norm(v*t) for v in V)
        if m>best: best=m
    return best

def enum_tight(n, R, K=120):
    delta=Fr(1,n+1); df=1.0/(n+1); out=[]
    for s in itertools.combinations(range(1,R+1), n):
        if reduce(gcd,s)!=1: continue
        if not_tight_fast(s, K, df): continue
        if gap_exact(list(s))==delta: out.append(s)
    return out

if __name__=="__main__":
    print("EXACTLY-TIGHT primitive LRC instances  (gamma = 1/(n+1) = forbidden arcs cover circle)")
    print("="*78)
    ranges={3:14,4:16,5:18,6:21,7:24}
    seq=[]
    for n in range(3,8):
        R=ranges[n]; t0=time.time(); found=enum_tight(n,R)
        seq.append(len(found))
        print(f"\n n={n}  (range 1..{R}):  {len(found)} tight sets   [{time.time()-t0:.1f}s]")
        for s in found:
            resid=tuple(v%(n+1) for v in s)
            has1 = any(v%(n+1) in (1,n) for v in s)         # some speed = +-1 mod (n+1)
            zerofree = all(v%(n+1)!=0 for v in s)
            print(f"    {str(s):26s} resid mod {n+1}: {resid}  +-1?{has1} 0-free?{zerofree}")
        # finiteness probe: does a bigger range add any?
    print("\nCOUNT SEQUENCE n=3..7 (bounded range):", seq)
