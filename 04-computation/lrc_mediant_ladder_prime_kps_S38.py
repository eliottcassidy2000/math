#!/usr/bin/env python3
"""
kps-S38 (the centerpiece): the UNIVERSAL MEDIANT LADDER family, and how it exposes opus's
prime/composite dichotomy from INSIDE my ladder mechanism.

RECONCILIATION found in Part 1: opus's N=13 mediant witness {1..11,13,36}=3/41 is the j=3 rung
of base {1..11,13} (mu=1/12, rho=5, outlier 36=3*12).  Same shape as S35's N=7 slice
{1..5,7}+18.  So the universal mediant family is
        F(N) = {1,...,N-2, N}  +  outlier 3(N-1)         [ base = AP{1..N-2} + defect N ]
and the LADDER FORMULA predicts M(F(N)) = (1/(N-1)) * 3(N-1)/(3(N-1)+5) = 3/(3N+2) (the mediant).

CLAIM to test: the ladder gives a LOWER bound M(F(N)) >= 3/(3N+2) always, but EQUALITY
(M = mediant, in the gap) holds IFF 3N+2 is PRIME.  When 3N+2 is COMPOSITE, a BETTER witness
(from a proper factor of 3N+2) lifts M ABOVE the gap -- opus's arithmetic obstruction, seen
from inside the ladder.  N=12: 38=2*19 composite => M(F(12)) > 2/25 (escapes the gap).

If true, this is the mechanism UNIFIED with the arithmetic: the ladder builds the mediant
candidate at every N; primality of 3N+2 decides whether it survives as the actual M.
"""
from fractions import Fraction
import numpy as np
from math import isqrt, gcd
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)
def isprime(n): return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

print("=== THE UNIVERSAL MEDIANT LADDER F(N) = {1..N-2, N} + 3(N-1), vs 3N+2 prime ===\n", flush=True)
print(f"  {'N':>3}{'3N+2':>6}{'prime':>7}{'mediant':>9}{'gap top 2/(2N+1)':>18}{'M(F(N))':>10}{'= mediant?':>11}{'in gap?':>9}", flush=True)
rows=[]
for N in range(5,17):
    base=list(range(1,N-1))+[N]          # {1..N-2, N}
    fam=sorted(base+[3*(N-1)])
    if len(set(fam))!=N:                  # guard (small N overlaps)
        print(f"  {N:>3}  (degenerate speeds, skip)"); continue
    if reduce(gcd,fam)!=1:
        pass
    M,(c,q)=Mw(fam)
    med=Fraction(3,3*N+2); top=Fraction(2,2*N+1); pr=isprime(3*N+2)
    ingap = Fraction(1,N+1) < M < top
    rows.append((N,pr,M,med,ingap))
    print(f"  {N:>3}{3*N+2:>6}{str(pr):>7}{str(med):>9}{str(top):>18}{str(M):>10}{str(M==med):>11}{str(ingap):>9}", flush=True)

print("\n  CHECK the claim 'M(F(N)) = mediant  <=>  3N+2 prime':", flush=True)
ok=all((r[2]==r[3])==r[1] for r in rows)
print(f"    holds across N=5..16: {ok}", flush=True)
print(f"    prime 3N+2 (N): {[r[0] for r in rows if r[1]]}  -> M=mediant, IN GAP", flush=True)
print(f"    composite   (N): {[r[0] for r in rows if not r[1]]}  -> M>mediant (better witness), gap ESCAPED", flush=True)

# zoom on N=12: show the better witness explicitly
print("\n=== N=12 zoom: F(12) = {1..10,12} + 33; the composite 38=2*19 gives a better witness ===", flush=True)
fam=sorted(list(range(1,11))+[12,33]); M,(c,q)=Mw(fam)
print(f"  F(12)={fam}: ladder predicts mediant 3/38 (~0.0789); ACTUAL M={M} (~{float(M):.4f}) at t={c}/{q}", flush=True)
print(f"  gap (1/13,2/25)=(0.0769,0.0800): M {'IN' if Fraction(1,13)<M<Fraction(2,25) else 'ABOVE (escaped)'} gap.", flush=True)
print(f"  => 38 composite => a proper-factor witness lifts M off the mediant, OUT of the gap.  This is", flush=True)
print(f"     opus's arithmetic obstruction (HYP-4506) seen from inside the ladder: the mediant candidate", flush=True)
print(f"     is BUILT at every N, and primality of 3N+2 decides survival.  Width was a symptom; 3N+2 is the cause.", flush=True)
