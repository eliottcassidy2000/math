#!/usr/bin/env python3
"""
kps-S44: is the case-2 covering system UNIFORMLY bounded (every non-AP pair-blocker clears at
q <= Q0, all heights) or does the clearing modulus grow with height (=> crux needs the height
bound)?  This is the open residual after S43 + the fleet's convergence (opus S124: (C) = one
residue-covering rigidity; the whole crux).

CLEAN SUB-RESULT: for q <= 12 the clearance band 2q/25 < 1, so clearing at q <=> avoid residue 0
<=> (for the relevant q) NO MULTIPLE OF q.  So:
    no multiple of q (q in {7,9,10,11,12}) => clears at q => M >= 1/q > 2/25.
This is a clean rational_point_margin certificate (t = 1/q, c=1) closing every family without a
small-prime multiple.  The ADVERSARIAL residual: 'super-blockers' divisible by many small q AND
blocking mod 25 -- do THOSE have a bounded clearing modulus?

Part 1: verify the clean sub-result.
Part 2: adversarial high-height super-blockers (CRT: block mod 25 + carry multiples of 7,9,10,11,
        12,13,...) -- find the minimal clearing modulus vs height.  Bounded => uniform covering.
"""
from fractions import Fraction
import numpy as np, random
from math import gcd
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd)
PAIRS=[{u,25-u} for u in [1,2,3,4,6,7,8,9,11,12]]
def is_blocker(v):
    if any(x%25==0 for x in v): return False
    R=set(x%25 for x in v); return all(p&R for p in PAIRS)
def clears_at(v,q):
    thr=Fraction(2,25)
    for c in range(1,q):
        if gcd(c,q)!=1: continue
        if all(Fraction(min((x*c)%q,q-(x*c)%q),q)>=thr for x in v): return True
    return False
def min_clear_mod(v,QMAX=200):
    for q in range(6,QMAX+1):
        if q==25: continue
        if clears_at(v,q): return q
    return None

print("=== Part 1: clean sub-result -- no multiple of q => clears at q (q<=12), M>=1/q>2/25 ===", flush=True)
for q in [7,9,10,11,12]:
    band=2*q/25
    print(f"  q={q}: clearance band 2q/25={band:.3f}<1 => avoid only residue 0; nonzero residue has margin>=1/{q}={1/q:.4f}>2/25", flush=True)
print("  => any family with no multiple of {7,9,10,11,12} clears (M>=1/12). Residual = has such multiples.\n", flush=True)

print("=== Part 2: ADVERSARIAL high-height super-blockers -- minimal clearing modulus vs height ===", flush=True)
# build 12 speeds: block mod 25 (residues = AP-like {1..12} mod 25) but each carries multiples of
# small primes via CRT to resist small-q clearing; scale up heights.
random.seed(7)
def crt(res_mod):  # res_mod: list of (mod, res); return smallest nonneg solution
    M=1; x=0
    for m,r in res_mod:
        g=gcd(M,m);
        if (r-x)%g!=0: return None
        lcm=M//g*m;
        # solve x + M*t ≡ r mod m
        t=((r-x)//g * pow(M//g, -1, m//g))%(m//g)
        x=x+M*t; M=lcm; x%=M
    return x
worst=0; examples=[]
smallprimes=[7,11,13,17,19,23,29,31]
for trial in range(3000):
    # target residues mod 25: a blocker set (one per pair, +2 non-units)
    base=list(range(1,13))
    random.shuffle(base)
    speeds=[]; ok=True
    used_primes=random.sample(smallprimes, random.randint(2,5))
    for i in range(12):
        cons=[(25, base[i]%25)]
        # make some speeds multiples of a chosen small prime (resist that q)
        if i < len(used_primes):
            cons.append((used_primes[i],0))
        # random large lift
        x=crt(cons)
        if x is None: ok=False; break
        L=random.choice([1,2,3,5,8])
        x=x+ L*(25*used_primes[0] if used_primes else 25)*random.randint(1,40)
        speeds.append(x if x>0 else x+25*7*11)
    if not ok: continue
    v=sorted(set(speeds))
    if len(v)!=12: continue
    if reduce(gcd,v)!=1:
        g=reduce(gcd,v); v=[x//g for x in v]
        if len(set(v))!=12: continue
    if not is_blocker(v): continue
    qm=min_clear_mod(v)
    if qm is None:
        examples.append(('NO CLEAR <=200', v, Mw(v))); continue
    if qm>worst:
        worst=qm; examples.append((qm, v, Mw(v), max(v)))
print(f"  worst (largest) minimal clearing modulus found among adversarial super-blockers: {worst}", flush=True)
for e in examples[-8:]:
    if e[0]=='NO CLEAR <=200': print(f"    !! no clear q<=200: {e[1]} M={e[2]}", flush=True)
    else: print(f"    q_min={e[0]:>3}  height={e[3]:>7}  M={e[2]}  v={e[1]}", flush=True)
print("\n  READING: if q_min stays small (<=~40) even for large-height super-blockers, the covering is", flush=True)
print("  UNIFORMLY BOUNDED (crux closes by finite covering). If q_min grows with height, the covering", flush=True)
print("  is NOT bounded and the rigidity needs the height/order bound (mac-mini/opus lane).", flush=True)
