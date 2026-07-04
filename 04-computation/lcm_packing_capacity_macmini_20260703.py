#!/usr/bin/env python3
"""
LCM-PACKING capacity (mac-mini-2026-07-03-S29): the ELEMENTARY bound on q*.
To block moduli {15..Q} by divisibility, each q must divide some speed, so lcm(15..Q) | prod(speeds) <= M^13.
=> Q <= Q_pack(M) where lcm(15..Q) fits in 13 factors each <= M. This is ELEMENTARY (lcm growth ~ e^Q).
Pack prime-POWERS p^{floor(log_p Q)} into 13 bins (each <= M), covering q=2..14, blocking {15..Q}. Measure the
actual first witness q* and confirm q* ~ Q_pack = O(log M) (the true adversarial capacity), and that the
'tiling correction' (free moduli that are still no-witness) is bounded.
"""
from math import gcd, log
from functools import reduce
import random
def gcd_all(xs): return reduce(gcd,xs)
def nd(x): x=x%1; return min(x,1-x)
def is_prime(n):
    if n<2: return False
    for d in range(2,int(n**0.5)+1):
        if n%d==0: return False
    return True
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def has_witness(sp,q):
    for a in range(1,q):
        if gcd(a,q)!=1: continue
        if all(nd(v*a/q)>=1/14 for v in sp): return True
    return False
def first_witness(sp,qmax=800):
    q=2
    while q<=qmax and not has_witness(sp,q): q+=1
    return q if q<=qmax else None

def prime_powers_upto(Q):
    """the prime powers p^k that appear in lcm(15..Q): for each prime p<=Q, the largest p^k<=Q."""
    pps=[]
    for p in range(2,Q+1):
        if is_prime(p):
            pk=p
            while pk*p<=Q: pk*=p
            pps.append(pk)
    return sorted(pps, reverse=True)

def pack(Q, M, nbins=13):
    """bin-pack the prime-powers of lcm(2..Q) into nbins bins each <= M (first-fit-decreasing). Return bins or None."""
    pps=prime_powers_upto(Q)
    bins=[1]*nbins
    for pp in pps:
        if pp>M: return None
        # place in the bin with smallest product that stays <= M
        cand=sorted(range(nbins), key=lambda i: bins[i])
        placed=False
        for i in cand:
            if bins[i]*pp<=M: bins[i]*=pp; placed=True; break
        if not placed: return None
    return bins

if __name__=="__main__":
    rng=random.Random(29)
    print("LCM-packing capacity: max Q blockable (lcm(2..Q) in 13 bins<=M) vs actual first witness q*.")
    print("="*82)
    print(f"{'M':>16} {'log10 M':>8} {'Q_pack (max blockable)':>22} {'q* (first witness)':>18} {'q*-Q_pack':>10}")
    for M in [10**4, 10**6, 10**9, 10**12, 10**15, 10**18]:
        # find max Q with a valid packing
        Qpack=14
        for Q in range(15, 800):
            if pack(Q, M) is not None: Qpack=Q
            else: break
        bins=pack(Qpack, M)
        # the packed speeds (13 bins) block {2..Qpack}; make 13 distinct nonzero, covering
        sp=sorted(set(b if b>1 else rng.randint(2,22) for b in bins))
        while len(sp)<13: sp.append(rng.randint(2,22))
        sp=sorted(set(sp))[:13]
        cov=is_covering(sp)
        qstar=first_witness(sp,800) if len(sp)==13 else None
        gap = (qstar-Qpack) if qstar else None
        print(f"{M:>16} {int(log(M,10)):>8} {Qpack:>22} {str(qstar):>18} {str(gap):>10}  cov={cov} len={len(sp)}")
    print("\n=> Q_pack = O(log M) (lcm(2..Q)~e^Q <= M^13 => Q<=13 log M) is the ELEMENTARY capacity.")
    print("   if q* ~ Q_pack (small bounded gap): the crux is ELEMENTARY lcm-packing + bounded tiling correction.")
    print("   q* CANNOT exceed Q_pack by much: beyond Q_pack some modulus divides no speed (lcm exhausted) =>")
    print("   witness there unless a bounded 2-term tiling, which is a p<=~95 correction (circle method, bounded).")
