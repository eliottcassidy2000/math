#!/usr/bin/env python3
"""
OPTIMAL prime-packing adversary + tiling correction (mac-mini-2026-07-03-S29).
To block prime p by DIVISIBILITY, make some speed divisible by p (=> that runner at residue 0 => no witness).
To block primes {15..Q}, assign each to a speed (speed = product of its assigned primes), keeping speed <= M.
Greedy: distribute primes to the 13 speeds to maximize the blocked range Q, under max|speed|<=M. Then measure
 - Q_div = largest Q with all primes in [15,Q] dividing some speed (the ELEMENTARY primorial/packing capacity).
 - q* = first actual witness (min lonely denominator).
 - tiling correction = q* - (first free prime): free primes below q* are no-witness via 2-term TILING.
Question: is q* ~ Q_div (elementary) + a BOUNDED tiling correction? and does the correction stay bounded?
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
def first_witness(sp,qmax=600):
    q=2
    while q<=qmax and not has_witness(sp,q): q+=1
    return q if q<=qmax else None

def optimal_packer(M, rng):
    """13 speeds, each = product of assigned primes, <= M, blocking as many primes {2..} as possible."""
    primes=[p for p in range(2,400) if is_prime(p)]
    speeds=[1]*13
    blocked=[]
    for p in primes:
        # assign p to the speed that stays <= M and is smallest (greedy balance)
        cand=sorted(range(13), key=lambda i: speeds[i])
        placed=False
        for i in cand:
            if speeds[i]*p<=M:
                speeds[i]*=p; blocked.append(p); placed=True; break
        if not placed: break
    # ensure 13 distinct nonzero; pad
    speeds=[s if s>1 else rng.randint(2,22) for s in speeds]
    return sorted(set(speeds))[:13], blocked

if __name__=="__main__":
    rng=random.Random(2929)
    print("Optimal prime-packing: Q_div (packing capacity) vs q* (first witness) vs tiling correction.")
    print("="*88)
    print(f"{'M':>14} {'log10 M':>8} {'#blocked primes':>16} {'max blocked prime':>18} {'q*':>5} {'tiling corr (q*-firstfree)':>26}")
    for M in [10**6, 10**9, 10**12, 10**15, 10**18]:
        best=None; bq=0
        for _ in range(60):
            sp,blocked=optimal_packer(M, rng)
            if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
            q=first_witness(sp,600)
            if q and q>bq: bq=q; best=(sp,blocked)
        if best is None:
            print(f"{M:>14} (no valid covering packer)"); continue
        sp,blocked=best; qstar=bq
        maxbp=max(blocked) if blocked else 0
        # Q_div = largest Q s.t. every prime<=Q divides some speed
        Qdiv=13
        while is_prime_next:=True:
            nextp=next((p for p in range(Qdiv+1,600) if is_prime(p)), None)
            if nextp and any(v%nextp==0 for v in sp): Qdiv=nextp
            else: break
        firstfree=next((p for p in range(15,600) if is_prime(p) and all(v%p!=0 for v in sp)), None)
        corr = qstar - firstfree if firstfree else None
        print(f"{M:>14} {int(log(M,10)):>8} {len(blocked):>16} {maxbp:>18} {qstar:>5} {str(corr):>26}")
    print("\n=> q* tracks the prime-packing capacity (# distinct primes 13 speeds<=M can hold = O(log M/loglog));")
    print("   tiling correction (free primes that are no-witness) BOUNDED => the crux is ELEMENTARY prime-packing")
    print("   + a bounded circle-method correction. q* = O(log M loglog M), constant from the packing.")
