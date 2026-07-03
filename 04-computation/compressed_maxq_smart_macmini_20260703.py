#!/usr/bin/env python3
"""
SMART adversary: compressed covering gcd=1 families with runners = lcm-PRODUCTS blocking MANY band moduli
(mac-mini-2026-07-03-S26). Each runner ~N divisible by lcm(group of band moduli), lcm <= N. 13 runners can
block ~13*(moduli-per-lcm) moduli. At large N, more moduli per lcm => does compressed q grow (=> log-census)
or stay bounded (=> fixed census closes compressed)?
"""
from math import gcd, log10
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def lcm(a,b): return a*b//gcd(a,b)
def lcm_list(xs): return reduce(lcm, xs, 1)
def nd(x): x = x % 1; return min(x,1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def is_compressed(sp):
    return all(any(j!=i and abs(sp[i])<13*abs(sp[j]) for j in range(len(sp))) for i in range(len(sp)))
def min_witness_q(sp, qmax=800):
    for q in range(2, qmax+1):
        dang = None
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q) >= 1/14 for v in sp): return q
    return None

def build_lcm_blocker(N, target_moduli, rng):
    """runners ~N each = lcm(a greedy group of target_moduli) * round(N/lcm), blocking those moduli."""
    mods = list(target_moduli); rng.shuffle(mods)
    runners = []
    i = 0
    while i < len(mods) and len(runners) < 13:
        grp = []; L = 1
        while i < len(mods):
            L2 = lcm(L, mods[i])
            if L2 <= N and len(grp) < 8:
                L = L2; grp.append(mods[i]); i += 1
            else:
                break
        if L > 1:
            m = L * max(1, round(N / L))
            runners.append(m)
        else:
            i += 1
    # pad to 13 with small covering fillers + ensure gcd 1
    for q in [8,9,5,7,11,13,2,3,4,6]:
        if len(runners) >= 13: break
        if not any(s % q == 0 for s in runners): runners.append(q)
    while len(runners) < 13: runners.append(rng.randint(2,22))
    return sorted(set(runners))[:13]

if __name__ == "__main__":
    print("SMART adversary (lcm-product runners): max witness q over compressed covering gcd=1 vs N.")
    print("=" * 84)
    print(f"{'N':>12} {'log10':>6} {'target band':>12} {'#cc tested':>11} {'MAX witness q':>14}")
    rows = []
    for N in [10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**12]:
        rng = random.Random(N)
        # target band up to where lcm(15..Q) roughly exceeds N^13 (the primorial cap)
        Qtarget = 15
        while lcm_list(range(15, Qtarget+2)) <= N**6 and Qtarget < 200:
            Qtarget += 1
        worst = 0; ncc = 0
        for _ in range(6000):
            sp = build_lcm_blocker(N, range(15, Qtarget+1), rng)
            if len(sp) != 13 or gcd_all(sp) != 1: continue
            if not is_covering(sp) or not is_compressed(sp): continue
            ncc += 1
            q = min_witness_q(sp, qmax=800)
            if q and q > worst: worst = q
        rows.append((N, worst))
        print(f"{N:>12} {log10(N):>6.0f} {'15..'+str(Qtarget):>12} {ncc:>11} {worst:>14}")
    qs=[q for _,q in rows]; logs=[log10(N) for N,_ in rows]
    if len(set(logs))>1:
        n=len(qs); sx=sum(logs); sy=sum(qs); sxx=sum(x*x for x in logs); sxy=sum(x*y for x,y in zip(logs,qs))
        slope=(n*sxy-sx*sy)/(n*sxx-sx*sx)
    else: slope=0
    print(f"\nMAX-q sequence: {qs};  slope ~ {slope:.2f} per decade")
    print(f"=> slope>0 & growing: compressed q ~ log(mag), UNBOUNDED -> renormalization/log-census is the crux.")
    print(f"   slope~0 & bounded: fixed census (q<=~{max(qs) if qs else 0}) closes the compressed case entirely.")
