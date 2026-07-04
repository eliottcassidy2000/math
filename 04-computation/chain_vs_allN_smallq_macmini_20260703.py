#!/usr/bin/env python3
"""
Chain vs all-~N band blocker: which needs large q? (mac-mini-2026-07-03-S27, per owner hint)
The GAP (opus-S54: ratio>13, not all-comparable) has two shapes:
 (A) GEOMETRIC CHAIN: runners at r_i ~ 12^i * base (each within 13x of the next, spanning widely). Distinct
     SCALES => residues mod small q ~ independent => likely a SMALL-q census.
 (B) ALL-~N lcm-blocker: 13 runners all ~N (near-equal MAGNITUDE) + small fillers, each blocking a band
     modulus. Same scale => the small-q census is defeated (my S26 q~128).
Determine the small-q census (min witness q) for each. If the CHAIN is small-q but all-~N is large-q, the
real crux is the all-~N near-equal-MAGNITUDE blocker, and the chain is census-able (owner's small-q census).
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def lcm(a,b): return a*b//gcd(a,b)
def lcm_list(xs): return reduce(lcm, xs, 1)
def nd(x): x=x%1; return min(x,1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def ratio(sp): return max(abs(v) for v in sp)/min(abs(v) for v in sp)
def min_witness_q(sp, qmax=400):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q)>=1/14 for v in sp): return q
    return None

def geometric_chain(base, ratio_step, rng):
    """13 runners r_i ~ ratio_step^i * base, covering, gcd=1."""
    sp=[]
    r=base
    for i in range(13):
        # perturb slightly, keep within 13x of neighbor
        sp.append(int(r) + rng.randint(0, max(1, int(r)//20)))
        r *= ratio_step
    return sorted(set(sp))

if __name__ == "__main__":
    rng = random.Random(27)
    print("(A) GEOMETRIC CHAIN band-blockers: min witness q vs base scale.")
    print("=" * 76)
    print(f"{'base':>8} {'ratio-step':>10} {'max/min':>10} {'#cov tested':>12} {'MAX witness q':>14}")
    for base in [1, 3, 5, 10]:
        for rstep in [3, 8, 12]:
            worst=0; nc=0
            for _ in range(3000):
                sp = geometric_chain(base + rng.randint(0,4), rstep, rng)
                if len(sp)!=13: continue
                # make covering by adjusting: only keep covering gcd=1
                if gcd_all(sp)!=1 or not is_covering(sp): continue
                nc+=1
                q=min_witness_q(sp, 400)
                if q and q>worst: worst=q
            if nc>0:
                print(f"{base:>8} {rstep:>10} {'~'+str(rstep**12):>10} {nc:>12} {worst:>14}")

    print("\n(B) ALL-~N lcm-blockers (near-equal MAGNITUDE + fillers): min witness q vs N.")
    print("=" * 76)
    def build_lcm_blocker(N, rng):
        mods=list(range(15, 80)); rng.shuffle(mods); runners=[]; i=0
        while i<len(mods) and len(runners)<13:
            grp=[]; L=1
            while i<len(mods):
                L2=lcm(L,mods[i])
                if L2<=N and len(grp)<8: L=L2; grp.append(mods[i]); i+=1
                else: break
            if L>1: runners.append(L*max(1,round(N/L)))
            else: i+=1
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(runners)>=13: break
            if not any(s%q==0 for s in runners): runners.append(q)
        while len(runners)<13: runners.append(rng.randint(2,22))
        return sorted(set(runners))[:13]
    print(f"{'N':>12} {'max/min':>12} {'#cov tested':>12} {'MAX witness q':>14}")
    for N in [10**5, 10**7, 10**9, 10**12]:
        worst=0; nc=0
        for _ in range(3000):
            sp=build_lcm_blocker(N, rng)
            if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
            nc+=1
            q=min_witness_q(sp, 400)
            if q and q>worst: worst=q
        print(f"{N:>12} {'~'+str(int(ratio(build_lcm_blocker(N,rng)))):>12} {nc:>12} {worst:>14}")
    print("\n=> compare: if CHAIN stays small-q while all-~N grows, the crux is the near-equal-MAGNITUDE")
    print("   blocker (all-~N), and the chain is census-able. The 13/7>1 wall bites only when all 13")
    print("   danger arcs are at ONE scale (all-~N) so they cannot be dodged by a small common denominator.")
