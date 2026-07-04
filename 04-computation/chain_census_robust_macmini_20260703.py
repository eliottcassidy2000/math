#!/usr/bin/env python3
"""
Robust confirmation: the compressed CHAIN band-blocker is UNIFORMLY small-q census-able (mac-mini-S27).
Adversarial chains: distinct-scale runners r_i (each within 13x of a neighbor), covering, gcd=1, engineered
to block small moduli. Is min-witness q uniformly bounded (unlike the all-~N cluster q~log N)?
WHY (hypothesis): distinct scales => residues r_i mod q are 'generic' => the 13 danger-arc bad-sets
{a: r_i a in danger} overlap generically => a free a exists at small q. The 13/7>1 union bound is defeated
by the SCALE DIVERSITY (each runner constrains a different residue direction).
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x=x%1; return min(x,1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def min_witness_q(sp, qmax=300):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q)>=1/14 for v in sp): return q
    return None

if __name__ == "__main__":
    rng = random.Random(271)
    print("Robust chain census: min-witness q over adversarial distinct-scale chains (ratio-step <=13).")
    print("=" * 84)
    print(f"{'ratio-step':>10} {'top scale ~':>14} {'#cov gcd1 tested':>16} {'MAX witness q':>14} {'#q>40':>6}")
    overall = 0
    for rstep in [2, 3, 5, 8, 12, 13]:
        for topexp in [6, 10, 15, 20]:
            base = rng.randint(1, 5)
            worst=0; nc=0; nover=0
            for _ in range(4000):
                # build a chain: r_0=base.., each ~rstep*prev + jitter, blocking-biased (make some divisible by band moduli)
                sp=[]; r=base + rng.randint(0,3)
                for i in range(13):
                    val=int(r)
                    # blocking bias: snap some to a multiple of a band modulus
                    if rng.random()<0.5:
                        m=rng.choice(range(15,50)); val=(val//m)*m or m
                    sp.append(max(1,val))
                    r *= rstep; r += rng.randint(0, max(1,int(r)//15))
                sp=sorted(set(sp))
                if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
                nc+=1
                q=min_witness_q(sp, 300)
                if q:
                    if q>worst: worst=q
                    if q>40: nover+=1
            overall=max(overall,worst)
            if nc>0:
                print(f"{rstep:>10} {'~'+str(rstep**12*base):>14} {nc:>16} {worst:>14} {nover:>6}")
    print(f"\nOVERALL max witness q over ALL adversarial chains = {overall}")
    print(f"=> if bounded (~<=40): the compressed CHAIN is UNIFORMLY census-able at small q (owner's small-q")
    print(f"   census); the crux is ONLY the all-~N one-scale wide cluster. Scale diversity defeats 13/7>1.")
