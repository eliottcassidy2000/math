#!/usr/bin/env python3
"""
S642 — LRC attack: the ±-transversal reduction of C'(n) for 2n-1 prime.
THEOREM (elementary): p=2n-1 prime, S a speed set (n-1 distinct), no v_i≡0 mod p. The bad multipliers
at shell p (a with some a*v_i in {±1}) are ∪_i {±v_i^{-1}}, ≤ 2(n-1)=p-1, with EQUALITY iff the
residues {v_i mod p} are a ±-TRANSVERSAL (one per ±-pair). So non-transversal => a good multiplier
exists => M >= 2/p = 2/(2n-1) > 1/n => LOOSE.  Hence C'(n) [2n-1 prime] reduces to the ±-transversals.
We verify the theorem and characterize the residual (transversal + must band-block EVERY shell).
"""
import itertools, random
from math import gcd
from fractions import Fraction as Fr

def good_mult_fast(S, m):
    """a multiplier a in (Z/m)* with all a*v off the central band {dist<=1} (i.e. min(r,m-r)>=2); or None."""
    for a in range(1, m):
        if gcd(a, m) != 1: continue
        if all(min((a*v) % m, m-((a*v) % m)) >= 2 for v in S):
            return a
    return None

def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True

def is_pm_transversal(residues, p):
    """residues (mod p, all nonzero): a ±-transversal iff all in distinct ±-pairs (one per pair)."""
    classes=set()
    for r in residues:
        if r%p==0: return False
        c=min(r%p, (-r)%p)
        if c in classes: return False
        classes.add(c)
    return len(classes)==(p-1)//2   # hits every pair

def banded_at(S, m):
    """does shell m FAIL to dodge? (no multiplier a with all a*v off band {0,±1}) = band-blocking."""
    return good_mult_fast(list(S), m) is None if isprime(m) else (
        all(any(min((a*v)%m,m-((a*v)%m))<2 for v in S) for a in range(1,m) if gcd(a,m)==1))

if __name__=="__main__":
    for n in (15,19):
        p=2*n-1
        if not isprime(p):
            print(f"n={n}: 2n-1={p} not prime, skip"); continue
        print(f"\n=== n={n}, 2n-1={p} (prime, unramified) ===")
        rng=random.Random(0)
        N=6000; nontrans_loose=0; nontrans=0; trans=0; trans_configs=[]
        for _ in range(N):
            S={n*rng.randint(1,2)}
            while len(S)<n-1: S.add(rng.randint(1,3*n))
            S=tuple(sorted(S))
            if len(S)!=n-1: continue
            g=0
            for x in S: g=gcd(g,x)
            if g!=1: continue
            res=[v%p for v in S]
            if any(r==0 for r in res): continue   # divisible-by-p handled separately (dominant)
            if is_pm_transversal(res,p):
                trans+=1
                if len(trans_configs)<2000: trans_configs.append(S)
            else:
                nontrans+=1
                # THEOREM: should be loose via shell p. verify good multiplier exists & M>=2/p
                a=good_mult_fast(list(S),p)
                if a is not None: nontrans_loose+=1
        print(f"  non-transversal configs: {nontrans}; of these LOOSE via shell-{p} dodge (M>=2/{p}): {nontrans_loose}  (theorem: ALL)")
        print(f"  ±-transversal configs (the RESIDUAL): {trans}  ({100*trans/max(trans+nontrans,1):.1f}% of configs)")
        # characterize residual: are transversals band-blocking at OTHER shells too? and loose how?
        if trans_configs:
            stillhard=0; loose_via_lower=0
            for S in trans_configs[:400]:
                # try dodge at any shell m<=2n-1 (m != p, or lower 1-clocks)
                escaped=False
                for m in range(2,2*n):
                    if any(v%m==0 for v in S): continue
                    if good_mult_fast(list(S),m) is not None: escaped=True; break
                if escaped: loose_via_lower+=1
                else: stillhard+=1
            print(f"  of {min(400,len(trans_configs))} transversal residual configs: loose via SOME shell m<=2n-1: {loose_via_lower}; still-hard (need B'): {stillhard}")
