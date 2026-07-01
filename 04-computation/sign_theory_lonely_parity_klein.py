#!/usr/bin/env python3
"""
sign_theory_lonely_parity_klein.py  --  klein-2026-06-30-S55

SIGN THEORY attack on the OPEN CRUX (the multi-metric danger sheaf COVERS the modulus site =
no lonely point survives all metrics).

The sign structure = the ANTIPODAL involution  iota: a -> -a  (on rotations mod D),
equivalently t -> 1-t (the complement map, THM-584).  Since ||x|| is EVEN (||-x||=||x||),
the loneliness landscape is iota-invariant, so:

  (P) PARITY LEMMA.  At an ODD prime p, radius r, the set of LONELY rotations
      L(p,r) = { a in (Z/p)* : sa mod p avoids {+-1,...,+-r} for all runners s }
      is iota-invariant with NO fixed point (a=-a => a=0, excluded), hence #L(p,r) is EVEN.
      => lonely rotations come in +- PAIRS; a witness is an iota-SYMMETRIC hole.

  (LEVER) "WITHIN-1 SUFFICES".  To prove coverage (#L=0) it is enough to show #L <= 1:
      parity then forces #L=0.  So any counting/measure bound only needs to reach
      p-2 covered rotations (all but one), and the sign-parity closes the last gap.
      This HALVES the strength required of the danger-count.

  (WITNESS as iota-pair) the (n+q)-witness hole is exactly the iota-pair {+1,-1} mod (n+q).

This script VERIFIES (P), exhibits the AP as the unique zero-lonely set with near-AP deviations
creating iota-PAIRS, and quantifies the "within-1 suffices" lever.
"""
from fractions import Fraction as F

def dist0(x,D):
    x%=D; return min(x,D-x)
def isprime(p): return p>1 and all(p%d for d in range(2,int(p**.5)+1))

def lonely_rotations(S, p, r):
    """rotations a in 1..p-1 with all runners avoiding {+-1..+-r} mod p (radius-r lonely)."""
    out=[]
    for a in range(1,p):
        if all(dist0(s*a,p)>r for s in S):
            out.append(a)
    return out

def check_parity(sets, n):
    print("(P) PARITY: #lonely rotations is EVEN at every odd prime (radius r=floor(p/n)):")
    ok=True
    for name,S in sets:
        line=[]
        for p in [pp for pp in range(3,4*n) if isprime(pp)]:
            r=p//n
            L=lonely_rotations(S,p,r)
            # iota-pairing check: a lonely <=> p-a lonely
            paired=all(((p-a) in L) for a in L)
            if len(L)%2!=0 or not paired: ok=False
            if L: line.append(f"p={p}(r={r}):#{len(L)}")
        print(f"   {name:14s}: {'all even & iota-paired' if ok else 'PARITY VIOLATION'};  lonely: {line[:8]}")
    return ok

if __name__=="__main__":
    n=14
    AP=list(range(1,n))                       # tight AP {1..13}, M=1/n
    GW=[1,2,3,4,5,6,7,8,9,10,11,13,24]        # the doubling (tight)
    swap11=[x for x in AP if x!=11]+[22]      # drop prime 11, add 22 (not tight -> (n+q)-witness)
    swap13=[x for x in AP if x!=13]+[26]      # drop prime 13, add 26
    sets=[("AP",AP),("GW",GW),("swap11(+22)",swap11),("swap13(+26)",swap13)]

    check_parity(sets,n)

    print()
    print("(WITNESS = iota-pair) the (n+q)-witness hole is the pair {+1,-1} mod (n+q):")
    for q,g in [(11,22),(13,26)]:
        D=n+q; a=pow(q,-1,D); S=[x for x in AP if x!=q]+[g]
        # residues that map to +-1 under a: speeds q and n; both absent
        holes=[res for res in range(1,D) if all(dist0(s*a,D)!=res for s in S) and dist0(res,D)==1]
        print(f"   q={q}: D=n+q={D}, a=q^-1={a}; uncovered radius-1 residues (the iota-pair) = {sorted(set(dist0(h,D) for h in holes))} -> {holes}")

    print()
    print("(LEVER) 'within-1 suffices': AP is zero-lonely EVERYWHERE; a single swap opens an EVEN hole")
    print("   (>=2, an iota-pair), never exactly 1 -- so a danger-count reaching p-2 (all but one) => full")
    print("   coverage by parity. Max #lonely a single swap opens (the pair count) at n=14:")
    for name,S in [("swap11(+22)",swap11),("swap13(+26)",swap13),("GW",GW)]:
        worst=0; where=None
        for p in [pp for pp in range(3,4*n) if isprime(pp)]:
            r=p//n; L=lonely_rotations(S,p,r)
            if len(L)>worst: worst=len(L); where=(p,r)
        print(f"   {name:14s}: max #lonely = {worst} at {where} (even: {worst%2==0})")
