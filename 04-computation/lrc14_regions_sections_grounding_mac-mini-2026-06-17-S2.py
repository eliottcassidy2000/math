#!/usr/bin/env python3
"""
lrc14_regions_sections_grounding — mac-mini-2026-06-17-S2

REGIONS-OF-THE-LOOP reframe (user's idea): track the n=14 SECTIONS of the circle,
not the 13 runners. At a grid time tau=a/14 (a coprime to 14), runner i sits in
section r_i = v_i*a mod 14 (sections 0..13, each an arc of length 1/14). The
observer (0) is LONELY at a/14  <=>  no r_i = 0  (then ||v_i a/14|| >= 1/14).
The user's "each runner uses its OWN section" = the residues {r_i} are DISTINCT
(a System of Distinct Representatives / perfect spreading), leaving section 0 clear.

This grounds the picture:
 (1) {1..13} is a PERFECT SDR at every unit a (runners biject onto sections 1..13);
 (2) covering-set hard core (contains 14m) FAILS the grid (a runner sits in section 0)
     -> the lonely time is OFF the 1/14 grid; we show the section occupancy there;
 (3) the (Z/14)* action permutes the section-assignment (the "switches");
 (4) first look at the overtaking / cyclic-order (tournament) structure.
"""
from fractions import Fraction as F
from math import gcd

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def Mexact(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2))
    best=F(0); at=None
    for t in C:
        v=g(S,t)
        if v>best: best=v; at=t
    return best, at
units14=[a for a in range(1,14) if gcd(a,14)==1]

def section_profile(S, a, N=14):
    """residues r_i = v_i*a mod N (the section each runner occupies at tau=a/N)."""
    return [ (v*a) % N for v in S ]

def analyze_grid(S, name):
    print(f"\n--- {name}: S={S} ---")
    print(f"  speeds mod 14: {[v%14 for v in S]}  (distinct & nonzero => SDR-able)")
    lonely_a=[]
    for a in units14:
        r=section_profile(S,a)
        zero = 0 in r
        sdr = len(set(r))==len(r)
        if not zero: lonely_a.append(a)
        tag = "OBSERVER HIT (sec 0)" if zero else ("PERFECT SDR" if sdr else f"clumped (only {len(set(r))} sections)")
        if a in (1,3,5): print(f"  a={a:2d}: sections {sorted(r)}  -> {tag}")
    print(f"  grid-14 lonely at a in {lonely_a}  (empty => off-grid lonely time needed)")

print("="*74)
print("(1)-(2) GRID-14 SECTION OCCUPANCY  r_i = v_i*a mod 14")
print("="*74)
analyze_grid(list(range(1,14)), "TIGHT AP {1..13} (no mult of 14)")
analyze_grid([1,2,3,4,5,6,7,8,9,10,11,13,84], "COVERING hard core {1..11,13,84} (84=6*14)")
analyze_grid([1,2,3,4,5,7,8,9,10,11,12,13,98], "interior-drop j=6 u 98")

print("\n"+"="*74)
print("(3) OFF-GRID lonely time for the hard core: where do runners sit?")
print("="*74)
for S,name in ([ [1,2,3,4,5,6,7,8,9,10,11,13,84], "hardcore 84"],
               [ [1,2,3,4,5,7,8,9,10,11,12,13,98], "drop6 u98"]):
    M,tau=Mexact(S)
    N=tau.denominator
    print(f"\n  {name}: M={M}={float(M):.5f}, lonely tau*={tau} (denom {N} = #fine sections)")
    pos=sorted((v*tau)%1 for v in S)
    # section index in 14 equal arcs
    secs=[int((float((v*tau)%1))*14) for v in S]
    from collections import Counter
    occ=Counter(secs)
    print(f"  14-section occupancy {{sec:count}}: {dict(sorted(occ.items()))}")
    print(f"  sections EMPTY (no runner): {sorted(set(range(14))-set(secs))}  (observer's gap lives in an empty/clear band)")
    print(f"  is each runner in its own 14-section? {len(set(secs))==len(secs)}")

print("\n"+"="*74)
print("(4) The (Z/14)* action permutes the SDR (the 'switches' for {1..13})")
print("="*74)
S=list(range(1,14))
print("  a -> section-permutation of runners (sigma_a(i)=v_i*a mod 14):")
for a in units14:
    perm=[(i*a)%14 for i in range(1,14)]
    print(f"   a={a:2d}: {perm}")
print("  (each a in (Z/14)* gives a DISTINCT bijection runners<->sections: the symmetry group of the nice case)")

print("\n"+"="*74)
print("(5) OVERTAKING structure (first look at the tournament): crossings & cyclic order")
print("="*74)
S=list(range(1,14))
M,tau=Mexact(S)
print(f"  at lonely tau*={tau}, cyclic order of runners by position:")
order=sorted(range(len(S)), key=lambda i: float((S[i]*tau)%1))
print(f"   {[S[i] for i in order]}  (positions {[float((S[i]*tau)%1) for i in order]})")
print("  pairwise relative laps |v_i-v_j| (crossings per period) = the wiring-diagram swap counts:")
print("   sum of |v_i-v_j| over pairs =", sum(abs(S[i]-S[j]) for i in range(len(S)) for j in range(i+1,len(S))))
