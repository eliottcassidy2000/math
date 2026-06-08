#!/usr/bin/env python3
"""
S724 — Erdos problem 2 (covering systems, min modulus) and its bridge to the LRC covering-depth.

ERDOS 2 (DISPROVED, $1000): a covering system is finitely many congruences x = a_i (mod n_i) with
DISTINCT moduli covering every integer. Erdos asked: can min(n_i) be arbitrarily large? Conjectured YES
(his favourite problem). Hough (2015) DISPROVED it: min modulus <= 10^16. Balister-Bollobas-Morris-
Sahasrabudhe-Tiba (2022, "distortion method"): min modulus <= 616000; and NO covering with all moduli odd
& squarefree. Best known CONSTRUCTION (lower bound): min modulus 42 (Owens, BYU thesis; Nielsen had 40).

THE BRIDGE TO THE REPO. The LRC is the CONTINUOUS covering problem and Erdos 2 is the DISCRETE one, with
the SAME inclusion-exclusion skeleton:
  - covering system: uncovered density = Sum_{I subseteq C} (-1)^|I| [I compatible] / lcm(moduli in I)
    (CRT; this is exactly S561's sieve rho(S) = Sum_T (-1)^|T| / lcm(T)).
  - LRC (THM-406 cor.b): lonely measure p_0 = Sum_{j} (-1)^j S_j,  S_j = sum of j-fold danger-arc overlaps.
  Same alternating-sum form; the covering system is the rational/CRT SKELETON of the LRC danger arcs.
  "Covers everything" <=> p_0 = 0; "leaves a gap" <=> p_0 > 0 (a lonely time / an uncovered residue).

THE METHOD BRIDGE. Hough/BBMST proved "spread-out (large-modulus) pieces always leave a positive
uncovered density" (a distortion / density-increment). That is the discrete sibling of the LRC's "the
cover always leaks" (p_0 > 0 for multiple-of-n configs; my S643). The distortion method is the tool the
LRC wants. And both share S561's lesson: density Sum 1/n_i can exceed 1 yet NOT cover -- the obstruction
is the inclusion-exclusion OVERLAP structure, not the density.

42 IN THE REPO: incidental (a(42) index in A038375; a Paley-table entry). The COVERING structure is what
is central (THM-406 covering-depth = the LRC master object; S561 rho = the covering density).

No numpy/sympy.
"""
from math import gcd
from itertools import combinations
from fractions import Fraction as Fr

def lcm(a,b): return a*b//gcd(a,b)
def lcm_list(L):
    r=1
    for x in L: r=lcm(r,x)
    return r

def covers(congs):
    """congs = list of (a,n). Does the union cover Z (i.e. Z/L)? Return (covers?, uncovered_density)."""
    L=lcm_list([n for _,n in congs])
    covered=[False]*L
    for a,n in congs:
        x=a%n
        while x<L: covered[x]=True; x+=n
    unc=covered.count(False)
    return unc==0, Fr(unc,L)

def uncovered_IE(congs):
    """uncovered density via inclusion-exclusion: Sum_I (-1)^|I| [compatible]/lcm(I) (= S561 rho form)."""
    C=congs; total=Fr(0)
    for r in range(0,len(C)+1):
        for I in combinations(range(len(C)),r):
            # compatible? all pairwise a_i = a_j mod gcd(n_i,n_j)
            ok=True
            for x in range(len(I)):
                for y in range(x+1,len(I)):
                    ai,ni=C[I[x]]; aj,nj=C[I[y]]
                    if (ai-aj)%gcd(ni,nj)!=0: ok=False;break
                if not ok: break
            if not ok: continue
            l=lcm_list([C[i][1] for i in I]) if I else 1
            total+=Fr((-1)**r, l)
    return total

def naive_density(congs): return sum(Fr(1,n) for _,n in congs)

if __name__=="__main__":
    print("="*88)
    print("S724 — Erdos 2 (covering systems, min modulus) <-> LRC covering-depth (THM-406)")
    print("="*88)

    # (1) the classic min-modulus-2 covering system; uncovered density 0 = covers
    print("\n(1) classic min-modulus-2 covering: {0(2),0(3),1(4),5(6),7(12)}")
    S=[(0,2),(0,3),(1,4),(5,6),(7,12)]
    cov,unc=covers(S); ie=uncovered_IE(S)
    print(f"   covers Z: {cov}; uncovered density (direct)={unc}; inclusion-exclusion rho={ie}; match={unc==ie}")
    print(f"   naive density Sum 1/n_i = {naive_density(S)} (>1: overlap {naive_density(S)-1} wasted)")
    print("   => uncovered density = S561's rho = THM-406's p_0-form (Sum_I (-1)^|I|/lcm). Covers <=> rho=0.")

    # (2) DENSITY IS NOT SUFFICIENT: Sum 1/n_i > 1 yet leaks (the S561 / Hough lesson)
    print("\n(2) density > 1 but does NOT cover (the obstruction is OVERLAP structure, not density)")
    # min modulus 3, distinct moduli, residues chosen poorly: density>1 but leaks
    bad=[(0,3),(1,3+0)]  # placeholder; build a real leaking example
    bad=[(0,3),(1,4),(2,5),(0,6),(3,8),(4,9),(5,10),(7,12)]
    cov,unc=covers(bad); dens=naive_density(bad)
    print(f"   moduli {[n for _,n in bad]} (min 3): Sum 1/n_i = {dens} (~{float(dens):.3f}); covers={cov}; "
          f"uncovered density={unc} (~{float(unc):.4f})")
    print("   => even with density >~1, large-ish min modulus leaves a positive uncovered density: the LEAK.")

    # (3) the structural growth: min modulus m needs SMOOTH moduli (lcm with many small primes)
    print("\n(3) STRUCTURE of large-min-modulus coverings: moduli must be SMOOTH (lcm has many small primes)")
    print("   record min moduli: 2 (Erdos) -> 3,4,... -> 40 (Nielsen 2009, ~10^50 congruences) ->")
    print("   42 (Owens, BYU) = 2*3*7; the moduli are divisors of a highly composite L; HOUGH(2015): min")
    print("   modulus <= 10^16; BBMST(2022): <= 616000 (distortion method) + no all-odd-squarefree covering.")
    print("   => Erdos 2 is DISPROVED: the minimum modulus is BOUNDED. 'Spread-out pieces always leak.'")

    # (4) the LRC bridge: same alternating-sum; 'cover always leaks' = Hough's bound
    print("\n(4) LRC BRIDGE: continuous covering. THM-406 p_0 = Sum_j (-1)^j S_j (danger-arc overlaps) has the")
    print("   SAME form as the covering-system rho = Sum_I (-1)^|I|/lcm (S561). The covering system is the")
    print("   rational/CRT SKELETON of the LRC danger arcs. 'Multiple-of-n LRC config is loose' (p_0>0, S643)")
    print("   is the continuous sibling of Hough's 'large-min-modulus covering leaks' -- both: spread pieces")
    print("   leave a positive gap. The DISTORTION / density-increment method (Hough/BBMST) is exactly the")
    print("   structural tool the LRC needs (S561: density alone can't close it; cut by structure).")

    print("\n(5) 42 in the repo: incidental (a(42) in A038375; a Paley-table entry). The COVERING is central:")
    print("    THM-406 covering-depth = LRC master object; S561 rho = covering density; THM-439 witness tower.")
    print("="*88)
