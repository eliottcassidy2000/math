#!/usr/bin/env python3
"""
lrc_treebolic_adelic_geometry_s547.py    oracle-2026-06-01-S547o

WILD STRUCTURE for the parity divide & the +2/x2 dual recursion: the natural numbers
are the shadow of the ADELIC / TREEBOLIC geometry.

THE PROPOSAL. Q has ONE archimedean completion R (a LINE, where + is natural) and ONE
per prime Q_p (a TREE, where /p is natural). So Z sits diagonally in R x prod_p Q_p
(the adeles). Two recursions:
  +2  = motion on the LINE (archimedean R); preserves the p-adic 'level' generically.
  x2  = descent in the 2-adic TREE (Q_2); the parity divide = the tree's FIRST branch.
The horocyclic product R x T_2 (the BS(1,2) treebolic space / Diestel-Leader DL(2,2))
is the geometry: + is the horocycle (horizontal), x2 is the tree (vertical).

THE COMPUTABLE CORE: every n = 2^{v2(n)} * odd(n) = (TREE-DEPTH v2, LINE-POSITION odd).
  x2: (a,m) -> (a+1, m)   -- clean VERTICAL (tree descent).
  +2: (a,m) -> erratic     -- HOROCYCLIC (weaves through tree levels).
The parity divide = a=0 (odd, tree top) vs a>=1 (even, descended). Doubled prime
n=2p = (a=1, m=prime). And the LRC channel cleanliness is read off (v2(n), odd(n)):
n=2p (a=1, m prime) -> clean; n=2*p^k -> filtered (S534/S546).

We illustrate the two recursions in (v2, odd) coords, the parity/doubled-prime
positions, and the LRC tie. Plus a light sketch of the adeles (other primes = other
trees) and why 2 is special (smallest prime + the archimedean sign).
"""
from sympy import isprime, factorint

def v2(n):
    a=0
    while n%2==0: n//=2; a+=1
    return a
def odd_part(n):
    while n%2==0: n//=2
    return n
def coord(n): return (v2(n), odd_part(n))

def main():
    print("="*72)
    print("THE TREEBOLIC DECOMPOSITION: n = 2^{v2} * odd = (tree-depth, line-position)")
    print("="*72)
    print("  n  -> (v2=tree depth, odd=line position)   parity   note")
    for n in range(1,21):
        a,m=coord(n)
        par="ODD (tree top)" if a==0 else f"even (depth {a})"
        note=""
        if a==1 and isprime(m): note=f"DOUBLED PRIME 2*{m}"
        if a>=1 and m==1: note=f"pure power 2^{a}"
        print(f"  {n:2d} -> (v2={a}, odd={m:2d})   {par:16s} {note}")
    print()

    print("="*72)
    print("THE TWO RECURSIONS as the two directions of the geometry")
    print("="*72)
    print("  x2 (TREE DESCENT, clean vertical): (a,m)->(a+1,m)")
    n=3
    for _ in range(5):
        a,m=coord(n); print(f"     {n:3d} = (v2={a}, odd={m})"); n*=2
    print("  +2 (HOROCYCLIC, erratic in tree depth): weaves through levels")
    n=2
    traj=[]
    for _ in range(10):
        a,m=coord(n); traj.append((n,a,m)); n+=2
    print("     "+"  ".join(f"{nn}:(d{a})" for nn,a,m in traj))
    print(f"     tree-depths visited by +2 from 2: {[a for _,a,_ in traj]} (ERRATIC = horocyclic)")
    print("  => x2 moves straight down ONE tree (clean); +2 is a horocycle weaving across")
    print("     tree levels. The two recursions = the two perpendicular directions of R x T_2.")
    print()

    print("="*72)
    print("PARITY DIVIDE = the 2-adic tree's FIRST branch; bridge x2 = descent")
    print("="*72)
    print("  odd numbers (v2=0) = the tree TOP (the 'odd line'); even = descended via x2.")
    print("  ODD = atomic/free (tree top); EVEN = doubled (one+ levels down). The bridge that")
    print("  carries an odd number n across the divide is x2: n -> 2n (down one level).")
    print()

    print("="*72)
    print("LRC TIE: channel cleanliness from (v2(n), odd(n)); the ADELIC interface")
    print("="*72)
    print("  n   (v2,odd)     n*=n/2^[even]   odd-part type        LRC channels")
    for n in (13,14,15,16,17,18,22):
        a,m=coord(n)
        ns = n//2 if n%2==0 else n
        f=factorint(m)
        if m==1: mt="(unit, pure 2-power)"
        elif len(f)==1 and list(f.values())[0]==1: mt="PRIME"
        elif len(f)==1: mt=f"prime-power {list(f.keys())[0]}^{list(f.values())[0]}"
        else: mt="composite"
        clean = "CLEAN" if (a==1 and "PRIME" in mt) or (a==0 and "PRIME" in mt) else "filtered/messy"
        print(f"  {n:2d}  (v2={a},odd={m:2d})   n*={ns:2d}          {mt:20s} {clean}")
    print("  => n=14=(v2=1,odd=7 PRIME) -> CLEAN (doubled prime); n=16=(v2=4,odd=1) pure 2-power,")
    print("     n=18=(v2=1,odd=9=3^2) prime-power -> filtered. The (v2, odd) treebolic coordinates")
    print("     of n predict its LRC channel structure. The ARCHIMEDEAN runner condition (||v t||")
    print("     >=1/n on R/Z) is obstructed by the p-ADIC channels (the trees) -- LRC is ADELIC.")
    print()
    print("="*72)
    print("WILD STRUCTURES (proposed): (1) BS(1,2) treebolic R x T_2 / Diestel-Leader DL(2,2)")
    print("(+ = horocycle, x2 = tree). (2) ADELES R x prod Q_p (+ at R, /p at the p-trees;")
    print("parity = Q_2's first branch). (3) the dyadic SOLENOID (inverse limit under x2).")
    print("(4) exp/log bridging + and x (primes = multiplicative frequencies). (5) WITT vectors")
    print("(+ and x unified via ghost components). (6) CAYLEY-DICKSON dim-doubling 1,2,4,8 (repo).")
    print("(7) STERN-BROCOT/Calkin-Wilf binary tree of Q+. WHY 2 is special: smallest prime AND")
    print("the archimedean place carries the SIGN (the order), so the +/x tension centers on 2.")
    print("="*72)

if __name__=="__main__":
    main()
