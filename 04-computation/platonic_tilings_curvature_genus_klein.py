#!/usr/bin/env python3
"""
platonic_tilings_curvature_genus_klein.py  --  klein-2026-06-30-S57

The 5 Platonic solids, the 3 regular plane tilings, and the 5 planar Bravais lattices: where the
counts come from (Schlafli {p,q}), and how they connect to the repo's proofs (hexagonal covering-min,
the antipodal duality of S55/S56 sign theory, and the CURVATURE = GENUS trichotomy of S56).

(1) SCHLAFLI CLASSIFICATION: a regular {p,q} (p-gons, q around each vertex) has combinatorial
    curvature kappa(p,q) = 1/p + 1/q - 1/2. Sign of kappa <=> sign of 4 - (p-2)(q-2):
      kappa>0  (p-2)(q-2)<4  -> SPHERE, chi=2, genus 0 : the 5 PLATONIC solids
      kappa=0  (p-2)(q-2)=4  -> PLANE,  chi=0, genus 1 : the 3 regular EUCLIDEAN tilings
      kappa<0  (p-2)(q-2)>4  -> HYPERBOLIC, chi<0, genus>=2 : infinitely many
    The count "5" (Platonic) and "3" (plane) are the integer solutions of one Diophantine inequality.

(2) DUALITY = the antipodal involution: (p,q) <-> (q,p). Self-dual = iota-FIXED; the rest in dual
    PAIRS. This mirrors the tournament complement / the covering-min iota (THM-584, HYP-3767) and the
    SC-spine (self-complementary = self-dual).

(3) CURVATURE = GENUS (the S56 bridge): the tiling trichotomy chi=2/0/<0 (g=0/1/>=2) IS the genus
    trichotomy of X_0(2p): genus 0,0,1,2,2 for p=3,5,7,11,13 (n=6,10,14,22,26). n=14 => genus 1 =>
    the FLAT / EUCLIDEAN / PLANE-TILING case -- and the covering-min lives on the hexagonal A_2 lattice
    (Phi_6=Eisenstein, HYP-3715). LRC-14 is the plane-tiling case; the Platonic solids are the
    genus-0 small-n analogues, hyperbolic the genus->2 large-n ones.

(4) CRYSTALLOGRAPHIC RESTRICTION: 2D lattice rotation orders are {1,2,3,4,6} -- 5-fold FORBIDDEN. The
    covering-min is 6-fold (hexagonal Phi_6), NOT 5-fold; 5-fold (icosahedral, golden ratio) is pushed
    to the Fibonacci/quasicrystal thread.
"""
from fractions import Fraction as F

def classify(p,q):
    kappa = F(1,p)+F(1,q)-F(1,2)
    disc = (p-2)*(q-2)   # <4 sphere, =4 plane, >4 hyperbolic
    if disc<4: return "sphere/Platonic (chi=2,g=0)", kappa
    if disc==4: return "PLANE tiling (chi=0,g=1)", kappa
    return "hyperbolic (chi<0,g>=2)", kappa

NAMES={(3,3):"tetrahedron",(3,4):"octahedron",(4,3):"cube",(3,5):"icosahedron",(5,3):"dodecahedron",
       (3,6):"triangular tiling",(6,3):"HEXAGONAL tiling",(4,4):"square tiling"}

if __name__=="__main__":
    print("="*74)
    print("(1)+(2) SCHLAFLI {p,q}: the counts 5 (Platonic) and 3 (plane) from (p-2)(q-2) vs 4")
    print("="*74)
    plat=[]; plane=[]
    for p in range(3,8):
        for q in range(3,8):
            cls,kappa=classify(p,q)
            if "Platonic" in cls: plat.append((p,q))
            elif "PLANE" in cls: plane.append((p,q))
    print(f"   PLATONIC (5): {[(pq,NAMES.get(pq,'')) for pq in plat]}")
    print(f"   PLANE    (3): {[(pq,NAMES.get(pq,'')) for pq in plane]}")
    print(f"   count check: Platonic={len(plat)} (=5?), plane={len(plane)} (=3?)")
    print()
    print("   DUALITY (p,q)<->(q,p) = the antipodal iota (THM-584/HYP-3767):")
    for grp,lst in [("Platonic",plat),("plane",plane)]:
        seen=set(); pairs=[]
        for (p,q) in lst:
            if (p,q) in seen: continue
            if p==q: pairs.append(f"{NAMES.get((p,q),(p,q))} [SELF-DUAL=iota-fixed]"); seen.add((p,q))
            else:
                pairs.append(f"{NAMES.get((p,q),(p,q))}<->{NAMES.get((q,p),(q,p))}"); seen.add((p,q)); seen.add((q,p))
        print(f"     {grp}: {pairs}")

    print()
    print("="*74)
    print("(3) CURVATURE = GENUS: the tiling trichotomy IS the X_0(2p) genus trichotomy (S56)")
    print("="*74)
    genus={3:0,5:0,7:1,11:2,13:2}   # genus X_0(2p), HYP-3587/3041
    print(f"   {'p':>3} {'n=2p':>5} {'genus X_0(2p)':>13} {'curvature regime':>22} {'tiling analogue':>18}")
    for p in [3,5,7,11,13]:
        g=genus[p]; n=2*p
        regime = "sphere (g=0)" if g==0 else ("PLANE (g=1)" if g==1 else "hyperbolic (g>=2)")
        til = "Platonic" if g==0 else ("hexagonal/A_2" if g==1 else "hyperbolic {p,q}")
        star = "  <== LRC-14, the flat case" if n==14 else ""
        print(f"   {p:>3} {n:>5} {g:>13} {regime:>22} {til:>18}{star}")
    print()
    print("   => LRC-14 (n=14, genus 1) = the EUCLIDEAN/PLANE-TILING case; covering-min = hexagonal A_2")
    print("      (Phi_6 = Eisenstein, HYP-3715/3706). Platonic = genus-0 small-n; hyperbolic = genus>=2 large-n.")

    print()
    print("="*74)
    print("(4) CRYSTALLOGRAPHIC RESTRICTION: 2D rotation orders {1,2,3,4,6}; 5-fold FORBIDDEN")
    print("="*74)
    # allowed n-fold: 2cos(2pi/n) integer <=> n in {1,2,3,4,6}
    import math
    allowed=[k for k in range(1,13) if abs(round(2*math.cos(2*math.pi/k))-2*math.cos(2*math.pi/k))<1e-9]
    print(f"   n-fold with 2cos(2pi/n) in Z (lattice-compatible): {allowed}  -- note 5 and 7+ EXCLUDED")
    print("   the 5 BRAVAIS lattices: oblique(2), rectangular(2), rhombic(2), square(4), HEXAGONAL(6)")
    print("   => the covering-min is 6-fold (hexagonal, Phi_6), NOT 5-fold. Icosahedral 5-fold (golden")
    print("      ratio) is lattice-FORBIDDEN -> lives in the Fibonacci/Zeckendorf/quasicrystal thread.")
