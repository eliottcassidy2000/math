#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S66
=======================
MORE crystallographic connections to LRC14 (extending S65's (2,3,n) angle-defect spine).

CENTERPIECE -- THE TRIPLE 6.  The generalized crystallographic restriction: an n-fold rotation
needs an ambient lattice of dimension >= phi(n) (the cyclotomic degree; Z[zeta_n] is the minimal
lattice, and a planar n-fold quasicrystal is a cut-and-project from R^phi(n)).  For LRC14:

    phi(7) = phi(14) = 6

so the 7-fold / 14-fold symmetry of the apex first lives in SIX dimensions (Z[zeta_14]).  And
    |(Z/14)^*| = 6 = phi(14)   with (Z/14)^* = {1,3,5,9,11,13}
is EXACTLY the LRC14 extremal lonely set (the units / touch-points at the cusp {1..13}, HYP-3571/3615).
So the "6" of LRC14 is three things at once:
    phi(14) (cyclotomic degree)  =  quasicrystal embedding dimension  =  #lonely units.
The lonely runners ARE the 6 independent directions of the 14-fold cyclotomic star; the 3 antipodal
unit-pairs (phi(14)/2 = 3) are the "physical" half of the 6D = 3+3 cut-and-project splitting.

OTHER CONNECTIONS:
- 3D crystallography has 7 crystal systems and 14 Bravais lattices (= 2*7 = LRC14); the hexagonal
  system (order 6) is the covering-min side.
- Coxeter/Weyl: crystallographic groups have branch labels m in {2,3,4,6}. Rank-2 crystallographic:
  A1xA1(2), A2(3, hexagonal Weyl S3), B2=C2(4, square), G2(6). The apex is I2(7) (dihedral heptagon),
  the smallest NON-crystallographic dihedral after H2=I2(5); LRC covering-min = A2/G2 (crystallographic),
  apex = I2(7) (non-crystallographic) = the 4cos^2(3pi/7) floor gap.
"""
from math import gcd
from collections import defaultdict

def phi(n): return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)

print("="*74)
print("(1) GENERALIZED CRYSTALLOGRAPHIC RESTRICTION: n-fold needs dim >= phi(n)")
print("="*74)
byd = defaultdict(list)
for n in range(1, 25): byd[phi(n)].append(n)
tag = {1: "trivial", 2: "2D CRYSTAL (5 Bravais lattices)", 4: "4D quasicrystal (5,8,10,12-fold: Penrose/Ammann/decagonal/dodecagonal)",
       6: "6D quasicrystal (7,9,14,18-fold)", 8: "8D", 10: "10D", 12: "12D"}
for d in sorted(byd):
    print(f"  dim phi={d:2d}: n in {byd[d]:}   {tag.get(d,str(d)+'D')}")

print()
print("="*74)
print("(2) THE TRIPLE 6 of LRC14")
print("="*74)
units = [k for k in range(1, 14) if gcd(k, 14) == 1]
print(f"  phi(7) = {phi(7)},  phi(14) = {phi(14)}   -> 7-fold & 14-fold symmetry live in dim 6 (Z[zeta_14])")
print(f"  (Z/14)^* = {units}   |.| = {len(units)} = phi(14) = the LRC14 extremal LONELY SET (units; HYP-3571/3615)")
print(f"  antipodal pairs = phi(14)/2 = {phi(14)//2}  (the 'physical' half of the 6D = 3+3 cut-and-project)")
print("  => 6 = cyclotomic degree = quasicrystal embed dim = #lonely units.  The runners ARE the cyclotomic star.")

print()
print("="*74)
print("(3) CRYSTAL SYSTEMS / BRAVAIS LATTICES by dimension (the 7 & 14 resonance)")
print("="*74)
tbl = [("1D", 1, 1), ("2D", 4, 5), ("3D", 7, 14), ("4D", 33, 64), ("5D", 59, 189), ("6D", 251, 841)]
for dim, systems, bravais in tbl:
    flag = "   <== 7 systems, 14 Bravais = 2*7 = LRC14" if dim == "3D" else ("   <== 5 lattices (S65)" if dim == "2D" else "")
    print(f"  {dim}: {systems:3d} crystal systems, {bravais:3d} Bravais lattices{flag}")

print()
print("="*74)
print("(4) COXETER / WEYL: crystallographic <=> labels m in {2,3,4,6}; apex = I2(7) non-crystallographic")
print("="*74)
print("  rank-2 crystallographic (Weyl/lattice-preserving): A1xA1 m=2 | A2 m=3 (hexagonal, |W|=6) |")
print("                                                     B2=C2 m=4 (square, |W|=8) | G2 m=6 (|W|=12)")
print("  NON-crystallographic dihedral: H2=I2(5) pentagon | I2(7) HEPTAGON = the apex | H3,H4 (5-fold)")
def cos_deg(m):  # degree of 2cos(pi/m) over Q (label m); crystallographic iff this is 1
    return phi(2*m)//2 if m > 2 else 1
for m in [2, 3, 4, 5, 6, 7]:
    crys = "crystallographic" if m in (2, 3, 4, 6) else "NON-crys (quasicrystal)"
    nm = {2: "A1xA1", 3: "A2", 4: "B2/C2", 5: "H2=I2(5)", 6: "G2", 7: "I2(7)=apex"}[m]
    print(f"    I2({m}) = {nm:11}: deg 2cos(pi/{m}) = {cos_deg(m)}   {crys}")
print("  => LRC covering-min = order-6 A2/G2 (crystallographic); apex = I2(7) (non-crys) = 4cos^2(3pi/7) floor gap.")
print("\nDONE.")
