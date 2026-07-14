#!/usr/bin/env python3
"""mac-mini-S92: FIND the global object that breaks arc-transitivity (S91: S_n-transitivity forbids
the 2*7 Atkin-Lehner factorization). Two candidates:
 (A) the TILING model (staircase, base-path fixed): does complement become FIXED-POINT-FREE here
     (resolving the S91 NEG-3 tension where iso-class complement fixes SC)?
 (B) the PALEY/QR tournament on Z_7 (the apex prime 7): does 14=2*7 appear combinatorially (odd
     cycles), and does its broken (Frobenius) symmetry carry the 2*7 structure?"""
from itertools import combinations
from math import comb

print("="*70); print("(B) PALEY / QR TOURNAMENT ON Z_7  (the apex '7' of 14=2*7)")
print("="*70)
p=7; QR={ (x*x)%p for x in range(1,p) }   # {1,2,4}
print(f"Z_{p}: QR = {sorted(QR)} (Paley connection set); nonQR = {sorted(set(range(1,p))-QR)}")
adj=[[1 if (j-i)%p in QR else 0 for j in range(p)] for i in range(p)]
scores=[sum(adj[i]) for i in range(p)]
c3=comb(p,3)-sum(comb(s,2) for s in scores)   # cyclic triangles = odd 3-cycles = OCF atoms
print(f"scores = {scores} (regular tournament, all (p-1)/2=3)")
print(f"DIRECTED 3-CYCLES c3 = {c3}  <-- {'== 14 = 2*7 !!' if c3==14 else c3}   (= the OCF odd-cycle count)")
# self-complementary? complement = multiply differences by a nonQR (e.g. 3): i->j in comp iff j-i in nonQR
def relabel_mult(a):  # vertex map x->a*x mod p is an iso from T to (a*T)
    return [[adj[(pow(a,p-2,p)*i)%p][(pow(a,p-2,p)*j)%p] for j in range(p)] for i in range(p)]
comp=[[adj[j][i] for i in range(p)] for j in range(p)]   # complement (reverse arcs)
# is comp isomorphic to adj? Paley is SC iff -1 is a nonQR, i.e. p=3 mod 4
sc = (p%4==3)
print(f"self-complementary (Paley SC iff p=3 mod4): {sc}  (p=7=3 mod4 => SC)")
# automorphism group: Frobenius Z_p ⋊ Z_{(p-1)/2} order p*(p-1)/2 = 21 (rotations + QR-mult)
print(f"Aut(Paley_7) = Frobenius Z_7 rtimes Z_3, order {p*(p-1)//2} = 21 -- breaks S_7 (order 5040).")
print(f"  => arcs are NOT all S_7-equivalent; Aut has orbits = the QR-difference classes. TRANSITIVITY BROKEN.")
print(f"  The '2' of 14: the Legendre Z_2 = (Z_7)*/QR (QR vs nonQR = complement). The '7': additive Z_7.")
print(f"  So Z_14 = Z_2(Legendre) x Z_7(additive) lives INSIDE the apex-7 Paley structure (CRT).")

print("\n"+"="*70); print("(A) TILING MODEL: is complement FIXED-POINT-FREE (unlike iso-classes)?")
print("="*70)
for n in [4,5,6,7]:
    m=comb(n-1,2)  # tiles = cycle space; base path = cut space (breaks arc-transitivity!)
    # complement = flip ALL m tiles (antipode of Q_m): fixed iff tiling = its flip => impossible for m>=1
    comp_fixed = 0 if m>=1 else 1
    print(f"  n={n}: m={m} tiles; complement=flip-all-tiles => fixed tilings = {comp_fixed} (FIXED-POINT-FREE)")
print("  => On the TILING model complement is fixed-point-FREE (matches the AL regular/fixed-point-free")
print("     action on cusps), UNLIKE iso-classes where it fixes SC (the S91 NEG-3 tension). The base path")
print("     splits arcs into CUT (base path, n-1) + CYCLE (tiles, m) => arc-transitivity BROKEN.")
print("\nLEAD: the apex-7 Paley tournament realizes 14=2*7 as its odd-cycle (OCF) count, with the 2 and 7")
print("both inside Z_7's arithmetic (Legendre Z_2 x additive Z_7 = Z_14). This is arithmetic, not S_n-symmetric.")
