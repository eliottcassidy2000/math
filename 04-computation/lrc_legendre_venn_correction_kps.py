"""
The Legendre recursion is a 3-SET INCLUSION-EXCLUSION (owner's correction, kps-S31s).
Sizes: A=B=h(n-1); C=D=h(n-2); E=F=h(n-3); G=h(n-4).  3 sub-tilings A,B (depth n-1) and D (depth n-2);
pairwise intersections A&B=C, A&D=E, B&D=F; triple A&B&D=G.
  corners (the sets)     : A, D, B
  edges  (pairwise unions): A+B-C=|A u B|, A+D-E=|A u D|, B+D-F=|B u D|
  center (full union)     : A+B+D-C-E-F+G = |A u B u D| = h(n)   [ODD]
The EVEN mode A+B-C is exactly the A&B EDGE (the pronic, conjugate pair, no D corner / no triple).
The n-2 appears as corner D(+) AND overlap C(-), CANCELLING in the net coefficient (my earlier 'skips n-2' was WRONG).
"""
def h(n): return ((n-1)**2)//4 if n>=1 else 0
print("n   h(n)  A+B+D-C-E-F+G (full Venn center)  A+B-C (A&B edge)   parity-correct form")
for n in range(5,16):
    A=B=h(n-1); C=D=h(n-2); E=F=h(n-3); G=h(n-4)
    center = A+B+D-C-E-F+G
    abedge = A+B-C
    if n%2==1:  # ODD -> Legendre full Venn
        ok = (center==h(n)); form=f"ODD/Legendre: center={center}={'h(n) OK' if ok else 'X'}"
    else:       # EVEN -> Eisenstein A&B edge
        ok = (abedge==h(n)); form=f"EVEN/Eisenstein: A&B edge={abedge}={'h(n) OK' if ok else 'X'}"
    print(f"{n:2d}  {h(n):3d}   center={center:3d} {'=h' if center==h(n) else '  '}            edge={abedge:3d} {'=h' if abedge==h(n) else '  '}    {form}")
print("\nVerify the 7 Venn PARTITION regions sum to h(n) (odd n=7):")
n=7; A=B=h(n-1); C=D=h(n-2); E=F=h(n-3); G=h(n-4)
# pairwise ints A&B=C, A&D=E, B&D=F, triple=G; D here is the THIRD set (call its size dD=h(n-2))
dD=h(n-2)
regions={"A only":A-C-E+G,"B only":B-C-F+G,"D only":dD-E-F+G,"A&B only":C-G,"A&D only":E-G,"B&D only":F-G,"center A&B&D":G}
print(f"  sizes A={A},B={B},D={dD} (corners); A&B=C={C},A&D=E={E},B&D=F={F} (edge ints); triple G={G}")
for r,v in regions.items(): print(f"   |{r}| = {v}")
print(f"   SUM = {sum(regions.values())} =?= h(7) = {h(7)}  {'PARTITION OK' if sum(regions.values())==h(7) else 'X'}")
print(f"\n  n-2 CANCELLATION: corner D (+h(n-2)) + overlap -C (-h(n-2)) = 0 net => coefficients (2,0,-2,1).")
print(f"  EVEN = degenerate (only A,B conjugate pair + overlap C = the A&B edge; D corner & triple G FOLDED away by complement).")
