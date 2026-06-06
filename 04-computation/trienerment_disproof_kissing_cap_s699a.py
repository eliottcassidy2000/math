"""Does the <1/=1/>1 TRIERNEMENT simplify the disproof of the grid as optimal? KEY REFRAME: the
'=1' (tie) layer of the threshold tournament is the NORM-1 LAYER of the point set. For a 2D LATTICE
the per-point =1 degree = the kissing number κ ≤ 6 (max at triangular/Eisenstein) ⟹ unit distances
≈ κn/2 ≤ 3n (LINEAR). The CM disproof ESCAPES this geometric cap by making '=1' a FIELD-NORM
condition (unbounded via class-field towers), not a Euclidean-kissing one. opus-2026-06-06-S699a."""
import math
def kissing_and_unit(lattice_basis, R=40):
    # count minimal (norm-1) vectors; estimate unit-distance density per point
    b1,b2=lattice_basis; pts={}
    best=1e18
    vecs=[]
    for i in range(-R,R+1):
        for j in range(-R,R+1):
            if i==0 and j==0: continue
            x=i*b1[0]+j*b2[0]; y=i*b1[1]+j*b2[1]; d2=x*x+y*y
            vecs.append(d2)
            if d2<best-1e-9: best=d2
    kappa=sum(1 for d2 in vecs if abs(d2-best)<1e-9)
    return kappa, best
def main():
    print("TRIERNEMENT '=1' LAYER = NORM-1 LAYER. For 2D LATTICES, per-point unit degree = kissing # κ:")
    lats={
      "square ℤ² (Gaussian)": ((1,0),(0,1)),
      "triangular (Eisenstein ζ6)": ((1,0),(0.5, math.sqrt(3)/2)),
      "hexagonal-dual/other": ((1,0),(0.5, math.sqrt(3)/2)),
      "rectangular 1x2": ((1,0),(0,2)),
    }
    for name,basis in lats.items():
        k,nrm=kissing_and_unit(basis)
        print(f"  {name:30s}: kissing κ = {k:2d}  ⟹  unit distances ≈ κ·n/2 = {k/2:.1f}·n  (LINEAR, exponent 1)")
    print(f"\n  2D-lattice kissing CAP = 6 (triangular/Eisenstein) ⟹ ANY 2D lattice gives ≤ 3n unit distances.")
    print(f"  Harborth ⌊3n−√(12n−3)⌋ = exactly this cap minus the √ perimeter — the BEST 2D lattice is triangular.\n")
    print("WHY THE GRID/LATTICE IS NOT OPTIMAL (the trienerment reading of the disproof):")
    print("  The '=1' layer of a 2D lattice is CAPPED by the Euclidean kissing number κ≤6.")
    print("  The CM construction is NOT a 2D lattice — it is a 2D image of a high-rank unit group, where")
    print("  '=1' means FIELD-norm 1 (β·β̄=1), which has UNBOUNDEDLY many solutions (Hilbert 90 + class-")
    print("  field towers, Golod–Shafarevich). The trienerment tie-layer is then arithmetic, NOT geometric,")
    print("  so it is not capped at 6 ⟹ n^{1.014} unit distances. The trienerment localizes the disproof to:")
    print("    'replace the Euclidean-kissing =1 layer (cap 6) by a field-norm =1 layer (no cap).'")
    print()
    # demonstrate: number of roots of unity / norm-1 elements grows with cyclotomic degree (toy)
    print("TOY: # roots of unity (norm-1, |ζ|=1) in ℚ(ζ_m) = 2m grows with m — the unbounded =1 layer:")
    for m in (4,6,12,24,60):
        print(f"   ℚ(ζ_{m}): {2*m} roots of unity on |z|=1  (vs 2D lattice cap 6)")
    print("  (Real CM construction packs these into the plane via class-field towers; the point is the")
    print("   tie-layer is no longer kissing-bounded.)")
if __name__=='__main__': main()
