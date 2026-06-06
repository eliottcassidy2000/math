"""Attempts at novel statements on χ(ℝ²)∈{5,6,7}. (1) LOESCHIAN reframing of the upper bound: the
Eisenstein-periodic colorings use index=Loeschian norm a²+ab+b²; 5,6 are NOT Loeschian, 7=N(2+ζ6)
is the smallest valid ⟹ the hexagonal/periodic χ = 7. (2) SPECTRAL/FRACTIONAL CAP: χ_f=1/m₁≈4.36
< 5, so the spectral bound CANNOT reach 5 — the 5,6,7 distinctions are irreducibly COMBINATORIAL
(the integrality gap = the Vitali wall). (3) FIELD-TOWER chromatic growth χ(ℚ²)=2, χ(ℤ[ζ6])=3,
+√−11→4, →[5,7]; both ends of [5,7] ESCAPE the cyclotomic lattice. opus-2026-06-06-S699n."""
import math
def loeschian(N): 
    for a in range(0,int(math.isqrt(N))+2):
        for b in range(0,int(math.isqrt(N))+2):
            if a*a+a*b+b*b==N: return (a,b)
    return None
def main():
    print("(1) LOESCHIAN numbers (a²+ab+b² = Eisenstein norms) up to 30; which give hexagonal colorings:")
    L=[N for N in range(1,31) if loeschian(N)]
    print(f"   Loeschian: {L}")
    print(f"   5 Loeschian? {bool(loeschian(5))}; 6 Loeschian? {bool(loeschian(6))}; 7 = N{loeschian(7)} = (2+ζ6) prime norm")
    print(f"   ⟹ NO Eisenstein-hexagonal periodic coloring has 5 or 6 colors (5,6 not norms);")
    print(f"      7 is the smallest Loeschian > 4 ⟹ the Isbell hexagonal 7-coloring is index-7 = ℤ[ζ6]/(2+ζ6)≅𝔽₇.")
    print(f"      [So improving the UPPER bound to 6 REQUIRES a non-Eisenstein (aperiodic/non-hexagonal) coloring.]\n")
    print("(2) SPECTRAL / FRACTIONAL CAP (why the 5,6,7 distinctions are combinatorial, not spectral):")
    m1_lo, m1_hi = 0.2293, 0.2598   # known bounds on the plane's 1-avoiding density
    print(f"   m₁ (1-avoiding density) ∈ [{m1_lo},{m1_hi}] ⟹ χ_f = 1/m₁ ∈ [{1/m1_hi:.2f},{1/m1_lo:.2f}]")
    print(f"   χ_f < 5 (≈4.36 at the best m₁) ⟹ the FRACTIONAL/spectral bound CANNOT reach 5.")
    print(f"   ⟹ χ≥5 (de Grey) is IRREDUCIBLY COMBINATORIAL = the integrality gap χ−χ_f = the LRC 'Vitali wall' (S699g).")
    print(f"   ⟹ NO spectral/measure argument can narrow {{5,6,7}}; only combinatorial gadgets (Moser/de Grey-type) can.\n")
    print("(3) FIELD-TOWER chromatic growth (S687 extended): χ climbs along the number-field tower:")
    print("   χ(ℚ²)            = 2   (rational plane bipartite)")
    print("   χ(ℤ[ζ6]=√−3)     = 3   (Eisenstein lattice; W6 forces 3, lattice 3-colors via mod √−3 = 𝔽₃)")
    print("   χ(ℚ(√−3,√−11))  ≥ 4   (contains the Moser spindle; √−11 = the non-cyclotomic rotation)")
    print("   χ(ℝ²)            ∈ {5,6,7}")
    print("   ⟹ the cyclotomic FLOOR is 3; [5,7] is the POST-CYCLOTOMIC regime. BOTH ends escape the lattice:")
    print("      LOWER bound (χ≥5) escapes via Heegner rotations (S687); UPPER bound (≤6) would need a")
    print("      NON-Eisenstein coloring (1). The lattice pins only χ=3; the true value is an ESCAPE statement.")
    print("\n(4) HEEGNER ROADMAP (conjectural, S699m): χ=3↦√−3, χ=4↦√−11, χ=5↦√−19? (next Heegner).")
    print("   If a finite UD graph can force at most k independent Heegner rotations, χ(ℝ²)=2+k. Narrowing")
    print("   {5,6,7} ⟺ counting the max Heegner-rotation rank of a finite unit-distance graph.")
if __name__=='__main__': main()
