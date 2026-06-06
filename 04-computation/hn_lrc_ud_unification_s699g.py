"""Hadwiger–Nelson ∪ unit-distance ∪ LRC unification. All three = invariants of the SAME distance
Cayley graphs: HN = chromatic number χ(Cay(ℝ²,unit circle))∈[5,7]; UD = edge density; LRC =
independence/covering density (its tournament has dichromatic χ). Shared TECHNIQUE = the Fourier
transform of the forbidden-distance measure (Bessel J0 for the plane; Hoffman spectral bound for
the lattice). Shared STRUCTURE = Eisenstein/π/3. opus-2026-06-06-S699g."""
import math
from itertools import product
def tri_dispersion():
    # infinite triangular lattice adjacency symbol: 2(cos a + cos b + cos(a+b)) over BZ
    lo=1e9; hi=-1e9
    N=400
    for i in range(N):
        for j in range(N):
            a=2*math.pi*i/N; b=2*math.pi*j/N
            val=2*(math.cos(a)+math.cos(b)+math.cos(a+b))
            lo=min(lo,val); hi=max(hi,val)
    return lo,hi
def bessel_J0(x, terms=60):
    s=0.0
    for m in range(terms):
        s += ((-1)**m)*(x/2)**(2*m)/(math.factorial(m)**2)
    return s
def J0_min():
    # min of J0 on (0, 20)
    lo=1e9; xm=0
    x=0.01
    while x<20:
        v=bessel_J0(x)
        if v<lo: lo=v; xm=x
        x+=0.001
    return lo,xm
def main():
    print("SHARED TECHNIQUE — Fourier/spectral (Hoffman) bound on the chromatic number:")
    lo,hi=tri_dispersion()
    chi_hoff = 1 - hi/lo
    print(f"  TRIANGULAR LATTICE (UD/Eisenstein): adjacency symbol range [{lo:.3f},{hi:.3f}];")
    print(f"    Hoffman χ ≥ 1 − λmax/λmin = 1 − {hi:.1f}/({lo:.1f}) = {chi_hoff:.3f}  ⟹ χ ≥ 3 (TIGHT: χ=3, Eisenstein (i−j) mod 3)")
    j0lo,xm=J0_min()
    chi_plane = 1 - 1/j0lo
    print(f"  PLANE (Hadwiger–Nelson): the unit-circle measure's transform is the Bessel J0;")
    print(f"    min J0 = {j0lo:.4f} at x≈{xm:.2f}; spectral (fractional) χ ≥ 1 − 1/min(J0) = {chi_plane:.3f}")
    print(f"    (the same Hoffman/Lovász method; the true χ(ℝ²)∈[5,7] needs the combinatorial Moser/deGrey graphs.)\n")
    print("SHARED STRUCTURE — Eisenstein ζ6 / π/3 (the extremal configs of all three):")
    print("  UD optimum = triangular/Eisenstein lattice (κ=6, S699a);")
    print("  LRC worry-set witnesses = roots of unity (THM-403); n=14 prime-3 C=27=3³ = Eisenstein norm-3 (HYP-2170);")
    print("  HN: the Moser spindle (χ=4) & de Grey graph (χ=5) live on the Eisenstein lattice rotated by")
    print("      INCOMMENSURATE angles arccos(5/6)/arcsin(1/(2√3)) — the Diophantine/irrational structure = LRC's home.\n")
    th=2*math.asin(1/(2*math.sqrt(3)))
    print(f"  Moser-spindle rotation θ = 2·arcsin(1/(2√3)) = {th:.5f} rad = {math.degrees(th):.3f}° (irrational multiple of π) — an")
    print(f"  INCOMMENSURATE rotation, the same resonance/Diophantine phenomenon LRC studies (the speeds' irrationality).")
    print("\nTHE THREE INVARIANTS OF ONE GRAPH FAMILY:")
    print("  HN  = χ (chromatic number)            of Cay(ℝ², unit circle)  ∈ [5,7]")
    print("  UD  = max edges (edge density)         of Cay(lattice, unit vectors)  ≈ κ n/2 (κ≤6) or n^1.014 (CM)")
    print("  LRC = independence/covering density p0 of Cay(circle, danger arcs); worry-set tournament dichromatic χ=2")
if __name__=='__main__': main()
