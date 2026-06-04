"""The shared object behind the 1.014 exponent: the angle π/3 (Eisenstein ζ6 = triangular UD
lattice = tournament Lee-Yang zeros at ±2π/3) and the constant Cl2(π/3)=1.0149. Compute it;
verify 7=Φ3(2),21=3Φ3(2); the SC shape parameter s=α1/√α2 growth; the tropical constant
κ_trop=Cl2(π/3)/(π logφ). CORRECTS s599s (forbidden H={7,21} only; u22∈[60,61] not 49).
opus-2026-06-03-S599u."""
import math
def Cl2(theta, N=400000):
    # Cl2(θ) = Σ_{k≥1} sin(kθ)/k²
    return sum(math.sin(k*theta)/(k*k) for k in range(1,N+1))
def main():
    th=math.pi/3
    cl=Cl2(th)
    print(f"Cl2(π/3) = {cl:.7f}   (max of the Clausen function; vol of ideal regular tetrahedron / 3... )")
    print(f"  Sawin/OpenAI UD disproof exponent (stated) = 1.014;  Cl2(π/3) = {cl:.6f};  diff = {abs(1.014-cl):.5f}")
    print(f"  UD EXCESS exponent δ = 0.014  vs  Cl2(π/3)-1 = {cl-1:.6f}")
    phi=(1+math.sqrt(5))/2
    kappa=cl/(math.pi*math.log(phi))
    print(f"  tropical dominance κ_trop = Cl2(π/3)/(π·logφ) = {kappa:.6f}   (HYP-707)")
    print()
    print("CYCLOTOMIC identity of the forbidden H-values (prime-3 angle):")
    Phi3=lambda x: x*x+x+1
    print(f"  Φ3(x)=x²+x+1;  Φ3(2)={Phi3(2)}  => 7 = Φ3(2);  21 = 3·Φ3(2) = {3*Phi3(2)}")
    print(f"  roots of Φ3 = primitive cube roots of unity = e^{{±2πi/3}} = the Eisenstein ζ3 (angle 2π/3)")
    print(f"  => the forbidden tournament H-values live at the SAME angle 2π/3 as the UD Eisenstein lattice ζ6 (π/3)")
    print()
    print("SHAPE PARAMETER s=α1/√α2 of SC H-maximizers (α2=1 'norm-1/CM' family):")
    # known SC H-maximizers (from repo): n=5 H=9, n=6 H=45 (α1=20,α2=1)
    data=[(5,9,4,1),(6,45,20,1)]  # (n, H, α1, α2) — α1=#3-cycles-ish
    for (n,H,a1,a2) in data:
        s=a1/math.sqrt(a2)
        print(f"  n={n}: H={H}, α1={a1}, α2={a2}; shape s=α1/√α2={s}; roots product ρ1ρ2=1/α2={1/a2} (CM norm-1)")
    print("  leading shape exponent = 3 (s~α1~n³/24, the 3-cycle count); the UD-analog excess 0.014")
    print("  is a SECONDARY/tropical correction (Cl2(π/3) in the log-sin limit), NOT the leading term.")
    print()
    print("HONEST corrections to S599s/p-r:")
    print("  * forbidden H = {7,21} ONLY (sporadic, prime-7=Φ3(2)); 63,189 ACHIEVABLE at n=8 (my 7·3^k was a sampling artifact)")
    print("  * u(22) ∈ [60,61] (Alexeev–Tikhonov); my triangular-lattice 49 = Harborth LATTICE optimum, NOT the true max")
    print("  * 1.014 ≈ Cl2(π/3) is numerical; Sawin's is a LOWER bound (>1.014), so equality UNVERIFIED")
if __name__=='__main__': main()
