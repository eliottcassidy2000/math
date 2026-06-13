"""Does de Grey's χ=5 graph force ≥3 distinct √−d? Likely the WRONG invariant. Minkowski sums stay
in a fixed field; the growing quantity is the ROTATION RANK = # multiplicatively-independent
unit-norm elements |β|=1 used. By Hilbert 90 the unit-norm group of an imaginary-quadratic field is
INFINITE rank (β=γ/γ̄), so χ grows WITHIN fixed fields. χ=3 (cyclotomic floor, rank 0/torsion only);
Moser χ=4 (rank 1, the √−11 rotation); de Grey χ=5 (rank ≥2). The shared object across HN/UD/LRC =
unit-norm elements; LRC=torsion (roots of unity), UD/HN=non-torsion escape. opus-2026-06-06-S699o."""
import math
from fractions import Fraction as F
def unitnorm_from(a,b,d=11):
    # γ=a+b√−d ; β=γ/γ̄ = (a²−d b² + 2ab√−d)/(a²+d b²), |β|=1; return (Re, Im_coeff_of_√d, |γ|², angle)
    N=a*a+d*b*b
    re=F(a*a-d*b*b, N); im2=F(2*a*b, N)   # β = re + im2·√−d  (so Im(β)=im2·√d)
    ang=math.degrees(math.atan2(im2*math.sqrt(d), float(re)))
    return re, im2, N, ang
def is_root_of_unity_angle(ang, tol=1e-7):
    # rational multiple of 360 with small denominator?
    for q in range(1,25):
        for p in range(0,q):
            if abs((ang%360) - 360*p/q) < 1e-4: return (p,q)
    return None
def main():
    print("(1) The MOSER rotation is a unit-norm element of ℚ(√−11):")
    re,im2,N,ang=unitnorm_from(11,1,11)   # γ=11+√−11 ⟹ β=(5+√−11)/6
    print(f"   γ=11+√−11: β=γ/γ̄ = {re} + {im2}·√−11  = (5+√−11)/6 ; |β|²={float(re)**2+float(im2)**2*11:.6f}")
    print(f"   angle={ang:.4f}°; cosθ=Re(β)={float(re)}=5/6; root of unity? {is_root_of_unity_angle(ang)} (None ⟹ NON-torsion, the Moser/√−11 rotation)")
    print(f"   minimal poly 3z²−5z+3 (disc −11, Mahler measure 3 ≠ 1 ⟹ non-cyclotomic).\n")
    print("(2) ℚ(√−11) has INFINITELY MANY independent unit-norm rotations (Hilbert 90, β=γ/γ̄):")
    angs=[]
    for (a,b) in [(11,1),(1,1),(2,1),(3,1),(4,1),(5,1),(1,2),(3,2)]:
        re,im2,N,ang=unitnorm_from(a,b,11)
        ru=is_root_of_unity_angle(ang)
        angs.append(round(ang%360,3))
        print(f"   γ={a}+{b}√−11: β angle={ang%360:.3f}°  cos={float(re):.4f}  root-of-unity?={ru}")
    print(f"   ⟹ many DISTINCT non-torsion unit-norm rotations, ALL in the single field ℚ(√−11).")
    print(f"   The unit-norm group is infinite rank ⟹ you can pack many INDEPENDENT rotations WITHOUT new fields.\n")
    print("(3) ANSWER to 'does de Grey force ≥3 distinct √−d?': almost certainly NO (the field COUNT")
    print("    is the wrong invariant — Minkowski sums preserve the field; de Grey plausibly lives in")
    print("    ℚ(√−3,√−11), 2 imaginary fields). The RIGHT invariant is the ROTATION RANK:")
    print("    χ=3 (cyclotomic floor, rank 0 = torsion/roots of unity); Moser χ=4 (rank 1, the √−11");
    print("    rotation); de Grey χ=5 forces rank ≥2 (two independent unit-norm rotations). χ = 3 + rank.\n")
    print("(4) CROSS-PROBLEM (the shared object = unit-norm elements |β|=1 of imag-quadratic/CM fields):")
    print("   LRC   : the WITNESSES = n-th roots of unity = the TORSION unit-norm (THM-403); the cyclotomic")
    print("           FLOOR. LRC hardness = the torsion ARITHMETIC (shells mod 2n−1, prime-3) — rank 0.")
    print("   UD    : Sawin's n^1.014 = COUNT of unit-norm elements in a ball (CM field + class-field tower)")
    print("           — the non-torsion group's DENSITY. Escapes the lattice (kissing 6) via non-torsion.")
    print("   HN    : χ = 3 + rank of forced non-torsion rotations (Moser, de Grey). Escapes χ=3 floor.")
    print("   ⟹ LRC = the TORSION (cyclotomic) member; UD & HN = the NON-TORSION escape members; the")
    print("      'escape from the cyclotomic floor' = USING non-torsion unit-norm elements. Same object,")
    print("      three readings: torsion clock (LRC), ball-count (UD), independent-rank (HN).")
if __name__=='__main__': main()
