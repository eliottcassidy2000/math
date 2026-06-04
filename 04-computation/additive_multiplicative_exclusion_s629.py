from fractions import Fraction as F
# The Vieta/symmetric-function exclusion: (sum, product) = (e1, e2) of a pair; both rational <=> the monic quadratic
# is rational <=> the pair is algebraic. This is the additive (sum) vs multiplicative (product) complementarity.
print("=== the additive-multiplicative exclusion (Vieta) and its parallels across the arc ===")
rows=[
("pi, e (transcendence)", "sum pi+e (ADDITIVE) vs product pi*e (MULTIPLICATIVE): NOT both rational (else pi,e algebraic)"),
("LRC two faces (HYP-2150)", "additive resonance Sigma m_i v_i=0 (mod 2n-1) vs multiplicative doubling x->2x (mod n): can't both be trivial"),
("Collatz (HYP-2175)", "additive (+1 in 3n+1) vs multiplicative (x3 vs 2-adic) ; the 2^K=3^L resonance"),
("independence polynomial", "alpha_1 (#3-cycles, the SUM/e1) vs alpha_2 (disjoint pairs, the PRODUCT/e2); norm-1 = alpha_2=1 = product collapses"),
("CM / norm-1 (HYP-2235)", "roots rho_1+rho_2 (sum) vs rho_1*rho_2=1 (product=|beta|^2): the multiplicative face fixed by conjugation"),
("rational vs irrational", "~ even vs odd ~ the 2-adic (H odd, Redei) ~ addition vs multiplication"),
("cyclotomic 7=Phi_3(2)", "roots of unity: sum=-1, product=1 BOTH rational (algebraic) -- the ALIGNED/algebraic extreme, opposite of pi,e"),
]
for a,b in rows: print(f"  {a:28}: {b}")
print()
# The cyclotomic extreme: for Phi_3 = X^2+X+1, roots w,w^2 (primitive cube roots); sum = -1 (rational), product = +1 (rational).
# Both symmetric functions rational => algebraic (roots of unity). The OPPOSITE of pi,e where at least one is irrational.
import cmath
w=cmath.exp(2j*cmath.pi/3)
print(f"cyclotomic Phi_3 roots w, w^2: sum = {(w+w*w).real:+.3f} (rational -1), product = {(w*w*w).real:+.3f} (rational +1) -> BOTH rational = the algebraic/aligned extreme")
print("pi,e: at least one of {sum, product} irrational = the transcendental/anti-aligned extreme.")
print("=> a SPECTRUM from cyclotomic (both rational = forbidden/resonant/collapse) to pi,e (one irrational = generic/lonely).")
