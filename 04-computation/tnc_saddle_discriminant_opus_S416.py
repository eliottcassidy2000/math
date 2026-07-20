"""opus-2026-07-20-S416 -- A SADDLE/DISCRIMINANT NECESSARY CONDITION FOR TNC.

STATE (not duplicating): THM-1530(B) proves TNC when min exponent = -1 via Lagrange-Burmann;
boxeph THM-1595 proves N=1 (all M), (2,2), (2,3) and closes (2,4), (3,3) by elimination;
klein THM-1550 reduces TNC to Pi(t) = ct.  Open: general M, N >= 2.

THE REFORMULATION.  Write Lambda = u^{-N} R(u), deg R = M+N, R(0) != 0.  Then
      CT(Lambda^m) = [u^{Nm}] R(u)^m ,
a DIAGONAL OF A POWER.  The rising-factorial content is already split off: THM-1530(A) has
E[P^m] ~ Gamma(Dm+1) * CT(Lambda^m) -- the RADIAL part contributes the factorial (Gamma), the
ANGULAR part contributes the toral constant term.  So TNC is exactly the factorial-free core.

THE SADDLE CONDITION.  By Cauchy,
      [u^{Nm}] R^m = (1/2 pi i) * contour_integral R(u)^m / u^{Nm+1} du ,
whose exponential rate is governed by the saddle of  log R(u) - N log u , i.e. by
      u R'(u) / R(u)  =  N .                                   (SADDLE EQUATION)
If u* is a nondegenerate saddle with R(u*) != 0, the coefficient grows like
(R(u*)/u*^N)^m / sqrt(2 pi m sigma) -- an m-th power times a nonzero algebraic prefactor,
which CANNOT vanish for all m.  So a TNC violator needs the saddle value to VANISH:
      R(u*) = 0  with  u* R'(u*) = N R(u*) = 0 ,  so (u* != 0)  R'(u*) = 0 .
      => u* is a MULTIPLE ROOT of R  =>  disc(R) = 0 .

**NECESSARY CONDITION (to test): any Lambda with CT(Lambda^m) = 0 for all m, other than a
single monomial, must have R with a REPEATED ROOT, i.e. disc(R) = 0.**

That is a codimension-1 algebraic condition on the (M+N+1) coefficients -- it would cut the
TNC search space from a full coefficient space to a hypersurface, at EVERY bidegree at once,
which is exactly the kind of reduction the bidegree-by-bidegree ladder lacks.
"""
import sympy as sp
from itertools import product

u, t = sp.symbols('u t')

def CT(Rpoly, N, m):
    """CT((u^-N R)^m) = [u^{Nm}] R^m"""
    e = sp.expand(Rpoly**m)
    return sp.Poly(e, u).coeff_monomial(u**(N*m)) if e != 0 else 0

print("="*78)
print("(1) SANITY: the saddle prediction on cases already settled")
print("="*78)
print("   Lambda = u^-N R(u).  Nonzero CT(Lambda^m) for some m  <=>  TNC holds there.")
print("   Prediction: disc(R) != 0  ==>  some CT(Lambda^m) != 0.")
tests = [
    ("N=1, R=1+u",            1, 1 + u),
    ("N=1, R=(1+u)^2",        1, (1 + u)**2),
    ("N=2, R=1+u^2",          2, 1 + u**2),
    ("N=2, R=(1+u)^4",        2, (1 + u)**4),
    ("N=2, R=1+u+u^2+u^3+u^4",2, 1 + u + u**2 + u**3 + u**4),
    ("N=3, R=(1+u)^6",        3, (1 + u)**6),
    ("N=2, R=(1-u)^2(1+u)^2", 2, (1 - u)**2*(1 + u)**2),
]
for lab, N, R in tests:
    Rp = sp.Poly(sp.expand(R), u)
    d = sp.discriminant(Rp)
    cts = [CT(sp.expand(R), N, m) for m in range(1, 7)]
    firstnz = next((i+1 for i, v in enumerate(cts) if v != 0), None)
    print(f"   {lab:26s} disc={str(d):>8s}  CT m=1..6 = {cts}  first nonzero m={firstnz}")

print()
print("="*78)
print("(2) THE SADDLE EQUATION u R'(u) = N R(u), and whether its roots hit R")
print("="*78)
for lab, N, R in tests:
    Rx = sp.expand(R)
    sad = sp.expand(u*sp.diff(Rx, u) - N*Rx)
    common = sp.gcd(sp.Poly(Rx, u), sp.Poly(sad, u))
    print(f"   {lab:26s} gcd(R, uR'-NR) = {sp.factor(common.as_expr())}")
print("   nontrivial gcd  <=>  the saddle sits on a root of R  <=>  candidate TNC violator")

print()
print("="*78)
print("(3) SEARCH: any Laurent with min<=-2, max>=2 and CT(Lambda^m)=0 for m=1..8?")
print("="*78)
print("   (THM-1530(C) tested 3456 three-term cases, zero hits.  Widen to FOUR terms,")
print("    and record whether every survivor has disc(R) = 0 as the saddle argument predicts.)")
hits = 0; scanned = 0; discs = []
N = 2
for c in product([-2,-1,0,1,2], repeat=4):
    if c[0] == 0 or c[-1] == 0: continue          # keep deg and the u^0 coefficient honest
    R = c[0] + c[1]*u + c[2]*u**2 + c[3]*u**3 + u**4
    scanned += 1
    cts = [CT(sp.expand(R), N, m) for m in range(1, 6)]
    if all(v == 0 for v in cts):
        hits += 1
        discs.append(sp.discriminant(sp.Poly(sp.expand(R), u)))
        if hits <= 5: print(f"      HIT: R={R}  disc={discs[-1]}")
print(f"   scanned {scanned} quartics at N=2; CT-vanishing through m=5: {hits}")
if hits == 0:
    print("   none -- consistent with TNC, and the saddle condition is not yet stressed.")
else:
    print(f"   all hits have disc = 0? {all(d == 0 for d in discs)}")
