"""opus-2026-07-20-S412 -- CLOSING THE CHARGE-0 CASE OVER C, AND A GAP IN MY OWN THM-1535.

SELF-AUDIT FIRST.  THM-1535 s2 argued: with only nonnegative charges, E[P^2] = E[P0^2] =
c^T H c with H_{ab} = (a+b)! positive definite, so c = 0.  That is valid for REAL c.  But
GMC is stated over C, and for COMPLEX c the form c^T H c (NO conjugate) can vanish with
c != 0 -- e.g. H = I, c = (1, i) gives c^T c = 1 + i^2 = 0.  My exhaustive sweep used
coefficients in {-1,0,1}, i.e. REAL, so it did not test this.  **The complex case was
assumed, not proved.**  This script closes it.

THE RIGHT REDUCTION.  In polar form z = r e^{i theta}, a monomial z^a zbar^b has modulus
r^{a+b} and phase e^{i(a-b)theta} -- so THE CHARGE GRADING IS EXACTLY THE FOURIER GRADING IN
theta, and E is (angular average) then (radial average).  Setting s = r^2, the radial law is
the EXPONENTIAL distribution e^{-s} ds, because E[(z zbar)^k] = k!.

Hence the charge-0 part P0 = sum_a c_a (z zbar)^a is just g(s) = sum_a c_a s^a, and

    E[P0^m]  =  int_0^infty g(s)^m e^{-s} ds .

So the whole charge-0 question becomes ONE-DIMENSIONAL:
    ** for which complex polynomials g is  int_0^infty g(s)^m e^{-s} ds = 0  for all m>=1 ? **

CLAIM: only g = 0.  Reason (asymptotic): if deg g = d >= 1 with leading coefficient c_d,
then as s -> infinity, |g(s)| -> infinity while arg g(s) -> arg c_d.  The PHASE STABILISES,
so for large m there is no oscillation left to cancel the integral: it is dominated by the
tail and behaves like c_d^m (dm)! .  No cancellation => nonzero.  If d = 0 then g = c_0 and
the integral is c_0^m, forcing c_0 = 0.
This script tests the claim numerically and by exact solving.
"""
import sympy as sp
from math import factorial
from itertools import product

s = sp.symbols('s')

def moment(gpoly, m):
    """int_0^infty g(s)^m e^{-s} ds, exactly: int s^k e^{-s} = k!"""
    e = sp.expand(gpoly**m)
    p = sp.Poly(e, s)
    return sp.expand(sum(c*factorial(k) for (k,), c in zip(p.monoms(), p.coeffs())))

print("="*78)
print("(0) THE GAP IS REAL: c^T H c can vanish for COMPLEX c with H positive definite")
print("="*78)
H2 = sp.Matrix([[1, 1], [1, 2]])                    # (a+b)! for a,b in {0,1}
c = sp.Matrix([sp.symbols('c0'), sp.symbols('c1')])
print("   H (a,b in {0,1}) =", H2.tolist())
sol = sp.solve([(c.T*H2*c)[0, 0]], sp.symbols('c0'))
print("   c^T H c = 0 solved for c0 :", sol, "  -> complex solutions exist")
print("   so the S411 argument as written covers REAL coefficients only.")

print()
print("="*78)
print("(1) THE REDUCTION: charge-0 case  <=>  int_0^inf g(s)^m e^{-s} ds = 0 for all m")
print("="*78)
print("   because s = |z|^2 is EXPONENTIAL(1) (E[(z zbar)^k] = k!).")
z, zb = sp.symbols('z zb')
def E1(e):
    e = sp.expand(e)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z, zb); t = 0
    for (a, b), co in zip(p.monoms(), p.coeffs()):
        if a == b: t += co*factorial(a)
    return sp.expand(t)
for g in [1 - s, s**2 - 4*s + 2, sp.I*s - 1]:
    P0 = g.subs(s, z*zb)
    lhs = [E1(sp.expand(P0**m)) for m in (1, 2, 3)]
    rhs = [moment(g, m) for m in (1, 2, 3)]
    print(f"   g={str(g):16s} E[P0^m] m=1..3 = {lhs}   int g^m e^-s = {rhs}   "
          f"{'MATCH' if lhs == rhs else 'mismatch'}")

print()
print("="*78)
print("(2) CAN A NONZERO COMPLEX g KILL ALL MOMENTS?  Exact solve, degree by degree.")
print("="*78)
for d in (1, 2, 3):
    cs = sp.symbols(f'a0:{d+1}')
    g = sum(cs[k]*s**k for k in range(d+1))
    eqs = [moment(g, m) for m in range(1, 2*d + 4)]
    sols = sp.solve(eqs, cs, dict=True)
    nz = [S for S in sols if any(sp.simplify(S.get(x, x)) != 0 for x in cs)]
    print(f"   deg {d}: solving E[g^m]=0 for m=1..{2*d+3}  ->  {len(sols)} solution(s); "
          f"NONZERO ones: {len(nz)}")
    for S in nz[:3]:
        print(f"      *** nonzero solution: {S}")
    if not nz:
        print(f"      only g = 0  =>  charge-0 part must VANISH at degree {d}, over C")

print()
print("="*78)
print("(3) THE ASYMPTOTIC REASON (numerical illustration of phase stabilisation)")
print("="*78)
print("   for g of degree d>=1, arg g(s) -> arg(leading coeff) as s -> inf, so the")
print("   integrand's phase stops oscillating and the tail cannot cancel.")
for gg in [s - 1, sp.I*s - 1, (1 + sp.I)*s**2 - 3*s + 1]:
    vals = [sp.nsimplify(moment(gg, m)) for m in range(1, 7)]
    print(f"   g={str(gg):22s} moments m=1..6: {[sp.simplify(v) for v in vals]}")
    print(f"       -> all zero? {all(sp.simplify(v) == 0 for v in vals)}")
