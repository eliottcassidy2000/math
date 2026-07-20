"""opus-2026-07-20-S405 -- THE NAMED NEXT STEP FROM THM-1375:
"do the two reflections coincide? compute the monodromy of a loop around the Jelonek
quartic and compare with the lambda = -1 action -- if they agree the counterexample is
'reflection = torus'."

Setting (all from canon, re-verified here):
  F is the owner's dim-3 Keller counterexample, det JF = -2, degree 3, etale everywhere.
  Its Jelonek set (asymptotic variety, where sheets escape to infinity) is Zariski's
  three-cuspidal quartic -- zero nodes, the worst rational quartic.
  The C*-action has weights (1,-1,-2) -> (-2,-1,1); at lambda = -1 this is exactly the
  sigma/tau involution of THM-1350: sigma = (-x,-y,z), tau = (a,-b,-c).
  Over a tau-fixed target (a,0,0) the fibre is 1 sigma-fixed point + 1 free sigma-orbit,
  so sigma acts on the 3-element fibre as a TRANSPOSITION.

PLAN.  Reduce the fibre to the roots of a UNIVARIATE CUBIC with coefficients polynomial
in the target (a,b,c).  Then three things fall out of one computation:
  (i)  leading coefficient = 0  <=>  a root escapes to infinity  =  the Jelonek set;
  (ii) discriminant of the cubic = the S_3 vs A_3 character (the repo's "Redei sign =
       discriminant character" reading, and Campbell's non-square criterion at d=3);
  (iii) monodromy = permutation of the three roots along a loop, trackable numerically.
"""
import sympy as sp

x, y, z, a, b, c, t = sp.symbols('x y z a b c t')

u = 1 + x*y
F1 = sp.expand(u**3*z + y**2*u*(4 + 3*x*y))
F2 = sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y))
F3 = sp.expand(2*x - 3*x**2*y - x**3*z)

print("="*76)
print("(0) SANITY: det JF and the sigma/tau equivariance at lambda = -1")
print("="*76)
J = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in (F1, F2, F3)])
print("  det JF =", sp.simplify(sp.expand(J.det())))
sub = {x: -x, y: -y}                      # sigma = (-x,-y,z)
S1 = sp.expand(F1.subs(sub, simultaneous=True))
S2 = sp.expand(F2.subs(sub, simultaneous=True))
S3 = sp.expand(F3.subs(sub, simultaneous=True))
print("  F(sigma p) - tau(F p) = ",
      [sp.simplify(S1 - F1), sp.simplify(S2 + F2), sp.simplify(S3 + F3)],
      " (tau = (a,-b,-c))")
print("  F(0,0,z) =", (F1.subs({x: 0, y: 0}), F2.subs({x: 0, y: 0}), F3.subs({x: 0, y: 0})))

print()
print("="*76)
print("(1) THE FIBRE CUBIC: eliminate to one variable")
print("="*76)
# F3 = c is linear in z when x != 0:  z = (2x - 3x^2 y - c)/x^3
zsol = sp.together((2*x - 3*x**2*y - c)/x**3)
G1 = sp.simplify(sp.expand(sp.numer(sp.together(F1.subs(z, zsol) - a))))
G2 = sp.simplify(sp.expand(sp.numer(sp.together(F2.subs(z, zsol) - b))))
print("  after substituting z, two equations in (x,y):")
print("    deg_x G1 =", sp.degree(G1, x), " deg_y G1 =", sp.degree(G1, y))
print("    deg_x G2 =", sp.degree(G2, x), " deg_y G2 =", sp.degree(G2, y))
R = sp.resultant(sp.Poly(G1, y), sp.Poly(G2, y))
R = sp.factor(sp.expand(R))
print("  resultant eliminating y: factored form has", len(sp.Mul.make_args(R)), "factors")
Rp = sp.Poly(sp.expand(R), x)
print("  total degree in x:", Rp.degree())

# strip the spurious x^k factor introduced by clearing x^3 denominators
core = sp.expand(sp.cancel(R / x**sp.Poly(R, x).monoms()[-1][0])) if R.has(x) else R
core = sp.factor(core)
print()
print("  factors and their x-degrees:")
for f in sp.Mul.make_args(core):
    base = f.base if f.is_Pow else f
    if base.has(x):
        print(f"    deg_x = {int(sp.degree(base, x)):2d}   {str(base)[:170]}")

print()
print("="*76)
print("(2) WHICH FACTOR IS THE FIBRE? -- test on a known fibre")
print("="*76)
print("  canon: F^-1(1,0,0) = {(0,0,1)} u {(+-i/2, +-3i, -26)}  (3 points, from THM-1300)")
for pt in [(0, 0, 1), (sp.I/2, 3*sp.I, -26), (-sp.I/2, -3*sp.I, -26)]:
    val = (sp.simplify(F1.subs({x: pt[0], y: pt[1], z: pt[2]})),
           sp.simplify(F2.subs({x: pt[0], y: pt[1], z: pt[2]})),
           sp.simplify(F3.subs({x: pt[0], y: pt[1], z: pt[2]})))
    print(f"    F{pt} = {val}")

print()
print("="*76)
print("(3) THE FIBRE CUBIC IS DEPRESSED, AND L IS THE JELONEK POLYNOMIAL")
print("="*76)
L  = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
P  = L*x**3 + (4 - 3*b*c)*x - 2*c
print("  P(x) = L*x^3 + (4-3bc)*x - 2c      (NO x^2 term: depressed)")
print("  L    =", L)
print("  cross-check vs canon ('-L at g=0 = b^2-16a'):  -L|_{c=0} =",
      sp.expand(-L.subs(c, 0)))
print("  tau-invariance of L:  L(a,-b,-c) - L(a,b,c) =",
      sp.simplify(sp.expand(L.subs({b: -b, c: -c}, simultaneous=True) - L)))
print("  tau acts on P by:  P_tau(x) + P(-x) =",
      sp.simplify(sp.expand(P.subs({b: -b, c: -c}, simultaneous=True) + P.subs(x, -x))))
print("  => sigma acts on the fibre EXACTLY as x -> -x on the roots of P.")

print()
print("  fibre over a tau-fixed target (b=c=0):  P = 16a x^3 + 4x = 4x(4a x^2 + 1)")
print("    roots: x = 0  and  x = +-1/(2 sqrt(-a)).")
print("    at a=1: x = 0, +-i/2  -- matches canon F^-1(1,0,0) = {(0,0,1)},{(+-i/2,...)}")

print()
print("="*76)
print("(4) MONODROMY OF A MERIDIAN AROUND THE JELONEK QUARTIC")
print("="*76)
print("  The a-axis {b=c=0} lies in the tau-FIXED locus, and meets {L=0} at a=0.")
print("  dL/da at the origin =", sp.diff(L, a).subs({a: 0, b: 0, c: 0}),
      "!= 0  => a=0 is a SMOOTH point of the quartic and the a-axis is TRANSVERSE.")
print("  So a = eps*exp(i theta) is a genuine meridian loop, and it stays inside the")
print("  tau-fixed locus -- so monodromy and sigma act on the SAME three roots.")
import numpy as np
def roots_at(av):
    return np.sort_complex(np.roots([16*av, 0.0, 4.0, 0.0]))
eps = 0.30
N = 4000
th = np.linspace(0, 2*np.pi, N)
cur = np.roots([16*(eps+0j), 0.0, 4.0, 0.0])
start = cur.copy()
for k in range(1, N):
    av = eps*np.exp(1j*th[k])
    nxt = np.roots([16*av, 0.0, 4.0, 0.0])
    matched = []
    pool = list(nxt)
    for r in cur:                      # continuity tracking: nearest unused root
        j = min(range(len(pool)), key=lambda i: abs(pool[i]-r))
        matched.append(pool.pop(j))
    cur = np.array(matched)
perm = []
for r in cur:
    j = min(range(3), key=lambda i: abs(start[i]-r))
    perm.append(j)
print(f"\n  start roots (a={eps}):      {np.round(start,6)}")
print(f"  after one full loop:       {np.round(cur,6)}")
print(f"  induced permutation of the start labels: {perm}")
fixed = [i for i in range(3) if perm[i]==i]
print(f"  fixed labels: {fixed}   moved: {[i for i in range(3) if perm[i]!=i]}")
print(f"  start root at the fixed label: {np.round(start[fixed[0]],6) if fixed else None}"
      "   (x = 0 is the sigma-FIXED root)")

print()
print("  sigma acts as x -> -x.  On the start fibre that is:")
for i, r in enumerate(start):
    j = min(range(3), key=lambda k: abs(start[k]-(-r)))
    print(f"    label {i}  ({np.round(r,6)})  -> label {j}  ({np.round(-r,6)})")
print()
print("  VERDICT: if the meridian permutation equals the sigma permutation, then")
print("  'reflection = torus' -- THM-1375's named open question resolves AFFIRMATIVELY.")
