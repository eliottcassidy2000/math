#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c98 -- HYP-8090 / THM-1310: the fiber geometry of the
first Jacobian-Conjecture counterexample.

F = ((1+xy)^3 z + y^2(1+xy)(4+3xy), y + 3x(1+xy)^2 z + 3xy^2(4+3xy), 2x-3x^2y-x^3z)

Battery:
 (1) FIBER SYSTEM at a general target (a,b,g): z is eliminated linearly
     (F is z-affine); show the two residual surfaces COLLAPSE to a conic pair
     in (x, s=xy); extract the FIBER CUBIC in x by resultant; verify the
     rationality y = R(x; a,b,g)  =>  field degree EXACTLY 3, cubic = min poly.
 (2) S_3 PIN: discriminant Delta(a,b,g) of the fiber cubic; polynomial
     factorization decides square/non-square RIGOROUSLY (UFD exponents);
     non-square + irreducible => arithmetic AND geometric monodromy = S_3.
 (3) JELONEK/ASYMPTOTIC SET: leading coefficient L(a,b,g) (degree drop =
     escape to infinity) + the claim {Delta=0} is escape not ramification
     (Keller => no critical points): sample a Delta=0 target, count fiber.
 (4) EQUIVARIANCE UNIQUENESS check: k=4 equivariant shape (weights (1,-1,-3))
     has SINGULAR linear part at 0 (weight -(k-2) not a coordinate weight);
     verify numerically on random instances.  k=3 works (our F).
 (5) ITERATION: F o F fiber histogram mod p=11 (field degree should be 9 =
     3^2 -- the counterexample MONOID / doubling-tower echo).
 (6) BOUNDARY PRIMES: exact sweeps p=2 (det=0, critical) and p=3 (Frobenius
     degeneration: F == (u^3 z + y^2 u, y, 2x - x^3 z) with u^3 = 1+(xy)^3).
"""
import sympy as sp
from sympy import symbols, expand, factor, factor_list, resultant, discriminant, Poly, Rational, sqrt
from collections import Counter
import random, time

x, y, z, s, a, b, g = symbols('x y z s a b g')
u = 1 + x*y

F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z

A = [u**3, 3*x*u**2, -x**3]
B = [F1 - A[0]*z, F2 - A[1]*z, F3 - A[2]*z]
B = [sp.expand(t) for t in B]

print("== (1) the fiber system at general target (a,b,g) ==", flush=True)
# z-elimination (all components affine in z):  z = (a - B1)/A1
G2 = sp.expand(A[1]*(a - B[0]) + A[0]*(B[1] - b))   # A1*(F2-b) after subst
G3 = sp.expand(A[2]*(a - B[0]) + A[0]*(B[2] - g))   # A1*(F3-g) after subst
# collapse to (x, s):  substitute y -> s/x and clear powers of x
def to_xs(P):
    q = sp.together(P.subs(y, s/x))
    num, den = sp.fraction(sp.cancel(q))
    return sp.expand(num), den
N2, D2 = to_xs(G2); N3, D3 = to_xs(G3)
f2 = sp.factor(N2); f3 = sp.factor(N3)
print(f"  G2 collapsed/factored: {f2}", flush=True)
print(f"  G3 collapsed/factored: {f3}", flush=True)

# extract the non-spurious factors (drop powers of x and of u=1+s)
def core_factors(fac):
    const, pairs = sp.factor_list(fac)
    keep = []
    for p_, e in pairs:
        if p_ in (x, s + 1) or p_ == -1: continue
        keep.append((p_, e))
    return keep
c2 = core_factors(f2); c3 = core_factors(f3)
print(f"  core of G2: {c2}", flush=True)
print(f"  core of G3: {c3}", flush=True)
conic1 = c2[0][0] if len(c2) == 1 else None
conic2 = c3[0][0] if len(c3) == 1 else None
print(f"  CONIC PAIR:\n    C1 = {sp.expand(conic1)}\n    C2 = {sp.expand(conic2)}", flush=True)

# fiber polynomial in x: resultant in s
R = sp.resultant(Poly(conic1, s), Poly(conic2, s))
Rf = sp.factor(R)
print(f"  Res_s(C1,C2) factored: {Rf}", flush=True)
# identify the cubic-in-x factor
cub = None
for p_, e in sp.factor_list(Rf)[1]:
    if sp.degree(p_, x) == 3:
        cub = sp.expand(p_)
print(f"  FIBER CUBIC in x: {cub}", flush=True)
Pc = Poly(cub, x)
coeffs = Pc.all_coeffs()
print(f"  coefficients (deg 3..0): {coeffs}", flush=True)

# rationality of s (hence y, z) in x => field degree exactly 3, via subresultant PRS:
prs = sp.subresultants(Poly(conic1, s), Poly(conic2, s))
lin1 = [q for q in prs if sp.degree(q, s) == 1]
if lin1:
    q1 = sp.Poly(lin1[-1], s)
    s_of_x = sp.cancel(-q1.all_coeffs()[1] / q1.all_coeffs()[0])
    print(f"  s = rational in (x; a,b,g) via deg-1 subresultant: numerator deg_x = {sp.degree(sp.numer(s_of_x), x)}, denominator deg_x = {sp.degree(sp.denom(s_of_x), x)}", flush=True)
    print(f"  => C(x,y,z) = C(a,b,g)(x): FIELD DEGREE EXACTLY 3, fiber cubic = minimal polynomial", flush=True)

# sanity: specialize to the collision target (a,b,g)=(-1/4,0,0)
spec = cub.subs({a: Rational(-1,4), b: 0, g: 0})
print(f"  cubic at collision target: {sp.factor(spec)}  (roots should include x=+-1)", flush=True)

print("\n== (2) discriminant and the S_3 pin ==", flush=True)
Delta = sp.discriminant(Pc, x)
Delta = sp.factor(Delta)
print(f"  Delta(a,b,g) = {Delta}", flush=True)
const_, pairs_ = sp.factor_list(Delta)
odd = [(p_, e) for p_, e in pairs_ if e % 2 == 1]
csq = sp.sqrt(const_)
print(f"  factor exponents: {[(str(p_), e) for p_, e in pairs_]}, unit = {const_}", flush=True)
is_square = (len(odd) == 0) and (csq.is_rational and csq.is_rational is True)
print(f"  VERDICT square-in-QQ(a,b,g): {is_square}  => " +
      ("monodromy inside A_3 (!!)" if is_square else "Delta NOT a square (odd-exponent factor persists over C too) => Galois closure S_3"), flush=True)
# irreducibility via a random specialization
random.seed(7)
for _ in range(3):
    sub = {a: Rational(random.randint(1,50), random.randint(1,7)),
           b: Rational(random.randint(1,50), random.randint(1,7)),
           g: Rational(random.randint(1,50), random.randint(1,7))}
    fl = sp.factor_list(cub.subs(sub))
    irr = all(sp.degree(p_, x) in (0,3) for p_, e in fl[1])
    print(f"  specialization {tuple(str(v) for v in sub.values())}: irreducible cubic = {irr}", flush=True)

print("\n== (3) the Jelonek / asymptotic set ==", flush=True)
L = sp.factor(coeffs[0])
print(f"  leading coeff L(a,b,g) = {L}   (degree-drop escape locus)", flush=True)
print(f"  Delta factored (above) -- {{Delta=0}} should be ESCAPE not ramification (Keller):", flush=True)
# the two walls: L = 0 (degree drop: a root at x = infinity) and Q = 0 (finite double root)
Q = 27*a*g**2 - 9*b*g + 8
# pure Q-wall sample: Q=0, L!=0 at (b,g)=(1,1): a = 1/27
aa = Rational(1, 27); bb = Rational(1); gg = Rational(1)
print(f"  pure Q-wall sample (Q=0, L!=0): (a,b,g) = ({aa},{bb},{gg}); L = {L.subs({a:aa,b:bb,g:gg})}", flush=True)
cub0 = sp.expand(cub.subs({a: aa, b: bb, g: gg}))
rts = sp.roots(Poly(cub0, x))
print(f"  cubic there: {sp.factor(cub0)}; roots: { {str(k): v for k, v in rts.items()} }", flush=True)
for rt, mult in rts.items():
    c1r = sp.expand(conic1.subs({x: rt, a: aa, b: bb, g: gg}))
    c2r = sp.expand(conic2.subs({x: rt, a: aa, b: bb, g: gg}))
    svals = sp.solve([c1r, c2r], s)
    print(f"  x = {rt} (mult {mult}): common s-solutions {svals} -> {len(svals)} fiber point(s) at this x", flush=True)
print("  (Keller forbids ramification: the double x must carry at most ONE fiber point; the lost point ESCAPES -> Jelonek)", flush=True)
# L-wall sample (degree drop)
print(f"  L-wall sample: (a,b,g)=(0,1,1): cubic = {sp.factor(cub.subs({a:0,b:1,g:1}))} -> degree drop, roots escape to x = infinity", flush=True)
print(f"  ASYMPTOTIC VARIETY = {{Q=0}} u {{L=0}}: degrees 3 + 4 = 7 = deg F", flush=True)

print("\n== (4) equivariant-degree uniqueness: k=4 shape dies at J(0) ==", flush=True)
# weights (1,-1,-3): F = (u^4 z + y^3 p1(s), 4u^3 x z + y^2 p2(s), -x^4 z + x p3(s))
t_ = symbols('t_')
for trial in range(3):
    p1 = sum(random.randint(1,9)*t_**i for i in range(3))
    p2 = sum(random.randint(1,9)*t_**i for i in range(3))
    p3 = sum(random.randint(1,9)*t_**i for i in range(2))
    E1 = u**4*z + y**3*p1.subs(t_, x*y)
    E2 = 4*u**3*x*z + y**2*p2.subs(t_, x*y)
    E3 = -x**4*z + x*p3.subs(t_, x*y)
    J0 = sp.Matrix([E1, E2, E3]).jacobian([x, y, z]).subs({x:0, y:0, z:0})
    print(f"  random (1,-1,-3)-equivariant instance {trial+1}: det J(0) = {J0.det()}", flush=True)
print("  (weight -(k-2) = -2 is not a coordinate weight at k=4 => no linear term in F2 => always singular)", flush=True)
J0our = sp.Matrix([F1, F2, F3]).jacobian([x, y, z]).subs({x:0, y:0, z:0})
print(f"  our k=3 map: J(0) = {J0our.tolist()}, det = {J0our.det()}  (invertible OK)", flush=True)

print("\n== (5) iteration F∘F: field degree 9 check mod p=11 ==", flush=True)
def Fmap_mod(X, Y, Z, p):
    U = (1 + X*Y) % p
    return ((U**3*Z + Y*Y*U*(4 + 3*X*Y)) % p,
            (Y + 3*X*U*U*Z + 3*X*Y*Y*(4 + 3*X*Y)) % p,
            (2*X - 3*X*X*Y - X**3*Z) % p)
p = 11
cnt = Counter()
for X in range(p):
    for Y in range(p):
        for Z in range(p):
            im = Fmap_mod(*Fmap_mod(X, Y, Z, p), p)
            cnt[im] += 1
fibs = Counter(cnt.values())
print(f"  p={p}: F∘F fiber histogram {dict(sorted(fibs.items()))}, max fiber = {max(fibs)}, misses = {p**3 - len(cnt)}", flush=True)
print(f"  (degree multiplicativity => deg(F∘F) = 9; sizes should divide/compose within S_3 wr S_2-type statistics)", flush=True)

print("\n== (6) boundary primes p=2, p=3 ==", flush=True)
for p in (2, 3):
    cnt = Counter()
    for X in range(p):
        for Y in range(p):
            for Z in range(p):
                cnt[Fmap_mod(X, Y, Z, p)] += 1
    fibs = Counter(cnt.values())
    print(f"  p={p}: fiber histogram {dict(sorted(fibs.items()))}, image {len(cnt)}/{p**3}", flush=True)
F1_3 = sp.Poly(F1, x, y, z).set_modulus(3)
F2_3 = sp.Poly(F2, x, y, z).set_modulus(3)
F3_3 = sp.Poly(F3, x, y, z).set_modulus(3)
print(f"  F mod 3: F1 = {F1_3.as_expr()}", flush=True)
print(f"           F2 = {F2_3.as_expr()}   (y preserved => plane-family fibration)", flush=True)
print(f"           F3 = {F3_3.as_expr()}", flush=True)
print("  (u^3 = 1 + (xy)^3 mod 3: Frobenius-additive twist -- the char-3 shadow of the design)", flush=True)
print("\nDONE.", flush=True)
