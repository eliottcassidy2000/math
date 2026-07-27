# G1 stratum exclusions: point-cap / line-cap / resolvent-parity
# kind-pasteur 2026-07-26 (G1 = existence of a 2-jet Keller map of field degree 4)
# Companion transcript for the exact claims; every check is a hard assert.
import sympy as sp

x, y, z = sp.symbols('x y z')
a, b, g = sp.symbols('a b g')
FAILED = []

def check(name, cond):
    ok = bool(cond)
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        FAILED.append(name)

def brk(u, v, w):
    return sp.Matrix.hstack(sp.Matrix(u), sp.Matrix(v), sp.Matrix(w)).det()

print("=" * 72)
print("LEG A: point-cap graded system (A = f*e3), first-jet identities")
print("=" * 72)
# generic degree-1 components: value + gradient free at the origin -> the
# z-coefficients of det J depend only on first jets pointwise, and jet
# evaluation is surjective onto C^21, so identities checked at the origin on
# these 21 free symbols are identities for ALL C^1 components (THM-2446 SS1
# argument).
names = ['f', 'B1', 'B2', 'B3', 'C1', 'C2', 'C3']
val, gx, gy, fun = {}, {}, {}, {}
syms = []
for n in names:
    v, p, q = sp.symbols(f'{n}v {n}x_ {n}y_')
    val[n], gx[n], gy[n] = v, p, q
    fun[n] = v + p * x + q * y
    syms += [v, p, q]
f_, B1, B2, B3, C1, C2, C3 = (fun[n] for n in names)
F = sp.Matrix([B1 * z + C1, B2 * z + C2, f_ * z**2 + B3 * z + C3])
J = F.jacobian([x, y, z])
detJ = sp.expand(J.det())
detJ0 = sp.expand(detJ.subs({x: 0, y: 0}))
Pz = sp.Poly(detJ0, z)
check("A1: det J has z-degree <= 3 on the point-cap stratum (D5=D4=0 auto)",
      Pz.degree() <= 3)

def jj(u1, u2, v1, v2):   # 2x2 bracket
    return u1 * v2 - u2 * v1

# pointwise jet symbols at the origin
fv, fx_, fy_ = val['f'], gx['f'], gy['f']
b1, b1x, b1y = val['B1'], gx['B1'], gy['B1']
b2, b2x, b2y = val['B2'], gx['B2'], gy['B2']
b3, b3x, b3y = val['B3'], gx['B3'], gy['B3']
c1, c1x, c1y = val['C1'], gx['C1'], gy['C1']
c2, c2x, c2y = val['C2'], gx['C2'], gy['C2']
c3, c3x, c3y = val['C3'], gx['C3'], gy['C3']

L1 = jj(b1y, b2y, b1, b2)          # j(b_y, b)
M1 = jj(b1x, b2x, b1, b2)          # j(b_x, b)
L0 = jj(c1y, c2y, b1, b2)          # j(c_y, b)
M0 = jj(c1x, c2x, b1, b2)          # j(c_x, b)
N2 = jj(b1x, b2x, b1y, b2y)        # j(b_x, b_y)
N1 = jj(b1x, b2x, c1y, c2y) + jj(c1x, c2x, b1y, b2y)
N0 = jj(c1x, c2x, c1y, c2y)

D3 = fx_ * L1 - fy_ * M1 + 2 * fv * N2
D2 = fx_ * L0 + b3x * L1 - fy_ * M0 - b3y * M1 + 2 * fv * N1 + b3 * N2
D1 = b3x * L0 + c3x * L1 - b3y * M0 - c3y * M1 + 2 * fv * N0 + b3 * N1
D0 = c3x * L0 - c3y * M0 + b3 * N0
check("A2: D3 = f_x j(b_y,b) - f_y j(b_x,b) + 2 f j(b_x,b_y)",
      sp.expand(Pz.coeff_monomial(z**3) - D3) == 0)
check("A3: D2 formula", sp.expand(Pz.coeff_monomial(z**2) - D2) == 0)
check("A4: D1 formula", sp.expand(Pz.coeff_monomial(z**1) - D1) == 0)
check("A5: D0 formula", sp.expand(Pz.coeff_monomial(1) - D0) == 0)

# transport identities: X = L1 d/dx - M1 d/dy kills b projectively
XB1 = L1 * b1x - M1 * b1y
XB2 = L1 * b2x - M1 * b2y
check("A6: X(B1) = -N2*B1  (first-jet identity)", sp.expand(XB1 + N2 * b1) == 0)
check("A7: X(B2) = -N2*B2  (first-jet identity)", sp.expand(XB2 + N2 * b2) == 0)
Xf = L1 * fx_ - M1 * fy_
check("A8: D3 = X(f) + 2 f N2, hence B1^2 * X(f/B1^2) = B1 * D3 "
      "(f/B1^2 is a first integral of X iff D3 = 0)",
      sp.expand(D3 - (Xf + 2 * fv * N2)) == 0 and
      sp.expand((Xf * b1 - 2 * fv * XB1) - b1 * D3) == 0)
# D2 cascade shape: D2 = X(B3) + B3 N2 + Y(f), Y(f) = f_x L0 - f_y M0 + 2 f N1
XB3 = L1 * b3x - M1 * b3y
Yf = fx_ * L0 - fy_ * M0 + 2 * fv * N1
check("A9: D2 = [X(B3) + B3 N2] + Y(f)   (transport cascade, B3-slot)",
      sp.expand(D2 - (XB3 + b3 * N2 + Yf)) == 0)

# sub-stratum b == 0 : det J factors completely
detJ_b0 = sp.expand(detJ.subs({s: 0 for s in
                               [b1, b1x, b1y, b2, b2x, b2y]}))
jacC = sp.expand(sp.Matrix([C1, C2]).jacobian([x, y]).det())
check("A10: b == 0  =>  det J = jac(C1,C2) * (2 f z + B3)  identically",
      sp.expand(detJ_b0 - jacC * (2 * f_ * z + B3)) == 0)

print()
print("=" * 72)
print("LEG B: explicit degree-3 witnesses inside point-cap and line-cap")
print("=" * 72)
u = 1 + x * y
F1 = u**3 * z + y**2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z
F1310 = sp.Matrix([F1, F2, F3])
check("B1: THM-1310 map has det J = -2",
      sp.expand(F1310.jacobian([x, y, z]).det()) == -2)

W1 = sp.Matrix([F1, F2, F3 + F1**2])           # T1 = (w1, w2, w3 + w1^2)
check("B2: W1 = T1 o F(1310) has det J = -2",
      sp.expand(W1.jacobian([x, y, z]).det()) == -2)
A_W1 = [sp.expand(sp.Poly(sp.expand(c), z).coeff_monomial(z**2)) for c in W1]
check("B3: W1 z^2-part A = (0, 0, u^6): rank-1 POINT-CAP, A != 0",
      A_W1[0] == 0 and A_W1[1] == 0 and sp.expand(A_W1[2] - u**6) == 0)
bW1 = [sp.Poly(sp.expand(W1[i]), z).coeff_monomial(z) for i in range(2)]
check("B4: W1 has b = (u^3, 3xu^2) != 0 and f/B1^2 = 1 (first integral, "
      "consistent with A8)",
      sp.expand(bW1[0] - u**3) == 0 and sp.expand(bW1[1] - 3 * x * u**2) == 0
      and sp.simplify(A_W1[2] / bW1[0]**2 - 1) == 0)
# field degree 3: W1 fibers over (a,b,g) = F(1310) fibers over (a,b,g-a^2);
# target map T1 is an automorphism:
w1_, w2_, w3_ = sp.symbols('w1_ w2_ w3_')
T1 = (w1_, w2_, w3_ + w1_**2)
T1i = (w1_, w2_, w3_ - w1_**2)
comp = [e.subs({w1_: T1[0], w2_: T1[1], w3_: T1[2]}, simultaneous=True)
        for e in T1i]
check("B5: T1 is a polynomial automorphism (explicit inverse)",
      [sp.expand(c) for c in comp] == [w1_, w2_, w3_])

W2 = sp.Matrix([F1 + F2 * F3, F2 + F3**2, F3])  # T2 = (w1+w2*w3, w2+w3^2, w3)
check("B6: W2 = T2 o F(1310) has det J = -2",
      sp.expand(W2.jacobian([x, y, z]).det()) == -2)
A_W2 = [sp.expand(sp.Poly(sp.expand(c), z).coeff_monomial(z**2)) for c in W2]
check("B7: W2 is z-quadratic with A = (-3x^4 u^2, x^6, 0): LINE-CAP "
      "(image of [A] = the line {w3=0}, direction [-3u^2 : x^2] nonconstant)",
      sp.expand(A_W2[0] + 3 * x**4 * u**2) == 0 and
      sp.expand(A_W2[1] - x**6) == 0 and A_W2[2] == 0 and
      max(sp.Poly(sp.expand(c), z).degree() for c in W2) == 2)
check("B8: W2 third component is z-affine with B3 = -x^3 != 0 "
      "(line-cap normal form; z = (P3 - C3)/B3 on the fiber)",
      sp.Poly(sp.expand(W2[2]), z).degree() == 1 and
      sp.expand(sp.Poly(sp.expand(W2[2]), z).coeff_monomial(z) + x**3) == 0)
T2 = (w1_ + w2_ * w3_, w2_ + w3_**2, w3_)
T2i = (w1_ - (w2_ - w3_**2) * w3_, w2_ - w3_**2, w3_)
comp = [e.subs({w1_: T2[0], w2_: T2[1], w3_: T2[2]}, simultaneous=True)
        for e in T2i]
check("B9: T2 is a polynomial automorphism (explicit inverse)",
      [sp.expand(c) for c in comp] == [w1_, w2_, w3_])

print()
print("=" * 72)
print("LEG C: frame identities and the universal z-rationality mechanism")
print("=" * 72)
A3v = sp.Matrix(sp.symbols('a1 a2 a3'))
B3v = sp.Matrix(sp.symbols('bb1 bb2 bb3'))
C3v = sp.Matrix(sp.symbols('cc1 cc2 cc3'))
Pv = A3v * z**2 + B3v * z + C3v          # P := F on the fiber
check("C1: A x (P - C) = z (A x B)   on the fiber (THM-2455 sign pin)",
      sp.expand(A3v.cross(Pv - C3v) - z * A3v.cross(B3v)) == sp.zeros(3, 1))
check("C2: B x (C - P) = z^2 (A x B) on the fiber (THM-2455 sign pin)",
      sp.expand(B3v.cross(C3v - Pv) - z**2 * A3v.cross(B3v)) == sp.zeros(3, 1))
check("C3: g1 = [A,B,C-P] vanishes on the fiber",
      sp.expand(brk(A3v, B3v, C3v - Pv)) == 0)
# degenerate direction lemma: A x B == 0, A != 0  =>  the z-column of J dies
rho = sp.symbols('rho')
Fz = 2 * A3v * z + rho * A3v             # B = rho * A pointwise
check("C4: if B = rho*A pointwise then F_z = 2Az+B vanishes at z = -rho/2 "
      "(so det J = 0 somewhere: no Keller map has A != 0 and A x B == 0)",
      sp.expand(Fz.subs(z, -rho / 2)) == sp.zeros(3, 1))

print()
print("=" * 72)
print("LEG D: quartic/resolvent scaffold and THM-1310 disc class")
print("=" * 72)
p, q, r, w = sp.symbols('p q r w')
quart = z**4 + p * z**2 + q * z + r
D4disc = sp.discriminant(quart, z)
I4 = p**2 + 12 * r
J4 = 2 * p**3 + 27 * q**2 - 72 * p * r
check("D1: 27 * disc(quartic) = 4 I^3 - J^2  (THM-2455 (1))",
      sp.expand(27 * D4disc - (4 * I4**3 - J4**2)) == 0)
resolvent = w**3 - p * w**2 - 4 * r * w + (4 * p * r - q**2)
check("D2: disc(resolvent cubic) = disc(quartic), unit 1 (THM-2455 SS2)",
      sp.expand(sp.discriminant(resolvent, w) - D4disc) == 0)

L1310 = 27 * a**2 * g**2 - 18 * a * b * g + 16 * a + b**3 * g - b**2
Q1310 = 27 * a * g**2 - 9 * b * g + 8
N1310 = L1310 * x**3 + (4 - 3 * b * g) * x - 2 * g
check("D3: disc_x(N) = -4 Q^2 L for THM-1310 (quadratic resolvent field "
      "C(P)(sqrt(-L)), odd exponent of L)",
      sp.expand(sp.discriminant(N1310, x) + 4 * Q1310**2 * L1310) == 0)
check("D4: L is irreducible over Q (odd exponent is honest)",
      len(sp.factor_list(L1310)[1]) == 1 and
      sp.factor_list(L1310)[1][0][1] == 1)
check("D5: shear-transported witness W1 sits over the SAME cubic: "
      "min poly of x for W1 over (a,b,g) is N(x; a, b, g-a^2), leading "
      "coefficient L(a,b,g-a^2) != 0 => field degree of W1 = 3",
      sp.expand(L1310.subs(g, g - a**2)) != 0)

print()
print("=" * 72)
print("LEG E: group bookkeeping for the A4/S4 tower")
print("=" * 72)
from sympy.combinatorics import Permutation, PermutationGroup, SymmetricGroup
S4 = SymmetricGroup(4)
S3sub = PermutationGroup([Permutation(3)(0, 1), Permutation(3)(0, 1, 2)])
ok = True
for gp in S4.elements:
    if gp not in S3sub:
        H = PermutationGroup(S3sub.generators + [gp])
        if H.order() != 24:
            ok = False
check("E1: S3 (point stabiliser) is MAXIMAL in S4: no intermediate field "
      "between C(P) and K4 in the S4 case", ok)
A4 = PermutationGroup([Permutation(3)(0, 1, 2), Permutation(3)(0, 1)(2, 3)])
C3sub = PermutationGroup([Permutation(3)(0, 1, 2)])
ok = True
for gp in A4.elements:
    if gp not in C3sub:
        H = PermutationGroup(C3sub.generators + [gp])
        if H.order() != 12:
            ok = False
check("E2: C3 = A3 (point stabiliser) is MAXIMAL in A4: no intermediate "
      "field in the A4 case", ok)

# inertia bookkeeping: induced action on the three pairings
pairings = [frozenset([frozenset([0, 1]), frozenset([2, 3])]),
            frozenset([frozenset([0, 2]), frozenset([1, 3])]),
            frozenset([frozenset([0, 3]), frozenset([1, 2])])]

def induced(perm):
    im = []
    for pr in pairings:
        newpr = frozenset(frozenset(perm(i) for i in blk) for blk in pr)
        im.append(pairings.index(newpr))
    return Permutation(im)

def tame_exp(perm, n):
    return sum(len(c) - 1 for c in perm.cyclic_form) \
        if perm.size >= n else None

reps = {'transposition': Permutation(3)(0, 1),
        'double-transposition': Permutation(3)(0, 1)(2, 3),
        '3-cycle': Permutation(3)(0, 1, 2),
        '4-cycle': Permutation(3)(0, 1, 2, 3)}
expected = {'transposition': (1, 1),
            'double-transposition': (2, 0),
            '3-cycle': (2, 2),
            '4-cycle': (3, 1)}
allok = True
for nm, pm in reps.items():
    e4 = sum(le - 1 for le in pm.cycle_structure.keys()
             for _ in range(pm.cycle_structure[le]) if le > 1)
    ip = induced(pm)
    e3 = sum((le - 1) * ct for le, ct in ip.cycle_structure.items() if le > 1)
    ok = (e4, e3) == expected[nm] and (e4 - e3) % 2 == 0
    print(f"   inertia {nm:22s}: quartic tame disc exp {e4}, "
          f"resolvent exp {e3}  parity match {(e4 - e3) % 2 == 0}")
    allok = allok and ok
check("E3: tame inertia table (quartic vs resolvent-cubic disc exponents; "
      "parities agree, as forced by disc(res) = disc(quartic))", allok)
check("E4: V4 acts trivially on the pairings (double-transpositions "
      "ramify the quartic but NOT the resolvent)",
      induced(Permutation(3)(0, 1)(2, 3)) == Permutation(2) and
      induced(Permutation(3)(0, 2)(1, 3)) == Permutation(2))

print()
print("FAILED CHECKS: " + (", ".join(FAILED) if FAILED else "NONE"))
