#!/usr/bin/env python3
r"""
jacobian_twojet_zgraded_identity_kps_S131.py
(kind-pasteur-2026-07-26-S131; companion to THM-2445)

Exact companion for the z-graded Jacobian decomposition of 2-jet
(z-quadratic) polynomial maps

    F(x,y,z) = A(x,y) z^2 + B(x,y) z + C(x,y),   A,B,C : C^2 -> C^3,

with J = [F_x | F_y | F_z] (columns) and bracket [u,v,w] = det(u|v|w):

    det J = sum_{k=0}^{5} D_k z^k,
    D5 = 2[A_x,A_y,A]
    D4 = [A_x,A_y,B] + 2[A_x,B_y,A] + 2[B_x,A_y,A]
    D3 = [A_x,B_y,B] + [B_x,A_y,B]
         + 2([A_x,C_y,A] + [C_x,A_y,A] + [B_x,B_y,A])
    D2 = [B_x,B_y,B] + [A_x,C_y,B] + [C_x,A_y,B]
         + 2([B_x,C_y,A] + [C_x,B_y,A])
    D1 = [B_x,C_y,B] + [C_x,B_y,B] + 2[C_x,C_y,A]
    D0 = [C_x,C_y,B].

Keller <=> D5=D4=D3=D2=D1=0 in C[x,y] and D0 = nonzero constant.

Verification legs, all exact symbolic (sympy), no sampling:
  1. Fully generic degree-1 components (27 free coefficients) --
     sufficient for a proof because both sides depend only on 1-jets
     and the jet evaluation is surjective onto C^27; any wrong integer
     bracket multiplicity leaves a nonzero universal trilinear form.
  2. Fully generic degree-2 components (54 free coefficients) --
     redundant confirmation.
  3. Canon positive control: the THM-1310 dim-3 Keller map embedded as
     A=0; requires D2=0, D1 = -6x^2(xy+1) + 6x^2(xy+1) = 0 (exact
     cancellation of the two brackets), D0 = -2.
  4. Hostile control: C1 -> C1 + x*y^3 must break Keller with exactly
     D1' = -6x^4y^2(xy+1)(xy+3),
     D0' = -2(9x^5y^5+30x^4y^4+5x^3y^3-9x^2y^2+1), D2' = 0.
  5. Euler slice reductions on the homogeneous leading map
     Phi = alpha z^2 + beta z + gamma (deg 2,3,4):
     E5 == 0 identically; E4 = -2[alpha_x,alpha_y,beta];
     E3 = (4/3)[beta_x,beta_y,alpha] - 4[alpha_x,alpha_y,gamma];
     and [Phi_x,Phi_y,Phi] == (z/4) det J(Phi).
  6. Tame witness: F = (x+z^2, y, z) satisfies the full graded system
     with D0 = 1 -- so system-consistency is NOT a wildness
     certificate in the 2-jet box (unlike the z-affine d=3 box).
  7. d=3 ruling recheck: the counterexample's B = (u^3, 3xu^2, -x^3)
     satisfies the cuspidal-cone identity B2^3 + 27 B1^2 B3 = 0.

Reproduction: python 04-computation/jacobian_twojet_zgraded_identity_kps_S131.py
"""
import sympy as sp

x, y, z = sp.symbols("x y z")


def fail(msg):
    raise SystemExit("CHECK FAILED: " + msg)


def bracket(u, v, w):
    return sp.Matrix.hstack(sp.Matrix(u), sp.Matrix(v), sp.Matrix(w)).det()


def d_coeffs(A, B, C):
    Ax = [sp.diff(a, x) for a in A]; Ay = [sp.diff(a, y) for a in A]
    Bx = [sp.diff(b, x) for b in B]; By = [sp.diff(b, y) for b in B]
    Cx = [sp.diff(c, x) for c in C]; Cy = [sp.diff(c, y) for c in C]
    D5 = 2 * bracket(Ax, Ay, A)
    D4 = bracket(Ax, Ay, B) + 2 * bracket(Ax, By, A) + 2 * bracket(Bx, Ay, A)
    D3 = (bracket(Ax, By, B) + bracket(Bx, Ay, B)
          + 2 * (bracket(Ax, Cy, A) + bracket(Cx, Ay, A) + bracket(Bx, By, A)))
    D2 = (bracket(Bx, By, B) + bracket(Ax, Cy, B) + bracket(Cx, Ay, B)
          + 2 * (bracket(Bx, Cy, A) + bracket(Cx, By, A)))
    D1 = bracket(Bx, Cy, B) + bracket(Cx, By, B) + 2 * bracket(Cx, Cy, A)
    D0 = bracket(Cx, Cy, B)
    return [D0, D1, D2, D3, D4, D5]


def detJ(A, B, C):
    F = [A[i] * z ** 2 + B[i] * z + C[i] for i in range(3)]
    J = sp.Matrix([[sp.diff(F[i], v) for v in (x, y, z)] for i in range(3)])
    return sp.expand(J.det())


def check_identity(A, B, C, tag):
    lhs = detJ(A, B, C)
    Ds = d_coeffs(A, B, C)
    rhs = sp.expand(sum(Ds[k] * z ** k for k in range(6)))
    if sp.expand(lhs - rhs) != 0:
        fail(f"identity residual nonzero [{tag}]")
    print(f"identity exact [{tag}]: PASS")
    return Ds


def generic_vec(name, deg):
    coeffs = []
    comp = []
    for i in range(3):
        e = sp.Integer(0)
        for dx in range(deg + 1):
            for dy in range(deg + 1 - dx):
                s = sp.Symbol(f"{name}{i}_{dx}{dy}")
                coeffs.append(s)
                e += s * x ** dx * y ** dy
        comp.append(e)
    return comp, coeffs


# ---- leg 1: generic degree 1 (proof-grade by 1-jet surjectivity) ----
A1, ca = generic_vec("a", 1)
B1, cb = generic_vec("b", 1)
C1v, cc = generic_vec("c", 1)
assert len(ca) + len(cb) + len(cc) == 27
check_identity(A1, B1, C1v, "generic deg 1, 27 symbols")

# ---- leg 2: generic degree 2 ----
A2, _ = generic_vec("p", 2)
B2, _ = generic_vec("q", 2)
C2, _ = generic_vec("r", 2)
check_identity(A2, B2, C2, "generic deg 2, 54 symbols")

# ---- leg 3: THM-1310 canon positive control (A = 0) ----
u = 1 + x * y
Bc = [u ** 3, 3 * x * u ** 2, -x ** 3]
Cc = [y ** 2 * u * (4 + 3 * x * y),
      y + 3 * x * y ** 2 * (4 + 3 * x * y),
      2 * x - 3 * x ** 2 * y]
Ds = check_identity([sp.Integer(0)] * 3, Bc, Cc, "THM-1310 embed")
if sp.expand(Ds[2]) != 0 or sp.expand(Ds[3]) != 0:
    fail("D2/D3 nonzero on canon map")
Bx = [sp.diff(b, x) for b in Bc]; By = [sp.diff(b, y) for b in Bc]
Cx = [sp.diff(c, x) for c in Cc]; Cy = [sp.diff(c, y) for c in Cc]
t1 = sp.expand(bracket(Bx, Cy, Bc))
t2 = sp.expand(bracket(Cx, By, Bc))
if t1 != sp.expand(-6 * x ** 2 * (x * y + 1)) or t2 != sp.expand(6 * x ** 2 * (x * y + 1)):
    fail(f"D1 cancellation parts wrong: {t1}, {t2}")
if sp.expand(Ds[1]) != 0 or sp.expand(Ds[0]) != -2:
    fail("canon D1/D0 wrong")
print("canon control: D1 = -6x^2(xy+1) + 6x^2(xy+1) = 0, D0 = -2: PASS")

# ---- leg 4: hostile control ----
Ch = [Cc[0] + x * y ** 3, Cc[1], Cc[2]]
Dh = check_identity([sp.Integer(0)] * 3, Bc, Ch, "hostile xy^3")
if sp.expand(Dh[2]) != 0:
    fail("hostile D2 changed (C never enters D2 when A=0)")
if sp.factor(Dh[1]) != sp.factor(-6 * x ** 4 * y ** 2 * (x * y + 1) * (x * y + 3)):
    fail(f"hostile D1' wrong: {sp.factor(Dh[1])}")
if sp.factor(Dh[0] + 2 * (9 * x**5 * y**5 + 30 * x**4 * y**4 + 5 * x**3 * y**3
                          - 9 * x**2 * y**2 + 1)) != 0:
    fail(f"hostile D0' wrong: {sp.factor(Dh[0])}")
print("hostile control: failure localized to D1'/D0', exact forms: PASS")

# ---- leg 5: Euler slice reductions on homogeneous leading maps ----
def generic_hom_vec(name, deg):
    comp = []
    for i in range(3):
        e = sp.Integer(0)
        for dx in range(deg + 1):
            s = sp.Symbol(f"{name}{i}_{dx}")
            e += s * x ** dx * y ** (deg - dx)
        comp.append(e)
    return comp

al = generic_hom_vec("al", 2)
be = generic_hom_vec("be", 3)
ga = generic_hom_vec("ga", 4)
Dh5 = d_coeffs(al, be, ga)
alx = [sp.diff(a, x) for a in al]; aly = [sp.diff(a, y) for a in al]
bex = [sp.diff(b, x) for b in be]; bey = [sp.diff(b, y) for b in be]
gax = [sp.diff(g, x) for g in ga]; gay = [sp.diff(g, y) for g in ga]
if sp.expand(Dh5[5]) != 0:
    fail("E5 not identically zero on homogeneous alpha")
E4 = sp.expand(Dh5[4] - (-2) * bracket(alx, aly, be))
E3 = sp.expand(Dh5[3] - (sp.Rational(4, 3) * bracket(bex, bey, al)
                         - 4 * bracket(alx, aly, ga)))
if E4 != 0 or E3 != 0:
    fail("E4/E3 slice reductions wrong")
Phi = [al[i] * z ** 2 + be[i] * z + ga[i] for i in range(3)]
Jp = sp.Matrix([[sp.diff(Phi[i], v) for v in (x, y, z)] for i in range(3)])
lhs = sp.expand(4 * bracket([sp.diff(p, x) for p in Phi],
                            [sp.diff(p, y) for p in Phi], Phi))
if sp.expand(lhs - z * sp.expand(Jp.det())) != 0:
    fail("Euler cone identity 4[Phi_x,Phi_y,Phi] = z det J(Phi) wrong")
print("Euler reductions: E5 vacuous, E4 = -2[ax,ay,b], "
      "E3 = (4/3)[bx,by,a] - 4[ax,ay,g], cone identity: PASS")

# ---- leg 6: tame witness ----
Dt = d_coeffs([sp.Integer(1), sp.Integer(0), sp.Integer(0)],
              [sp.Integer(0), sp.Integer(0), sp.Integer(1)],
              [x, y, sp.Integer(0)])
if any(sp.expand(Dt[k]) != 0 for k in range(1, 6)) or Dt[0] != 1:
    fail("tame witness (x+z^2, y, z) fails graded system")
print("tame witness F=(x+z^2,y,z): graded-consistent, D0=1, field degree 1: PASS")

# ---- leg 7: d=3 ruling cusp identity ----
if sp.expand(Bc[1] ** 3 + 27 * Bc[0] ** 2 * Bc[2]) != 0:
    fail("d=3 cuspidal ruling identity broken")
print("d=3 ruling: B2^3 + 27 B1^2 B3 = 0: PASS")

print("ALL CHECKS PASSED")
