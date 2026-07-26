# P2 of THM-2446 on the CONIC-CAP stratum: A = nu = (x^2, x*y, y^2)
# (nondegenerate conical A: values span C^3, direction map dominant onto the
# smooth conic XZ = Y^2).
#
#   STEP 1  structural identities on this stratum, arbitrary B:
#             D4 = 2*(3n - theta n),  n := N.B,  N = (y^2,-2xy,x^2),
#                                     theta = x d/dx + y d/dy
#             E := [w_x,w_y,w] = n * [nu, B, theta B]
#           so D4 == 0 iff n is homogeneous cubic, and homogeneous B give
#           E = 0 automatically.
#   STEP 2  the deg<=3 box: D4 == 0 has solution dimension exactly 16 =
#           { B1 linear free (6) } + { B2 = beta*nu_p + gamma*nu_q,
#           beta,gamma linear (4) } + { B3 tangential cubic (6) }; B^(0)=0.
#   STEP 3  weaker form: D5 = D4 = 0 alone does not force E = 0
#           (B = (2x^2, xy, x), S3 != 0 but not imposed).
#   STEP 4  COUNTEREXAMPLE killing P2 on this stratum:
#             B = (0,0,x) + s*x*nu_q - 2*s^2*x*nu
#           satisfies D5 = D4 = S3 = 0 identically and E = -s*x^8 != 0.
#           (The B3-part is the z-gauge shift B -> B + 2*phi*A with
#           phi = -s^2 x; the gauge-essential part is B1 + s*x*nu_q.)
#   STEP 5  structure notes: S3(B1,B1) == 0 identically for linear B1;
#           on the slice B3 = 0 the S3 system forces B2 proportional to nu
#           (hence E = 0): the counterexample NEEDS the cubic tail.
#
# All checks hard asserts; final line ALL CHECKS PASSED.

import sympy as sp

x, y, s = sp.symbols('x y s')

def br(u, v, w):
    return sp.Matrix.hstack(u, v, w).det()

def cross(u, v):
    return sp.Matrix(3, 1, [u[1]*v[2] - u[2]*v[1],
                            u[2]*v[0] - u[0]*v[2],
                            u[0]*v[1] - u[1]*v[0]])

nu  = sp.Matrix([x**2, x*y, y**2])
nup = sp.Matrix([2*x, y, 0])
nuq = sp.Matrix([0, x, 2*y])
N   = sp.Matrix([y**2, -2*x*y, x**2])

def D4f(A, B):
    Ax, Ay, Bx, By = A.diff(x), A.diff(y), B.diff(x), B.diff(y)
    return br(Ax, Ay, B) + 2*br(Ax, By, A) + 2*br(Bx, Ay, A)

def S3f(A, B):
    Ax, Ay, Bx, By = A.diff(x), A.diff(y), B.diff(x), B.diff(y)
    return br(Ax, By, B) + br(Bx, Ay, B) + 2*br(Bx, By, A)

def Ef(A, B):
    w = cross(A, B)
    return br(w.diff(x), w.diff(y), w)

# ---------------------------------------------------------------- STEP 1
Bg = sp.Matrix(3, 1, [sp.Function('b%d' % i)(x, y) for i in (1, 2, 3)])
n  = (N.T*Bg)[0]
th = lambda f: x*sp.diff(f, x) + y*sp.diff(f, y)
assert sp.simplify(D4f(nu, Bg) - 2*(3*n - th(n))) == 0
thB = x*Bg.diff(x) + y*Bg.diff(y)
assert sp.simplify(Ef(nu, Bg) - n*br(nu, Bg, thB)) == 0
print("STEP 1 (A = nu, arbitrary B):")
print("  D4 = 2*(3n - theta n),  n = N.B,  N = (y^2, -2xy, x^2)")
print("  E  = n * [nu, B, theta B]")
print("  => D4 == 0 iff n homogeneous cubic; homogeneous B => E = 0.")

# ---------------------------------------------------------------- STEP 2
mons = [x**i*y**j for i in range(4) for j in range(4) if i + j <= 3]
cs, Bgen = [], sp.zeros(3, 1)
for comp in range(3):
    for k, m in enumerate(mons):
        c = sp.Symbol('c_%d_%d' % (comp, k))
        cs.append(c)
        Bgen[comp] += c*m
eqs = sp.Poly(sp.expand(D4f(nu, Bgen)), x, y).coeffs()
Mmat, rhs = sp.linear_eq_to_matrix(eqs, cs)
nullity = len(cs) - Mmat.rank()
t1 = sp.symbols('t1_0:6'); t2 = sp.symbols('t2_0:4'); t3 = sp.symbols('t3_0:6')
B1 = sp.Matrix([t1[0]*x + t1[1]*y, t1[2]*x + t1[3]*y, t1[4]*x + t1[5]*y])
B2 = (t2[0]*x + t2[1]*y)*nup + (t2[2]*x + t2[3]*y)*nuq
B3 = ((t3[0]*x**2 + t3[1]*x*y + t3[2]*y**2)*nup
      + (t3[3]*x**2 + t3[4]*x*y + t3[5]*y**2)*nuq)
Bpar = sp.expand(B1 + B2 + B3)
assert sp.expand(D4f(nu, Bpar)) == 0
assert nullity == 16
print("STEP 2: deg<=3 box: D4==0 solution space has dim 16, matched by")
print("        B = B1(free linear) + tangential quadratic + tangential cubic;")
print("        the constant part B^(0) is forced to 0.")

# closed form of E on the D4-solved family
nn   = sp.expand((N.T*Bpar)[0])
Ecls = sp.expand(nn*br(nu, B1, B2 + 2*B3))
assert sp.expand(Ef(nu, Bpar) - Ecls) == 0
print("        on this family: E = n * [nu, B1, B2 + 2*B3],  n = N.B1.")

# ---------------------------------------------------------------- STEP 3
Bw = sp.Matrix([0, 0, x]) + x*nup
assert sp.expand(D4f(nu, Bw)) == 0
Ew  = sp.factor(sp.expand(Ef(nu, Bw)))
S3w = sp.factor(sp.expand(S3f(nu, Bw)))
assert Ew != 0
print("STEP 3: weaker form on the stratum: B = (2x^2, xy, x):")
print("        D4 = 0,  S3 =", S3w, " (not imposed),  E =", Ew, "!= 0")
print("        => D5 = D4 = 0 alone does NOT force E = 0 here either.")

# ---------------------------------------------------------------- STEP 4
Bc = sp.expand(sp.Matrix([0, 0, x]) + s*x*nuq - 2*s**2*x*nu)
d4 = sp.expand(D4f(nu, Bc))
s3 = sp.expand(S3f(nu, Bc))
E  = sp.factor(sp.expand(Ef(nu, Bc)))
assert d4 == 0 and s3 == 0 and E == -s*x**8
print("STEP 4: CONIC-STRATUM COUNTEREXAMPLE (one-parameter family):")
print("        A =", list(nu))
print("        B =", list(Bc))
print("        D5/2 = 0, D4 =", d4, ", S3 =", s3, " identically;  E =", E)
w = cross(nu, Bc)
print("        w = A x B =", [sp.factor(c) for c in w])
Bc1 = sp.expand(Bc.subs(s, 1))
E1  = sp.factor(sp.expand(Ef(nu, Bc1)))
assert E1 == -x**8
print("        s = 1 witness: B =", list(Bc1), ";  E =", E1)
print("        => P2 REFUTED on the nondegenerate conic-cap stratum.")

# gauge decomposition note: B3-part is a z-translation gauge shift
phi = -s**2*x
assert sp.expand(Bc - (sp.Matrix([0, 0, x]) + s*x*nuq + 2*phi*nu)) == sp.zeros(3, 1)
print("        gauge note: B = [(0,0,x) + s*x*nu_q] + 2*phi*A, phi = -s^2*x;")
print("        w and E are invariant under B -> B + 2*phi*A.")

# ---------------------------------------------------------------- STEP 5
s3_lin = sp.expand(S3f(nu, B1))
assert s3_lin == 0
print("STEP 5a: S3(B1) == 0 identically for EVERY linear B1 (no condition).")

# slice B3 = 0: the pure-B2 part of S3 (degree-4 coefficients) forces
# B2 proportional to nu, hence E = 0 on that slice.
s3_2 = sp.expand(S3f(nu, B2))
eqs2 = sp.Poly(s3_2, x, y).coeffs()
sol2 = sp.solve(eqs2, list(t2), dict=True)
assert len(sol2) == 1 and sp.simplify(sol2[0][t2[1]]) == 0 \
       and sp.simplify(sol2[0][t2[2]]) == 0 and sol2[0][t2[0]] == t2[3]
B2forced = sp.expand(B2.subs(sol2[0]))
assert sp.expand(B2forced - 2*t2[3]*nu) == sp.zeros(3, 1)
print("STEP 5b: on the slice B3 = 0, S3 == 0 forces B2 = 2*t*nu (gauge")
print("         shift), hence E = 0: the counterexample NEEDS the cubic")
print("         tangential tail B3 to unlock a non-gauge B2.")

print("ALL CHECKS PASSED")
