#!/usr/bin/env python3
"""symmetric_transport_verify_kps_S128c112.py -- kind-pasteur-2026-07-20-S128c112

INDEPENDENT VERIFICATION OF A HEADLINE CLAIM BEFORE IT ENTERS CANON.

A research agent reports that applying the de Bondt-van den Essen SYMMETRIC
reduction directly to the THM-1300 Jacobian counterexample F -- skipping the
Bass-Connell-Wright cubic reduction, which the symmetric step does not require --
yields an explicit GRADIENT / symmetric-Jacobian non-injective Keller map on C^6.
The repo's backlog records this transport as "OPEN, unclaimed", so it is checked
rather than trusted.

THE CONSTRUCTION (de Bondt-van den Essen 2005; independently Meng).

    Phi(x,y) = ( F(x), JF(x)^T y )                        on C^{2n}
    T(u,v)   = ( u + i v, (u - i v)/2 )                   linear
    G        = T^{-1} . Phi . T

The twist T is the whole content: the naive potential <y, H(x)> has NON-nilpotent
Hessian, which is exactly the negative opus-S422 recorded; conjugating by T fixes it.

DET IS A PROOF, NOT A COMPUTATION.  JPhi is BLOCK LOWER TRIANGULAR,

    JPhi = [[ JF , 0 ], [ * , JF^T ]]

so det JPhi = det JF * det JF^T = (det JF)^2, and since T is linear,
det JG = det(T)^{-1} * det JPhi * det(T) = (det JF)^2 = 4.  Constant, non-zero:
G is Keller.  No 6x6 symbolic determinant is needed, and none is attempted --
an earlier run of this script died inside simplify() on exactly that.
The identity is spot-checked numerically below rather than assumed.

CHECKED HERE
  1. T^{-1} . T = id.
  2. JG is SYMMETRIC (polynomial identity, via expand -- no simplify).
     On C^6 (simply connected) a polynomial field with symmetric Jacobian IS a
     gradient, so this is the gradient property; P is not needed to establish it.
  3. det JG = (det JF)^2 = 4, spot-checked at random rational points.
  4. The three collision points of F transport by z_a = (a/2, -i a/2) and collide.
  5. HONEST NEGATIVE: Hess P nilpotent?  Needed for Zhao's HN normal form.
"""
import random
from sympy import symbols, I, Matrix, expand, Rational, diff, factor

random.seed(11)
x1, x2, x3 = symbols('x1 x2 x3')
u1, u2, u3, v1, v2, v3 = symbols('u1 u2 u3 v1 v2 v3')
Zv = [u1, u2, u3, v1, v2, v3]

uu = 1 + x1 * x2
F = Matrix([
    uu**3 * x3 + x2**2 * uu * (4 + 3 * x1 * x2),
    x2 + 3 * x1 * uu**2 * x3 + 3 * x1 * x2**2 * (4 + 3 * x1 * x2),
    2 * x1 - 3 * x1**2 * x2 - x1**3 * x3,
])
X = Matrix([x1, x2, x3])
JF = F.jacobian(X)
detF = expand(JF.det())

print("=" * 78)
print("INPUT: the THM-1300 counterexample")
print("=" * 78)
print("  det JF = %s   constant: %s" % (detF, detF.free_symbols == set()))
pts = [(0, 0, Rational(-1, 4)), (1, Rational(-3, 2), Rational(13, 2)),
       (-1, Rational(3, 2), Rational(13, 2))]
ims = [tuple(expand(F[i].subs({x1: p[0], x2: p[1], x3: p[2]})) for i in range(3))
       for p in pts]
print("  three points share the image %s : %s" % (ims[0], len(set(ims)) == 1))

# ---------------------------------------------------------------- (1)
U = Matrix([u1, u2, u3])
V = Matrix([v1, v2, v3])
Tx = U + I * V
Ty = (U - I * V) / 2


def Tinv(xv, yv):
    return ((xv + 2 * yv) / 2, -I * (xv - 2 * yv) / 2)


bu, bv = Tinv(Tx, Ty)
res = [expand(e) for e in list(bu - U) + list(bv - V)]
print()
print("(1) T^{-1}(T(u,v)) - (u,v) = %s   -> inverse correct: %s"
      % (res, all(e == 0 for e in res)))

# ---------------------------------------------------------------- build G
sub = {x1: Tx[0], x2: Tx[1], x3: Tx[2]}
Phi_x = F.subs(sub)
Phi_y = JF.subs(sub).T * Ty
Gu, Gv = Tinv(Phi_x, Phi_y)
G = Matrix([expand(Gu[0]), expand(Gu[1]), expand(Gu[2]),
            expand(Gv[0]), expand(Gv[1]), expand(Gv[2])])
JG = G.jacobian(Matrix(Zv))

# ---------------------------------------------------------------- (2)
print()
print("=" * 78)
print("(2) IS JG SYMMETRIC?  (the gradient property)")
print("=" * 78)
bad = []
for i in range(6):
    for j in range(i + 1, 6):
        if expand(JG[i, j] - JG[j, i]) != 0:
            bad.append((i, j))
print("  asymmetric entries : %d   -> JG SYMMETRIC : %s" % (len(bad), not bad))
print("  (on C^6, simply connected, symmetric Jacobian => G - z is a GRADIENT)")

# ---------------------------------------------------------------- (3)
print()
print("=" * 78)
print("(3) IS G KELLER?  det JG = (det JF)^2 by block triangularity")
print("=" * 78)
print("  JPhi = [[JF, 0], [*, JF^T]] is block lower triangular, so")
print("     det JPhi = det JF * det JF^T = (%s)^2 = %s" % (detF, expand(detF**2)))
print("  and T linear gives det JG = det JPhi.  Spot-check at random points:")
okdet = True
for t in range(4):
    ptn = {vv: Rational(random.randint(-4, 4), random.randint(1, 3)) for vv in Zv}
    dv = expand(JG.subs(ptn).det())
    ok = (dv == expand(detF**2))
    okdet = okdet and ok
    print("     point %d : det JG = %-8s  equals (det JF)^2 : %s" % (t + 1, dv, ok))
print("  all spot-checks agree : %s" % okdet)

# ---------------------------------------------------------------- (4)
print()
print("=" * 78)
print("(4) DOES THE COLLISION TRANSPORT?   z_a = (a/2, -i a/2)")
print("=" * 78)
zs, imgs = [], []
for p in pts:
    a = Matrix(list(p))
    z = list(a / 2) + list(-I * a / 2)
    zs.append(tuple(z))
    s = dict(zip(Zv, z))
    imgs.append(tuple(expand(G[i].subs(s)) for i in range(6)))
for z, im in zip(zs, imgs):
    print("  z = %-42s ->  G(z) = %s" % (str(z), im))
print("  three z distinct        : %s" % (len(set(zs)) == 3))
print("  all three images equal  : %s" % (len(set(imgs)) == 1))
print()
print("  VERDICT: G is an explicit NON-INJECTIVE KELLER map on C^6 whose Jacobian")
print("  is SYMMETRIC -- i.e. a counterexample in the SYMMETRIC/GRADIENT category.")

# ---------------------------------------------------------------- (5)
print()
print("=" * 78)
print("(5) HONEST NEGATIVE -- is Hess P nilpotent?  (Zhao's HN normal form)")
print("=" * 78)
H = [expand(G[i] - Zv[i]) for i in range(6)]
HJ = Matrix(6, 6, lambda i, j: diff(H[i], Zv[j]))
ptn = {vv: Rational(random.randint(-3, 3), random.randint(1, 3)) for vv in Zv}
Hn = HJ.subs(ptn)
pw = Hn
nil = False
for k in range(2, 8):
    pw = expand(pw * Hn)
    if all(expand(e) == 0 for e in pw):
        nil = True
        print("  Hess(P)^%d = 0 at the sample point -> nilpotent there" % k)
        break
if not nil:
    print("  Hess(P)^k != 0 for k <= 7 at the sample point -> NOT nilpotent")
    print("  char poly : %s" % factor(Hn.charpoly().as_expr()))
print()
print("  Non-nilpotent is EXPECTED here: skipping BCW leaves P inhomogeneous, and")
print("  Zhao's vanishing conjecture needs a HOMOGENEOUS QUARTIC P with nilpotent")
print("  Hessian.  So this C^6 object is a symmetric-case counterexample but NOT")
print("  yet a VC witness; that needs the BCW cubic reduction first, which is the")
print("  much larger object (the agent's estimate: N = 79 cubic, then 158).")
