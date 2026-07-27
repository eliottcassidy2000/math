"""Exact certificate for THM-2546.

This script proves, by polynomial identities over QQ, that both transverse
coordinates of the sporadic Keller map are integral, resolves the full
L = 0 slice, and distinguishes genuine affine fibre points from roots of
coordinate eliminants that are only shadows of the two escaping x-sheets.

The source coordinates are (x,y,z), the target coordinates are (a,b,c), and
r = b-y.  In particular the sign convention is

    P_r = 27*a^2*c - 18*a*r + 3*b*r^2 - 2*r^3,
    P_r' = -6*D_r,  D_r = r^2-b*r+3*a.

opus/codex 2026-07-27.  Companion to THM-2473 and THM-2546.
"""

import sympy as sp


x, y, z, a, b, c, r, W = sp.symbols("x y z a b c r W")
u = 1 + x * y
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
p = 4 - 3 * b * c
q = b**2 - 12 * a

Pr = 27 * a**2 * c - 18 * a * r + 3 * b * r**2 - 2 * r**3
Dr = r**2 - b * r + 3 * a


def zero(expr):
    """Assert an exact rational expression is zero."""
    assert sp.cancel(expr) == 0


print("=" * 78)
print("[1] Global r-coordinate identities and the sign convention")
print("=" * 78)
zero(sp.diff(Pr, r) + 6 * Dr)
target_sub = {a: F1, b: F2, c: F3, r: F2 - y}
zero(Pr.subs(target_sub, simultaneous=True))
zero((x * Dr - r).subs(target_sub, simultaneous=True))
disc_r = sp.factor(sp.discriminant(Pr, r))
assert disc_r == -2916 * a**2 * L
print("  P_r(r) =", Pr)
print("  D_r(r) =", Dr)
print("  P_r'(r) = -6 D_r(r)  [PASS; the sign is negative]")
print("  on every fibre: P_r(r)=0 and x D_r(r)=r  [PASS]")
print("  disc_r(P_r) =", disc_r)
print("  Tr(r)=3b/2 and Tr(y)=3b/2, whereas Tr(x)=0 (depressed x-core).")
print("  Thus a=0 is a coordinate-collision component independent of L=0.")
print()

print("=" * 78)
print("[2] Global constant-lead z-cubic: z is integral (PROVED)")
print("=" * 78)
Az = (
    324 * a**2 * c**2
    - 216 * a * b * c
    + 408 * a
    - 15 * b**3 * c
    + 6 * b**2
)
Bz = 27 * a**2 * c**2 - 18 * a * b * c + 52 * a + b**3 * c + 14 * b**2
Cz = (
    729 * a**4 * c**4
    - 972 * a**3 * b * c**3
    + 2322 * a**3 * c**2
    + 54 * a**2 * b**3 * c**3
    + 270 * a**2 * b**2 * c**2
    - 3735 * a**2 * b * c
    - 338 * a**2
    - 36 * a * b**4 * c**2
    + 122 * a * b**3 * c
    + 1372 * a * b**2
    + b**6 * c**2
    - 2 * b**5 * c
    - 80 * b**4
)
Qz = 8 * W**3 + Az * W**2 + 6 * L * Bz * W + L * Cz
zero(Qz.subs({a: F1, b: F2, c: F3, W: z}, simultaneous=True))
zero(Az - 12 * L - 9 * (b**2 * p - 2 * q))
print("  Q_z(W) = 8 W^3 + A_z W^2 + 6 L B_z W + L C_z")
print("  A_z =", Az)
print("  B_z =", Bz)
print("  C_z =", Cz)
print("  Q_z(F(x,y,z); z) = 0 identically  [PASS]")
print("  lead(Q_z)=8 is constant, hence z is integral over QQ[a,b,c].")
print("  A_z - 12L = 9(b^2 p - 2q)  [PASS]")
print("  Therefore on L=0: Q_z(W)=W^2(8W+A_z).")
print()

print("=" * 78)
print("[3] Exact L=0 factorization in r and the escaping double shadow")
print("=" * 78)
r0 = 3 * a * (b - 9 * a * c) / q
r1 = 3 * (b**3 - 16 * a * b + 36 * a**2 * c) / (2 * q)
zero(
    q**3 * (Pr + 2 * (r - r0) ** 2 * (r - r1))
    - 27 * a**2 * (b**3 - 108 * a**2 * c - 6 * q * r) * L
)
zero(q**2 * Dr.subs(r, r0) - 27 * a**2 * L)
print("  r0 =", r0)
print("  r1 =", r1)
print("  modulo L, with q != 0: P_r=-2(r-r0)^2(r-r1)  [PASS]")
print("  q^2 D_r(r0)=27a^2 L, so D_r(r0)=0 on L  [PASS]")
print("  For a != 0, x D_r=r excludes r0 from the affine fibre:")
print("  the double r-root is the shadow of two x-sheets at infinity.")
print()

print("=" * 78)
print("[4] The finite survivor is a global rational section modulo L")
print("=" * 78)
xs = 2 * c / p
rs = 3 * c * q / (2 * p)
ys = b - rs
zs = -sp.Rational(9, 8) * (b**2 * p - 2 * q)
zero(p * q * (rs - r1) - 6 * b * L)
zero(p**2 * (Dr.subs(r, rs) - 3 * q / 4) - 12 * L)

sub_survivor = {x: xs, y: ys, z: zs}
H = 54 * a * c**2 - 27 * b * c + 28
d1 = sp.cancel(F1.subs(sub_survivor, simultaneous=True) - a)
d2 = sp.cancel(F2.subs(sub_survivor, simultaneous=True) - b)
d3 = sp.cancel(F3.subs(sub_survivor, simultaneous=True) - c)
zero((3 * b * c - 4) ** 3 * d1 - 4 * H * L)
zero((3 * b * c - 4) ** 2 * d2 + 36 * c * L)
zero(d3)
print("  (x_s,y_s,z_s) =")
print("   ", (xs, ys, zs))
print("  p q (r_s-r1)=6bL and p^2(D_r(r_s)-3q/4)=12L  [PASS]")
print("  global identities:")
print("    (3bc-4)^3(F1(x_s,y_s,z_s)-a)=4(54ac^2-27bc+28)L")
print("    (3bc-4)^2(F2(x_s,y_s,z_s)-b)=-36cL")
print("    F3(x_s,y_s,z_s)-c=0")
print("  Hence this is an affine fibre point at every L=0, p!=0 target.")
print("  On L, z_s=-A_z/8; the double W=0 roots are the two escape shadows.")
print()

print("=" * 78)
print("[5] Exceptional-locus audit")
print("=" * 78)
# p=0 is exactly the omitted curve once L=0.
zero(L.subs(c, 4 / (3 * b)) - q**2 / (3 * b**2))
zero(
    Pr.subs({a: b**2 / 12, c: 4 / (3 * b)}, simultaneous=True)
    + 2 * (r - b / 2) ** 3
)
print("  p=0: L=q^2/(3b^2), so L=p=0 is")
print("    (a,b,c)=(4/(27t^2),4/(3t),t), t!=0.")
print("    The x-core is the nonzero constant -2c: the fibre is empty.")

# q=0 on L splits into the z-axis and omitted curve.
zero(L.subs(a, b**2 / 12) - b**2 * p**2 / 48)
print("  q=0: L=b^2 p^2/48, hence the L-slice splits into")
print("    b=0 (the z-axis, one point (c/2,0,0)) or p=0 (empty).")
print("    P_r is triple on both pieces; a triple coordinate root is not")
print("    a three-point affine fibre.")

# a=0 is an independent coordinate collision.
zero(Pr.subs(a, 0) - r**2 * (3 * b - 2 * r))
zero(L.subs(a, 0) - b**2 * (b * c - 1))
print("  a=0: P_r=r^2(3b-2r), while L=b^2(bc-1).")
print("    Off L, r=0 is two finite points with the same y=b.")
print("    On bc=1, b!=0, those two sheets escape and the sole point is")
print("    (2/b,-b/2,9b^2/8).")
good_a0 = {x: 2 / b, y: -b / 2, z: 9 * b**2 / 8}
bad_a0 = {x: 2 / b, y: b, z: -9 * b**2 / 8}
for actual, expected in zip(
    (F1.subs(good_a0), F2.subs(good_a0), F3.subs(good_a0)),
    (0, b, 1 / b),
):
    zero(actual - expected)
for actual, expected in zip(
    (F1.subs(bad_a0), F2.subs(bad_a0), F3.subs(bad_a0)),
    (-3 * b**2 / 8, b / 4, 1 / b),
):
    zero(actual - expected)

# The z-coordinate triple-shadow locus on L.
a_ztriple = b**2 * (3 * b * c - 2) / 24
zero(
    L.subs(a, a_ztriple)
    - b**2 * p**2 * (9 * b**2 * c**2 + 12 * b * c - 28) / 192
)
print("  A_z=0 on L (so Q_z=8W^3) consists of the z-axis, the omitted")
print("    curve, and two curves 9(bc)^2+12bc-28=0 (b p != 0).")
print("    On the last two curves one finite point has z_s=0 and the other")
print("    two W=0 roots are escape shadows: again, not three fibre points.")

# Exact boundary and mechanism controls.
assert [sp.cancel(v) for v in (
    F1.subs({x: 2, y: sp.Rational(5, 6), z: -sp.Rational(7, 8)}),
    F2.subs({x: 2, y: sp.Rational(5, 6), z: -sp.Rational(7, 8)}),
    F3.subs({x: 2, y: sp.Rational(5, 6), z: -sp.Rational(7, 8)}),
)] == [sp.Rational(2, 27), 1, 1]
zero(
    Pr.subs({a: sp.Rational(2, 27), b: 1, c: 1})
    + 2 * (r - sp.Rational(2, 3)) ** 2 * (r - sp.Rational(1, 6))
)
zero(Qz.subs({a: sp.Rational(2, 27), b: 1, c: 1}) - W**2 * (8 * W + 7))
print("  generic control (2/27,1,1):")
print("    survivor (2,5/6,-7/8);")
print("    P_r=-2(r-2/3)^2(r-1/6), Q_z=W^2(8W+7).  [PASS]")

zero(L.subs({c: 0, a: b**2 / 16}, simultaneous=True))
boundary = (
    xs.subs({c: 0, a: b**2 / 16}, simultaneous=True),
    ys.subs({c: 0, a: b**2 / 16}, simultaneous=True),
    zs.subs({c: 0, a: b**2 / 16}, simultaneous=True),
)
for actual, expected in zip(boundary, (0, b, -sp.Rational(63, 16) * b**2)):
    zero(actual - expected)
print("  c=0 boundary control: a=b^2/16 gives (0,b,-63b^2/16). [PASS]")
print()

print("=" * 78)
print("[6] Why the double (r,z) root is at infinity")
print("=" * 78)
print("  The x-core is E=Lx^3+px-2c.  At generic L=0,p!=0 two roots")
print("  satisfy x ~ +/-sqrt(-p/L), whereas y stays bounded by P_r.")
print("  From F3=c, z=2/x^2-3y/x-c/x^3 -> 0.")
print("  Thus (r,z)->(r0,0) records two compactified sheets, not two")
print("  finite affine points.  Only (r1,-A_z/8) is the survivor.")
print()

print("=" * 78)
print("[7] Pointwise parity rigidity (logical consequence, unchanged)")
print("=" * 78)
print("  For a real etale polynomial map, conjugation pairs every nonreal")
print("  point of a full complex fibre.  Hence #real = complex degree (mod 2).")
print("  An even-degree wild map cannot have an odd-real-count FULL fibre;")
print("  HYP-9030 test (ii) remains released as stated.")
print()
print("all exact checks passed")
