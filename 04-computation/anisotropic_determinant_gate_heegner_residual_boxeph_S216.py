#!/usr/bin/env python3
"""Audit S216's determinant/Heegner comparison without transferring structures.

What is exact and useful:

* In the chosen saturated basis ``u,z`` of a THM-2053 plane, the certificate

      max_i |a*z_i-b*u_i| <= (a^2+b^2)/91

  has a determinant on the left and the Euclidean coordinate norm on the
  right.  Failure is exactly membership in a union of 2n open tangent disks.
* Classical reduced-form enumeration gives the displayed imaginary-quadratic
  class numbers.

What this script deliberately does *not* infer: the Euclidean form has
discriminant -4, not -7; the 2n disks do not form one quadratic-form class;
Paley forms have no supplied map to the THM-2053 plane; local anisotropy does
not select a rank/Euler branch or imply loneliness; class number one does not
compress the LRC(14) atlas.
"""

from math import gcd, isqrt


def sep(title):
    print("\n" + "=" * 72 + "\n" + title + "\n" + "=" * 72)


# Thirteen columns c_i=(u_i,z_i).  The first two columns have determinant 1,
# so the displayed u,z generate a saturated rank-two lattice.  Also u+z has
# full positive support, as required by the LRC application of THM-2053.
COLS = [
    (1, 0), (0, 1), (1, 1), (2, 1), (1, 2), (3, 1), (1, 3),
    (2, 3), (3, 2), (1, 4), (4, 1), (5, 3), (2, 5),
]


def det_value(a, b, u, z):
    """det((a,b),(u,z)) in the convention used by THM-2053."""
    return a * z - b * u


def gate_fails(a, b):
    norm = a * a + b * b
    return any(91 * abs(det_value(a, b, u, z)) > norm for u, z in COLS)


def in_tangent_disk_union(a, b):
    """Integer-only version of THM-2053 equation (13)."""
    for u, z in COLS:
        radius4 = 91 * 91 * (u * u + z * z)
        for sigma in (-1, 1):
            # Four times the squared distance from
            # (a,b) to (91*sigma/2)*(z,-u).
            distance4 = (2 * a - 91 * sigma * z) ** 2 + (2 * b + 91 * sigma * u) ** 2
            if distance4 < radius4:
                return True
    return False


sep("A  THM-2053: a determinant certificate in chosen Euclidean coordinates")
print("  n=13; det(c_1,c_2)=1, so the displayed basis is saturated; u+z has positive support.")
print("  certificate: 91*max_i|a*z_i-b*u_i| <= a^2+b^2  ==>  M(a*u+b*z) >= 1/14")
print("  a^2+b^2 is the chosen-coordinate Euclidean form [1,0,1], discriminant -4.")
print("  It is not a plane-intrinsic form of discriminant -7 and changes under a GL(2,Z) basis change.")
for a, b in ((30, 1), (200, 3), (700, 700), (2000, 1)):
    det_max = max(abs(det_value(a, b, u, z)) for u, z in COLS)
    norm = a * a + b * b
    print(
        f"  d=({a:4d},{b:4d}): max|det|={det_max:5d}; norm/91={norm/91:9.1f}; "
        f"certified={91*det_max <= norm}"
    )


sep("B  Exact failure carrier: a union of 26 open tangent disks")
max_col_norm_sq = max(u * u + z * z for u, z in COLS)
envelope_sq = 91 * 91 * max_col_norm_sq
bound = isqrt(envelope_sq - 1)
fail_count = 0
tested = 0
for a in range(-bound, bound + 1):
    for b in range(-bound, bound + 1):
        norm = a * a + b * b
        if (a, b) == (0, 0) or gcd(abs(a), abs(b)) != 1 or norm >= envelope_sq:
            continue
        tested += 1
        fails = gate_fails(a, b)
        in_disks = in_tangent_disk_union(a, b)
        assert fails == in_disks
        fail_count += int(fails)
print(f"  columns={len(COLS)} give 2n={2*len(COLS)} disks; every disk has the origin on its boundary.")
print(f"  exact round envelope: a^2+b^2 < 91^2*{max_col_norm_sq} = {envelope_sq}.")
print(f"  checked all {tested} primitive nonzero lattice points in that envelope.")
print(f"  uncertified points={fail_count}, or {fail_count//2} directions modulo d~-d; disk/gate tests agree exactly.")
print("  'Uncertified' is not 'unsafe': THM-2053 gives a sufficient certificate, not an equivalence.")


def is_reduced(a, b, c):
    return (
        abs(b) <= a <= c
        and not (abs(b) == a and b < 0)
        and not (a == c and b < 0)
    )


def class_number(discriminant):
    """Count reduced primitive positive-definite forms of a negative discriminant."""
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    count = 0
    a = 1
    while 3 * a * a <= abs(discriminant):
        for b in range(-a, a + 1):
            numerator = b * b - discriminant
            if numerator % (4 * a):
                continue
            c = numerator // (4 * a)
            if gcd(gcd(a, abs(b)), c) == 1 and is_reduced(a, b, c):
                count += 1
        a += 1
    return count


sep("C  Independent classical check: selected imaginary-quadratic class numbers")
for discriminant in (-3, -4, -7, -8, -11, -15, -19, -20, -23, -43, -67, -163):
    print(f"  discriminant {discriminant:4d}: class number {class_number(discriminant)}")
print("  In particular h(-7)=1.  This classifies forms of discriminant -7 only.")
print("  The THM-2053 coordinate norm above has discriminant -4, and its 26 translated disks are not")
print("  ideal classes or equivalence classes of binary quadratic forms.  No atlas compression follows.")


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


sep("D  Legendre-symbol audit: (-p/p)=0, while (-1/p) is a different symbol")
for p in (3, 5, 7, 11, 13):
    print(f"  p={p:2d}: (-p/p)={legendre(-p,p):+d}; (-1/p)={legendre(-1,p):+d}")
print("  Because p divides the discriminant -p, (disc/p)=0, not +/-1.")
print("  The separate identity (-1/p)=(-1)^((p-1)/2) records p mod 4; it is not (disc/p).")


sep("E  Scope of the surviving observation")
print("  SURVIVES: determinant identity; quadratic-vs-linear safe tail; exact 26-disk failure carrier;")
print("            finite primitive residual for each chosen basis; classical class-number table.")
print("  DOES NOT FOLLOW: an LRC14 discriminant -7; a Paley-to-plane identification; rank= isotropic;")
print("                   Euler=anisotropic; anisotropic=>lonely; or class-number-one atlas compression.")
print("  Remaining LRC14 work is still to construct/compress the actual plane atlas and discharge each")
print("  plane's exact tangent-disk lattice points by a valid loneliness argument.")
