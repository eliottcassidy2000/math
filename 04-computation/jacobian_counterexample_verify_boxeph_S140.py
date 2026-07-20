#!/usr/bin/env python3
"""
jacobian_counterexample_verify_boxeph_S140.py  (HYP-8070)

INDEPENDENT EXACT VERIFICATION of the reported Jacobian-conjecture counterexample
(owner communication, 2026-07-19), plus its mod-p spectrum through this project's
instruments.

THE MAP  F : C^3 -> C^3:
  f1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
  f2 = y + 3x(1+xy)^2 z + 3x y^2 (4+3xy)
  f3 = 2x - 3x^2 y - x^3 z

CLAIMS TO VERIFY (each sufficient data is fully checkable in exact arithmetic):
 (1) det JF = -2 identically (symbolic, exact rational polynomial arithmetic);
 (2) F(0,0,-1/4) = F(1,-3/2,13/2) = F(-1,3/2,13/2) = (-1/4, 0, 0);
 (1)+(2) together: a Keller map that is not injective => JC_3 is FALSE (and by
 stabilization JC_n false for all n >= 3; by DC_n => JC_n, DC_n false for n >= 3).
 (3) the Z/2 equivariance F o sigma = tau o F with sigma = diag(-1,-1,1),
     tau = diag(1,-1,-1) -- the collision fiber is sigma-organized as
     {fixed point} u {2-orbit}: 3 = 1 + 2.
 (4) mod-p spectrum (our gate/rung instrument aimed at the map): for small p,
     full enumeration of F : F_p^3 -> F_p^3 -- image deficiency, fiber-size
     histogram, generic degree signature.  det = -2 => the only bad prime is 2.

Pure Python, exact integers/Fractions.  boxeph-2026-07-19-S140.
"""

from fractions import Fraction as Fr
from itertools import product

# ---------------- exact polynomial arithmetic in Q[x,y,z] ----------------------

def pmul(a, b):
    out = {}
    for (i, j, k), c in a.items():
        for (l, m, n), d in b.items():
            key = (i + l, j + m, k + n)
            out[key] = out.get(key, 0) + c * d
    return {k: v for k, v in out.items() if v != 0}

def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}

def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c != 0}

def pdiff(p, var):
    out = {}
    for (i, j, k), c in p.items():
        e = (i, j, k)[var]
        if e:
            key = list((i, j, k)); key[var] -= 1
            out[tuple(key)] = out.get(tuple(key), 0) + c * e
    return out

def peval(p, x, y, z):
    return sum(c * x**i * y**j * z**k for (i, j, k), c in p.items())

X, Y, Z, ONE = {(1,0,0): 1}, {(0,1,0): 1}, {(0,0,1): 1}, {(0,0,0): 1}
XY1 = padd(ONE, pmul(X, Y))                     # 1+xy
A43 = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))   # 4+3xy

f1 = padd(pmul(pmul(pmul(XY1, XY1), XY1), Z),
          pmul(pmul(pmul(Y, Y), XY1), A43))
f2 = padd(Y, pscale(pmul(pmul(X, pmul(XY1, XY1)), Z), 3),
          pscale(pmul(pmul(X, pmul(Y, Y)), A43), 3))
f3 = padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3),
          pscale(pmul(pmul(pmul(X, X), X), Z), -1))
F = [f1, f2, f3]

# ---------------- (1) det JF symbolically ---------------------------------------
J = [[pdiff(f, v) for v in range(3)] for f in F]
det = padd(
    pmul(J[0][0], padd(pmul(J[1][1], J[2][2]), pscale(pmul(J[1][2], J[2][1]), -1))),
    pscale(pmul(J[0][1], padd(pmul(J[1][0], J[2][2]), pscale(pmul(J[1][2], J[2][0]), -1))), -1),
    pmul(J[0][2], padd(pmul(J[1][0], J[2][1]), pscale(pmul(J[1][1], J[2][0]), -1))))
print("(1) det JF =", det, "  -> CONSTANT -2:", det == {(0, 0, 0): -2})
assert det == {(0, 0, 0): -2}, "det JF is NOT identically -2!"

# ---------------- (2) the three-point collision ----------------------------------
pts = [(Fr(0), Fr(0), Fr(-1, 4)), (Fr(1), Fr(-3, 2), Fr(13, 2)),
       (Fr(-1), Fr(3, 2), Fr(13, 2))]
imgs = [tuple(peval(f, *p) for f in F) for p in pts]
print("(2) images:", imgs)
assert imgs[0] == imgs[1] == imgs[2] == (Fr(-1, 4), Fr(0), Fr(0)), "collision FAILS"
print("    THREE distinct preimages of (-1/4, 0, 0)  => F NOT injective.")
print("    (1)+(2): a non-injective Keller map exists  =>  JC_3 IS FALSE.")

# ---------------- (3) equivariance -----------------------------------------------
def sflip(p):   # substitute (x,y,z) -> (-x,-y,z)
    return {k: (v if (k[0] + k[1]) % 2 == 0 else -v) for k, v in p.items()}
eq = (sflip(f1) == f1) and (sflip(f2) == pscale(f2, -1)) and (sflip(f3) == pscale(f3, -1))
print("(3) F o diag(-1,-1,1) = diag(1,-1,-1) o F :", eq)
assert eq
print("    collision fiber = sigma-fixed point (0,0,-1/4) + sigma-2-orbit: 3 = 1 + 2.")

# ---------------- (4) mod-p spectrum ----------------------------------------------
print("\n(4) mod-p spectrum (full enumeration of F on F_p^3):")
print("%-5s %-10s %-12s %-10s %s" % ("p", "|image|", "deficiency", "def/p^2", "fiber-size histogram"))
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29):
    # integer-coefficient versions of f1..f3 mod p
    fs = [{k: int(v) % p for k, v in f.items()} for f in F]
    cnt = {}
    for x, y, z in product(range(p), repeat=3):
        img = tuple(sum(c * pow(x, i, p) * pow(y, j, p) * pow(z, k, p) for (i, j, k), c in f.items()) % p
                    for f in fs)
        cnt[img] = cnt.get(img, 0) + 1
    hist = {}
    for v in cnt.values():
        hist[v] = hist.get(v, 0) + 1
    imgsz = len(cnt)
    defc = p**3 - imgsz
    print("%-5d %-10d %-12d %-10.3f %s" % (p, imgsz, defc, defc / p**2,
                                           dict(sorted(hist.items()))))
print("\nreading: constant det=-2 => etale away from p=2; fiber sizes bounded by the")
print("geometric degree; the collision point forces generic degree >= 3.  The")
print("deficiency growth and the dominant fiber size are the map's arithmetic")
print("fingerprint -- the same instrument as our improperness-mod-p ladders.")
print("DONE.")
