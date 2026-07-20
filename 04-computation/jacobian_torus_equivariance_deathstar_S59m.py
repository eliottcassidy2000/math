#!/usr/bin/env python3
"""
death-star-2026-07-19-S59m (HYP-8075) -- the C*-equivariant structure of the
JC counterexample.

Claim 1 (weighted homogeneity): under weights w(x,y,z) = (1,-1,-2), the
components are homogeneous: w(F1) = -2, w(F2) = -1, w(F3) = +1.  Hence
F(lam.x, lam^-1 y, lam^-2 z) = (lam^-2 F1, lam^-1 F2, lam F3): C*-equivariant,
and kind-pasteur's sigma = (-x,-y,z) with tau = (a,-b,-c) is the lam = -1
element (lam^-2 = 1, lam^-1 = -1, lam = -1 ... note: tau acts as (+,-,-)).

Claim 2 (the collision is a torus phenomenon): the fiber over the a-axis
consists of the FIXED branch {x = y = 0} (F(0,0,z) = (z,0,0), bijective onto
the axis) plus the ORBIT branch C* . (1, -3/2, 13/2), which maps 2:1 onto the
punctured axis: F(lam, -3/(2 lam), 13/(2 lam^2)) = (-1/(4 lam^2), 0, 0).
At lam = +-1 the two orbit points land on a = -1/4 together with the fixed
point (0,0,-1/4): the verified triple collision = 1 (fixed) + 2 (orbit).
Exact Laurent arithmetic in lam, no dependencies.
"""
from fractions import Fraction

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(x+y for x, y in zip(ka, kb))
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}

def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c != 0}

def pscale(p, s):
    return {k: c*s for k, c in p.items() if c*s != 0}

X = {(1,0,0): 1}; Y = {(0,1,0): 1}; Z = {(0,0,1): 1}; ONE = {(0,0,0): 1}
U = padd(ONE, pmul(X, Y)); U2 = pmul(U, U); U3 = pmul(U2, U)
W = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))
F1 = padd(pmul(U3, Z), pmul(pmul(pmul(Y, Y), U), W))
F2 = padd(Y, pscale(pmul(pmul(X, U2), Z), 3), pscale(pmul(pmul(X, pmul(Y, Y)), W), 3))
F3 = padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3), pscale(pmul(pmul(pmul(X, X), X), Z), -1))

WT = (1, -1, -2)
def weights(p):
    return sorted(set(sum(w*e for w, e in zip(WT, k)) for k in p))
print("Claim 1 -- weight sets under w(x,y,z) = (1,-1,-2):")
for nm, p, want in [("F1", F1, -2), ("F2", F2, -1), ("F3", F3, 1)]:
    ws = weights(p)
    print(f"  {nm}: weights {ws}  homogeneous of weight {want}: {ws == [want]}")

# Claim 2: Laurent evaluation at (lam, -3/(2 lam), 13/(2 lam^2)); key = lam-exponent
def eval_laurent(p):
    r = {}
    for (i, j, k), c in p.items():
        coef = Fraction(c) * Fraction(-3, 2)**j * Fraction(13, 2)**k
        e = i - j - 2*k
        r[e] = r.get(e, Fraction(0)) + coef
    return {e: c for e, c in r.items() if c != 0}

img = [eval_laurent(p) for p in (F1, F2, F3)]
print("\nClaim 2 -- F(lam, -3/(2 lam), 13/(2 lam^2)) as Laurent polynomials in lam:")
print("  F1 =", img[0], " (claim: {-2: -1/4})")
print("  F2 =", img[1], " (claim: 0)")
print("  F3 =", img[2], " (claim: 0)")
ok = img[0] == {-2: Fraction(-1, 4)} and img[1] == {} and img[2] == {}
print("  orbit-branch law F = (-1/(4 lam^2), 0, 0):", ok)
print("\nfixed branch: F(0,0,z) = (z, 0, 0) -- read off (x=y=0 => u=1).")
print("VERDICT: collision = 1 fixed-branch + 2 orbit-branch points; the fiber")
print("over every point of the punctured a-axis contains exactly this 1+2")
print("pattern from these two branches; geometric degree >= 3.")
