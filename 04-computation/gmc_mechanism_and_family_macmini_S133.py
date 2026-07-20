#!/usr/bin/env python3
"""
The GMC(4) counterexample IS the alternating binomial identity   (mac-mini-S133)
================================================================================
Follow-up to gmc_counterexample_verify_macmini_S133.py, which verified the claim
exactly.  Two things went wrong / came right in that run and both are pursued here.

CORRECTION TO MY OWN TEST.  That script flagged the degree-1 part of P as "still a
counterexample" because E[P^m] = 0 and E[QP^m] was nonzero SOMEWHERE.  That is the
WRONG criterion.  A Mathieu-Zhao subspace requires E[QP^m] = 0 only for m >> 0, so a
genuine counterexample needs E[QP^m] != 0 for INFINITELY MANY m.  The degree-1 part
gives [1,0,0,0,0] -- it is eventually zero, hence NOT a counterexample; it satisfies
the conjecture.  Everything below uses the correct "not eventually zero" criterion.

THE MECHANISM.  Writing the surviving reading as
    P = (1+W)(Zbar(1-Z) + Wbar),   Q = W
and pairing Gaussian-style (E[Z^a Zbar^b] = delta_ab a!), the two claims reduce to
ONE elementary identity.  Using C(m,k)k!(m-k)! = m! throughout:

    E[P^m]     = m! * sum_{j=0}^{m}   (-1)^j C(m,j)         = m! * 0    = 0
    E[Q P^m]   = m! * sum_{j=1}^{m} (-1)^{j-1} C(m,j)       = m! * C(m,0) = m!

So BOTH sides are the SAME alternating binomial sum, and multiplying by Q merely
SHIFTS THE SUMMATION INDEX BY ONE.  The shift exposes the boundary term C(m,0) = 1
that the complete alternating sum had cancelled.  That is the entire counterexample:
  * the factor (1-Z) manufactures the alternating signs,
  * the Gaussian pairing manufactures the factorials that cancel the binomial
    denominators exactly,
  * Q shifts the index and cashes in the boundary term.

CONSEQUENCES, all derived and then verified below:
  (F1)  Q = Z2^r  gives  E[Q P^m] = (-1)^{r+1} m! sum_{j=0}^{r-1} (-1)^j C(m,j)
        -- an infinite FAMILY, growth m! times a degree-(r-1) polynomial in m.
        r=1 is the published case; r=2 gives m!(m-1).
  (F2)  Deforming (1-Z) -> (1-lambda Z) gives E[P^m] = m! (1-lambda)^m, so
        lambda = 1 is FORCED EXACTLY.  The construction is rigid.
  (F3)  minimality, re-tested with the correct criterion.
"""
from math import factorial, comb
from itertools import product
from fractions import Fraction

def pmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2], k1[3]+k2[3])
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def padd(*ps):
    out = {}
    for p in ps:
        for k, c in p.items(): out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}

def pneg(p): return {k: -c for k, c in p.items()}
def pscal(p, s): return {k: c*s for k, c in p.items() if c*s}

def expect(p):
    return sum(c * factorial(a) * factorial(cc)
               for (a, b, cc, d), c in p.items() if a == b and cc == d)

Z    = {(1,0,0,0): 1}; ZB = {(0,1,0,0): 1}
W    = {(0,0,1,0): 1}; WB = {(0,0,0,1): 1}
ONE  = {(0,0,0,0): 1}

def makeP(lam=1):
    """P = (1+W)( Zbar(1 - lam Z) + Wbar )"""
    return pmul(padd(ONE, W), padd(pmul(ZB, padd(ONE, pscal(Z, -lam))), WB))

MM = 9
print("=" * 78)
print("PART 1 -- the closed forms, derived by hand and checked against brute force")
print("=" * 78)
P = makeP()
print(f"{'m':>3} {'E[P^m]':>10} {'m!*sum(-1)^j C(m,j)':>22} {'E[QP^m]':>12} {'m!':>10} {'ok':>5}")
Pm = ONE; allok = True
for m in range(1, MM + 1):
    Pm = pmul(Pm, P)
    lhs = expect(Pm)
    pred_l = factorial(m) * sum((-1)**j * comb(m, j) for j in range(m + 1))
    rhs = expect(pmul(W, Pm))
    ok = (lhs == pred_l == 0) and (rhs == factorial(m))
    allok &= ok
    print(f"{m:>3} {lhs:>10} {pred_l:>22} {rhs:>12} {factorial(m):>10} {str(ok):>5}")
print(f"  closed forms hold for m = 1..{MM}: {allok}")

print()
print("=" * 78)
print("PART 2 -- (F1) the family Q = Z2^r :  m! times a degree-(r-1) polynomial")
print("=" * 78)
print(f"{'r':>3} " + " ".join(f"{'m='+str(m):>12}" for m in range(1, 7)) + "   closed form")
for r in range(1, 5):
    Qr = ONE
    for _ in range(r): Qr = pmul(Qr, W)
    row = []; ok = True
    Pm = ONE
    for m in range(1, 7):
        Pm = pmul(Pm, P)
        got = expect(pmul(Qr, Pm))
        pred = (-1)**(r+1) * factorial(m) * sum((-1)**j * comb(m, j) for j in range(r))
        if got != pred: ok = False
        row.append(got)
    label = {1: "m!", 2: "m!(m-1)", 3: "m!(1-m+C(m,2))", 4: "m!(...deg 3...)"}[r]
    print(f"{r:>3} " + " ".join(f"{v:>12}" for v in row) + f"   {label}  match={ok}")
print("  So the published counterexample is the r=1 member of an INFINITE FAMILY;")
print("  every r >= 1 breaks GMC, with growth m! * (degree r-1 polynomial in m).")

print()
print("=" * 78)
print("PART 3 -- (F2) lambda-rigidity: E[P^m] = m! (1-lambda)^m, so lambda=1 is FORCED")
print("=" * 78)
print(f"{'lambda':>7} " + " ".join(f"{'m='+str(m):>10}" for m in range(1, 6)))
for lam in (0, 1, 2, 3, -1):
    Pl = makeP(lam); Pm = ONE; row = []
    for m in range(1, 6):
        Pm = pmul(Pm, Pl); row.append(expect(Pm))
    pred = [factorial(m) * (1 - lam)**m for m in range(1, 6)]
    print(f"{lam:>7} " + " ".join(f"{v:>10}" for v in row) +
          f"   predicted {pred}  match={row == pred}")
print("  Only lambda = 1 kills the left-hand side.  The construction is RIGID -- it sits")
print("  at the unique point where the alternating binomial sum vanishes.")

print()
print("=" * 78)
print("PART 4 -- (F3) minimality, with the CORRECT 'not eventually zero' criterion")
print("=" * 78)
def is_counterexample(Pp, Qq, M=10):
    """E[P^m]=0 for all tested m, AND E[QP^m] != 0 at the LARGEST tested m."""
    Pm = ONE; tail = []
    for m in range(1, M + 1):
        Pm = pmul(Pm, Pp)
        if expect(Pm) != 0: return False, None
        tail.append(expect(pmul(Qq, Pm)))
    return (tail[-1] != 0 and tail[-2] != 0), tail

ok, tail = is_counterexample(P, W)
print(f"  full P:  counterexample = {ok},  E[QP^m] tail = {tail}")

deg = {}
for k, c in P.items(): deg.setdefault(sum(k), {})[k] = c
for d in sorted(deg):
    okd, td = is_counterexample(deg[d], W)
    print(f"  degree-{d} part alone: counterexample = {okd}   E[QP^m] = {td}"
          f"{'   <-- eventually ZERO, so it SATISFIES GMC' if td and not okd else ''}")

drops = 0
for k in list(P):
    Pd = {kk: c for kk, c in P.items() if kk != k}
    okk, _ = is_counterexample(Pd, W)
    if okk: drops += 1
print(f"  dropping any ONE of the 6 terms still a counterexample in {drops} cases"
      f"  => {'NOT term-minimal' if drops else 'TERM-MINIMAL'}")

print()
print("=" * 78)
print("PART 5 -- how few REAL Gaussians are needed?  exhaustive small search")
print("=" * 78)
def dfac2(a):
    """E[X^a] for X ~ N(0,1): (a-1)!! if a even else 0."""
    if a % 2: return 0
    return factorial(a) // (2**(a//2) * factorial(a//2))

def rmul(p, q, n):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(k1[i]+k2[i] for i in range(n))
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def rexpect(p):
    tot = 0
    for k, c in p.items():
        v = 1
        for a in k:
            v *= dfac2(a)
            if v == 0: break
        tot += c*v
    return tot

def monos(n, dmax):
    out = []
    for d in range(1, dmax + 1):
        for k in product(range(d + 1), repeat=n):
            if sum(k) == d: out.append(k)
    return out

for n, dmax in ((1, 3), (2, 3), (3, 2)):
    ms = monos(n, dmax)
    RONE = {tuple([0]*n): 1}
    found = []
    for coeffs in product((-1, 0, 1), repeat=len(ms)):
        Pp = {m: c for m, c in zip(ms, coeffs) if c}
        if not Pp: continue
        # E[P^m] = 0 for m = 1..6 ?
        Pm = RONE; good = True
        pw = []
        for m in range(1, 7):
            Pm = rmul(Pm, Pp, n)
            if rexpect(Pm) != 0: good = False; break
            pw.append(dict(Pm))
        if not good: continue
        # some monomial Q with E[QP^m] != 0 at the top two m ?
        for Q in ms:
            vals = [rexpect(rmul({Q: 1}, x, n)) for x in pw]
            if vals[-1] != 0 and vals[-2] != 0:
                found.append((Pp, Q, vals)); break
        if len(found) >= 3: break
    print(f"  n={n} real Gaussians, deg<={dmax}, coeffs in {{-1,0,1}} "
          f"({3**len(ms)} polys): counterexamples found = {len(found)}")
    for Pp, Q, vals in found[:3]:
        print(f"      P={Pp}  Q={Q}  E[QP^m]={vals}")

print()
print("SUMMARY")
print("  The counterexample is the alternating binomial identity in disguise.  Q shifts")
print("  the summation index by one and cashes the boundary term C(m,0)=1, turning a sum")
print("  that vanishes into one that equals m!.  It is rigid in lambda, it is the r=1")
print("  member of an infinite family, and the degree-1 piece alone is NOT a")
print("  counterexample (it is eventually zero, so it satisfies the conjecture).")
