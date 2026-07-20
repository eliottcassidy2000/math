#!/usr/bin/env python3
"""
GMC(3) verified, and the MASTER THEOREM behind both counterexamples  (mac-mini-S134)
====================================================================================
Owner supplied a claimed GMC(3) counterexample (5-term quartic, 3 real Gaussians):

    Z = (X+iY)/sqrt2,  W = Zbar,  U = T^2/2,
    P = (1+Z)(W - (2+Z)U),   Q = Z,     E[P^m] = 0,  E[Q P^m] = m!.

PART A verifies it exactly, twice (complex formalism + real coordinates in a,b,T).

PARTS B-D then answer the question I asked at the end of S133 -- "is the alternating
binomial collapse the ONLY mechanism?" -- with a flat NO, and replace it by something
much better.  Set the construction up in general:

    P = (1+Z)(W - g(Z) U),  Q = Z,   U >= 0 independent of the complex pair (Z,W).

Using E[Z^a W^b] = a! delta_ab, i.e. E[W^r F(Z)] = r! [s^r] F(s), and writing
Phi(x) = sum_k E[U^k] x^k / k!  for the MOMENT GENERATING FUNCTION of U:

    E[P^m] = m! [s^m] (1+s)^m * Phi(-s g(s))          <-- ONE line, all m at once

So everything is controlled by the single series  f(s) := Phi(-s g(s)).

  (B) f is FORCED.  Requiring [s^m](1+s)^m f(s) = 0 for every m >= 1 pins f uniquely:
      by Lagrange inversion, sum_m t^m [s^m] f(s)(1+s)^m = f(t/(1-t))/(1-t), so the
      requirement is f(t/(1-t)) = f(0)(1-t), i.e.  f(s) = f(0)/(1+s).  Nothing else works.
      And then E[Q P^m] = m! [s^m] s (1+s)^{m-1} = m! AUTOMATICALLY -- the m! is FORCED,
      not a coincidence of either example.

  (C) THE MASTER EQUATION is therefore   Phi(-s g(s)) = 1/(1+s).
      For U = chi^2_d / 2 (i.e. d real Gaussians), Phi(x) = (1-x)^{-d/2}, so
          (1 + s g(s))^{-d/2} = (1+s)^{-1}   =>   1 + s g(s) = (1+s)^{2/d}
          =>   g(s) = ((1+s)^{2/d} - 1)/s.
      d=1 gives g = 2+s   -- EXACTLY the owner's (2+Z), now derived rather than guessed.
      d=2 gives g = 1     -- EXACTLY the S133 four-term example P = (1+W)(Wbar - |Z|^2).
      The two known counterexamples are ONE construction at d = 1 and d = 2.

  (D) MINIMALITY.  g is a POLYNOMIAL iff 2/d is a positive integer, i.e. d in {1,2}.
      Total real dimension is 2 (for Z,W) + d, so n in {3,4} and n=3 is MINIMAL for the
      family.  n = 2 would need d = 0, i.e. U constant, giving Phi(x) = e^{cx} and
      -c s g(s) = -log(1+s): NOT a polynomial.  So THIS FAMILY CANNOT REACH GMC(2).
"""
from fractions import Fraction as F
from math import factorial

# ================================================================= exact helpers
def poch(a, k):
    """rising factorial (a)_k with a rational."""
    r = F(1)
    for i in range(k): r *= (a + i)
    return r

# ---- complex formalism: monomials Z^i W^j U^k, coefficients rational -----------
def cmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2])
            out[k] = out.get(k, F(0)) + c1*c2
    return {k: c for k, c in out.items() if c}

def cadd(*ps):
    out = {}
    for p in ps:
        for k, c in p.items(): out[k] = out.get(k, F(0)) + c
    return {k: c for k, c in out.items() if c}

def cexp(p, dU):
    """E[Z^i W^j U^k] = i! * delta_ij * (d/2)_k   for U = chi^2_d / 2."""
    tot = F(0)
    for (i, j, k), c in p.items():
        if i == j: tot += c * factorial(i) * poch(F(dU, 2), k)
    return tot

Z = {(1,0,0): F(1)}; W = {(0,1,0): F(1)}; U = {(0,0,1): F(1)}
ONE = {(0,0,0): F(1)}
def scal(p, s): return {k: c*F(s) for k, c in p.items() if c*F(s)}

def build_P(g):
    """P = (1+Z)(W - g(Z) U);  g given as a list of coeffs [g0,g1,...] in Z."""
    gz = {}
    for i, gi in enumerate(g):
        if gi: gz[(i,0,0)] = F(gi)
    return cmul(cadd(ONE, Z), cadd(W, scal(cmul(gz, U), -1)))

print("=" * 78)
print("PART A -- exact verification of the GMC(3) counterexample")
print("=" * 78)
P3 = build_P([2, 1])                       # g(Z) = 2 + Z,  U = T^2/2  (d = 1)
print(f"  P has {len(P3)} terms; monomials (i,j,k) = Z^i W^j U^k:")
for k, c in sorted(P3.items()): print(f"      Z^{k[0]} W^{k[1]} U^{k[2]}  coeff {c}")
print()
print(f"{'m':>3} {'E[P^m]':>10} {'E[QP^m]':>12} {'m!':>12} {'ok':>5}")
Pm = ONE; allok = True
for m in range(1, 11):
    Pm = cmul(Pm, P3)
    lhs = cexp(Pm, 1); rhs = cexp(cmul(Z, Pm), 1)
    ok = (lhs == 0 and rhs == factorial(m)); allok &= ok
    print(f"{m:>3} {str(lhs):>10} {str(rhs):>12} {factorial(m):>12} {str(ok):>5}")
print(f"  GMC(3) counterexample CONFIRMED for m = 1..10: {allok}")

# ---- real-coordinate cross-check: a,b ~ N(0,1/2), T ~ N(0,1) -------------------
print()
print("  Real-coordinate cross-check.  Z = a+ib, W = a-ib with a,b ~ N(0,1/2)")
print("  (so a = X/sqrt2, b = Y/sqrt2 and ZW = a^2+b^2 ~ Exp(1)), U = T^2/2.")
def rmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2])
            re = c1[0]*c2[0] - c1[1]*c2[1]; im = c1[0]*c2[1] + c1[1]*c2[0]
            o = out.get(k, (F(0), F(0))); out[k] = (o[0]+re, o[1]+im)
    return {k: v for k, v in out.items() if v[0] or v[1]}

def radd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            o = out.get(k, (F(0), F(0))); out[k] = (o[0]+v[0], o[1]+v[1])
    return {k: v for k, v in out.items() if v[0] or v[1]}

def dfac(a):                       # (a-1)!! for even a, else 0
    if a % 2: return 0
    return factorial(a) // (2**(a//2) * factorial(a//2))

def rexp(p):
    """a,b ~ N(0,1/2): E[a^k] = (k-1)!! (1/2)^{k/2};  T ~ N(0,1): E[T^k] = (k-1)!!"""
    tre = tim = F(0)
    for (ka, kb, kt), (re, im) in p.items():
        if ka % 2 or kb % 2 or kt % 2: continue
        v = F(dfac(ka), 2**(ka//2)) * F(dfac(kb), 2**(kb//2)) * dfac(kt)
        tre += re*v; tim += im*v
    return tre, tim

RZ = {(1,0,0): (F(1), F(0)), (0,1,0): (F(0), F(1))}      # a + i b
RW = {(1,0,0): (F(1), F(0)), (0,1,0): (F(0), F(-1))}     # a - i b
RU = {(0,0,2): (F(1,2), F(0))}                            # T^2 / 2
RONE = {(0,0,0): (F(1), F(0))}
rneg = lambda p: {k: (-v[0], -v[1]) for k, v in p.items()}
# P = (1+Z)(W - (2+Z)U)
RP = rmul(radd(RONE, RZ), radd(RW, rneg(rmul(radd({k:(v[0]*2,v[1]*2) for k,v in RONE.items()}, RZ), RU))))
print(f"  expanded in (a,b,T): {len(RP)} monomials")
Pm = RONE; ok2 = True
for m in range(1, 7):
    Pm = rmul(Pm, RP)
    lre, lim = rexp(Pm); rre, rim = rexp(rmul(RZ, Pm))
    good = (lre == 0 and lim == 0 and rre == factorial(m) and rim == 0)
    ok2 &= good
    print(f"    m={m}: E[P^m] = {lre}+{lim}i   E[QP^m] = {rre}+{rim}i   (m! = {factorial(m)})  {good}")
print(f"  real-coordinate cross-check agrees: {ok2}")
print("  => genuinely THREE independent real standard Gaussians.")

# ================================================================= PART B
print()
print("=" * 78)
print("PART B -- f(s) = 1/(1+s) is FORCED, and so is the m!")
print("=" * 78)
print("  Requirement: [s^m](1+s)^m f(s) = 0 for every m >= 1.")
print("  Lagrange inversion with phi(s)=1+s (so sigma = t/(1-t), phi'=1):")
print("      sum_m t^m [s^m] f(s)(1+s)^m = f(t/(1-t))/(1-t).")
print("  Vanishing for all m>=1 forces f(t/(1-t))/(1-t) = f(0), i.e. with u=t/(1-t),")
print("      f(u) = f(0)/(1+u).    UNIQUE.")
print()
def coeff_series(fcoef, m, r=None):
    """[s^r] (1+s)^m * f(s);  r defaults to m.  (Independent exponents -- the earlier
    version tied r to m and mis-reported the Q column.)"""
    from math import comb
    if r is None: r = m
    if r < 0: return F(0)
    return sum(fcoef[j] * comb(m, r-j) for j in range(min(len(fcoef), r+1)))
NF = 14
f_inv = [F((-1)**j) for j in range(NF)]                  # 1/(1+s)
print(f"{'m':>3} {'[s^m](1+s)^m/(1+s)':>20} {'[s^m]s(1+s)^m/(1+s)':>22}")
for m in range(1, 9):
    a = coeff_series(f_inv, m)
    b = coeff_series(f_inv, m, m-1)   # [s^m] s*(1+s)^m f = [s^{m-1}](1+s)^m f
    print(f"{m:>3} {str(a):>20} {str(b):>22}")
print("  Left column identically 0, right column identically 1  =>  E[QP^m] = m! is")
print("  FORCED by the same uniqueness, not a coincidence of either example.")
# uniqueness check: perturb f and watch the left column break
print()
print("  Uniqueness spot-check -- perturb one Taylor coefficient of f:")
for j in (0, 1, 2, 3):
    g = list(f_inv); g[j] += F(1)
    bad = [coeff_series(g, m) for m in range(1, 6)]
    print(f"    f + s^{j}:  [s^m](1+s)^m f = {bad}   {'still 0' if not any(bad) else 'BREAKS'}")

# ================================================================= PART C
print()
print("=" * 78)
print("PART C -- the MASTER EQUATION  Phi(-s g(s)) = 1/(1+s),  and it RECOVERS both g")
print("=" * 78)
print("  For U = chi^2_d / 2 :  E[U^k] = (d/2)_k,  Phi(x) = (1-x)^{-d/2}.")
print("  Master equation  =>  1 + s g(s) = (1+s)^{2/d},  g(s) = ((1+s)^{2/d} - 1)/s.")
print()
from math import comb
for d in (1, 2, 3, 4):
    e = F(2, d)
    # Taylor coefficients of (1+s)^{2/d}
    co = [poch(e - (i - 1), 0) for i in range(0)]        # placeholder
    co = []
    for i in range(7):
        c = F(1)
        for t in range(i): c *= (e - t)
        c /= factorial(i)
        co.append(c)
    g = [co[i+1] for i in range(6)]                      # ((1+s)^{2/d}-1)/s
    ispoly = all(x == 0 for x in g[2:])
    print(f"  d={d} (n = 2+d = {2+d}):  2/d = {e},  g(s) = {[str(x) for x in g[:4]]}..."
          f"   polynomial? {ispoly}")
print()
print("  d=1 -> g = 2 + s   : the owner's (2+Z).  DERIVED, not guessed.")
print("  d=2 -> g = 1       : the S133 four-term example P = (1+W)(Wbar - |Z|^2).")
print("  The two counterexamples are ONE construction at d = 1 and d = 2.")

# verify d=2 reproduces the S133 example numerically in this formalism
P4 = build_P([1])
Pm = ONE; ok4 = True
for m in range(1, 9):
    Pm = cmul(Pm, P4)
    if cexp(Pm, 2) != 0 or cexp(cmul(Z, Pm), 2) != factorial(m): ok4 = False
print(f"  check: d=2, g=1 gives E[P^m]=0 and E[QP^m]=m! for m=1..8:  {ok4}")

# ================================================================= PART D
print()
print("=" * 78)
print("PART D -- MINIMALITY: n = 3 is optimal for this family, n = 2 is UNREACHABLE")
print("=" * 78)
print("  g(s) = ((1+s)^{2/d} - 1)/s is a POLYNOMIAL iff 2/d is a positive integer,")
print("  i.e. iff d in {1, 2}.  Total real dimension = 2 (for Z,W) + d, so n in {3,4}.")
print("  Hence n = 3 is MINIMAL within the family, attained uniquely at d = 1.")
print()
print("  n = 2 would need d = 0, i.e. U constant = c.  Then Phi(x) = e^{cx} and the")
print("  master equation reads  e^{-c s g(s)} = 1/(1+s),  i.e.  c s g(s) = log(1+s).")
print("  log(1+s)/s = 1 - s/2 + s^2/3 - ... is NOT a polynomial, so no such g exists.")
print("  => THIS FAMILY CANNOT REACH GMC(2).  General GMC(2) stays open.")
print()
print("  For completeness, the same obstruction in the Gamma family: U ~ Gamma(alpha)")
print("  has Phi(x) = (1-x)^{-alpha}, forcing 1+sg = (1+s)^{1/alpha}, polynomial iff")
print("  1/alpha is a positive integer.  alpha = 1/2 (d=1) and alpha = 1 (d=2) are the")
print("  only ones realizable as a sum of squares of Gaussians.")
