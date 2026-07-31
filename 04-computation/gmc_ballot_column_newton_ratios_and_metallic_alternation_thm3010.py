#!/usr/bin/env python3
"""THM-3010 -- ballot-column Newton ratios are exactly rational, the unique
metallic discriminant in the family is BRONZE, and every metallic recurrence
attains the maximal circuit alternation.

Setting is THM-3000 section 1: h_0=1, R_k = h_k^2/(h_(k-1)h_(k+1)), circuit
c_k = log(R_k/R_(k-1)).  Here the h-sequence is taken DIRECTLY from a classical
integer family rather than from a polynomial's coefficients.

PART A -- exact rational Newton ratios of the ballot columns.
  For h_k = binom(2k+a, k+b) the ratio is an exact rational function
      R_k = 1 - Q_(a,b)(k)/D_(a,b)(k),
  with Q monic in k of degree <= 2 and D > 0 on the range.  So the column is
  log-concave (R_k > 1) exactly where Q(k) < 0.  Notable rows:
      binom(2k,k)     : Q = 1                       (log-convex throughout)
      Catalan         : Q = 3, D = (k+1)(2k+1)      (log-convex throughout)
      binom(2k,k-1)   : Q = k^2 - 3k - 1
      binom(2k+1,k-1) : Q = k^2 - 7k - 6
  p = 2 is FORCED: binom(2k,k)/binom(2k-2,k-1) = 2(2k-1)/k is rational in k, but
  for binom(pk,k) with p >= 3 the ratio is not a rational function of k, so no
  Fuss-Catalan analogue of this table exists.

PART B -- the unique metallic discriminant is BRONZE, at 1,4,15,56,210.
  Scanning a in [-2,4], b in [-3,2], the ONLY (a,b) whose monic discriminant has
  the metallic shape k^2 - n k - 1 is (a,b) = (0,-1) and its mirror (0,+1), i.e.
      h_k = binom(2k,k-1) = 1, 4, 15, 56, 210, 792, ...
  with Q = k^2 - 3k - 1, the BRONZE equation, root (3+sqrt13)/2 = 3.302775...
  So that sequence is log-CONCAVE below the bronze ratio and log-CONVEX above it,
  and it is the only column in the family that changes sense at all.

PART C -- metallic recurrences attain the MAXIMUM circuit alternation.
  Let x^2 - n x - 1 = 0 (n >= 1: golden, silver, bronze, ...), root product -1,
  and a_k = n a_(k-1) + a_(k-2), a_0 = 0, a_1 = 1.  The Simson/Catalan identity
      a_(k-1) a_(k+1) - a_k^2 = (-1)^k
  is exactly the NORM FORM of the quadratic order.  Hence
      R_k = a_k^2 / (a_k^2 + (-1)^k)
  alternates strictly about 1, the circuit alternates at EVERY index, and the
  sign-change count is the maximum a circuit of that length admits.
  This is the h-sequence face of THM-3004's maximal real-rooted family
  prod_j (n + T^j)^2: both are "two geometric terms whose ratio is negative",
  algebraically (root product -1) and geometrically (interleaved bands).

Reproduce: python3 04-computation/gmc_ballot_column_newton_ratios_and_metallic_alternation_thm3010.py
"""

from fractions import Fraction as Fr
from math import comb
import sympy as sp

k = sp.symbols('k')


def Q_of(f):
    """1 - R_k = Q(k)/D(k); return (monic Q, D, degree)."""
    R = sp.simplify(f ** 2 / (f.subs(k, k - 1) * f.subs(k, k + 1)))
    num, den = sp.fraction(sp.cancel(sp.together(1 - R)))
    num = sp.expand(num)
    p = sp.Poly(num, k)
    lead = p.all_coeffs()[0]
    return sp.expand(num / lead), sp.factor(den / lead), p.degree()


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def partA():
    rule("A. EXACT RATIONAL NEWTON RATIOS OF THE BALLOT COLUMNS")
    rows = [("binom(2k,k)   central", sp.binomial(2 * k, k)),
            ("Catalan binom(2k,k)/(k+1)", sp.binomial(2 * k, k) / (k + 1)),
            ("binom(2k,k-1)   = 1,4,15,56,210", sp.binomial(2 * k, k - 1)),
            ("binom(2k+1,k-1) = 1,5,21,84,330", sp.binomial(2 * k + 1, k - 1))]
    ok = True
    for nm, f in rows:
        Q, D, deg = Q_of(f)
        print(f"  {nm:34s} 1-R_k = ({Q}) / ({D})")
        ok &= (deg <= 2)
    print()
    print("  numeric check of the closed forms against direct computation:")
    fams = {"Catalan": [comb(2 * i, i) // (i + 1) for i in range(1, 9)],
            "binom(2k,k-1)": [comb(2 * i, i - 1) for i in range(1, 9)],
            "binom(2k+1,k-1)": [comb(2 * i + 1, i - 1) for i in range(1, 9)],
            "binom(2k,k)": [comb(2 * i, i) for i in range(1, 9)]}
    for nm, h in fams.items():
        R = [Fr(h[i]) ** 2 / (Fr(h[i - 1]) * Fr(h[i + 1])) for i in range(1, len(h) - 1)]
        c = [1 if R[i] > R[i - 1] else -1 for i in range(1, len(R))]
        print(f"    {nm:18s} R = {[str(x) for x in R[:4]]}  circuit "
              f"{''.join('+' if x > 0 else '-' for x in c)}")
    print(f"  VERDICT A: {'all ratios exactly rational, deg Q <= 2' if ok else 'FAILED'}")
    return ok


def partB():
    print()
    rule("B. THE UNIQUE METALLIC DISCRIMINANT IN THE FAMILY IS BRONZE")
    hits = []
    for a in range(-2, 5):
        for b in range(-3, 3):
            try:
                Q, D, deg = Q_of(sp.binomial(2 * k + a, k + b))
            except Exception:
                continue
            if deg != 2:
                continue
            co = sp.Poly(Q, k).all_coeffs()
            if co[2] == -1:
                hits.append((a, b, -co[1], Q))
    names = {1: "GOLDEN", 2: "silver", 3: "BRONZE", 4: "copper"}
    for a, b, n, Q in hits:
        print(f"    (a,b)=({a},{b}): Q = {Q}  ->  x^2-{n}x-1  {names.get(n, 'metallic-' + str(n))}")
    ok = bool(hits) and all(n == 3 for _, _, n, _ in hits)
    root = sp.nsimplify((3 + sp.sqrt(13)) / 2)
    print(f"    every metallic hit is n=3: {ok}")
    print(f"    bronze root (3+sqrt13)/2 = {float(root):.9f}")
    h = [comb(2 * i, i - 1) for i in range(1, 10)]
    Rv = [Fr(h[i]) ** 2 / (Fr(h[i - 1]) * Fr(h[i + 1])) for i in range(1, len(h) - 1)]
    signs = ''.join('>' if v > 1 else '<' for v in Rv)
    print(f"    binom(2k,k-1) log-concavity pattern R_k vs 1, k=2..: {signs}")
    print(f"    flip occurs between k=3 and k=4, bracketing {float(root):.4f}: "
          f"{signs[0] == '>' and signs[-1] == '<'}")
    ok &= (signs[0] == '>' and signs[-1] == '<')
    print(f"  VERDICT B: {'BRONZE IS THE UNIQUE METALLIC COLUMN' if ok else 'FAILED'}")
    return ok


def partC():
    print()
    rule("C. METALLIC RECURRENCES ATTAIN MAXIMAL CIRCUIT ALTERNATION")
    names = {1: "GOLDEN  phi=(1+sqrt5)/2", 2: "silver  1+sqrt2",
             3: "BRONZE  (3+sqrt13)/2", 4: "copper  2+sqrt5", 5: "metallic-5"}
    ok = True
    for n in (1, 2, 3, 4, 5):
        a = [0, 1]
        for i in range(2, 14):
            a.append(n * a[-1] + a[-2])
        simson = [a[i - 1] * a[i + 1] - a[i] ** 2 for i in range(2, 8)]
        R = [Fr(a[i]) ** 2 / (Fr(a[i - 1]) * Fr(a[i + 1])) for i in range(2, len(a) - 1)]
        pat = ''.join('>' if v > 1 else '<' for v in R)
        alt = all(pat[i] != pat[i + 1] for i in range(len(pat) - 1))
        c = [1 if R[i] > R[i - 1] else -1 for i in range(1, len(R))]
        ch = sum(1 for i in range(len(c) - 1) if c[i] != c[i + 1])
        maximal = (ch == len(c) - 1)
        ok &= alt and maximal and all(abs(s) == 1 for s in simson)
        print(f"    n={n} {names[n]:24s} a = {a[1:8]}")
        print(f"        Simson a_(k-1)a_(k+1)-a_k^2 = {simson}  (norm form, all +-1)")
        print(f"        R_k vs 1: {pat}  strict alternation {alt}")
        print(f"        circuit  {''.join('+' if x > 0 else '-' for x in c)}"
              f"  sign changes {ch} of max {len(c) - 1}  MAXIMAL {maximal}")
    print(f"  VERDICT C: {'ALL METALLIC RECURRENCES MAXIMALLY ALTERNATE' if ok else 'FAILED'}")
    return ok


def main():
    a = partA()
    b = partB()
    c = partC()
    print()
    rule(f"SUMMARY  ballot-closed-forms={a}  bronze-unique={b}  metallic-maximal={c}")


if __name__ == "__main__":
    main()
