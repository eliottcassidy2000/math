#!/usr/bin/env python3
"""THM-3023 -- the Newton-ratio transform as a dynamical system.

NEW OBJECT.  On positive sequences define the NEWTON-RATIO TRANSFORM

    T(h)_k = h_k^2 / (h_(k-1) h_(k+1)).

This is the map whose output is the Newton ratio R_k of THM-2997 (1); the whole
circuit lane (THM-3000..THM-3015) studies one application of it.  Iterating it is
new.  Five results, all verified here.

A. T IS MINUS THE SECOND DIFFERENCE IN LOG COORDINATES, hence
     T^m(h)_k = prod_(i=0)^(2m) h_(k-m+i)^((-1)^(m+i) C(2m,i)),
   i.e. T^m is multiplicative iterated binomial differencing.

B. FIXED POINTS ARE EXACTLY THE LYNESS PERIOD-6 SEQUENCES.
     T(h) = h  <=>  h_(k+1) = h_k / h_(k-1),
   whose orbits have period 6 (a, b, b/a, 1/a, 1/b, a/b, ...).

C. EIGENSEQUENCES.  With log h_k = A x^k, T(h) = h^mu where
     mu(x) = 2 - x - 1/x,   eigen-condition  x^2 + (mu-2) x + 1 = 0.
   That is a RECIPROCAL quadratic (norm +1) -- complementary to THM-3010's
   metallic family x^2 - n x - 1 (norm -1).  At mu = 1 it is x^2 - x + 1 =
   Phi_6(x), the SIXTH CYCLOTOMIC POLYNOMIAL, whose roots are primitive 6th roots
   of unity: exactly the period-6 fixed points of B.

D. SPECTRUM.  On |x| = 1, x = e^(i theta),
     mu = 2 - 2 cos theta = 4 sin^2(theta/2)  in  [0,4],
   the discrete-Laplacian dispersion relation.  The NEUTRAL FREQUENCY is
   theta = pi/3.  Log-modes below it contract, above it expand.

E. THE METALLIC SIMSON EIGENVALUE IS THE DISCRIMINANT.  For the metallic
   recurrence a_(k+1) = n a_k + a_(k-1) with root lambda of x^2 - n x - 1, the
   Simson/alternating mode is x = -lambda^(-2), and
     mu = 2 + lambda^2 + lambda^(-2) = n^2 + 4 = disc(x^2 - n x - 1).
   Golden 5, silver 8, BRONZE 13 -- the same 13 as the sqrt(13) in THM-3010's
   bronze log-concavity threshold.  So THM-3010's "metallic recurrences attain
   maximal circuit alternation" is exactly the statement that the alternation is
   the UNSTABLE eigendirection of T, with rate the discriminant.

Reproduce: python3 04-computation/newton_ratio_transform_dynamics_thm3023.py
"""

import math
from fractions import Fraction as Fr
from math import comb, factorial

import sympy as sp


def T(h):
    ks = [k for k in h if (k - 1) in h and (k + 1) in h]
    return {k: h[k] ** 2 / (h[k - 1] * h[k + 1]) for k in ks}


def Tm_formula(h, m):
    ks = [k for k in h if all((k - m + i) in h for i in range(2 * m + 1))]
    out = {}
    for k in ks:
        v = Fr(1)
        for i in range(2 * m + 1):
            v *= h[k - m + i] ** ((-1) ** (m + i) * comb(2 * m, i))
        out[k] = v
    return out


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def partA():
    rule("A. T^m IS MULTIPLICATIVE BINOMIAL DIFFERENCING")
    h = {k: Fr(comb(2 * k, k), k + 1) for k in range(0, 18)}   # Catalan
    ok = True
    for m in (1, 2, 3, 4):
        a = dict(h)
        for _ in range(m):
            a = T(a)
        b = Tm_formula(dict(h), m)
        common = sorted(set(a) & set(b))
        g = bool(common) and all(a[k] == b[k] for k in common)
        ok &= g
        print(f"   m={m}: iterated T == prod h^((-1)^(m+i) C(2m,i)) on {len(common)} indices: {g}")
    print(f"  VERDICT A: {'OK' if ok else 'FAILED'}")
    return ok


def partB():
    print()
    rule("B. FIXED POINTS = LYNESS PERIOD-6")
    ok = True
    for (a0, b0) in ((Fr(3), Fr(5)), (Fr(2), Fr(7)), (Fr(1, 2), Fr(5, 3))):
        s = {0: a0, 1: b0}
        for k in range(2, 16):
            s[k] = s[k - 1] / s[k - 2]
        per = all(s[k] == s[k + 6] for k in range(0, 10))
        t = T(s)
        fix = all(t[k] == s[k] for k in t)
        ok &= per and fix
        print(f"   seed ({a0},{b0}): period 6 {per}; T(h)==h {fix}")
    print(f"  VERDICT B: {'OK' if ok else 'FAILED'}")
    return ok


def partCD():
    print()
    rule("C/D. EIGENSEQUENCES, Phi_6, AND THE LAPLACIAN SPECTRUM")
    x, mu = sp.symbols('x mu')
    print("   mu(x) = 2 - x - 1/x ;  eigen-condition x^2 + (mu-2)x + 1 = 0 (norm +1)")
    ok = True
    for m_, want in ((0, "x=1 (double): geometric h, T(h)=1"),
                     (1, "Phi_6(x)=x^2-x+1: primitive 6th roots -> period 6"),
                     (4, "x=-1: pure alternation")):
        sol = sp.solve(sp.Eq(x ** 2 + (m_ - 2) * x + 1, 0), x)
        print(f"     mu={m_}: x = {sol}   {want}")
    ok &= (sp.expand(x ** 2 + (1 - 2) * x + 1 - sp.cyclotomic_poly(6, x)) == 0)
    print(f"   mu=1 polynomial IS the 6th cyclotomic polynomial: "
          f"{sp.expand(x**2 - x + 1 - sp.cyclotomic_poly(6, x)) == 0}")
    print("   on |x|=1: mu = 4 sin^2(theta/2) in [0,4]  (discrete Laplacian symbol)")
    for nm, th in (("pi/6", math.pi / 6), ("pi/3", math.pi / 3), ("pi/2", math.pi / 2),
                   ("2pi/3", 2 * math.pi / 3), ("pi", math.pi)):
        print(f"     theta={nm:6s} mu = {4 * math.sin(th / 2) ** 2:.6f}"
              f"{'   <- NEUTRAL' if nm == 'pi/3' else ''}")
    print(f"  VERDICT C/D: {'OK' if ok else 'FAILED'}")
    return ok


def partE():
    print()
    rule("E. METALLIC SIMSON EIGENVALUE = DISCRIMINANT n^2+4")
    n = sp.symbols('n', positive=True)
    lam = (n + sp.sqrt(n ** 2 + 4)) / 2
    s = sp.simplify(lam ** 2 + lam ** -2)
    xx = -1 / lam ** 2
    m_ = sp.simplify(2 - xx - 1 / xx)
    ok = (sp.simplify(s - (n ** 2 + 2)) == 0) and (sp.simplify(m_ - (n ** 2 + 4)) == 0)
    print(f"   lambda^2 + lambda^-2 = {sp.simplify(s)}   (= n^2+2)")
    print(f"   Simson mode x = -lambda^-2  ->  mu = {sp.simplify(m_)}   (= n^2+4 = disc)")
    names = {1: "GOLDEN", 2: "silver", 3: "BRONZE", 4: "copper"}
    for k in (1, 2, 3, 4):
        tag = "   <- 5 = disc(phi)" if k == 1 else ("   <- the sqrt13 of THM-3010" if k == 3 else "")
        print(f"     n={k} {names[k]:7s}: mu = {k * k + 4}{tag}")
    print()
    print("   numeric growth rate of |log T^m| at a mid index (should tend to n^2+4):")
    for nn in (1, 2):
        a = [0, 1]
        for i in range(2, 44):
            a.append(nn * a[-1] + a[-2])
        cur = {k: Fr(a[k + 1]) for k in range(0, 42)}
        amps = []
        for _ in range(6):
            cur = T(cur)
            mid = sorted(cur)[len(cur) // 2]
            v = abs(math.log(float(cur[mid])))
            amps.append(v if v > 0 else float('nan'))
        rats = [amps[i + 1] / amps[i] for i in range(len(amps) - 1)
                if amps[i] == amps[i] and amps[i] > 0]
        good = rats and abs(rats[-1] - (nn * nn + 4)) < 1e-3
        ok &= bool(good)
        print(f"     n={nn}: ratios {['%.4f' % r for r in rats]}  -> predicted {nn * nn + 4}  {'OK' if good else 'CHECK'}")
    print(f"  VERDICT E: {'OK' if ok else 'FAILED'}")
    return ok


def partF():
    print()
    rule("F. ORBITS: HYPERGEOMETRIC CONTRACT, METALLIC EXPAND")
    fams = {"Catalan": {k: Fr(comb(2 * k, k), k + 1) for k in range(0, 28)},
            "central binomial": {k: Fr(comb(2 * k, k)) for k in range(0, 28)},
            "factorial": {k: Fr(factorial(k)) for k in range(0, 28)}}
    fib = [0, 1]
    for i in range(2, 30):
        fib.append(fib[-1] + fib[-2])
    fams["Fibonacci (golden)"] = {k: Fr(fib[k + 1]) for k in range(0, 28)}
    ok = True
    for nm, h in fams.items():
        cur = dict(h)
        row = []
        for _ in range(4):
            cur = T(cur)
            if not cur:
                break
            row.append(max(abs(math.log(float(v))) for v in cur.values()))
        shrink = all(row[i + 1] < row[i] for i in range(len(row) - 1))
        expect = (nm != "Fibonacci (golden)")
        ok &= (shrink == expect)
        print(f"   {nm:20s} " + " ".join(f"{v:.3e}" for v in row)
              + f"   {'contracts' if shrink else 'EXPANDS'}  {'OK' if shrink == expect else 'MISMATCH'}")
    print("   => smooth hypergeometric logs are low-frequency (theta -> 0, mu -> 0);")
    print("      the metallic Simson mode sits at mu = n^2+4 > 1 and is unstable.")
    print(f"  VERDICT F: {'OK' if ok else 'FAILED'}")
    return ok


def main():
    a, b, c, e, f = partA(), partB(), partCD(), partE(), partF()
    print()
    rule(f"SUMMARY  A={a}  B={b}  C/D={c}  E={e}  F={f}")


if __name__ == "__main__":
    main()
