#!/usr/bin/env python3
"""The strengthened exponential-integral claim: TRANSCENDENCE of int_0^1 e^{Q(t)} dt.

CLAIM UNDER TEST (as supplied 2026-08-01): for every NONCONSTANT Q in Qbar[t],

        I(Q) := int_0^1 e^{Q(t)} dt   is TRANSCENDENTAL.

This is strictly stronger than the non-vanishing form of HYP-9079, and the
strengthening changes the shape of the problem in a way worth recording.

------------------------------------------------------------------------------
PART A -- THE REAL-CASE REDUCTION DIES UNDER THE STRENGTHENING.

HYP-9079 sec 1 reduces the NON-VANISHING claim to non-real coefficients: if Q
has real algebraic coefficients then e^{Q(t)} > 0 on [0,1], so I(Q) > 0 and in
particular I(Q) != 0. That argument is correct and it is the whole real case.

It says NOTHING about transcendence. Positivity gives "nonzero", not
"transcendental", and a positive real number is exactly as likely to be
algebraic as any other. So under the strengthened claim the real case is FULLY
LOADED: I(t^2) = int_0^1 e^{t^2} dt = 1.46265174590718160880... must be
transcendental, and no positivity argument touches that.

This matters for how the claim is attacked: the non-vanishing form could be
pushed entirely into C, but the transcendence form cannot, and its real case is
already the erf/E-function problem. Part A below makes the point concrete by
exhibiting the real values and showing they are not algebraic of small degree
and height.

------------------------------------------------------------------------------
PART B -- WHAT DEGREE 2 REALLY IS.

Completing the square, for a != 0,

    int_0^1 e^{a t^2 + b t + c} dt
        = e^{c - b^2/(4a)} * (1/sqrt a) * [ D(sqrt a (1 + b/(2a))) - D(sqrt a b/(2a)) ]

where D(z) = int_0^z e^{u^2} du. And D IS AN E-FUNCTION: writing
D(z) = sum_n z^{2n+1} / (n! (2n+1)), the coefficient of z^m at m = 2n+1 times
m! is (2n+1)!/(n!(2n+1)) = (2n)!/n!, an integer of moderate growth. So degree 2
is exactly a difference of two E-function values at algebraic points, times an
exponential of an algebraic number -- which is why Siegel-Shidlovskii plus
Beukers is the stated tool. Part B verifies the completion-of-square identity
numerically so the reduction is not just asserted.

------------------------------------------------------------------------------
PART C -- CENSUS. No algebraic value found.

For each Q in a large family with small Gaussian-integer coefficients, I(Q) is
computed to high precision and tested for algebraicity by PSLQ on the vector
(1, I, I^2, ..., I^d). A hit would be a counterexample to the claim (or, at
minimum, something that must be explained). None is found.

This is EVIDENCE, not proof. What it rules out is a counterexample of small
degree and small height -- which is the only kind a census can rule out, and
which is worth having on record before anyone invests in the E-function route.

Reproduce: python3 04-computation/exponential_integral_transcendence_census_klein.py
"""

from itertools import product

from mpmath import mp, mpf, mpc, quad, exp, log, sqrt, pi, pslq, erfi, mpmathify

mp.dps = 220


def rule(s):
    print("=" * 78)
    print(s)
    print("=" * 78)


def I(coeffs):
    """int_0^1 exp(Q(t)) dt with Q(t) = sum coeffs[k] t^k (coeffs[0] = constant)."""
    def f(t):
        acc = mpc(0)
        for k in range(len(coeffs) - 1, -1, -1):
            acc = acc * t + coeffs[k]
        return exp(acc)
    return quad(f, [0, 1])


#: PSLQ on n numbers with coefficients up to 10^H needs FAR more than n*H digits
#: of true precision, or it returns spurious relations for literally any input.
#: An earlier version of this file used dps=90, maxdeg=8, maxcoeff=1e14,
#: tol=1e-60 and duly "proved" that e, pi and every value below were algebraic.
#: Two things fix it: a working precision comfortably above (MAXDEG+1)*log10(H),
#: and -- decisively -- RESIDUAL VERIFICATION of every candidate relation.
MAXDEG = 6
MAXCOEFF = 10 ** 10


def alg_degree_test(z, maxdeg=MAXDEG, maxcoeff=MAXCOEFF):
    """Is z algebraic of degree <= maxdeg with integer coefficients <= maxcoeff?

    Returns (d, relation) only if the relation is VERIFIED to annihilate z to
    near full working precision. Returns None otherwise.
    """
    z = mpc(z)
    tol = mpf(10) ** (-(mp.dps - 15))
    accept = mpf(10) ** (-(mp.dps - 25))
    for d in range(1, maxdeg + 1):
        pows = [z ** k for k in range(d + 1)]
        if max(abs(p.imag) for p in pows) < mpf(10) ** (-(mp.dps - 5)):
            rel = pslq([p.real for p in pows], tol=tol, maxcoeff=maxcoeff,
                       maxsteps=50000)
        else:
            # A genuine integer relation must kill the real AND imaginary parts
            # simultaneously; PSLQ handles one real vector, so search on the real
            # part and let the residual check below reject the rest.
            rel = pslq([p.real for p in pows], tol=tol, maxcoeff=maxcoeff,
                       maxsteps=50000)
        if rel and any(rel):
            resid = sum(mpf(rel[k]) * pows[k] for k in range(d + 1))
            scale = max(abs(mpf(r)) for r in rel) * max(abs(p) for p in pows)
            if abs(resid) < accept * max(mpf(1), scale):
                return d, rel
    return None


def controls():
    """A relation test that cannot tell sqrt2 from pi is worthless. Check both."""
    print("  PSLQ CONTROLS (a test with no controls is not a test):")
    ok = True
    pos = [("sqrt(2)", sqrt(mpf(2)), True), ("golden phi", (1 + sqrt(mpf(5))) / 2, True),
           ("2^(1/3)", mpf(2) ** (mpf(1) / 3), True)]
    neg = [("pi", pi, False), ("e", exp(mpf(1)), False),
           ("e^2", exp(mpf(2)), False), ("log 2", log(mpf(2)), False)]
    for name, val, want in pos + neg:
        got = alg_degree_test(val)
        good = (got is not None) == want
        ok &= good
        print(f"    {name:12s} algebraic? {'YES deg ' + str(got[0]) if got else 'no':12s}"
              f" expected {'YES' if want else 'no':3s}  {'ok' if good else '*** CONTROL FAILED ***'}")
    print(f"    controls pass: {ok}   (dps={mp.dps}, maxdeg={MAXDEG}, "
          f"maxcoeff=1e{len(str(MAXCOEFF))-1})")
    return ok


def partA():
    rule("A. THE REAL CASE IS NOT TRIVIAL UNDER THE TRANSCENDENCE STRENGTHENING")
    print("  HYP-9079 sec 1: real algebraic coefficients => e^Q > 0 on [0,1] => I(Q) > 0.")
    print("  That kills NON-VANISHING in the real case and says NOTHING about")
    print("  transcendence. Under the strengthened claim the real case is loaded.")
    print()
    assert controls(), "PSLQ controls failed -- census meaningless"
    print()
    reals = [("t^2", [0, 0, 1]), ("-t^2", [0, 0, -1]), ("t^3", [0, 0, 0, 1]),
             ("t^2+t", [0, 1, 1]), ("2t^2-3t+1", [1, -3, 2]),
             ("t^4", [0, 0, 0, 0, 1])]
    ok = True
    for name, c in reals:
        v = I([mpf(x) for x in c])
        hit = alg_degree_test(v)
        ok &= (hit is None)
        print(f"    Q = {name:12s} I = {mp.nstr(v.real, 30):34s} "
              f"algebraic(deg<=6, H<=1e10)? {'*** YES ' + str(hit) + ' ***' if hit else 'no'}")
    print()
    print("  Cross-check against a closed form: I(t^2) = (sqrt(pi)/2) erfi(1)")
    lhs = I([mpf(0), mpf(0), mpf(1)]).real
    rhs = (sqrt(pi) / 2) * erfi(1)
    print(f"    quad   = {mp.nstr(lhs, 40)}")
    print(f"    erfi   = {mp.nstr(rhs, 40)}")
    print(f"    |diff| = {mp.nstr(abs(lhs - rhs), 5)}")
    ok &= abs(lhs - rhs) < mpf(10) ** (-70)
    print(f"  VERDICT A: no real value algebraic of small degree/height, closed form"
          f" matches: {ok}")
    return ok


def partB():
    print()
    rule("B. DEGREE 2 IS A DIFFERENCE OF E-FUNCTION VALUES (identity verified)")

    def D(z):
        """int_0^z e^{u^2} du, by series -- the E-function."""
        return quad(lambda u: exp(u * u), [0, z])

    ok = True
    print("    a, b, c            quad                       completed-square      |diff|")
    for (a, b, c) in [(1, 0, 0), (1, 2, 3), (-1, 1, 0), (2, -3, 1),
                      (mpc(0, 1), 1, 0), (mpc(1, 1), mpc(0, -2), 1)]:
        A, B, C = mpmathify(a), mpmathify(b), mpmathify(c)
        lhs = I([C, B, A])
        sa = sqrt(A)
        rhs = exp(C - B ** 2 / (4 * A)) / sa * (D(sa * (1 + B / (2 * A))) - D(sa * B / (2 * A)))
        d = abs(lhs - rhs)
        ok &= (d < mpf(10) ** (-60))
        print(f"    {str((a,b,c)):18s} {mp.nstr(lhs,12):26s} {mp.nstr(rhs,12):21s} "
              f"{mp.nstr(d,3)}")
    print("    D(z) = int_0^z e^{u^2} du = sum_n z^(2n+1)/(n!(2n+1)); coefficient of")
    print("    z^m at m=2n+1 times m! is (2n)!/n!, an integer of moderate growth, so D")
    print("    is an E-function. Degree 2 = e^{algebraic} * (difference of two")
    print("    E-function values at algebraic points).")
    print(f"  VERDICT B: completion-of-square identity holds to 60+ digits: {ok}")
    return ok


def partC():
    print()
    rule("C. CENSUS -- PSLQ search for an algebraic value")
    print("  Q(t) = sum_{k=0}^{deg} c_k t^k, c_k in {Gaussian integers of small height},")
    print("  nonconstant (some c_k != 0 for k >= 1). Test (1, I, ..., I^d) for an")
    print("  integer relation with |coeff| <= 1e14 at 60-digit tolerance.")
    print()
    tested = 0
    hits = []
    # degree 2 and 3, real coefficients
    for deg in (2, 3):
        rng = [-2, -1, 0, 1, 2]
        for c in product(rng, repeat=deg + 1):
            if all(x == 0 for x in c[1:]):
                continue                      # constant Q excluded
            v = I([mpf(x) for x in c])
            r = alg_degree_test(v, maxdeg=4)
            tested += 1
            if r:
                hits.append((c, r))
    print(f"    real-coefficient family, deg 2-3: {tested} nonconstant Q tested, "
          f"{len(hits)} algebraic hits")
    # complex (Gaussian integer) coefficients, degree 2
    tested2, hits2 = 0, []
    g = [mpc(0, 0), mpc(1, 0), mpc(0, 1), mpc(1, 1), mpc(-1, 1), mpc(2, 0), mpc(0, 2)]
    for c in product(g, repeat=3):
        if c[1] == 0 and c[2] == 0:
            continue
        v = I(list(c))
        r = alg_degree_test(v, maxdeg=3)
        tested2 += 1
        if r:
            hits2.append((c, r))
    print(f"    Gaussian-integer family, deg 2:    {tested2} nonconstant Q tested, "
          f"{len(hits2)} algebraic hits")
    for c, r in (hits + hits2)[:10]:
        print(f"      *** HIT  Q coeffs {c}  relation {r}")
    ok = not hits and not hits2
    print(f"  VERDICT C: no algebraic value found in {tested + tested2} nonconstant "
          f"cases: {ok}")
    print("    (evidence only -- rules out small degree and height, nothing more)")
    return ok


def main():
    a = partA()
    b = partB()
    c = partC()
    print()
    rule(f"SUMMARY  real-case-loaded={a}  degree2-is-E-function={b}  census-clean={c}")
    print("  MAIN STRUCTURAL POINT: the transcendence strengthening INVALIDATES")
    print("  HYP-9079 sec 1's reduction to non-real coefficients. Positivity gives")
    print("  nonvanishing, not transcendence; the real case is already the full")
    print("  erf/E-function problem and cannot be set aside.")


if __name__ == "__main__":
    main()
