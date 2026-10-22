#!/usr/bin/env python3
"""The six requested Rogers-Ramanujan-type families, checked exactly to O(q^K),
and their three shared theta families (opus-2026-10-22-ssj).

Statements are transcribed from McLaughlin-Sills-Zimmer, "Rogers-Ramanujan-Slater
type identities", Electron. J. Combin. DS15 (2nd ed. 2022): sections 2.7 (S.31-33),
2.8 (S.34, S.36, S.38, S.39), 2.9 (S.40c-42c), 2.14 (S.59-61), 2.27 (S.90-93).
Notation: f(a,b) = sum_{n in Z} a^{n(n+1)/2} b^{n(n-1)/2} (Ramanujan), f(-q) = (q;q)_inf,
psi(q) = f(q,q^3).

Part 1 checks all 17 identities (plus the Rogers-Ramanujan pair as a positive control
and two deliberately wrong variants as hostile controls, which must FAIL).
Part 2 checks the three pairings of families through common theta functions
  theta_N(a; q) = f(-q^a, -q^(N-a)) = (q^a, q^(N-a), q^N; q^N)_inf:
    Rogers-Selberg (mod 7)   RS_a(q)   = theta_7(a; q)   / f(-q^2)
    Rogers (mod 14)          R14_2a(q) = theta_7(a; q^2) / f(-q)
    Bailey (mod 9)           B_a(q)    = theta_9(a; q)   / f(-q^3)
    Dyson (mod 27)           D_3a(q)   = theta_9(a; q^3) / f(-q)
    Goellnitz-Gordon (mod 8) GG_a(q)   = theta_8(a; q)   / psi(-q)
    Jackson-Slater (mod 8)   JS_a(q)   = theta_8(a; -q)  / f(-q^2)
  so that, as identities between the SUM sides,
    RS_a-sum(q^2) * f(-q^4)  = R14_2a-sum(q) * f(-q),
    B_a-sum(q^3) * f(-q^9)   = D_3a-sum(q) * f(-q)        (a = 1, 2, 4),
    GG_a-sum(-q) * psi(q)    = JS_a-sum(q) * f(-q^2)      (a = 1, 3).
  The identity index sets are (Z/N)^x/{+-1} = {1,2,3}, {1,3}, {1,2,4} for N = 7, 8, 9;
  Dyson adds the non-unit a = 3 (theta_9(3; q) = (q^3; q^3)_inf).
Run: python3 -B qseries_identities.py > qseries_identities.out
"""
from __future__ import annotations

import sys

K = 700  # truncation: all series are exact modulo q^(K+1)


def zero():
    return [0] * (K + 1)


def one():
    s = zero()
    s[0] = 1
    return s


def mul(a, b):
    out = [0] * (K + 1)
    for i, x in enumerate(a):
        if x:
            lim = K - i
            for j in range(lim + 1):
                y = b[j]
                if y:
                    out[i + j] += x * y
    return out


def inv(a):
    """1/a for a series with constant term +-1."""
    assert a[0] in (1, -1)
    out = [0] * (K + 1)
    c0 = a[0]
    out[0] = c0
    for n in range(1, K + 1):
        s = 0
        for k in range(1, n + 1):
            if a[k]:
                s += a[k] * out[n - k]
        out[n] = -s * c0
    return out


def subst(a, m, sign=1):
    """a(sign * q^m)."""
    out = zero()
    for i, x in enumerate(a):
        if x and i * m <= K:
            out[i * m] += x * (sign ** i)
    return out


def mono(c, e):
    s = zero()
    if e <= K:
        s[e] = c
    return s


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def poch(sign, s, step, n):
    """(sign * q^s; q^step)_n = prod_{k<n} (1 - sign * q^(s + k*step))."""
    out = one()
    for k in range(n):
        e = s + k * step
        if e > K:
            break
        fac = one()
        fac[e] -= sign
        out = mul(out, fac)
    return out


def poch_inf(sign, s, step):
    n = (K - s) // step + 1 if s <= K else 0
    return poch(sign, s, step, max(n, 0))


def ramanujan_f(sa, alpha, sb, beta):
    """f(sa q^alpha, sb q^beta) = sum_n sa^{n(n+1)/2} sb^{n(n-1)/2} q^{alpha n(n+1)/2 + beta n(n-1)/2}."""
    out = zero()
    n = 0
    while True:
        done = True
        for m in {n, -n}:
            e = alpha * m * (m + 1) // 2 + beta * m * (m - 1) // 2
            if e <= K:
                done = False
                out[e] += (sa ** ((m * (m + 1) // 2) % 2)) * (sb ** ((m * (m - 1) // 2) % 2))
        if done and n > 0:
            break
        n += 1
    return out


def fm(k):
    """f(-q^k) = (q^k; q^k)_inf."""
    return poch_inf(1, k, k)


def theta(N, a, sign=1, m=1):
    """theta_N(a; sign * q^m) = f(-x^a, -x^(N-a)) at x = sign * q^m."""
    base = ramanujan_f(-1, a, -1, N - a)
    return subst(base, m, sign)


def psi_minus():
    """psi(-q) = f(-q, -q^3)."""
    return ramanujan_f(-1, 1, -1, 3)


def psi_plus():
    return ramanujan_f(1, 1, 1, 3)


# --- sum sides -----------------------------------------------------------------

def sum_series(term, nmax=None):
    out = zero()
    n = 0
    while True:
        t = term(n)
        if t is None:
            break
        out = add(out, t)
        n += 1
    return out


def S_term(qexp, num_factors, den_factors):
    """q^qexp * prod(num) / prod(den), or None if qexp > K."""
    if qexp > K:
        return None
    s = mono(1, qexp)
    for f in num_factors:
        s = mul(s, f)
    for f in den_factors:
        s = mul(s, inv(f))
    return s


def rogers_selberg(a):
    if a == 1:
        return sum_series(lambda n: S_term(2 * n * (n + 1), [], [poch(1, 2, 2, n), poch(-1, 1, 1, 2 * n + 1)]))
    if a == 2:
        return sum_series(lambda n: S_term(2 * n * (n + 1), [], [poch(1, 2, 2, n), poch(-1, 1, 1, 2 * n)]))
    if a == 3:
        return sum_series(lambda n: S_term(2 * n * n, [], [poch(1, 2, 2, n), poch(-1, 1, 1, 2 * n)]))


def rogers14(a2):
    if a2 == 2:
        return sum_series(lambda n: S_term(n * (n + 2), [], [poch(1, 1, 2, n + 1), poch(1, 1, 1, n)]))
    if a2 == 4:
        return sum_series(lambda n: S_term(n * (n + 1), [], [poch(1, 1, 2, n + 1), poch(1, 1, 1, n)]))
    if a2 == 6:
        return sum_series(lambda n: S_term(n * n, [], [poch(1, 1, 2, n), poch(1, 1, 1, n)]))


def gollnitz_gordon(a):
    # S.36 (a = 3): sum q^{n^2} (-q;q^2)_n/(q^2;q^2)_n ; S.34 (a = 1): exponent n(n+2)
    e = (lambda n: n * n) if a == 3 else (lambda n: n * (n + 2))
    return sum_series(lambda n: S_term(e(n), [poch(-1, 1, 2, n)], [poch(1, 2, 2, n)]))


def jackson_slater(a):
    # S.39 (a = 3): sum q^{2n^2}/(q;q)_{2n} ; S.38 (a = 1): q^{2n(n+1)}/(q;q)_{2n+1}
    if a == 3:
        return sum_series(lambda n: S_term(2 * n * n, [], [poch(1, 1, 1, 2 * n)]))
    return sum_series(lambda n: S_term(2 * n * (n + 1), [], [poch(1, 1, 1, 2 * n + 1)]))


def bailey9(a):
    if a == 1:
        return sum_series(lambda n: S_term(3 * n * (n + 1), [poch(1, 1, 1, 3 * n + 1)],
                                           [poch(1, 3, 3, n), poch(1, 3, 3, 2 * n + 1)]))
    if a == 2:
        def t(n):
            if 3 * n * (n + 1) > K:
                return None
            fac = one()
            if 3 * n + 2 <= K:
                fac[3 * n + 2] -= 1
            return S_term(3 * n * (n + 1), [poch(1, 1, 1, 3 * n), fac], [poch(1, 3, 3, n), poch(1, 3, 3, 2 * n + 1)])
        return sum_series(t)
    if a == 4:
        return sum_series(lambda n: S_term(3 * n * n, [poch(1, 1, 1, 3 * n)],
                                           [poch(1, 3, 3, n), poch(1, 3, 3, 2 * n)]))


def dyson27(b):
    if b == 1:
        return sum_series(lambda n: S_term(n * (n + 3), [poch(1, 3, 3, n)], [poch(1, 1, 1, 2 * n + 2), poch(1, 1, 1, n)]))
    if b == 2:
        return sum_series(lambda n: S_term(n * (n + 2), [poch(1, 3, 3, n)], [poch(1, 1, 1, 2 * n + 2), poch(1, 1, 1, n)]))
    if b == 3:
        return sum_series(lambda n: S_term(n * (n + 1), [poch(1, 3, 3, n)], [poch(1, 1, 1, 2 * n + 1), poch(1, 1, 1, n)]))
    if b == 4:
        def t(n):
            if n == 0:
                return one()
            return S_term(n * n, [poch(1, 3, 3, n - 1)], [poch(1, 1, 1, 2 * n - 1), poch(1, 1, 1, n)])
        return sum_series(t)


def rr(a):
    e = (lambda n: n * n) if a == 1 else (lambda n: n * (n + 1))
    return sum_series(lambda n: S_term(e(n), [], [poch(1, 1, 1, n)]))


def check(name, lhs, rhs, expect=True):
    ok = lhs == rhs
    first_bad = next((i for i in range(K + 1) if lhs[i] != rhs[i]), None)
    status = "OK " if ok else "DIFF"
    verdict = "as expected" if ok == expect else "UNEXPECTED"
    print(f"  {status} {name:58s} {'' if ok else '(first difference at q^%d)' % first_bad} [{verdict}]")
    return ok == expect


def main():
    print(f"python {sys.version.split()[0]}; all series exact mod q^{K + 1}")
    allok = True
    f1, f2, f3, f4, f9 = fm(1), fm(2), fm(3), fm(4), fm(9)
    i1, i2, i3 = inv(f1), inv(f2), inv(f3)
    ipm = inv(psi_minus())
    print("== Part 1: the identities ==")
    # positive control: Rogers-Ramanujan
    allok &= check("control RR  G: f(-q^2,-q^3)/f(-q)", rr(1), mul(theta(5, 2), i1))
    allok &= check("control RR  H: f(-q,-q^4)/f(-q)", rr(2), mul(theta(5, 1), i1))
    # hostile controls: swapped products must fail
    allok &= check("hostile: G against f(-q,-q^4)/f(-q)", rr(1), mul(theta(5, 1), i1), expect=False)
    allok &= check("hostile: RS S.33 against f(-q,-q^6)/f(-q^2)", rogers_selberg(3), mul(theta(7, 1), i2), expect=False)
    for a, s in ((1, "S.31"), (2, "S.32"), (3, "S.33")):
        allok &= check(f"Rogers-Selberg {s}: f(-q^{a},-q^{7 - a})/f(-q^2)", rogers_selberg(a), mul(theta(7, a), i2))
    for a, s in ((1, "S.34"), (3, "S.36")):
        allok &= check(f"Goellnitz-Gordon {s}: f(-q^{a},-q^{8 - a})/psi(-q)", gollnitz_gordon(a), mul(theta(8, a), ipm))
    for a, s in ((1, "S.38"), (3, "S.39")):
        allok &= check(f"Jackson-Slater {s}: f(q^{a},q^{8 - a})/f(-q^2)", jackson_slater(a), mul(theta(8, a, sign=-1), i2))
    for a, s in ((1, "S.40c"), (2, "S.41c"), (4, "S.42c")):
        allok &= check(f"Bailey {s}: f(-q^{a},-q^{9 - a})/f(-q^3)", bailey9(a), mul(theta(9, a), i3))
    for a2, s in ((2, "S.59"), (4, "S.60"), (6, "S.61")):
        allok &= check(f"Rogers mod 14 {s}: f(-q^{a2},-q^{14 - a2})/f(-q)", rogers14(a2), mul(theta(14, a2), i1))
    for b, s in ((1, "S.90"), (2, "S.91"), (3, "S.92"), (4, "S.93")):
        allok &= check(f"Dyson mod 27 {s}: f(-q^{3 * b},-q^{27 - 3 * b})/f(-q)", dyson27(b), mul(theta(27, 3 * b), i1))

    print("\n== Part 2: three theta families shared by the six families ==")
    # theta-level identities (definitional, but checked)
    for a in (1, 2, 3):
        allok &= check(f"theta_14(2a={2 * a}; q) = theta_7(a; q^2)", theta(14, 2 * a), theta(7, a, m=2))
    for b in (1, 2, 3, 4):
        allok &= check(f"theta_27(3b={3 * b}; q) = theta_9(b; q^3)", theta(27, 3 * b), theta(9, b, m=3))
    allok &= check("theta_9(3; q) = (q^3;q^3)_inf  (the non-unit Dyson member)", theta(9, 3), f3)
    for a in (1, 3):
        allok &= check(f"theta_8(a={a}; -q) = f(q^{a},q^{8 - a})", theta(8, a, sign=-1), ramanujan_f(1, a, 1, 8 - a))
    allok &= check("psi(q)/f(-q^2) = 1/(q;q^2)_inf", mul(psi_plus(), i2), inv(poch_inf(1, 1, 2)))
    # sum-side corollaries
    for a, a2 in ((1, 2), (2, 4), (3, 6)):
        lhs = mul(subst(rogers_selberg(a), 2), f4)
        rhs = mul(rogers14(a2), f1)
        allok &= check(f"RS_{a}-sum(q^2) f(-q^4) = R14_{a2}-sum(q) f(-q)", lhs, rhs)
    for a, b in ((1, 1), (2, 2), (4, 4)):
        lhs = mul(subst(bailey9(a), 3), f9)
        rhs = mul(dyson27(b), f1)
        allok &= check(f"B_{a}-sum(q^3) f(-q^9) = D_{3 * b}-sum(q) f(-q)", lhs, rhs)
    for a in (1, 3):
        lhs = mul(subst(gollnitz_gordon(a), 1, sign=-1), psi_plus())
        rhs = mul(jackson_slater(a), f2)
        allok &= check(f"GG_{a}-sum(-q) psi(q) = JS_{a}-sum(q) f(-q^2)", lhs, rhs)
    # Dyson's non-unit member is an eta quotient: (q^9;q^9)/(q;q)
    allok &= check("Dyson S.92 product = (q^9;q^9)_inf/(q;q)_inf", mul(theta(27, 9), i1), mul(f9, i1))
    print("\nindex sets: (Z/7)^x/+-1 = {1,2,3} (RS, R14 = 2*{1,2,3}); (Z/8)^x/+-1 = {1,3} (GG, JS);")
    print("            (Z/9)^x/+-1 = {1,2,4} (Bailey; Dyson = 3*{1,2,4} plus the non-unit 3*3)")
    print("\nALL CHECKS PASSED" if allok else "\nSOME CHECK FAILED")
    if not allok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
