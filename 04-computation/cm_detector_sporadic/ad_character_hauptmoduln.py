#!/usr/bin/env python3
"""Argyres-Douglas (A_1,A_{2n}) Schur characters (= Virasoro (2,2n+3) characters) as Hauptmoduln of the odd-level
Apery-like families (exploratory; kind-pasteur-2026-10-18-S2).  Exact q-series checks (rational arithmetic):

  n=1 (Lee-Yang, level 5):  Zagier D mirror map  t_D(q) = (chi_0/chi_{-1/5})^5 = q prod (1-q^k)^{5(k/5)}       [Gamma_1(5)]
                             A274786 (weight 2):  x = H/(1+22H+125H^2), H = (chi_0 chi_{-1/5})^6 = (eta5/eta1)^6   [Gamma_0(5)+]
  n=2 (level 7):            Cooper s7:            x = h/(1+13h+49h^2),  h = (chi_1 chi_2 chi_3)^2 = (eta7/eta1)^4  [Gamma_0(7)+]
  n=3 (level 9):            Zagier B mirror map  t_B(q) = -P(-q),  P = chi_1 chi_2 chi_3 chi_4 = (eta9/eta1)^3   [Gamma_0(9), twisted]

The character side uses the Andrews-Gordon/Rocha-Caridi product forms with the q^{h-c/24} prefactors summed; the
mirror maps are computed from the Frobenius log-solution of each recurrence (no modular input).
"""
from fractions import Fraction

NT = 60


def mul(x, y):
    z = [Fraction(0)] * NT
    for i, xi in enumerate(x):
        if xi:
            for j in range(NT - i):
                z[i + j] += xi * y[j]
    return z


def inv(x):
    r = [Fraction(0)] * NT; r[0] = 1 / x[0]
    for n in range(1, NT):
        r[n] = -sum(x[i] * r[n - i] for i in range(1, n + 1)) / x[0]
    return r


def prod_excluding(mod, excluded):
    """prod_{k>=1, k mod `mod` not in excluded} (1-q^k)^{-1} as a power series."""
    s = [Fraction(0)] * NT; s[0] = Fraction(1)
    for k in range(1, NT):
        if k % mod in excluded:
            continue
        for i in range(k, NT):
            s[i] += s[i - k]
    return s


def vir2_char_series(K, i):
    """Andrews-Gordon: chi^{(2,2K+1)}_i (without q-prefactor) = prod_{k != 0, +-i mod 2K+1} (1-q^k)^{-1}, i = 1..K."""
    M = 2 * K + 1
    return prod_excluding(M, {0, i % M, (-i) % M})


def frob_mirror_2nd(a, b, c):
    u = [Fraction(1), Fraction(b)]; g = [Fraction(0)]
    for n in range(1, NT + 1):
        u.append(((a * n * n + a * n + b) * u[n] - c * n * n * u[n - 1]) / Fraction((n + 1) ** 2))
    for n in range(NT):
        gm = g[n - 1] if n else Fraction(0); um = u[n - 1] if n else Fraction(0)
        g.append(((a * n * n + a * n + b) * g[n] - c * n * n * gm - 2 * (n + 1) * u[n + 1] + (2 * a * n + a) * u[n]
                  - 2 * c * n * um) / Fraction((n + 1) ** 2))
    r = mul(g[:NT], inv(u[:NT]))
    e = [Fraction(0)] * NT; e[0] = Fraction(1)
    for n in range(1, NT):
        e[n] = sum(k * r[k] * e[n - k] for k in range(1, n + 1)) / n
    qser = [Fraction(0)] + e[:NT - 1]
    t = [Fraction(0)] * NT; t[1] = Fraction(1)
    for n in range(2, NT):
        comp = [Fraction(0)] * NT; pw = [Fraction(0)] * NT; pw[0] = Fraction(1)
        for k in range(1, n + 1):
            pw = mul(pw, t)
            if qser[k]:
                comp = [cc + qser[k] * pp for cc, pp in zip(comp, pw)]
        t[n] -= comp[n]
    return t


def shift_q(series, power):
    out = [Fraction(0)] * NT
    for i in range(NT - power):
        out[i + power] = series[i]
    return out


def main():
    ok = {}
    # Lee-Yang (2,5): the vacuum chi_0 has no level-1 state, so it is i = 1 (excludes k = +-1 mod 5: H(q), parts
    # +-2 mod 5); chi_{-1/5} is i = 2 (G(q), parts +-1 mod 5).  Prefactors: chi_0 ~ q^{11/60}, chi_{-1/5} ~ q^{-1/60};
    # ratio ~ q^{1/5}, product ~ q^{1/6}.
    c0 = vir2_char_series(2, 1); cm = vir2_char_series(2, 2)
    r = mul(c0, inv(cm)); r5 = r
    for _ in range(4):
        r5 = mul(r5, r)
    r5 = shift_q(r5, 1)                                   # (q^{1/5})^5 = q
    tD = frob_mirror_2nd(11, 3, -1)
    ok['Zagier D mirror map == (chi_0/chi_{-1/5})^5'] = (tD == r5)
    H = mul(c0, cm); H6 = H
    for _ in range(5):
        H6 = mul(H6, H)
    H6 = shift_q(H6, 1)                                   # (q^{1/6})^6
    # compare with (eta(5tau)/eta(tau))^6 = q prod (1-q^{5k})^6/(1-q^k)^6
    base = prod_excluding(1, set())                       # prod (1-q^k)^{-1}
    p5 = [Fraction(0)] * NT; p5[0] = Fraction(1)
    for k in range(1, NT):
        if 5 * k < NT:
            for i in range(NT - 1, 5 * k - 1, -1):
                p5[i] -= p5[i - 5 * k]
    eq = mul(base, p5)
    eq6 = eq
    for _ in range(5):
        eq6 = mul(eq6, eq)
    ok['(chi_0 chi_{-1/5})^6 == (eta5/eta1)^6'] = (H6 == shift_q(eq6, 1))
    # (2,7): chi_i, i = 1..3; product ~ q^{1/2}; squared ~ q
    P7 = mul(mul(vir2_char_series(3, 1), vir2_char_series(3, 2)), vir2_char_series(3, 3))
    h = shift_q(mul(P7, P7), 1)
    p7 = [Fraction(0)] * NT; p7[0] = Fraction(1)
    for k in range(1, NT):
        if 7 * k < NT:
            for i in range(NT - 1, 7 * k - 1, -1):
                p7[i] -= p7[i - 7 * k]
    eq7 = mul(base, p7); eq7_4 = mul(mul(eq7, eq7), mul(eq7, eq7))
    ok['(chi_1 chi_2 chi_3)^2 == (eta7/eta1)^4'] = (h == shift_q(eq7_4, 1))
    # (2,9): chi_i, i = 1..4; product ~ q^{1}
    P9 = mul(mul(vir2_char_series(4, 1), vir2_char_series(4, 2)), mul(vir2_char_series(4, 3), vir2_char_series(4, 4)))
    P9 = shift_q(P9, 1)
    p9 = [Fraction(0)] * NT; p9[0] = Fraction(1)
    for k in range(1, NT):
        if 9 * k < NT:
            for i in range(NT - 1, 9 * k - 1, -1):
                p9[i] -= p9[i - 9 * k]
    eq9 = mul(base, p9); eq9_3 = shift_q(mul(mul(eq9, eq9), eq9), 1)
    ok['chi_1 chi_2 chi_3 chi_4 (2,9) == (eta9/eta1)^3'] = (P9 == eq9_3)
    tB = frob_mirror_2nd(9, 3, 27)
    twisted = [-(x * (-1) ** i) for i, x in enumerate(P9)]
    ok['Zagier B mirror map == -P(-q)'] = (tB == twisted)
    for k, v in ok.items():
        print(f"{k}: {'TRUE' if v else 'FALSE'} (through q^{NT - 1})")


if __name__ == '__main__':
    main()
