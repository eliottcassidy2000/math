#!/usr/bin/env python3
"""EXACT classification of ALL rational CM points of the two Fricke-invariant weight-2 families with AD-character
Hauptmoduln (exploratory follow-up to P4/P8; kind-pasteur-2026-10-18-S2).

  level 7 (Cooper s7):    t = (eta(tau)/eta(7tau))^4 = (chi_1 chi_2 chi_3)^{-2}   [(A1,A4) characters]
                          m = 1/x = t + 13 + 49/t,   j(7tau) = (t^2+13t+49)(t^2+5t+1)^3/t = m((m-8)t-48)^3
  level 5 (A274786):      s = (eta(tau)/eta(5tau))^6 = (chi_0 chi_{-1/5})^{-6}  [(A1,A2) = Lee-Yang]
                          m = 1/x = s + 22 + 125/s,  j(5tau) = (s^2+10s+5)^3/s = ((m-12)s-120)^3/s

Argument.  For rational m, the variable t (resp. s) is quadratic over Q, so J = j(N tau) lies in a field of degree
<= 2; if it is a CM value its order has class number h(D) <= 2.  All such D are known (h = 1: 13 orders; h = 2:
29 orders; CITED: Heegner-Baker-Stark, Baker-Stark-Montgomery-Weinberger/Arno; recomputed here by counting reduced
forms for |D| <= 2000).  For each such D, R_D(m) = Res_J(P(m,J), H_D(J)) with P(m,J) = Res_t(J - formula, quadratic)
is an integer polynomial whose rational roots are exactly the rational m with a CM point of discriminant D.
H_D is computed from the reduced forms at 120 digits and rounded (coefficients are integers; the rounding margin is
printed).  Output: the complete list of rational (and integral) CM parameters m.
"""
import json
from math import gcd, isqrt
import mpmath as mp
from sympy import symbols, resultant, Poly, Rational, factor_list, ZZ, QQ, expand

mp.mp.dps = 120
m, J, t = symbols('m J t')


def reduced_forms(D):
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or gcd(gcd(a, abs(b)), c) != 1:
                continue
            if b < 0 and (a == c):
                continue
            out.append((a, b, c))
        a += 1
    return out


def hilbert_poly(D):
    forms = reduced_forms(D)
    js = [1728 * mp.kleinj((-b + mp.sqrt(-D) * 1j) / (2 * a)) for a, b, c in forms]
    coeffs = [mp.mpc(1)]
    for jv in js:
        new = [mp.mpc(0)] * (len(coeffs) + 1)
        for i, cf in enumerate(coeffs):
            new[i] += cf
            new[i + 1] -= cf * jv
        coeffs = new
    ints = [int(mp.nint(cf.real)) for cf in coeffs]
    margin = max(float(abs(cf - iv)) for cf, iv in zip(coeffs, ints))
    return Poly(sum(iv * J ** (len(ints) - 1 - k) for k, iv in enumerate(ints)), J), margin


def classify(formula_num, quad, label):
    P = resultant(formula_num, quad, t)          # polynomial in m, J
    P = Poly(expand(P), J, m)
    hlist = [D for D in range(-3, -2001, -1) if D % 4 in (0, 1) and len(reduced_forms(D)) <= 2]
    found = {}
    worst = 0.0
    for D in hlist:
        H, margin = hilbert_poly(D)
        worst = max(worst, margin)
        R = resultant(P.as_expr(), H.as_expr(), J)
        R = Poly(expand(R), m)
        if R.is_zero:
            found.setdefault('IDENTICALLY ZERO', []).append(D)
            continue
        _, facs = factor_list(R.as_expr(), m)
        for f, e in facs:
            fp = Poly(f, m)
            if fp.degree() == 1:
                a1, a0 = fp.all_coeffs()
                root = Rational(-a0, a1)
                found.setdefault(str(root), set()).add(D)
    print(f"{label}: class-number<=2 discriminants used: {len(hlist)} (max |D| = {-min(hlist)}); worst rounding margin of H_D: {worst:.2e}")
    rat = sorted(((Rational(k), sorted(v)) for k, v in found.items() if k != 'IDENTICALLY ZERO'), key=lambda kv: kv[0])
    for r, Ds in rat:
        print(f"   m = {r}{'   (integer)' if r.q == 1 else ''}   D in {Ds}")
    return {'discriminants_used': hlist, 'rounding_margin': worst,
            'rational_CM_m': [(str(r), Ds) for r, Ds in rat],
            'integral_CM_m': [int(r) for r, Ds in rat if r.q == 1]}


def main():
    out = {}
    # level 7:  t*J - (t^2+13t+49)(t^2+5t+1)^3 with t^2 + (13-m) t + 49 = 0
    out['level7_s7'] = classify(t * J - (t ** 2 + 13 * t + 49) * (t ** 2 + 5 * t + 1) ** 3, t ** 2 + (13 - m) * t + 49, 'level 7 (s7)')
    # level 5:  t*J - (t^2+10t+5)^3 with t^2 + (22-m) t + 125 = 0   (t plays the role of s)
    out['level5_A274786'] = classify(t * J - (t ** 2 + 10 * t + 5) ** 3, t ** 2 + (22 - m) * t + 125, 'level 5 (A274786)')
    surv = {'level7_s7': [-10648, -64, -8, -1, 2, 8, 24, 27, 125],
            'level5_A274786': [-15228, -828, -28, -3, 4, 18, 36, 72, 147]}
    for k in out:
        integral = sorted(set(out[k]['integral_CM_m']) - {0})
        out[k]['integral_nonzero_CM_m'] = integral
        out[k]['equals_detector_survivors'] = integral == sorted(surv[k])
        print(k, 'integral nonzero CM m:', integral, '== detector survivors:', integral == sorted(surv[k]))
    json.dump(out, open('integral_cm_classification.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
