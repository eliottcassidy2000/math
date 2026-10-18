#!/usr/bin/env python3
"""P8 of PREREGISTRATION.md (fourth pre-registration, pushed in 543bc0f3ba before this script was written).
kind-pasteur-2026-10-18-S2.

Lee-Yang weight-2 family on Gamma_0(5)+:  H = (eta(5tau)/eta(tau))^6 = (chi_0 chi_{-1/5})^6,  x = H/(1+22H+125H^2),
sum U_n x^n = (5 E_2(5tau) - E_2(tau))/4,   (n+1)^3 U_{n+1} = 2(2n+1)(11n^2+11n+3) U_n + 4n(4n^2-1) U_{n-1}.
"""
import json
from fractions import Fraction
import numpy as np
import mpmath as mp

import cm_archimedean as ca
from sporadic_cm_detector import primes_between
from ordinary_detector_and_mechanism import scan, mechanism
from sympy import Poly, symbols, gcd as pgcd, div

X = symbols('x')
NMAX = 520


def level5_U(N=NMAX):
    u = [Fraction(1), Fraction(6)]
    for n in range(1, N):
        u.append((2 * (2 * n + 1) * (11 * n * n + 11 * n + 3) * u[n] + 4 * n * (4 * n * n - 1) * u[n - 1]) / Fraction((n + 1) ** 3))
    assert all(v.denominator == 1 for v in u), 'non-integral'
    return [int(v) for v in u]


def E2(tau):
    q = ca.qq(tau)
    s = mp.mpf(0); n = 1
    while True:
        term = n * q ** n / (1 - q ** n)
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5) and n > 3:
            break
        n += 1
    return 1 - 24 * s


def step0(U):
    rep = {}
    for tau in (mp.mpc('0.1234', '0.71'), mp.mpc('-0.31', '0.93'), mp.mpc('0.05', '1.3')):
        H = (ca.eta(5 * tau) / ca.eta(tau)) ** 6
        s = 1 / H
        ly = [ca.vir_char(2, 5, 1, k, tau) for k in (1, 2)]
        rep.setdefault('(chi0 chi-1/5)^6 = H', []).append(float(abs((ly[0] * ly[1]) ** 6 / H - 1)))
        rep.setdefault('j(tau) = (s^2+250s+3125)^3/s^5', []).append(float(abs((s * s + 250 * s + 3125) ** 3 / s ** 5 / ca.jfun(tau) - 1)))
        rep.setdefault('j(5tau) = (s^2+10s+5)^3/s', []).append(float(abs((s * s + 10 * s + 5) ** 3 / s / ca.jfun(5 * tau) - 1)))
    tau = mp.mpc('0.1', '1.2')
    H = (ca.eta(5 * tau) / ca.eta(tau)) ** 6
    x = H / (1 + 22 * H + 125 * H * H)
    lhs = mp.fsum(U[n] * x ** n for n in range(len(U)))
    rhs = (5 * E2(5 * tau) - E2(tau)) / 4
    rep['sum U_n x^n = (5E2(5tau)-E2(tau))/4 at tau=0.1+1.2i'] = float(abs(lhs / rhs - 1))
    return rep


def density_set(u, M, thr=0.3):
    ps = primes_between(7, 500)
    ms = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    hits = np.zeros(len(ms)); tot = np.zeros(len(ms))
    for p in ps:
        coeffs = [x % p for x in u[:p]]
        y = np.arange(p, dtype=np.int64); val = np.zeros(p, dtype=np.int64)
        for c in coeffs[::-1]:
            val = (val * y + c) % p
        inv = np.zeros(p, dtype=np.int64)
        for r in range(1, p):
            inv[r] = pow(r, -1, p)
        zero_r = (val[inv] == 0); zero_r[0] = False
        res = np.mod(ms, p); ok = res != 0
        tot += ok; hits += ok & zero_r[res]
    d = hits / tot
    return {int(m): round(float(v), 3) for m, v in zip(ms[d >= thr], d[d >= thr])}


def level5_point(m):
    m = mp.mpf(m)
    disc = (22 - m) ** 2 - 500
    sq = mp.sqrt(mp.mpc(disc))
    out = []
    for sgn in (1, -1):
        s = ((m - 22) + sgn * sq) / 2
        j0 = (s * s + 10 * s + 5) ** 3 / s
        tau, res = ca.cm_test(j0)
        out.append({'s': mp.nstr(s, 25), 'j(5tau)': mp.nstr(j0, 30), 'tau_reduced': mp.nstr(tau, 25), 'matches': res})
    cm = all(e['matches'] and e['matches'][0]['pass'] for e in out)
    return {'m': int(m), 'CM': cm, 'D': sorted({e['matches'][0]['D'] for e in out if e['matches']}),
            'forms': [e['matches'][0]['form'] if e['matches'] else None for e in out], 'roots': out}


def main():
    U = level5_U()
    rep = {'U_head': U[:12]}
    rep['step0'] = step0(U[:400])
    print('STEP 0:', json.dumps(rep['step0'], indent=1), flush=True)
    M = 100000
    rep['B1prime'] = scan(U, M, False)
    rep['B1prime_ord'] = scan(U, M, True)
    rep['density>=0.3'] = density_set(U, M)
    print("P8a: B1' survivors", rep['B1prime'], "| B1'_ord", rep['B1prime_ord'], '| density>=0.3', rep['density>=0.3'], flush=True)
    rep['P8b'] = [level5_point(m) for m in rep['B1prime']]
    rep['P8b_controls'] = [level5_point(m) for m in (3, 5, 11, 100)]
    for e in rep['P8b'] + rep['P8b_controls']:
        print('P8b', e['m'], 'CM' if e['CM'] else 'no CM match', e['D'], e['forms'], flush=True)
    mech = {}
    S = Poly(1 - 44 * X - 16 * X ** 2, X)
    for p in primes_between(7, 199):
        mm = mechanism(U, p, True)
        F = Poly([U[k] % p for k in range(p)][::-1], X, modulus=p)
        g = pgcd(F, F.diff(X)); rad = div(F, g)[0]
        simple = div(rad, pgcd(rad, g))[0]
        mm['simple part | 1-44x-16x^2'] = div(Poly(S.as_expr(), X, modulus=p), simple)[1].is_zero if simple.degree() else True
        mm['deg simple part'] = simple.degree()
        mech[p] = mm
    ok = all(v['rad|V'] and v['simple part | 1-44x-16x^2'] for v in mech.values())
    rep['P8c'] = {'holds_all_p': ok, 'per_p': {str(k): v for k, v in mech.items()}}
    print('P8c: rad(F_p)|V and simple part | singular polynomial for all 7<=p<=199:', ok,
          '| primes with simple roots:', [p for p, v in mech.items() if v['deg simple part']], flush=True)

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (mp.mpf, mp.mpc)):
            return mp.nstr(o, 25)
        return o
    json.dump(clean(rep), open('level5_leeyang_family.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
