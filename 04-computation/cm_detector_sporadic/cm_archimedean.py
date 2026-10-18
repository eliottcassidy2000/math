#!/usr/bin/env python3
"""Archimedean CM checks P3', P4, P5 of PREREGISTRATION.md (second pre-registration, pushed in 2583387d46 before this
script was written).  kind-pasteur-2026-10-18-S2.

Step 0. Verify, by q-series at random tau, the identities the tests rely on:
  - (A1,A4) = Virasoro (2,7) Schur characters (Rocha-Caridi formula): (chi_1 chi_2 chi_3)^2 = h = (eta(7tau)/eta(tau))^4;
  - (A1,A2) = Lee-Yang (2,5):  (chi_0/chi_{-1/5})^5 = r5 = q prod (1-q^n)^{5(n/5)};
  - the level-7 j-maps, t = 1/h:  (t^2+13t+49)(t^2+5t+1)^3/t  and  (t^2+13t+49)(t^2+245t+2401)^3/t^7  (which one is
    j(tau) and which is j(7tau) is decided numerically, not assumed);
  - the level-5 icosahedral j-map in r5.
Step 1 (P4). s7 block-test survivors m (|m| <= 1e5, minus the singular fibres -1, 27): both roots h of
  49h^2+(13-m)h+1 = 0, the pre-registered j = (t^2+13t+49)(t^2+5t+1)^3/t, matched against j(tau_Q) for reduced
  primitive forms Q with -4000 <= D < 0 (tau found by Newton in the SL2(Z) fundamental domain, then the finite list of
  reduced forms with that tau is searched exactly).  Cross-check: h(gamma tau_Q) = h_0 for a coset representative.
Step 2 (P5). Controls m = 3, 5, 11, 100 for s7.
Step 3 (P3'). S_D = {0<|m|<=3000 : D-density >= 0.30}; j via the verified level-5 map; same CM match.
Usage: python cm_archimedean.py
"""
import json
from fractions import Fraction
from math import gcd, isqrt
import numpy as np
import mpmath as mp

from sporadic_cm_detector import zagier_rec, primes_between

mp.mp.dps = 80
TWO_PI_I = 2j * mp.pi
DMAX = 4000


# ----------------------------------------------------------------------------------------------- q-series helpers
def qq(tau):
    return mp.exp(TWO_PI_I * tau)


def eta(tau):
    q = qq(tau)
    return mp.exp(TWO_PI_I * tau / 24) * mp.qp(q)


def jfun(tau):
    return 1728 * mp.kleinj(tau)


def legendre5(n):
    return [0, 1, -1, -1, 1][n % 5]


def r5(tau, nmax=None):
    """q prod_{n>=1} (1-q^n)^{5 (n/5)} (Rogers-Ramanujan continued fraction to the fifth power)."""
    q = qq(tau)
    val = q
    n = 1
    tol = mp.mpf(10) ** (-mp.mp.dps - 5)
    while True:
        qn = q ** n
        e = 5 * legendre5(n)
        if e:
            val *= (1 - qn) ** e
        if abs(qn) < tol and n > 5:
            break
        n += 1
    return val


def vir_char(p, pp, r, s, tau, nterms=None):
    """Rocha-Caridi character chi_{r,s} of the Virasoro (p,p') minimal model, including q^{h - c/24}."""
    q = qq(tau)
    c = 1 - mp.mpf(6) * (p - pp) ** 2 / (p * pp)
    h = mp.mpf((pp * r - p * s) ** 2 - (p - pp) ** 2) / (4 * p * pp)
    tot = mp.mpf(0)
    k = 0
    tol = mp.mpf(10) ** (-mp.mp.dps - 5)
    while True:
        terms = []
        for kk in ([k] if k == 0 else [k, -k]):
            e1 = mp.mpf((2 * p * pp * kk + pp * r - p * s) ** 2 - (pp * r - p * s) ** 2) / (4 * p * pp)
            e2 = mp.mpf((2 * p * pp * kk + pp * r + p * s) ** 2 - (pp * r - p * s) ** 2) / (4 * p * pp)
            terms.append(q ** int(mp.nint(e1)) - q ** int(mp.nint(e2)))
        tot += sum(terms)
        if k > 3 and max(abs(t) for t in terms) < tol:
            break
        k += 1
    return mp.exp(TWO_PI_I * tau * (h - c / 24)) * tot / mp.qp(q)


# ----------------------------------------------------------------------------------------------- CM machinery
def reduce_tau(tau):
    for _ in range(1000):
        tau = mp.mpc(tau.real - mp.nint(tau.real), tau.imag)
        if abs(tau) < 1 - mp.mpf(10) ** (-mp.mp.dps + 10):
            tau = -1 / tau
        else:
            break
    if tau.real > mp.mpf(1) / 2:
        tau -= 1
    return tau


def tau_from_j(j0):
    """tau in the standard fundamental domain with j(tau) = j0 (secant/Newton from a q ~ 1/j start or a grid)."""
    j0 = mp.mpc(j0)
    starts = []
    if abs(j0) > 1e5:
        q0 = 1 / (j0 - 744)
        starts.append(mp.log(q0) / TWO_PI_I)
    else:
        with mp.workdps(15):
            best = []
            for a in np.linspace(-0.5, 0.5, 41):
                for b in np.linspace(0.866, 2.5, 41):
                    z = mp.mpc(a, b)
                    best.append((abs(jfun(z) - j0), z))
            best.sort(key=lambda t: t[0])
            starts.extend(z for _, z in best[:4])
    for z0 in starts:
        try:
            z = mp.findroot(lambda z: jfun(z) - j0, mp.mpc(z0), tol=mp.mpf(10) ** (-mp.mp.dps + 15), maxsteps=200)
        except (ValueError, ZeroDivisionError):
            continue
        z = reduce_tau(z)
        if abs(jfun(z) - j0) <= mp.mpf(10) ** (-40) * max(1, abs(j0)):
            return z
    raise RuntimeError(f"tau_from_j failed for j0 = {j0}")


def cm_match(tau, dmax=DMAX):
    """All reduced primitive forms (a,b,c), b^2-4ac = D in [-dmax, 0), with tau_Q = tau (to 1e-40)."""
    out = []
    for a in range(1, isqrt(dmax // 3) + 2):
        b = int(mp.nint(-2 * a * tau.real))
        if abs(-2 * a * tau.real - b) > mp.mpf(10) ** (-40):
            continue
        Dabs = int(mp.nint((2 * a * tau.imag) ** 2))
        if abs((2 * a * tau.imag) ** 2 - Dabs) > mp.mpf(10) ** (-40) or Dabs > dmax:
            continue
        D = -Dabs
        if (b * b - D) % (4 * a):
            continue
        c = (b * b - D) // (4 * a)
        if gcd(gcd(a, abs(b)), c) != 1 or not (abs(b) <= a <= c):
            continue
        tq = (-b + mp.sqrt(mp.mpf(Dabs)) * 1j) / (2 * a)
        out.append({'form': (a, b, c), 'D': D, 'tau_Q': tq})
    return out


def cm_test(j0, dmax=DMAX):
    tau = tau_from_j(j0)
    ms = cm_match(tau, dmax)
    res = []
    for mt in ms:
        jq = jfun(mt['tau_Q'])
        rel = abs(jq - j0) / max(1, abs(j0))
        res.append({'form': mt['form'], 'D': mt['D'], 'rel_err': float(rel), 'pass': bool(rel < 1e-20)})
    return tau, res


# ----------------------------------------------------------------------------------------------- level-7 formulas
def F7(t):     # pre-registered formula
    return (t * t + 13 * t + 49) * (t * t + 5 * t + 1) ** 3 / t


def G7(t):
    return (t * t + 13 * t + 49) * (t * t + 245 * t + 2401) ** 3 / t ** 7


def h7(tau):
    return (eta(7 * tau) / eta(tau)) ** 4


def ICOS(t):   # classical icosahedral map in r^5
    return -(t ** 4 - 228 * t ** 3 + 494 * t ** 2 + 228 * t + 1) ** 3 / (t * (t * t + 11 * t - 1) ** 5)


def ICOS5(t):  # candidate for j(5 tau)
    return -(t ** 4 + 12 * t ** 3 + 14 * t ** 2 - 12 * t + 1) ** 3 / (t ** 5 * (t * t + 11 * t - 1))


def step0():
    rep = {}
    for tau in (mp.mpc('0.1234', '0.71'), mp.mpc('-0.31', '0.93'), mp.mpc('0.05', '1.3')):
        h = h7(tau)
        chis = [vir_char(2, 7, 1, s, tau) for s in (1, 2, 3)]
        rep.setdefault('A1A4_product_identity', []).append(float(abs((chis[0] * chis[1] * chis[2]) ** 2 / h - 1)))
        t = 1 / h
        rep.setdefault('F7_vs_j(7tau)', []).append(float(abs(F7(t) / jfun(7 * tau) - 1)))
        rep.setdefault('F7_vs_j(tau)', []).append(float(abs(F7(t) / jfun(tau) - 1)))
        rep.setdefault('G7_vs_j(tau)', []).append(float(abs(G7(t) / jfun(tau) - 1)))
        r = r5(tau)
        ly = [vir_char(2, 5, 1, s, tau) for s in (1, 2)]   # s=1: h=0 (vacuum), s=2: h=-1/5
        rep.setdefault('LeeYang_ratio_identity', []).append(float(abs((ly[0] / ly[1]) ** 5 / r - 1)))
        rep.setdefault('LeeYang_product_identity', []).append(float(abs(ly[0] * ly[1] / (eta(5 * tau) / eta(tau)) - 1)))
        rep.setdefault('ICOS_vs_j(tau)', []).append(float(abs(ICOS(r) / jfun(tau) - 1)))
        rep.setdefault('ICOS5_vs_j(5tau)', []).append(float(abs(ICOS5(r) / jfun(5 * tau) - 1)))
    return rep


# ----------------------------------------------------------------------------------------------- tests
def level7_point(m, label):
    """Both roots h of 49h^2+(13-m)h+1 = 0, the pre-registered j, CM match, and the Hauptmodul cross-check."""
    m = mp.mpf(m)
    disc = (13 - m) ** 2 - 196
    sq = mp.sqrt(mp.mpc(disc))
    out = []
    for sgn in (1, -1):
        h0 = (-(13 - m) + sgn * sq) / 98
        j0 = F7(1 / h0)
        tau, res = cm_test(j0)
        entry = {'h0': mp.nstr(h0, 25), 'j_prereg': mp.nstr(j0, 30), 'tau_reduced': mp.nstr(tau, 25), 'matches': res}
        if res and res[0]['pass']:
            # Cross-check with the eta-quotient itself.  If F7(t(tau')) = j(7 tau') (step 0), then 7 tau' = gamma tau_Q,
            # and modulo Gamma_0(7) the 8 solutions are tau' = (tau_Q + k)/7 (k = 0..6) and tau' = -1/(7 tau_Q), where
            # h(-1/(7 tau)) = 1/(49 h(tau)) (Fricke).  One of the 8 values of h must equal h0.
            tq = mp.mpc(cm_match(tau)[0]['tau_Q'])
            cands = [('(tau_Q+%d)/7' % k, (tq + k) / 7, h7((tq + k) / 7)) for k in range(7)]
            cands.append(('-1/(7 tau_Q)', -1 / (7 * tq), 1 / (49 * h7(tq))))
            lab, tp, hv = min(cands, key=lambda c: abs(c[2] - h0))
            entry['hauptmodul_crosscheck'] = {'tau_prime': lab, 'tau_prime_value': mp.nstr(tp, 25),
                                              'abs_h_diff': float(abs(hv - h0)),
                                              'm_times_x(tau_prime)': mp.nstr(m * hv / (1 + 13 * hv + 49 * hv ** 2), 25)}
        out.append(entry)
    verdict = all(e['matches'] and e['matches'][0]['pass'] for e in out)
    return {'m': int(m), 'label': label, 'roots': out, 'CM': verdict,
            'D': sorted({e['matches'][0]['D'] for e in out if e['matches']})}


def d_density_table(M=3000):
    u = zagier_rec(11, 3, -1)
    ps = primes_between(7, 500)
    ms = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    hits = np.zeros(len(ms)); tot = np.zeros(len(ms))
    for p in ps:
        zero = np.zeros(p, dtype=bool)
        for r in range(1, p):
            ir = pow(r, -1, p); s = 0; w = 1
            for k in range(p):
                s = (s + u[k] * w) % p
                w = w * ir % p
            zero[r] = (s == 0)
        res = np.mod(ms, p)
        ok = res != 0
        tot += ok
        hits += ok & zero[res]
    dens = hits / tot
    return {int(m): float(d) for m, d in zip(ms, dens)}


def main():
    report = {}
    report['step0_identities'] = step0()
    print('STEP 0 (relative errors):', json.dumps(report['step0_identities'], indent=1), flush=True)

    s7_surv = [-10648, -64, -8, 2, 8, 24, 125]
    report['P4'] = [level7_point(m, 'survivor') for m in s7_surv]
    report['P4_singular_posthoc'] = [level7_point(m, 'singular fibre (post hoc)') for m in (-1, 27)]
    report['P5'] = [level7_point(m, 'control') for m in (3, 5, 11, 100)]
    for blk in ('P4', 'P4_singular_posthoc', 'P5'):
        for e in report[blk]:
            print(blk, e['m'], 'CM' if e['CM'] else 'no CM match', e['D'],
                  [(r['matches'][0]['form'] if r['matches'] else None) for r in e['roots']],
                  [r.get('hauptmodul_crosscheck', {}).get('abs_h_diff') for r in e['roots']], flush=True)

    dens = d_density_table(3000)
    SD = sorted(m for m, d in dens.items() if d >= 0.30)
    report['P3prime_SD'] = SD
    report['P3prime_density_hist'] = np.histogram(list(dens.values()), bins=[0, .05, .1, .2, .3, .4, .5, .6, 1])[0].tolist()
    print('S_D =', SD, ' density histogram', report['P3prime_density_hist'], flush=True)
    p3 = []
    for m in SD:
        t0 = Fraction(1, m)
        num = (t0 ** 4 - 228 * t0 ** 3 + 494 * t0 ** 2 + 228 * t0 + 1) ** 3
        den = t0 * (t0 * t0 + 11 * t0 - 1) ** 5
        jexact = -num / den
        j0 = mp.mpf(jexact.numerator) / jexact.denominator
        tau, res = cm_test(j0)
        p3.append({'m': m, 'density': dens[m], 'j_exact': str(jexact), 'tau_reduced': mp.nstr(tau, 25),
                   'matches': res, 'CM': bool(res and res[0]['pass'])})
        print('P3prime', m, round(dens[m], 3), 'j =', jexact, '->', [(r['form'], r['D'], r['pass']) for r in res], flush=True)
    report['P3prime'] = p3

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (mp.mpf, mp.mpc)):
            return mp.nstr(o, 25)
        return o
    json.dump(clean(report), open('cm_archimedean.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
