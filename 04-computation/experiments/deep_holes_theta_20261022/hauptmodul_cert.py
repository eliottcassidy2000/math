"""1/eta_{pi_X} is a Hauptmodul for every Niemeier root system X (session opus-2026-10-22-S53).

pi_X = product over the components of X of the Frame shape of the Coxeter element (all components share the Coxeter
number h). For f = 1/eta_{pi_X} = prod_d eta(d tau)^(r_d) on Gamma_0(h):
  1. Ligozat's formula gives the order of f at every cusp (an eta quotient has no zeros or poles in H);
  2. the Atkin-Lehner involutions W_e (e || h) fixing f are found numerically (mpmath, 30 digits);
  3. certificate: the poles of f are exactly the S-orbit of infinity, all simple, S acting freely on it.
Then f descends to X_0(h)/S with a single simple pole, so X_0(h)/S has genus 0 and f is a Hauptmodul.
The group Gamma_0(h)+S is compared with the lambency of X in Cheng-Duncan-Harvey's umbral moonshine.
Needs mpmath. Usage: python hauptmodul_cert.py"""
import math
import re
import mpmath as mp

mp.mp.dps = 30
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def coxeter_frame(t, r):
    if t == 'a':
        return {r + 1: 1, 1: -1}
    if t == 'd':
        out = {}
        for a, m in ((2 * r - 2, 1), (2, 1), (r - 1, -1), (1, -1)):
            out[a] = out.get(a, 0) + m
        return out
    return {6: {12: 1, 2: 1, 3: 1, 6: -1, 4: -1, 1: -1},
            7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
            8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[r]


TYPES = ['a1^24', 'a2^12', 'a3^8', 'a4^6', 'a5^4d4', 'd4^6', 'a6^4', 'a7^2d5^2', 'a8^3', 'a9^2d6', 'd6^4',
         'a11d7e6', 'e6^4', 'a12^2', 'd8^3', 'a15d9', 'a17e7', 'd10e7^2', 'd12^2', 'a24', 'd16e8', 'e8^3', 'd24']
# lambencies of Cheng-Duncan-Harvey, "Umbral moonshine and the Niemeier lattices" (2014), Table 2
LAMBENCY = {'d4^6': '6+3', 'd6^4': '10+5', 'd8^3': '14+7', 'd10e7^2': '18+9', 'd12^2': '22+11', 'd16e8': '30+15',
            'd24': '46+23', 'e6^4': '12+4', 'e8^3': '30+6,10,15'}


def eta(tau):
    q = mp.exp(2j * mp.pi * tau)
    s = mp.mpf(1)
    n = 1
    while True:
        t = q ** n
        s *= (1 - t)
        if abs(t) < mp.mpf(10) ** (-40):
            break
        n += 1
    return mp.exp(2j * mp.pi * tau / 24) * s


def f_val(r, tau):
    v = mp.mpc(1)
    for d, e in r.items():
        v *= eta(d * tau) ** e
    return v


def ligozat(r, N, d):
    g = math.gcd(d, N // d)
    s = sum(math.gcd(d, dl) ** 2 * e / dl for dl, e in r.items())
    return N * s / (24 * g * d)


def series(r, L):
    s = [0] * L; s[0] = 1
    for a, m in r.items():
        for k in range(1, L // a + 1):
            for _ in range(abs(m)):
                if m > 0:
                    s = [s[t] - (s[t - a * k] if t >= a * k else 0) for t in range(L)]
                else:
                    for t in range(a * k, L):
                        s[t] += s[t - a * k]
    return s


for name in TYPES:
    comps = []
    for t, rr, ex in re.findall(r'([ade])(\d+)(?:\^(\d+))?', name):
        comps += [(t, int(rr))] * (int(ex) if ex else 1)
    pi = {}
    for t, rr in comps:
        for a, m in coxeter_frame(t, rr).items():
            pi[a] = pi.get(a, 0) + m
    pi = {a: m for a, m in pi.items() if m}
    h = max(pi)
    k = len(comps)
    r = {a: -m for a, m in pi.items()}                     # f = 1/eta_pi
    assert sum(a * m for a, m in pi.items()) == 24 and sum(pi.values()) == 0
    divs = [d for d in range(1, h + 1) if h % d == 0]
    orders = {d: ligozat(r, h, d) for d in divs}
    poles = sorted(d for d in divs if orders[d] < -1e-12)
    hall = [e for e in divs if math.gcd(e, h // e) == 1]
    S = []
    for e in hall:
        # W_e = [[e x, y], [h z, e w]] with e x w - (h/e) y z = 1; take x = z = 1
        y = next(y for y in range(-10 * h, 10 * h) if (1 + (h // e) * y) % e == 0)
        w = (1 + (h // e) * y) // e
        t0 = mp.sqrt(e) / h
        tau = mp.mpc(-mp.mpf(e * w) / h + mp.mpf('0.0137'), t0 * mp.mpf('1.13'))
        tw = (e * tau + y) / (h * tau + e * w)
        a, b = f_val(r, tau), f_val(r, tw)
        if abs(a - b) < mp.mpf(10) ** (-15) * (1 + abs(a)):
            S.append(e)
    orbit = sorted(h // e for e in S)
    simple = all(abs(orders[d] + 1) < 1e-12 for d in poles)
    single_cusp = all(math.gcd(d, h // d) <= 2 for d in poles)       # phi(gcd) = 1: one cusp with denominator d
    cert = poles == orbit and simple and single_cusp and abs(orders[h] + 1) < 1e-12
    lam = str(h) + ('+' + ','.join(str(e) for e in S if e != 1) if len(S) > 1 else '')
    lam_ok = lam == LAMBENCY.get(name, str(h))
    ser = series(r, 8)
    Tser = [ser[0], ser[1] + k] + ser[2:]                      # T_X = f + k = q^-1 + 0 + ...
    pistr = ' '.join(f'{a}^{m}' for a, m in sorted(pi.items(), key=lambda x: -x[0]))
    check(cert and lam_ok and ser[1] == -k,
          f'{name:9s} pi_X = {pistr:28s} poles at cusps 1/d, d in {poles}; invariant under W_e, e in {S}; '
          f'Hauptmodul of Gamma_0({lam}) (CDH lambency {LAMBENCY.get(name, str(h))}); T_X = q^-1 + {Tser[2:6]}...')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
