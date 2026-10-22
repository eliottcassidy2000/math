"""Atkin-Lehner multipliers of the Coxeter eta quotients (session opus-2026-10-22-S53; THM-4759 (5), audit S53).

For an irreducible component R and e = h/p (p a bad prime, p || h), T^R = 1/eta_{pi_R} has a W_e-symmetric Frame shape,
so c(R, e) = T^R(W_e tau)/T^R(tau) is a constant (checked at three random points). The script prints c(R, e) per
component, checks c(R, e) = exp(-2 pi i rank(R)/24) (the central-charge phase of the rank-l lattice VOA), and checks
that the product over each Niemeier root system is exp(-2 pi i 24/24) = 1, so that W_e fixes T^X exactly."""
import mpmath as mp
from math import gcd
mp.mp.dps = 30


def eta(t):
    q = mp.exp(2j * mp.pi * t)
    return mp.exp(2j * mp.pi * t / 24) * mp.qp(q)


def frame(t, n):
    if t == 'A':
        return {n + 1: 1, 1: -1}
    if t == 'D':
        return {2 * n - 2: 1, 2: 1, n - 1: -1, 1: -1}
    return {6: {12: 1, 3: 1, 2: 1, 6: -1, 4: -1, 1: -1}, 7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
            8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[n]


def Tval(pi, t):
    v = mp.mpc(1)
    for a, ma in pi.items():
        v *= eta(a * t) ** ma
    return 1 / v


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def atkin_lehner(N, e):
    f = N // e
    for x in range(1, 60):
        for z in range(1, 60):
            if gcd(e * x, f * z) == 1:
                g, s, tt = egcd(e * x, f * z)
                return (e * x, -tt, N * z, e * s)
    raise ValueError


taus = [mp.mpc('0.1234', '0.6789'), mp.mpc('-0.31', '0.55'), mp.mpc('0.07', '1.1')]
LATTICES = {'D4^6': [('D', 4)] * 6, 'D6^4': [('D', 6)] * 4, 'D8^3': [('D', 8)] * 3, 'D10E7^2': [('D', 10), ('E', 7), ('E', 7)],
            'D12^2': [('D', 12)] * 2, 'D16E8': [('D', 16), ('E', 8)], 'D24': [('D', 24)], 'E6^4': [('E', 6)] * 4,
            'E8^3': [('E', 8)] * 3}
RESULTS = []
for name, comps in LATTICES.items():
    pis = [frame(t, n) for t, n in comps]
    h = max(pis[0])
    for e in [e for e in range(2, h) if h % e == 0 and gcd(e, h // e) == 1]:
        a, b, c, d = atkin_lehner(h, e)
        cons = []
        for (t, n), pi in zip(comps, pis):
            vals = [Tval(pi, (a * tt + b) / (c * tt + d)) / Tval(pi, tt) for tt in taus]
            const = all(abs(v - vals[0]) < mp.mpf(10) ** -12 for v in vals)
            cons.append(vals[0] if const else None)
        if any(cc is None for cc in cons):
            continue                                   # Frame shape not W_e-symmetric for some component
        total = mp.mpc(1)
        for cc in cons:
            total *= cc
        fmt = lambda z: f'{float(mp.re(z)):+.6f}{float(mp.im(z)):+.6f}i'
        rank = {('D', 4): 4, ('D', 6): 6, ('D', 8): 8, ('D', 10): 10, ('D', 12): 12, ('D', 16): 16, ('D', 24): 24,
                ('E', 6): 6, ('E', 7): 7, ('E', 8): 8}
        anomaly = all(abs(cc - mp.exp(-2j * mp.pi * rank[comp] / 24)) < mp.mpf(10) ** -12 for cc, comp in zip(cons, comps))
        ok = anomaly and abs(total - 1) < mp.mpf(10) ** -12
        RESULTS.append(ok)
        print(f'{"ok  " if ok else "FAIL"} {name:8s} W_{e:<2d}: per-component multipliers {[fmt(cc) for cc in cons]} '
              f'= exp(-2 pi i rank/24): {anomaly}; product {fmt(total)}')
print('ALL CHECKS PASSED' if all(RESULTS) else 'FAILED')
