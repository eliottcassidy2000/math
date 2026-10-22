"""Audit check (S53): W_e-invariance of the full T^X = 1/eta_{pi^X} for the 23 Niemeier root systems, with genuine
Atkin-Lehner matrices (a = e x, c = h z, d = e w, det e), at two random points; invariant set must be S_X of THM-4759 (4).
(A single component's eta quotient has degree l < 24 and only picks up a root of unity: e.g. D4 gets exp(-i pi/3).)"""
import mpmath as mp
from math import gcd
mp.mp.dps = 30
def eta(t):
    return mp.exp(2j * mp.pi * t / 24) * mp.qp(mp.exp(2j * mp.pi * t))
def frame(t, n):
    if t == 'A': return {n + 1: 1, 1: -1}
    if t == 'D': return {2 * n - 2: 1, 2: 1, n - 1: -1, 1: -1}
    return {6: {12: 1, 3: 1, 2: 1, 6: -1, 4: -1, 1: -1}, 7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
            8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[n]
X = {'A1^24': [('A', 1)] * 24, 'A2^12': [('A', 2)] * 12, 'A3^8': [('A', 3)] * 8, 'A4^6': [('A', 4)] * 6,
     'A5^4D4': [('A', 5)] * 4 + [('D', 4)], 'D4^6': [('D', 4)] * 6, 'A6^4': [('A', 6)] * 4,
     'A7^2D5^2': [('A', 7)] * 2 + [('D', 5)] * 2, 'A8^3': [('A', 8)] * 3, 'A9^2D6': [('A', 9)] * 2 + [('D', 6)],
     'D6^4': [('D', 6)] * 4, 'A11D7E6': [('A', 11), ('D', 7), ('E', 6)], 'E6^4': [('E', 6)] * 4, 'A12^2': [('A', 12)] * 2,
     'D8^3': [('D', 8)] * 3, 'A15D9': [('A', 15), ('D', 9)], 'A17E7': [('A', 17), ('E', 7)],
     'D10E7^2': [('D', 10), ('E', 7), ('E', 7)], 'D12^2': [('D', 12)] * 2, 'A24': [('A', 24)],
     'D16E8': [('D', 16), ('E', 8)], 'E8^3': [('E', 8)] * 3, 'D24': [('D', 24)]}
SX = {'D4^6': {3}, 'D6^4': {5}, 'D8^3': {7}, 'D10E7^2': {9}, 'D12^2': {11}, 'D16E8': {15}, 'D24': {23}, 'E6^4': {4},
      'E8^3': {6, 10, 15}}
def egcd(a, b):
    if b == 0: return a, 1, 0
    g, x, y = egcd(b, a % b); return g, y, x - (a // b) * y
def AL(N, e):
    f = N // e
    for x in range(1, 60):
        for z in range(1, 60):
            if gcd(e * x, f * z) == 1:
                g, s, t = egcd(e * x, f * z)
                return e * x, -t, N * z, e * s
ok = True
taus = [mp.mpc('0.1234', '0.6789'), mp.mpc('-0.2718', '0.5772')]
for name, comps in X.items():
    pi = {}
    for t, n in comps:
        for a, m in frame(t, n).items():
            pi[a] = pi.get(a, 0) + m
    pi = {a: m for a, m in pi.items() if m}
    assert sum(a * m for a, m in pi.items()) == 24 and sum(pi.values()) == 0
    h = max(pi)
    inv = set()
    ratios = {}
    for e in [e for e in range(2, h + 1) if h % e == 0 and gcd(e, h // e) == 1]:
        a, b, c, d = AL(h, e)
        rs = []
        for tau in taus:
            T = lambda t: 1 / mp.fprod([eta(k * t) ** m for k, m in pi.items()])
            rs.append(T((a * tau + b) / (c * tau + d)) / T(tau))
        if abs(rs[0] - 1) < mp.mpf(10) ** -12 and abs(rs[1] - 1) < mp.mpf(10) ** -12:
            inv.add(e)
        ratios[e] = mp.nstr(rs[0], 6)
    good = inv == SX.get(name, set())
    ok &= good
    print(f'{"ok  " if good else "FAIL"} {name:9s} h={h:2d}: W_e-invariant for e in {sorted(inv)} (expected {sorted(SX.get(name, set()))})')
print('ALL OK' if ok else 'FAIL')
