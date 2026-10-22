"""Audit check (S53): (a) the exponent law of THM-4759 (5') and the hand formulas in its proof, (b) the W_e component
rule including the multiplier (numerically, mpmath), (c) the Co_0 orbit-size arithmetic of the M24-point remark,
(d) the inner-product distribution of (2^24)/sqrt 8 with the 196560 minimal vectors from the Golay-code description."""
from fractions import Fraction
from math import gcd
import itertools
import mpmath as mp


def exps(t, n):
    if t == 'A':
        return list(range(1, n + 1)), n + 1
    if t == 'D':
        return sorted([2 * k + 1 for k in range(n - 1)] + [n - 1]), 2 * n - 2
    return {6: ([1, 4, 5, 7, 8, 11], 12), 7: ([1, 5, 7, 9, 11, 13, 17], 18), 8: ([1, 7, 11, 13, 17, 19, 23, 29], 30)}[n]


BAD = {'A': set(), 'D': {2}, 6: {2, 3}, 7: {2, 3}, 8: {2, 3, 5}}
ok_all = True
eq_cases = []
for t, rng in (('A', range(1, 301)), ('D', range(4, 301)), ('E', (6, 7, 8))):
    for n in rng:
        m, h = exps(t, n) if t != 'E' else exps('E', n)
        l = len(m)
        bad = BAD[t] if t != 'E' else BAD[n]
        for e in [e for e in range(1, h + 1) if h % e == 0]:
            S = sum((x % e) * (e - x % e) for x in m)
            rhs = Fraction(l * e * (e - 1), 6)
            p = h // e
            pred_eq = e == 1 or (p in bad and gcd(p, e) == 1)
            if S < rhs or (S == rhs) != pred_eq:
                ok_all = False
                print('EXPONENT LAW FAILS', t, n, e, S, rhs)
            if S == rhs and e > 1:
                eq_cases.append((t if t != 'E' else f'E{n}', n, e))
            # hand formulas of the proof
            if t == 'A' and e > 1:
                d = h // e
                assert S == d * e * (e * e - 1) // 6 and S - rhs == Fraction((e - 1) * (d + 1) * e, 6)
            if t == 'D' and e > 1:
                if e % 2:
                    assert (n - 1) % e == 0 and 6 * S == (n - 1) * (e * e - 1)
                    assert S - rhs == Fraction((e - 1) * (n - 1 - e), 6)
                else:
                    ep = e // 2; tt = (n - 1) // ep
                    assert (n - 1) % ep == 0
                    assert 3 * S == (n - 1) * (2 * ep * ep + 1) + 3 * ep * ep * (tt % 2)
                    br = tt + ep * tt + 1 - 2 * ep + 3 * (tt % 2) * ep
                    assert S - rhs == Fraction(ep * br, 3) and br > 0
            if e == h:
                assert 6 * sum(x * (h - x) for x in m) == l * h * (h + 1)          # strange formula
print('(a) exponent law and hand formulas:', 'OK' if ok_all else 'FAIL',
      '; equality cases e > 1:', sorted(set((c[0], c[2]) if c[0] != 'D' else ('D(n even)', 'n-1') for c in eq_cases)))
assert all(c[0] != 'D' or (c[1] % 2 == 0 and c[2] == c[1] - 1) for c in eq_cases)
assert not any(c[0] == 'A' for c in eq_cases)

# (b) component rule with multiplier: T = 1/eta_pi, test T(W_e tau) = T(tau) at a random point (W_e an Atkin-Lehner matrix)
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


def atkin_lehner(N, e):
    # W_e = [[e x, y], [N z, e w]] with det e: e x w - (N/e) y z = 1
    f = N // e
    for x in range(1, 50):
        for z in range(1, 50):
            # need e*x*w - f*y*z = 1 with w, y integers: solvable iff gcd(e x, f z) = 1
            if gcd(e * x, f * z) == 1:
                # solve e x w - f z y = 1
                g, s, tt = egcd(e * x, f * z)
                return (e * x, -tt, N * z, e * s)   # e x * s + f z * tt = 1 -> w = s, y = -tt
    raise ValueError


def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


tau = mp.mpc('0.1234', '0.6789')
rule_ok = True
cases = [('A', n) for n in range(1, 13)] + [('D', n) for n in range(4, 13)] + [('E', 6), ('E', 7), ('E', 8)]
for t, n in cases:
    pi = frame(t, n) if t != 'E' else frame('E', n)
    h = max(pi)
    for e in [e for e in range(2, h + 1) if h % e == 0 and gcd(e, h // e) == 1]:
        a, b, c, d = atkin_lehner(h, e)
        assert a * d - b * c == e
        tw = (a * tau + b) / (c * tau + d)
        inv = abs(Tval(pi, tw) - Tval(pi, tau)) < mp.mpf(10) ** -15
        bad = BAD[t] if t != 'E' else BAD[n]
        pred = (h // e) in bad
        if inv != pred:
            rule_ok = False
            print('COMPONENT RULE MISMATCH', t, n, e, inv, pred)
print('(b) component rule with exact multiplier at a random point: A1..A12, D4..D12, E6..E8:', 'OK' if rule_ok else 'FAIL')

# (c) orbit sizes
Co0 = 2 ** 22 * 3 ** 9 * 5 ** 4 * 7 ** 2 * 11 * 13 * 23
M24 = 2 ** 10 * 3 ** 3 * 5 * 7 * 11 * 23
U62 = 2 ** 15 * 3 ** 6 * 5 * 7 * 11
assert Co0 == 8315553613086720000 and M24 == 244823040 and U62 == 9196830720
o1, o2 = Co0 // M24, Co0 // (2 * U62)
print(f'(c) |Co0|/|M24| = {o1}, |Co0|/|U6(2):2| = {o2}, sum = {o1 + o2} (theta_Leech coefficient of q^6 is 34417656000: '
      f'{o1 + o2 == 34417656000})')

# (d) inner products of w = (2^24) with the minimal vectors (sqrt 8 scaling) via the Golay code
def golay():
    # extended binary Golay code from the generator of the [23,12] cyclic code, x^11 + x^10 + x^6 + x^5 + x^4 + x^2 + 1
    g = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]
    rows = []
    for s in range(12):
        r = [0] * 23
        for i, b in enumerate(g):
            r[(i + s) % 23] ^= b
        r.append(sum(r) % 2)
        rows.append(r)
    words = set()
    for mask in range(1 << 12):
        w = [0] * 24
        for i in range(12):
            if mask >> i & 1:
                w = [x ^ y for x, y in zip(w, rows[i])]
        words.add(tuple(w))
    return words


C = golay()
assert len(C) == 4096 and sorted(set(sum(w) for w in C)) == [0, 8, 12, 16, 24]
from collections import Counter
dist = Counter()
# shape (+-4, +-4, 0^22)
for i, j in itertools.combinations(range(24), 2):
    for si in (4, -4):
        for sj in (4, -4):
            dist[2 * (si + sj) // 8] += 1
# shape (+-2^8 on an octad, even number of minus signs)
for w in C:
    if sum(w) == 8:
        for k in range(0, 9, 2):
            from math import comb
            dist[2 * 2 * (8 - 2 * k) // 8] += comb(8, k)
# shape (-+3, +-1^23): x = (-3 at i0, +1 elsewhere) with signs flipped on a codeword S
for w in C:
    for i0 in range(24):
        x = [1] * 24; x[i0] = -3
        x = [-xi if w[i] else xi for i, xi in enumerate(x)]
        dist[2 * sum(x) // 8] += 1
assert sum(dist.values()) == 196560
print('(d) (2^24)/sqrt8 . lambda over the 196560 minimal vectors:', dict(sorted(dist.items())))
