"""Even critical points of Chamberland's C (session opus-2026-10-06-S19, NUMERICAL; note section 2.6 (f); reproduce: python3 <this> 1200): do they (a) stay left-near-integer to the end, (b) flip into a right
tube, or (c) escape (|xi| > XI_ESC)?  Compare the actual high-precision orbit with the prediction of the scaled maps
(exact O(xi,h), E(xi,h) evaluated in mp floats, valid for negative xi too); the agreement is a numerical consistency check only, since the scaled maps are exact."""
import sys
from multiprocessing import Pool
from mpmath import mp, mpf, cos, sin, pi, findroot

XI_ESC = mpf('1.3')
A_TUBE = mpf('0.8')


def C(x):
    return x + mpf(1) / 4 - (2 * x + 1) / 4 * cos(pi * x)


def Cp(x):
    return 1 - cos(pi * x) / 2 + (2 * x + 1) * pi / 4 * sin(pi * x)


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def sinc(u):
    return 1 if u == 0 else sin(u) / u


def O(xi, h):
    S = sinc(pi * xi * h / 2)
    return (3 + h) / 2 * (xi * (1 + cos(pi * xi * h) / 2) - pi ** 2 * xi ** 2 / 8 * S ** 2)


def E(xi, h):
    S = sinc(pi * xi * h / 2)
    return (1 + h) / 2 * (xi * (1 - cos(pi * xi * h) / 2) + pi ** 2 * xi ** 2 / 8 * S ** 2)


def classify_actual(m):
    mp.dps = 400
    c = findroot(Cp, mpf(m) - mpf('0.2') / (2 * m + 1))
    x, k = c, m
    while k != 1:
        xi = (2 * k + 1) * (x - k)
        if abs(xi) > XI_ESC:
            return 'escape'
        if 0 <= xi <= A_TUBE:
            return 'flip'          # entered a right tube: shadows forever (THM-4563)
        x, k = C(x), T(k)
    return 'left-to-1'


def classify_scaled(m):
    mp.dps = 50
    c = findroot(Cp, mpf(m) - mpf('0.2') / (2 * m + 1))
    xi = (2 * m + 1) * (c - m)
    k = m
    while k != 1:
        if abs(xi) > XI_ESC:
            return 'escape'
        if 0 <= xi <= A_TUBE:
            return 'flip'
        h = mpf(1) / (2 * k + 1)
        xi = O(xi, h) if k % 2 else E(xi, h)
        k = T(k)
    return 'left-to-1'


def both(m):
    return m, classify_actual(m), classify_scaled(m)


if __name__ == '__main__':
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    with Pool(12) as pool:
        res = pool.map(both, range(2, M + 1, 2), chunksize=8)
    from collections import Counter
    print('actual:', Counter(r[1] for r in res))
    print('scaled-map prediction agrees on', sum(1 for r in res if r[1] == r[2]), 'of', len(res))
    bad = [r for r in res if r[1] != r[2]]
    print('disagreements (first 10):', bad[:10])
