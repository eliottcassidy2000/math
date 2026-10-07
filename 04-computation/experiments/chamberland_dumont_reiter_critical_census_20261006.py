"""Census of critical-orbit fates for Chamberland's C and Dumont-Reiter's D (session opus-2026-10-06-S19, NUMERICAL).
Note: 05-knowledge/results/collatz_cycles_tubes_debt_walk_openai_20261006.md, section 2.6.  Reproduce: python3 <this> 2001 1500
Odd critical points: track scaled deviation xi_t = (2T^t(n)+1)(F^t(c_n) - T^t(n)) until T^t(n) = 1, record xi at entry to tau(1)
and the final attractor (C: A1 if xi_entry basin..., determined by further iteration).
Even critical points (left of even m): iterate with high precision; record fate (A1, A2, 0, UP, UNDEC) and whether the orbit
ever leaves the near-integer regime |xi| <= 1.2."""
import sys
from multiprocessing import Pool
from mpmath import mp, mpf, cos, sin, pi, findroot

MAXN = int(sys.argv[1]) if len(sys.argv) > 1 else 400
DPS_ODD = 60
DPS_EVEN = int(sys.argv[2]) if len(sys.argv) > 2 else 1500


def Cf(x):
    return x + mpf(1) / 4 - (2 * x + 1) / 4 * cos(pi * x)


def Cfp(x):
    return 1 - cos(pi * x) / 2 + (2 * x + 1) * pi / 4 * sin(pi * x)


def Df(x):
    s = sin(pi * x / 2) ** 2
    return (mpf(3) ** s * x + s) / 2


def Dfp(x):
    s = sin(pi * x / 2) ** 2
    sp = pi / 2 * sin(pi * x)
    return (mpf(3) ** s * (1 + x * mp.log(3) * sp) + sp) / 2


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


A2C = None


def fate(F, x, maxit=6000):
    # attractors: A1 = {1,2}; C has A2 = {1.19253190704664, 2.13865633551671}; 0
    A2 = (mpf('1.19253190704664'), mpf('2.13865633551671'))
    for it in range(maxit):
        if abs(x - 1) < mpf('1e-25') or abs(x - 2) < mpf('1e-25'):
            return 'A1', it
        if F is Cf and (abs(x - A2[0]) < mpf('1e-12') or abs(x - A2[1]) < mpf('1e-12')):
            return 'A2', it
        if x < mpf('0.25'):
            return '0', it
        if x > mpf(10) ** 30:
            return 'UP', it
        x = F(x)
    return 'UNDEC', maxit


def odd_case(args):
    name, n = args
    mp.dps = DPS_ODD
    F, Fp = (Cf, Cfp) if name == 'C' else (Df, Dfp)
    c = findroot(Fp, mpf(n) + mpf('0.3') / (2 * n + 1) * (2 if name == 'C' else 1.2))
    x, m = c, n
    xi_max = 0
    while m != 1:
        xi = (2 * m + 1) * (x - m)
        xi_max = max(xi_max, float(xi))
        if not (-1e-20 <= float(xi) <= 0.9):
            return (name, n, 'LEFT_TUBE', float(xi), None)
        x, m = F(x), T(m)
    xi_entry = float(3 * (x - 1))
    fa, it = fate(F, x)
    return (name, n, fa, xi_entry, xi_max)


def even_case(args):
    name, m0 = args
    mp.dps = DPS_EVEN
    F, Fp = (Cf, Cfp) if name == 'C' else (Df, Dfp)
    c = findroot(Fp, mpf(m0) - mpf('0.2') / (2 * m0 + 1) * (1 if name == 'C' else 1.8))
    x, m = c, m0
    escaped_at = None
    for t in range(400):
        if m == 1:
            break
        xi = (2 * m + 1) * (x - m)
        if escaped_at is None and abs(float(xi)) > 1.2:
            escaped_at = t
        x, m = F(x), T(m)
    fa, it = fate(F, x)
    return (name, m0, fa, escaped_at)


if __name__ == '__main__':
    with Pool(12) as pool:
        odd = pool.map(odd_case, [(nm, n) for nm in ('C', 'D') for n in range(1, MAXN + 1, 2)], chunksize=4)
        even = pool.map(even_case, [(nm, m) for nm in ('C', 'D') for m in range(2, min(MAXN, 300) + 1, 2)], chunksize=2)
    from collections import Counter
    for nm in ('C', 'D'):
        rows = [r for r in odd if r[0] == nm]
        cnt = Counter(r[2] for r in rows)
        print(f'{nm} odd critical points n <= {MAXN}: {dict(cnt)}')
        a2 = [r[1] for r in rows if r[2] == 'A2']
        print(f'   A2-captured odd n: {a2[:40]}')
        print(f'   max xi along orbits: {max(r[4] for r in rows if r[4] is not None):.4f}; entry xi range '
              f'[{min(r[3] for r in rows):.3g}, {max(r[3] for r in rows):.3g}]')
        erows = [r for r in even if r[0] == nm]
        ecnt = Counter(r[2] for r in erows)
        esc = sum(1 for r in erows if r[3] is not None)
        print(f'{nm} even critical points m <= {min(MAXN,300)}: fates {dict(ecnt)}; escaped near-integer regime: {esc}/{len(erows)}')
