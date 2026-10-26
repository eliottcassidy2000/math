#!/usr/bin/env python3
"""Pollack-Pomerance at finite truncation y, and the effective truncation y_eff(N) of the data (opus-2026-10-26-S57).

PP (3.1) at finite y:  Delta_y = sum over even a | A_y of dens(T_a) * exp(-a/s(a)),  dens(T_a) = phi(A_y/a)/A_y
(exact at finite y; PP's phi(A_y)/(A_y a) is its limit as y -> oo),
T_a = {n : gcd(n, A_y) = a}, A_y = lcm(1..y).  PP prove Delta_q P_q-type monotonicity: Delta_y increases to Delta.
We tabulate Delta_y for y = 2..Y, and match the empirical local untouchable densities:
  - from this session's exact data: windows [N/2, N] for N = 2^e (e = 12..29) (mult_even_1e9.npy);
  - from OEIS A057978 (counts to 10^12): decade windows [10^(k-1), 10^k].
Usage: python pp_truncated_20261026.py bigdir [Y]
"""
import sys, os, json, math
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
A057978 = {7: 1574973, 8: 16246940, 9: 165826606, 10: 1681871718, 11: 16988116409, 12: 171128671374}


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]


def delta_y(y):
    ps = primes_upto(y)
    exps = []
    for p in ps:
        e = 0
        while p ** (e + 1) <= y:
            e += 1
        exps.append(e)
    Ay = 1
    for p, e in zip(ps, exps):
        Ay *= p ** e
    phi_ratio = 1.0
    for p in ps:
        phi_ratio *= (1 - 1 / p)
    # enumerate divisors a of A_y with v_2(a) >= 1, tracking sigma(a)
    total = 0.0
    stack = [(0, 1, 1)]          # (prime index, a, sigma(a))
    # iterative DFS over primes
    divs = [(1, 1)]
    for p, e in zip(ps, exps):
        new = []
        for a, sg in divs:
            pk, sp = 1, 1
            for j in range(e + 1):
                if j > 0:
                    pk *= p
                    sp += pk
                if p == 2 and j == 0:
                    continue
                new.append((a * pk, sg * sp))
        divs = new
    # exact finite-y class density: dens(T_a) = phi(A_y/a)/A_y  (PP's phi(A_y)/(a A_y) is its y -> oo limit)
    def phi(n):
        r = n
        for p in ps:
            if n % p == 0:
                r = r // p * (p - 1)
        return r
    wsum = 0.0
    for a, sg in divs:
        w = phi(Ay // a) / Ay
        wsum += w
        total += w * math.exp(-a / (sg - a))
    assert abs(wsum - 0.5) < 1e-9, wsum
    return total, len(divs)


def local_density_from_mult(mult, lo, hi):
    """even untouchables in (lo, hi] divided by (hi - lo); odd untouchables (only 5) ignored for lo >= 5."""
    j0, j1 = lo // 2 + 1, hi // 2
    return np.count_nonzero(mult[j0:j1 + 1] == 0) / (hi - lo)


if __name__ == '__main__':
    bigdir = sys.argv[1]
    Y = int(sys.argv[2]) if len(sys.argv) > 2 else 47
    tab = {}
    for y in range(2, Y + 1):
        if y > 2 and all(y % p for p in range(2, y)) or y in (4, 8, 9, 16, 25, 27, 32):
            d, nd = delta_y(y)
            tab[y] = d
            print(f'y={y:>3}: Delta_y = {d:.6f}  ({nd} even divisors of A_y)', flush=True)
        elif y == 2:
            d, nd = delta_y(2)
            tab[2] = d
            print(f'y=  2: Delta_y = {d:.6f}', flush=True)
    mult = np.load(os.path.join(bigdir, 'mult_even_1000000000.npy'))
    emp = []
    for e in range(12, 30):
        N = 2 ** e
        dens = local_density_from_mult(mult, N // 2, N)
        emp.append(('2^%d' % e, N, dens))
    for k in range(8, 13):
        dens = (A057978[k] - A057978[k - 1]) / (9 * 10 ** (k - 1))
        emp.append(('10^%d (A057978 decade)' % k, 10 ** k, dens))
    ys = sorted(tab)
    print('\nempirical local density and the effective truncation y_eff (Delta_y interpolated in y):')
    out_rows = []
    for name, N, dens in emp:
        # find y with Delta_y = dens by linear interpolation on the tabulated (increasing) values
        yeff = None
        for i in range(len(ys) - 1):
            a, b = ys[i], ys[i + 1]
            if tab[a] <= dens <= tab[b]:
                yeff = a + (dens - tab[a]) / (tab[b] - tab[a]) * (b - a)
                break
        out_rows.append({'window_end': name, 'N': N, 'local_density': dens, 'y_eff': yeff})
        print(f'   {name:>24}: local density {dens:.6f}   y_eff = {yeff if yeff is None else round(yeff, 2)}   '
              f'(log N = {math.log(N):.2f}, log log N = {math.log(math.log(N)):.3f})', flush=True)
    json.dump({'Delta_y': tab, 'empirical': out_rows}, open(os.path.join(HERE, 'pp_truncated.json'), 'w'), indent=1)
