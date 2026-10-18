#!/usr/bin/env python3
"""Frontier-energy profile E_Y(N)/N as a function of the frontier height Y, at fixed N (S42).

Companion to root_transport_frontier_energy_20261017.py. For a few N it computes, for a geometric range of Y, the
frontier coefficients c_Y(u) (first seed <= Y on each source path) and E_Y(N) = sum_u c_Y(u)^2, together with the
Cauchy-Schwarz bound F_Y E_Y >= R(N)^2. As Y decreases from 4N to 1 the profile runs from the initial energy E_0 down
through the ideal regime; Y = 1 gives R(N)^2 itself. This is the scale-parametrised analogue of the depth profile E_d.
Usage: python root_transport_frontier_profile_20261017.py N1 N2 ...
Session opus-2026-10-17-S42.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_energy_profile_20261017 import initial
from root_transport_frontier_energy_20261017 import step


def frontier(S, C, Y):
    acc = {}
    for s, c in zip(S.tolist(), C.tolist()):
        x = s; sign = 1
        while x > Y:
            x, sg = step(x)
            sign *= sg
        acc[x] = acc.get(x, 0) + sign * c
    return acc


def parity(n):
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        c += 1
    return 1 - 2 * (c & 1)


def main():
    out = []
    for N in [int(a) for a in sys.argv[1:]]:
        S, C = initial(N)
        Ys = sorted(set([1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 10 ** 6, 3 * 10 ** 6, 4 * N]))
        Ys = [y for y in Ys if y <= 4 * N]
        R = None
        print('N = %d' % N)
        for Y in Ys:
            acc = frontier(S, C, Y)
            E = sum(v * v for v in acc.values())
            Rv = sum(v * parity(u) for u, v in acc.items())
            R = Rv if R is None else R
            assert Rv == R
            FY = sum(1 for s in range(1, Y + 1, 2) if s % 8 != 5)
            out.append(dict(N=N, Y=Y, E_over_N=E / N, F_Y=FY, R=Rv))
            print('   Y = %8d: E_Y/N = %.4f   F_Y = %6d   frontier seeds hit %6d   R^2/(F_Y E_Y) = %.4f'
                  % (Y, E / N, FY, len(acc), Rv * Rv / max(FY * E, 1)), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_frontier_profile_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
