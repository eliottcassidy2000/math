#!/usr/bin/env python3
"""Excesses of arithmetic orders of the Paley tournament P_p are class-number multiples (S40, THM-4666 companion).

Orders of Z/p (vertex 0 first, then the listed order of the nonzero residues):
  natural       1, 2, ..., p-1                      excess p h(-p)          (THM-4666 (1))
  discrete log  g^0, g^1, ..., g^(p-2)              excess 1 for every primitive root g (Jacobsthal: proof in the note)
  inverse       1^-1, 2^-1, ..., (p-1)^-1           excess 3 h(-p)                      (observed, all p <= 3000)
  square        sorted by x^2 mod p (x before p-x)  excess (1 - 2 (2/p)) h(-p)  (PROVED for every pair-block order: cross terms
                                                    between two {x, -x} pairs cancel; Dirichlet half sum)
(dlog_other_root uses g^e with e the least integer >= 3 prime to p - 1, another primitive root.)
For every prime p = 3 (mod 4) up to PMAX the script computes these excesses exactly and fits them against h(-p),
separately for p = 3 and p = 7 (mod 8).
Usage: python paley_arithmetic_orders_20261017.py PMAX
Session opus-2026-10-17-S40.
"""
import sys, os, json
import numpy as np
from numba import njit
from sympy import primitive_root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number, is_prime


@njit(cache=True)
def excess_of(order, chi, p):
    s = 0
    n = order.shape[0]
    for i in range(n):
        u = order[i]
        for j in range(i + 1, n):
            s += chi[(order[j] - u) % p]
    return s


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    stats = {}
    rows = []
    for p in range(7, PMAX + 1, 4):
        if not is_prime(p):
            continue
        chi = kronecker_table(-p)
        h = class_number(-p)
        g = primitive_root(p)
        orders = {
            'natural': [0] + list(range(1, p)),
            'dlog': [0] + [pow(g, k, p) for k in range(p - 1)],
            'dlog_other_root': [0] + [pow(pow(g, next(e for e in range(3, p) if np.gcd(e, p - 1) == 1), p), k, p) for k in range(p - 1)],
            'inverse': [0] + [pow(k, -1, p) for k in range(1, p)],
            'square': sorted(range(p), key=lambda x: ((x * x) % p, x)),
        }
        ex = {k: int(excess_of(np.array(v, dtype=np.int64), chi, p)) for k, v in orders.items()}
        rows.append(dict(p=p, h=h, cls=p % 8, **ex))
        for k, v in ex.items():
            key = (k, p % 8)
            stats.setdefault(key, set()).add(v / h if k != 'dlog' and k != 'dlog_other_root' else v)
    for key in sorted(stats):
        vals = sorted(stats[key])
        print('%-16s p = %d (mod 8): %s %s' % (key[0], key[1], 'excess / h in' if key[0] not in ('dlog', 'dlog_other_root') else 'excess in',
                                           vals if len(vals) <= 6 else (vals[:3], '...', vals[-3:], '%d distinct' % len(vals))))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_arithmetic_orders_20261017.json'), 'w'))


if __name__ == '__main__':
    main()
