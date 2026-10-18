#!/usr/bin/env python3
"""Does the class number decide whether the natural order is an optimal ranking of the Paley tournament? (S40)

The natural order 0 < 1 < ... < p-1 of P_p has forward-minus-backward excess p h(-p) (THM-4666). This script runs an
iterated insertion local search (best reinsertion of each vertex, random kicks, restarts from random orders and from the
natural order) on P_p for primes p = 3 (mod 4) up to PMAX, and reports, against L(1, chi_(-p)) = pi h(-p) / sqrt(p), the
best excess found and whether the natural order was beaten. The best excess found is a certified lower bound for the
maximum excess 2 MAS(P_p) - p(p-1)/2 (the order is printed for the record holders); it is not an upper bound.
Usage: python paley_ranking_local_search_20261017.py [PMIN:]PMAX [RESTARTS] [SWEEPS] [extra primes ...]
Runs: '131 30 150' (exploratory, p <= 131) and '132:500 12 80' (pre-registered, paley_ranking_preregistration_20261017.md).
Session opus-2026-10-17-S40.
"""
import sys, os, json, math
import numpy as np
from numba import njit
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


@njit(cache=True)
def improve(order, chi, p):
    """Best-improvement insertion sweeps until no vertex can be moved with positive gain. Returns the excess."""
    n = order.shape[0]
    improved = True
    while improved:
        improved = False
        for i in range(n):
            v = order[i]
            # gains of moving v from position i to position k
            best_gain = 0
            best_k = i
            g = 0
            for k in range(i + 1, n):          # v jumps over order[k], which it preceded
                g += -2 * chi[(order[k] - v) % p]
                if g > best_gain:
                    best_gain = g; best_k = k
            g = 0
            for k in range(i - 1, -1, -1):     # v jumps over order[k], which it followed
                g += -2 * chi[(v - order[k]) % p]
                if g > best_gain:
                    best_gain = g; best_k = k
            if best_gain > 0:
                if best_k > i:
                    for t in range(i, best_k):
                        order[t] = order[t + 1]
                else:
                    for t in range(i, best_k, -1):
                        order[t] = order[t - 1]
                order[best_k] = v
                improved = True
    return excess_of(order, chi, p)


@njit(cache=True)
def search(chi, p, restarts, sweeps, seed):
    np.random.seed(seed)
    best = -10 ** 9
    best_order = np.arange(p)
    nat = np.arange(p)
    nat_local = improve(nat.copy(), chi, p)
    for r in range(restarts):
        if r % 4 == 0:
            order = np.arange(p)
        else:
            order = np.random.permutation(p)
        cur = improve(order, chi, p)
        cur_order = order.copy()
        for s in range(sweeps):
            trial = cur_order.copy()
            # kick: a few random block moves
            for t in range(1 + p // 20):
                i = np.random.randint(p); k = np.random.randint(p)
                v = trial[i]
                if k > i:
                    for q in range(i, k):
                        trial[q] = trial[q + 1]
                else:
                    for q in range(i, k, -1):
                        trial[q] = trial[q - 1]
                trial[k] = v
            val = improve(trial, chi, p)
            if val >= cur:
                cur = val; cur_order = trial
        if cur > best:
            best = cur; best_order = cur_order.copy()
    return best, best_order, nat_local


def main():
    rng_arg = sys.argv[1]
    PMIN, PMAX = (int(rng_arg.split(':')[0]), int(rng_arg.split(':')[1])) if ':' in rng_arg else (3, int(rng_arg))
    RESTARTS = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    SWEEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    extra = [int(x) for x in sys.argv[4:]]
    primes = [p for p in range(3, PMAX + 1, 4) if is_prime(p) and p >= PMIN] + extra
    rows = []
    print('    p   h   L(1,chi)  natural  local-from-natural  best found  best/p^1.5  natural beaten')
    for p in primes:
        chi = kronecker_table(-p)
        h = class_number(-p)
        w = 6 if p == 3 else 2
        nat = 2 * p * h // w
        best, order, nat_local = search(chi, p, RESTARTS, SWEEPS, 12345 + p)
        assert excess_of(order, chi, p) == best
        L = math.pi * h / math.sqrt(p) if p > 3 else math.pi / (3 * math.sqrt(3))   # w = 6 at p = 3 (fixed after the S40 audit; earlier outputs print 1.8138)
        rows.append(dict(p=p, h=h, L=L, natural=nat, natural_local=int(nat_local), best=int(best),
                         best_scaled=best / p ** 1.5, beaten=bool(best > nat), order=[int(x) for x in order]))
        print('%5d %3d %9.4f %8d %12d %14d %10.4f   %s' % (p, h, L, nat, nat_local, best, best / p ** 1.5,
                                                         'YES' if best > nat else 'no'), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    tag = '' if PMIN <= 3 else '_p%d_%d' % (PMIN, PMAX)
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_ranking_local_search%s_20261017.json' % tag), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
