#!/usr/bin/env python3
"""All rotation orders of a Paley tournament in O(p^2): the Farey sweep (S40).

The rotation order R_alpha lists x in {0, ..., p-1} by frac(alpha x) (vertex 0 first). As alpha runs over (0, 1) the order
changes only at Farey fractions a/b with b <= p - 1. Crossing a/b:
  * each residue class c != 0 (mod b) is a contiguous block, ordered by decreasing x just before and increasing x just
    after: it is reversed, changing the excess by 2 chi(b) A(k_c), where k_c is its size and A(k) = sum_(t<k) S(t) is the
    area of the Legendre walk S(t) = sum_(j<=t) (j/p);
  * the class of 0 (minus the vertex 0) moves from the end of the order to just after 0 and is reversed: change
    2 chi(b) (S(k_0) + A(k_0)) (the cross terms sum to 2 sum_(m <= k_0) chi(m b)).
With p - 1 = q b + s (0 <= s < b): k_0 = q, s classes of size q + 1 and b - 1 - s of size q, so the jump is

    Delta_b = 2 chi(b) [ S(q) + (b - s) A(q) + s A(q + 1) ],

independent of a. Hence the excess of R_alpha on the Farey interval after the n-th fraction is p h(-p) + the sum of the
first n jumps, the best rotation order is the maximal prefix sum, and the sweep must end at the reversed natural order:

    sum_(b=2)^(p-1) phi(b) Delta_b = -2 p h(-p).                                   (Farey class number identity)

This script (1) checks the jump formula against direct excess computation on every Farey interval for small p, (2) checks
the identity, (3) returns the best rotation order (a certificate alpha = mediant of the interval) for given primes.
Usage: python paley_rotation_farey_sweep_20261017.py [check PMAX] [scan p1 p2 ...] [range PMIN PMAX LMIN]
Session opus-2026-10-17-S40.
"""
import sys, os, json, math
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number, is_prime


@njit(cache=True)
def walk_tables(chi, p):
    S = np.zeros(p + 1, dtype=np.int64)
    A = np.zeros(p + 2, dtype=np.int64)
    for t in range(1, p):
        S[t] = S[t - 1] + chi[t]
    for k in range(1, p + 1):
        A[k + 1 if k + 1 <= p + 1 else p + 1] = 0
    acc = 0
    A[1] = 0
    for k in range(2, p + 1):
        acc += S[k - 1]
        A[k] = acc
    return S, A


@njit(cache=True)
def jump(b, p, chi, S, A):
    q = (p - 1) // b
    s = (p - 1) - q * b
    cb = chi[b % p]
    return 2 * cb * (S[q] + (b - s) * A[q] + s * A[q + 1])


@njit(cache=True)
def sweep(chi, p, h):
    """Farey sequence of order N = p - 1 from 0/1 to 1/1 by the next-term recurrence; returns the maximal prefix excess,
    the fraction pair bounding its interval, the final excess and the identity sum."""
    S, A = walk_tables(chi, p)
    N = p - 1
    a, b, c, d = 0, 1, 1, N
    E = p * h
    best = E; ba, bb, bc, bd = 0, 1, 1, N
    total = 0
    while c <= N and not (c == 1 and d == 1):
        # cross the fraction c/d
        dl = jump(d, p, chi, S, A)
        E += dl
        total += dl
        k = (N + b) // d
        na, nb = c, d
        nc, nd = k * c - a, k * d - b
        a, b, c, d = na, nb, nc, nd
        if E > best:
            best = E; ba, bb, bc, bd = a, b, c, d       # interval (a/b, c/d) after the crossing
    return best, ba, bb, bc, bd, E, total


@njit(cache=True)
def excess_of(order, chi, p):
    s = 0
    n = order.shape[0]
    for i in range(n):
        u = order[i]
        for j in range(i + 1, n):
            s += chi[(order[j] - u) % p]
    return s


def rotation_order(p, num, den):
    keys = [((num * x) % den, x) for x in range(p)]
    return np.array([x for _, x in sorted(keys)], dtype=np.int64)


def farey(N):
    a, b, c, d = 0, 1, 1, N
    out = [(0, 1)]
    while c <= N:
        k = (N + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        out.append((a, b))
    return out


def check(PMAX):
    bad_jump = 0; bad_id = 0; n = 0
    for p in [q for q in range(7, PMAX + 1, 4) if is_prime(q)]:
        chi = kronecker_table(-p); h = class_number(-p)
        S, A = walk_tables(chi, p)
        F = farey(p - 1)
        E = p * h
        for i in range(len(F) - 1):
            (a1, b1), (a2, b2) = F[i], F[i + 1]
            num, den = a1 + a2, b1 + b2                  # mediant: strictly inside the interval
            if i > 0:
                E += int(jump(b1, p, chi, S, A))
            if p <= 47:
                e_direct = int(excess_of(rotation_order(p, num, den), chi, p))
                if e_direct != E:
                    bad_jump += 1
        best, ba, bb, bc, bd, Efin, total = sweep(chi, p, h)
        if total != -2 * p * h or Efin != -p * h:
            bad_id += 1
        n += 1
    print('check: %d primes p = 3 (mod 4) up to %d; jump-formula failures (direct excess on every Farey interval, p <= 47): %d; identity failures: %d'
          % (n, PMAX, bad_jump, bad_id))


def scan(primes):
    rows = []
    for p in primes:
        chi = kronecker_table(-p); h = class_number(-p)
        best, ba, bb, bc, bd, Efin, total = sweep(chi, p, h)
        num, den = ba + bc, bb + bd
        e_cert = int(excess_of(rotation_order(p, num, den), chi, p))
        L = math.pi * h / math.sqrt(p)
        rows.append(dict(p=p, h=h, L=L, natural=p * h, best_rotation=int(best), alpha=[int(num), int(den)], certificate_excess=e_cert,
                         beaten=bool(best > p * h)))
        print('p = %d  L = %.3f  natural %d  best rotation %d at alpha = %d/%d (direct check %d)  %s'
              % (p, L, p * h, best, num, den, e_cert, 'BEATEN' if best > p * h else 'natural is the best rotation order'), flush=True)
    return rows


def main():
    args = sys.argv[1:]
    out = {}
    if args and args[0] == 'check':
        check(int(args[1]))
    elif args and args[0] == 'scan':
        out['scan'] = scan([int(a) for a in args[1:]])
    elif args and args[0] == 'range':
        PMIN, PMAX, LMIN = int(args[1]), int(args[2]), float(args[3])
        ps = []
        for p in range(PMIN, PMAX + 1):
            if p % 4 == 3 and is_prime(p):
                if math.pi * class_number(-p) / math.sqrt(p) >= LMIN:
                    ps.append(p)
        out['range'] = scan(ps)
        here = os.path.dirname(os.path.abspath(__file__))
        json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                         'paley_rotation_farey_sweep_%d_%d_20261017.json' % (PMIN, PMAX)), 'w'), indent=1)


if __name__ == '__main__':
    main()
