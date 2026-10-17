#!/usr/bin/env python3
"""Faithful realizations of the Exoo-Ismailescu graph EI17 over finite fields, and what they say about its field (S39).

Direction structure (from the certified realization, as in ei17_placement_q_zeta55_20261012.py): 13 unit directions,
gauge U5 = 1, free U0, U1, U2, U3, U4, U6 and
    U7 = U1 + U2 - U3 - 1,          U8 = -U0 - U2 + U4 + 1,         U9 = -U0 + U1 - U2 + U3 + U6 + 1,
    U10 = 2U2 - U3 - U6 - 1,        U11 = U1 + U2 + U3 - U4 + U6 - 1,   U12 = -U1 - U2 - U6 + 2,
each dependent direction of norm one.  Vertex positions are the integer combinations printed by the structure step.

Split model over F_p: a direction is x in F_p^x with conjugate 1/x; a point is a pair (z, zbar); norm z zbar.  Reduction
modulo a prime of Q(zeta55) above p = 1 (mod 55) (residue field F_p, conjugation sigma_-1 acting as x -> 1/x on units)
sends a realization with all directions in Q(zeta55) to an F_p-point of this model with 17 distinct vertices, whenever the
prime is of good reduction for it.  Elimination: choose U1, U2; U3 from |U7| = 1 and U6 from |U12| = 1 (quadratics);
check |U10| = 1; U4 from |U11| = 1, U0 from |U8| = 1 (quadratics); check |U9| = 1.  Count solutions with 17 distinct
vertices ('faithful').
Usage: python ei17_finite_field_realizations_20261017.py PMAX [comma-separated prime list]
Session opus-2026-10-17-S39.
"""
import sys
import numpy as np
from numba import njit
from sympy import primerange

# vertex v = c0 + c1 U0 + c2 U1 + c3 U2 + c4 U3 + c5 U4 + c6 U6 (from the structure step)
VERT = {1: (-1, 0, -1, 0, -1, 0, -1), 2: (0, 0, -1, -1, -1, 0, -1), 3: (0, 0, -1, -1, 0, 0, -1), 4: (-1, 1, -1, 0, -1, 0, -1),
        5: (0, 0, -1, -1, 0, 0, 0), 6: (0, 0, -1, -1, -1, 1, -1), 7: (-1, 0, 0, 0, -1, 0, -1), 8: (-1, 0, -1, 1, -1, 0, -1),
        9: (-1, 0, 0, 0, -1, 0, 0), 10: (1, 0, -1, -1, -1, 0, -1), 11: (-1, 0, 0, 0, 0, 0, 0), 12: (1, 0, -1, -1, 0, 0, -1),
        13: (-1, 1, -1, 1, -1, 0, -1), 14: (0, 0, 0, -1, 0, 0, 0), 15: (-1, 0, 0, 1, -1, 0, -1), 16: (1, 0, -1, -1, -1, 1, -1),
        17: (0, 0, 0, 0, 0, 0, 0)}
VM = np.array([VERT[v] for v in range(1, 18)], dtype=np.int64)
PE = [(1, 4), (1, 7), (1, 8), (2, 3), (2, 6), (2, 10), (3, 5), (3, 7), (3, 12), (4, 6), (4, 13), (4, 14), (5, 8),
      (5, 9), (5, 14), (6, 11), (6, 16), (7, 9), (7, 15), (8, 13), (8, 15), (9, 10), (9, 11), (10, 12), (10, 16),
      (11, 12), (11, 17), (13, 17), (14, 15), (14, 17), (16, 17)]
EDGES = np.array([(i - 1, j - 1) for i, j in PE], dtype=np.int64)


@njit(cache=True)
def inv(a, p):
    return pow_mod(a % p, p - 2, p)


@njit(cache=True)
def pow_mod(a, e, p):
    r = 1
    a %= p
    while e > 0:
        if e & 1:
            r = r * a % p
        a = a * a % p
        e >>= 1
    return r


@njit(cache=True)
def unit_solutions(c, cb, p, sqrt_tab, out):
    """(see below); when the equation degenerates (c = cb = 0: every x works) all p - 1 values are written to out."""
    c2 = c % p; cb2 = cb % p
    if cb2 == 0 and c2 == 0:
        for x in range(1, p):
            out[x - 1] = x
        return p - 1
    return unit_solutions_nd(c, cb, p, sqrt_tab, out)


@njit(cache=True)
def unit_solutions_nd(c, cb, p, sqrt_tab, out):
    """x in F_p^x with (c - x)(cb - 1/x) = 1, i.e. cb x^2 - c cb x + c = 0; writes up to 2 roots into out, returns count."""
    c %= p; cb %= p
    if cb == 0:
        # (c - x)(-1/x) = 1  ->  c = 0 (handled above as degenerate); otherwise no solution
        return 0
    D = (c * cb % p) * ((c * cb - 4) % p) % p
    if D == 0:
        x = (c * cb % p) * inv(2 * cb, p) % p
        if x == 0:
            return 0
        out[0] = x
        return 1
    s = sqrt_tab[D]
    if s < 0:
        return 0
    k = 0
    i2 = inv(2 * cb, p)
    for sg in (1, -1):
        x = (c * cb + sg * s) % p * i2 % p
        if x != 0:
            out[k] = x
            k += 1
    return k


@njit(cache=True)
def is_unit(z, zb, p):
    return (z % p) * (zb % p) % p == 1


@njit(cache=True)
def count(p, sqrt_tab, VM, EDGES, edge_fail):
    faithful = 0
    total = 0
    degenerate = 0
    r3 = np.zeros(p, dtype=np.int64); r6 = np.zeros(p, dtype=np.int64)
    r4 = np.zeros(p, dtype=np.int64); r0 = np.zeros(p, dtype=np.int64)
    sols = np.zeros((64, 6), dtype=np.int64)
    zz = np.zeros(17, dtype=np.int64); zb = np.zeros(17, dtype=np.int64)
    for U1 in range(1, p):
        i1 = inv(U1, p)
        for U2 in range(1, p):
            i2 = inv(U2, p)
            # U7 = (U1 + U2 - 1) - U3
            n3 = unit_solutions(U1 + U2 - 1, i1 + i2 - 1, p, sqrt_tab, r3)
            n6 = unit_solutions(2 - U1 - U2, 2 - i1 - i2, p, sqrt_tab, r6)
            if n3 == p - 1 or n6 == p - 1:
                degenerate += 1
            for a in range(n3):
                U3 = r3[a]; i3 = inv(U3, p)
                for b in range(n6):
                    U6 = r6[b]; i6 = inv(U6, p)
                    if not is_unit(2 * U2 - U3 - U6 - 1, 2 * i2 - i3 - i6 - 1, p):
                        continue
                    # U11 = (U1 + U2 + U3 + U6 - 1) - U4
                    n4 = unit_solutions(U1 + U2 + U3 + U6 - 1, i1 + i2 + i3 + i6 - 1, p, sqrt_tab, r4)
                    if n4 == p - 1:
                        degenerate += 1
                    for c4 in range(n4):
                        U4 = r4[c4]; i4 = inv(U4, p)
                        # U8 = (1 - U2 + U4) - U0
                        n0 = unit_solutions(1 - U2 + U4, 1 - i2 + i4, p, sqrt_tab, r0)
                        if n0 == p - 1:
                            degenerate += 1
                        for d0 in range(n0):
                            U0 = r0[d0]; i0 = inv(U0, p)
                            if not is_unit(-U0 + U1 - U2 + U3 + U6 + 1, -i0 + i1 - i2 + i3 + i6 + 1, p):
                                continue
                            total += 1
                            vals = (1, U0, U1, U2, U3, U4, U6)
                            ivals = (1, i0, i1, i2, i3, i4, i6)
                            for v in range(17):
                                s = 0; sb = 0
                                for t in range(7):
                                    s += VM[v, t] * vals[t]
                                    sb += VM[v, t] * ivals[t]
                                zz[v] = s % p; zb[v] = sb % p
                            ok = True
                            for v in range(17):
                                for w in range(v + 1, 17):
                                    if zz[v] == zz[w] and zb[v] == zb[w]:
                                        ok = False
                                        break
                                if not ok:
                                    break
                            if ok:
                                for e in range(EDGES.shape[0]):
                                    a_, b_ = EDGES[e, 0], EDGES[e, 1]
                                    if (zz[a_] - zz[b_]) % p * ((zb[a_] - zb[b_]) % p) % p != 1:
                                        edge_fail[0] += 1
                                if faithful < 64:
                                    sols[faithful, 0] = U0; sols[faithful, 1] = U1; sols[faithful, 2] = U2
                                    sols[faithful, 3] = U3; sols[faithful, 4] = U4; sols[faithful, 5] = U6
                                faithful += 1
    return faithful, total, degenerate, sols


def sqrt_table(p):
    t = -np.ones(p, dtype=np.int64)
    for x in range(p):
        y = x * x % p
        if t[y] < 0:
            t[y] = x
    return t


if __name__ == '__main__':
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    res = {}
    print('  p  p mod 55  faithful  all-solutions  degenerate-branches   (every faithful solution checked on all 31 edges)')
    plist = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else list(primerange(7, PMAX))
    for p in plist:
        if p in (11, 5):
            continue
        ef = np.zeros(1, dtype=np.int64)
        f, tot, deg, sols = count(p, sqrt_table(p), VM, EDGES, ef)
        assert ef[0] == 0, ('edge check failed', p)
        res[p] = (f, tot, deg)
        print('%5d  %3d   %6d   %8d   %8d' % (p, p % 55, f, tot, deg), flush=True)
    ones = [p for p in res if p % 55 == 1]
    print('primes = 1 mod 55:', [(p, res[p][0]) for p in ones])
    mean = sum(res[p][0] for p in res) / len(res)
    print('mean number of faithful F_p-realizations over all primes: %.3f' % mean)
