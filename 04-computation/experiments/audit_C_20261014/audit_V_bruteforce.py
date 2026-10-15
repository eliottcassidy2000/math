#!/usr/bin/env python3
"""Audit C: independent brute force of V_q = I_{K6-P6}(q) with NO use of the audited reduction.

Definition (paley_next_cores_20261014i.md, eq. (4)): vertices 0..5, edges of K6 minus the path 0-1-2-3-4-5,
each oriented lower label -> higher label; V_q = sum over injective x: {0..5} -> F_q of prod chi(x_b - x_a).

Methods (injectivity enforced explicitly on all 15 pairs, not via chi(0) = 0):
  * 'full'  : enumerate every injective 6-tuple (no normalisation)                       q = 7, 11, 19, 23
  * 'trans' : fix x0 = 0 only (translation), multiply by q                               q = 31, 43, 47, 9, 25, 27
  * 'affine': fix x0 = 0, x1 = 1, multiply by q(q-1) (only to reach larger q)            q = 59, 67, 71, 49
Then compare with (i) the recorded values V_7 = 336, V_11 = -12320, V_19 = -76608, V_31 = 52080 and
(ii) THM-4653: V_q = q(q-1)(-eps M_q + 7 + 13 eps - 2 eps q), with M_q from this audit's own evaluators.
Also checks the reverse orientation gives the same value (10 edges: reversal multiplies by eps^10 = 1) and that
the source's 'canonical representative = -V' is not the natural orientation (reversing ONE edge flips the sign).

Usage: python3 -u audit_V_bruteforce.py > audit_V_bruteforce.out
"""
import sys
import time
import itertools
import numpy as np
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014')
from audit_common import chi_prime, M_exact_direct, GF, M_field, W_from_M

EDGES = [(a, b) for a in range(6) for b in range(a + 1, 6) if b != a + 1]
assert len(EDGES) == 10 and EDGES == [(0, 2), (0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 4), (2, 5), (3, 5)]
RECORDED = {7: 336, 11: -12320, 19: -76608, 31: 52080}


def brute(q, sub, chi, mode, edges=EDGES):
    """sub[a, b] = a - b (field indices); chi[...] quadratic character. Returns V_q (int)."""
    # loop over (x0, x1, x2) in python, vectorise (x3, x4, x5) as a q^3 grid
    X3 = np.arange(q)[:, None, None]
    X4 = np.arange(q)[None, :, None]
    X5 = np.arange(q)[None, None, :]
    distinct_345 = (X3 != X4) & (X3 != X5) & (X4 != X5)
    total = 0
    if mode == 'full':
        r0, r1 = range(q), None
    elif mode == 'trans':
        r0, r1 = [0], None
    elif mode == 'affine':
        r0, r1 = [0], [1]
    else:
        raise ValueError(mode)
    for x0 in r0:
        for x1 in (r1 if r1 is not None else range(q)):
            if x1 == x0:
                continue
            for x2 in range(q):
                if x2 in (x0, x1):
                    continue
                vals = {0: x0, 1: x1, 2: x2}
                # scalar part (edges among 0,1,2): only (0,2)
                w = 1
                for (a, b) in edges:
                    if a <= 2 and b <= 2:
                        w *= int(chi[sub[vals[b], vals[a]]])
                if w == 0:
                    continue
                prod = np.ones((q, q, q), dtype=np.int64) * w
                grids = {3: X3, 4: X4, 5: X5}
                for (a, b) in edges:
                    if b <= 2:
                        continue
                    if a <= 2:
                        prod = prod * chi[sub[grids[b], vals[a]]]
                    else:
                        prod = prod * chi[sub[grids[b], grids[a]]]
                mask = distinct_345.copy()
                for v in (x0, x1, x2):
                    mask &= (X3 != v) & (X4 != v) & (X5 != v)
                total += int((prod * mask).sum())
    if mode == 'trans':
        total *= q
    elif mode == 'affine':
        total *= q * (q - 1)
    return total


def prime_tables(q):
    x = np.arange(q)
    return (x[:, None] - x[None, :]) % q, chi_prime(q)


def field_tables(F):
    q = F.q
    sub = np.empty((q, q), dtype=np.int64)
    allx = np.arange(q)
    for a in range(q):
        sub[:, a] = F.sub_vec(allx, a)        # sub[x, a] = x - a
    return sub, F.chi


def main():
    t0 = time.time()
    rows = []
    jobs = [(7, 1, 'full'), (11, 1, 'full'), (19, 1, 'full'), (23, 1, 'full'),
            (31, 1, 'trans'), (43, 1, 'trans'), (47, 1, 'trans'),
            (3, 2, 'trans'), (5, 2, 'trans'), (3, 3, 'trans'), (7, 2, 'affine'),
            (59, 1, 'affine'), (67, 1, 'affine'), (71, 1, 'affine'),
            (5, 1, 'full'), (13, 1, 'full'), (17, 1, 'full')]
    for p, r, mode in jobs:
        q = p ** r
        if r == 1:
            sub, chi = prime_tables(q)
            eps = int(chi[q - 1])
            M = M_exact_direct(q, chi)[0]
        else:
            F = GF(p, r)
            sub, chi = field_tables(F)
            eps = F.eps
            M = M_field(F)[0]
        V = brute(q, sub, chi, mode)
        Vthm = q * (q - 1) * W_from_M(q, eps, M)
        ok = (V == Vthm)
        rec = RECORDED.get(q)
        line = (f"q = {p}^{r} = {q:3d} [{mode:6s}] eps = {eps:+d}  V_q(brute) = {V}  W_q = {V // (q * (q - 1))}  "
                f"M_q = {M}  THM-4653 V = {Vthm}  {'AGREE' if ok else 'DISAGREE'}")
        if rec is not None:
            line += f"  recorded {rec}: {'reproduced' if rec == V else 'MISMATCH'}"
        print(line, f"[{time.time() - t0:.1f}s]", flush=True)
        rows.append((q, ok, rec is None or rec == V))
    # orientation controls at q = 7, 11: reversing all ten edges keeps V; reversing one edge negates it
    for q in (7, 11):
        sub, chi = prime_tables(q)
        rev_all = [(b, a) for (a, b) in EDGES]
        rev_one = [(EDGES[0][1], EDGES[0][0])] + EDGES[1:]
        Vn = brute_general(q, sub, chi, EDGES)
        Va = brute_general(q, sub, chi, rev_all)
        Vo = brute_general(q, sub, chi, rev_one)
        print(f"orientation control q = {q}: natural {Vn}, all edges reversed {Va}, one edge reversed {Vo}")
    print(f"SUMMARY: {sum(ok for _, ok, _ in rows)}/{len(rows)} fields agree with THM-4653; "
          f"recorded values all reproduced: {all(r for _, _, r in rows)}  [{time.time() - t0:.1f}s]")


def brute_general(q, sub, chi, edges):
    """plain itertools enumeration of injective 6-tuples for arbitrary oriented edge list (q <= 11)."""
    tot = 0
    for x in itertools.permutations(range(q), 6):
        v = 1
        for (a, b) in edges:
            v *= int(chi[sub[x[b], x[a]]])
            if v == 0:
                break
        tot += v
    return tot


if __name__ == '__main__':
    main()
