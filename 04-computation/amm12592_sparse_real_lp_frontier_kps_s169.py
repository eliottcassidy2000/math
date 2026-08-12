#!/usr/bin/env python3
"""AMM 12592: sparse causal real-LP probe at the exact golden floor.

The full transportation relaxation is conjugated to junk-flow coordinates.
Pascal transport has width 2 or 3, so the matrix has O(R^2) variables and
nonzeros rather than a dense Bernstein-to-power representation.  This is a
NUMERICAL feasibility probe.  It does not certify a rational or integer point.

Reproduction:
  python 04-computation/amm12592_sparse_real_lp_frontier_kps_s169.py 32
  python 04-computation/amm12592_sparse_real_lp_frontier_kps_s169.py 128
  python 04-computation/amm12592_sparse_real_lp_frontier_kps_s169.py 256
  python 04-computation/amm12592_sparse_real_lp_frontier_kps_s169.py 512
"""

from __future__ import annotations

import argparse
import math
import time
from math import comb

import cvxpy as cp
import numpy as np
from scipy.sparse import csr_matrix


def fib_pair(n):
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def five_pow_le_phi2m(d, m):
    if d < 0:
        return True
    F, F1 = fib_pair(2 * m)
    L = 2 * F1 - F
    A = 2 * 5**d - L
    return A <= 0 or A * A < 5 * F * F


def floor_gamma_star(m):
    d = int(m * 0.5979874356654402)
    while five_pow_le_phi2m(d + 1, m):
        d += 1
    while not five_pow_le_phi2m(d, m):
        d -= 1
    return d


def two_G_coeffs(R):
    g = [2]
    b = 1
    for j in range(1, R):
        b = b * (R - j) // j
        g.append((-b if j & 1 else b) - 1)
    return g


def t4_row_load(R, d):
    return [((-1) ** (d - t)) * comb(R - 2 - t, d - t)
            - comb(d + 1, t + 1) + 2 * comb(d, t)
            for t in range(d + 1)]


def ratio(a, b):
    if a == 0:
        return 0.0
    return (-1.0 if a < 0 else 1.0) * math.exp(math.log(abs(a)) - math.log(b))


def build(R, D0):
    profile = [floor_gamma_star(R + i) + D0 for i in range(R)]
    ds = profile[:-1]
    offsets = [0]
    for d in ds:
        offsets.append(offsets[-1] + d)
    nvar = offsets[-1]

    def vid(i, t):
        return offsets[i] + t

    bounds = [(None, None)] * nvar
    scales = []
    d0 = ds[0]
    w0 = t4_row_load(R, d0)
    s0 = []
    for t in range(d0):
        C = comb(d0, t)
        S = max(1, 2 * C, abs(w0[t]))
        s0.append(S)
        lo = -2 * (comb(d0 - 1, t) if t <= d0 - 1 else 0)
        hi = 2 * (comb(d0 - 1, t - 1) if t >= 1 else 0)
        bounds[vid(0, t)] = (ratio(w0[t] - hi, S), ratio(w0[t] - lo, S))
    assert 0 <= w0[d0] <= 2
    scales.append(s0)

    g = two_G_coeffs(R)
    indptr, indices, data, rhs = [0], [], [], []

    def add(entries, b):
        for c, v in entries:
            if v:
                indices.append(c)
                data.append(v)
        indptr.append(len(indices))
        rhs.append(b)

    for i in range(1, len(ds)):
        d, dp = ds[i], ds[i - 1]
        delta = d - dp
        K = (1, 1) if delta == 0 else (1, 2, 1)
        feed0 = feed1 = 0
        if d + i <= R - 1:
            feed0 += g[d + i]
        if delta == 1 and d - 1 + i <= R - 1:
            z = g[d - 1 + i]
            feed0 += z
            feed1 += z
        si = []
        for t in range(d):
            f = feed0 if t == 0 else (feed1 if t == 1 else 0)
            incoming = abs(f) + sum(
                K[a] * scales[i - 1][t - a]
                for a in range(len(K)) if 0 <= t - a < dp)
            si.append(max(1, 2 * comb(d, t), incoming))
        scales.append(si)
        for t in range(d + 1):
            f = feed0 if t == 0 else (feed1 if t == 1 else 0)
            incoming = abs(f) + sum(
                K[a] * scales[i - 1][t - a]
                for a in range(len(K)) if 0 <= t - a < dp)
            S = si[t] if t < d else max(1, 2 * comb(d, t), incoming)
            prev = [(vid(i - 1, t - a), K[a] * ratio(scales[i - 1][t - a], S))
                    for a in range(len(K)) if 0 <= t - a < dp]
            lo = -2 * (comb(d - 1, t) if t <= d - 1 else 0)
            hi = 2 * (comb(d - 1, t - 1) if t >= 1 else 0)
            fn = ratio(f, S)
            neg = [(vid(i, t), -1.0)] if t < d else []
            pos = [(vid(i, t), 1.0)] if t < d else []
            add(prev + neg, ratio(hi, S) - fn)
            add([(c, -v) for c, v in prev] + pos, fn - ratio(lo, S))

    # Terminal row: bern(j)/x is the same width-2/3 Pascal transport.
    i = len(ds) - 1
    dp, dl = ds[i], profile[-1]
    K = (1, 1) if dl - dp == 0 else (1, 2, 1)
    for t in range(dl + 1):
        N = max(1, 2 * comb(dl, t), sum(
            K[a] * scales[i][t - a]
            for a in range(len(K)) if 0 <= t - a < dp))
        prev = [(vid(i, t - a), K[a] * ratio(scales[i][t - a], N))
                for a in range(len(K)) if 0 <= t - a < dp]
        add(prev, ratio(2 * comb(dl, t), N))
        add([(c, -v) for c, v in prev], 0.0)

    A = csr_matrix((np.asarray(data), np.asarray(indices, dtype=np.int32),
                    np.asarray(indptr, dtype=np.int64)),
                   shape=(len(rhs), nvar))
    return profile, A, np.asarray(rhs), bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("R", type=int)
    ap.add_argument("D0", type=int, nargs="?", default=0)
    ap.add_argument("--time-limit", type=float, default=900)
    args = ap.parse_args()
    t0 = time.time()
    profile, A, b, bounds = build(args.R, args.D0)
    print(f"R={args.R} D0={args.D0} vars={A.shape[1]} inequalities={A.shape[0]} "
          f"nnz={A.nnz} d0={profile[0]} dlast={profile[-1]} "
          f"build_s={time.time()-t0:.3f}", flush=True)
    print(f"A_range=[{A.data.min():.6e},{A.data.max():.6e}] "
          f"b_range=[{b.min():.6e},{b.max():.6e}]", flush=True)
    x = cp.Variable(A.shape[1])
    lo = np.asarray([z[0] if z[0] is not None else -np.inf for z in bounds])
    hi = np.asarray([z[1] if z[1] is not None else np.inf for z in bounds])
    problem = cp.Problem(cp.Minimize(0), [A @ x <= b, x >= lo, x <= hi])
    value = problem.solve(solver="CLARABEL", max_iter=1000, tol_feas=1e-9,
                          tol_gap_abs=1e-9, time_limit=args.time_limit)
    success = problem.status in ("optimal", "optimal_inaccurate")
    print(f"status={problem.status} success={success} value={value} "
          f"elapsed_s={time.time()-t0:.3f}", flush=True)
    if x.value is not None:
        viol = max(np.max(A @ x.value - b), np.max(lo - x.value),
                   np.max(x.value - hi))
        print(f"max_checked_violation={viol:.16e} "
              f"max_abs_scaled_variable={np.max(np.abs(x.value)):.16e}")
    print("scope=NUMERICAL_REAL_LP_FEASIBILITY_NOT_RATIONAL_CERTIFICATE_"
          "NOT_INTEGER_WITNESS_NOT_AMM_CLOSURE")


if __name__ == "__main__":
    main()
