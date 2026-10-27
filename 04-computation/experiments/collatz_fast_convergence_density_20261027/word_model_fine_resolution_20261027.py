#!/usr/bin/env python3
"""
word_model_fine_resolution_20261027.py  (opus-2026-10-27-S64) -- the pre-registered test of HYP-9425 (4).

Observed: floor counts #{n in B : tau(n) < 3.307 ln n} on the 36 quarter-dyadic cells B of [2^19, 2^28)
(k = 20..28, quarters i = 0..3), from the exact stopping times (tau_all of fast_convergence_density_20261027.py).

Model (independent 3-adic): an integer reaching 1 after k T-steps (one-division map) with d odd steps is
    n = (2^k - c_w) / 3^d,   c_w = sum over odd steps j (0-based) of 3^(d_j') 2^j,  d_j' = #odd steps after j,
and the model counts C(k, d) 3^(-d) P_(k,d)(n in B and k + d < 3.307 ln n), P over uniformly random words of shape
(k, d) whose last step is even (2 -> 1). Monte Carlo with NS words per shape; shapes with k + d <= 3.307 ln(2^28)
and 2^k/3^d in [2^18, 2^29] (others cannot land in the cells).

PRE-REGISTERED (HYP-9425 (4), committed before this run): Pearson correlation of log(model) vs log(observed)
over the 36 cells > 0.8, and sign agreement of consecutive-quarter changes (27 comparisons) >= 70%.
Refutation: correlation < 0.5 or sign agreement at 50%.
"""
import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else (__file__.rsplit("/", 1)[0] or "."))
from fast_convergence_density_20261027 import tau_all, FLOOR_C  # noqa: E402

KMAX = 28
NS = 20000
LN2, LN3 = math.log(2.0), math.log(3.0)
rng = np.random.default_rng(20261027)


def observed_cells(tau):
    cells = {}
    for k in range(20, KMAX + 1):
        base = 1 << (k - 1)
        for i in range(4):
            lo = base + (base * i) // 4
            hi = base + (base * (i + 1)) // 4
            n = np.arange(lo, hi, dtype=np.float64)
            t = tau[lo:hi].astype(np.float64)
            cells[(k, i)] = int(np.count_nonzero(t < FLOOR_C * np.log(n)))
    return cells


def cell_of(n):
    """(k, i) of an integer n in [2^19, 2^28), else None."""
    if n < (1 << 19) or n >= (1 << KMAX):
        return None
    k = n.bit_length()          # n in [2^(k-1), 2^k)
    base = 1 << (k - 1)
    i = ((n - base) * 4) // base
    return (k, int(i))


def model_cells():
    cells = {(k, i): 0.0 for k in range(20, KMAX + 1) for i in range(4)}
    Lmax = FLOOR_C * KMAX * LN2
    shapes = 0
    for k in range(20, 120):
        for d in range(1, k):
            if k + d > Lmax:
                break
            size = k * LN2 - d * LN3
            if size < 18 * LN2 or size > 29 * LN2:
                continue
            shapes += 1
            # the last step is even: choose the d odd positions among 0..k-2
            weight = math.comb(k, d) * (3.0 ** (-d))
            # Monte Carlo
            hits = {}
            for _ in range(NS):
                pos = np.sort(rng.choice(k - 1, size=d, replace=False))
                # c_w = sum_i 3^(d-1-i) 2^(pos_i)  (i-th odd step, 0-based, has d-1-i odd steps after it)
                c = 0
                for idx, j in enumerate(pos):
                    c += (3 ** (d - 1 - idx)) << int(j)
                num = (1 << k) - c
                if num <= 0:
                    continue
                n = num // (3 ** d)        # the model treats integrality as the 3^-d factor; use the real part
                if n <= 1:
                    continue
                if k + d >= FLOOR_C * math.log(n):
                    continue
                cell = cell_of(n)
                if cell is None:
                    continue
                hits[cell] = hits.get(cell, 0) + 1
            # words with an odd last step are excluded by construction; C(k,d) counts all words, so rescale by the
            # fraction of words with an even last step, (k-d)/k
            frac_even_last = (k - d) / k
            for cell, h in hits.items():
                cells[cell] += weight * frac_even_last * h / NS
    return cells, shapes


def main():
    X = 1 << KMAX
    tau, _ = tau_all(X)
    obs = observed_cells(tau)
    mod, shapes = model_cells()
    print("shapes sampled: %d, NS = %d per shape" % (shapes, NS))
    print("  k i   observed      model   ratio")
    xs, ys = [], []
    signs_obs, signs_mod = [], []
    for k in range(20, KMAX + 1):
        for i in range(4):
            o, m = obs[(k, i)], mod[(k, i)]
            print("%3d %d  %9d  %11.1f  %6.3f" % (k, i, o, m, (m / o) if o else float("nan")))
            xs.append(math.log(max(m, 1e-9)))
            ys.append(math.log(max(o, 1)))
            if i > 0:
                signs_obs.append(np.sign(obs[(k, i)] - obs[(k, i - 1)]))
                signs_mod.append(np.sign(mod[(k, i)] - mod[(k, i - 1)]))
    corr = float(np.corrcoef(xs, ys)[0, 1])
    agree = sum(1 for a, b in zip(signs_obs, signs_mod) if a == b) / len(signs_obs)
    total_o = sum(obs.values())
    total_m = sum(mod.values())
    print("total observed %d, total model %.1f (ratio %.3f)" % (total_o, total_m, total_m / total_o))
    print("Pearson correlation of log(model) vs log(observed) over %d cells: %.3f" % (len(xs), corr))
    print("sign agreement of consecutive-quarter changes: %.3f (%d comparisons)" % (agree, len(signs_obs)))
    ok = corr > 0.8 and agree >= 0.70
    refuted = corr < 0.5 or agree <= 0.5
    print(("PASS" if ok else ("REFUTED" if refuted else "INCONCLUSIVE")) +
          " HYP-9425 (4): correlation > 0.8 and sign agreement >= 0.70")


if __name__ == "__main__":
    main()
