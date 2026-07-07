#!/usr/bin/env python3
r"""
lrc_tent_fLP_k9_kps_S73.py   (kind-pasteur-2026-07-07-S73, HYP-5147 follow-through)

THE FULL f-LP: optimize the shifted-tent theorem's weight function.

For a nonneg f supported on (0, 1/7], the theorem machinery gives, for every k-family,
    mu_{1/7}(E) >= 1 - k(k-1) * int f / m(f),
    m(f) = min over gap vectors {g in [0,1/7]^k, sum g = 1} of RINGSUM_f(g),
    RINGSUM_f(g) = sum_{l=1}^{k-1} sum_{i=1}^{k} f(S_{i,l}),  S_{i,l} = g_i + ... + g_{i+l-1} (cyclic)
(only terms with S <= 1/7 pay, since supp f = (0,1/7]; every term is an ordered-pair
difference, and pair equidistribution gives E[sum over ALL k(k-1) ordered pairs] =
k(k-1) int f -- ring terms with S > 1/7 contribute >= 0 there, so the bound is valid).

The tent f = (s - beta)+ realizes: k=8: 3/4 (DISCHARGES the honest bar 0.6750);
k=9: 31/63 = 0.4921 vs bar 0.5622; k=10: 8/35 = 0.2286 vs bar 0.4521.
QUESTION: does the OPTIMAL f close k=9 or k=10?

METHOD: column-generation LP.  Variables: f on a grid of (0, 1/7] (N bins, piecewise
constant).  Constraints: RINGSUM_f(g) >= 1 for every g in a working set of gap configs.
Objective: minimize int f.  After each solve, search for violated configs (adversary:
minimize RINGSUM over the polytope by projected local descent + structured seeds);
add violations, repeat until clean.  Floor = 1 - k(k-1) * int f.

NOTE (honest): the certified floor is modulo the adversary search being exhaustive --
the tent value is a PROVED baseline; any improvement here is 'verified-LP' grade until
the binding configs are handled analytically.
"""
import random
import numpy as np
from scipy.optimize import linprog

THR = 1.0 / 7.0

def make_ringsum(k, N):
    """returns function: given f-vector (grid on (0,THR]) and gap vector, RINGSUM."""
    def bin_of(s):
        if s <= 0 or s > THR:
            return -1
        b = int(s / THR * N)
        return min(b, N - 1)
    def ringsum_vec(g):
        """coefficient vector: counts per bin over all ring sums <= THR."""
        v = np.zeros(N)
        kk = len(g)
        for i in range(kk):
            s = 0.0
            for l in range(1, kk):
                s += g[(i + l - 1) % kk]
                if s > THR:
                    break
                b = bin_of(s)
                if b >= 0:
                    v[b] += 1
        return v
    return ringsum_vec, bin_of

def adversary_min(fvals, k, N, rng, tries=400, steps=600):
    """minimize RINGSUM_f over the gap polytope; returns (value, g)."""
    ringsum_vec, bin_of = make_ringsum(k, N)
    def val(g):
        return float(np.dot(ringsum_vec(g), fvals))
    best_v, best_g = None, None
    seeds = []
    # structured two/three-value seeds
    for t in range(k):
        for lo_frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            lo = lo_frac * THR
            s_rest = 1 - t * lo
            if s_rest < 0 or k - t == 0:
                continue
            hi = s_rest / (k - t)
            if hi <= THR + 1e-12 and hi >= 0:
                seeds.append([lo] * t + [hi] * (k - t))
    for _ in range(tries):
        g = [rng.uniform(0, THR) for _ in range(k - 1)]
        s = sum(g)
        if s < 1 and 1 - s <= THR:
            seeds.append(g + [1 - s])
    for g0 in seeds:
        g = list(g0)
        if abs(sum(g) - 1) > 1e-9 or max(g) > THR + 1e-9:
            continue
        v = val(g)
        # projected pairwise-transfer descent
        for _ in range(steps):
            i, j = rng.randrange(k), rng.randrange(k)
            if i == j:
                continue
            eps = rng.uniform(-0.02, 0.02) * THR
            gi, gj = g[i] + eps, g[j] - eps
            if 0 <= gi <= THR and 0 <= gj <= THR:
                g2 = list(g)
                g2[i], g2[j] = gi, gj
                v2 = val(g2)
                if v2 < v - 1e-12:
                    g, v = g2, v2
        if best_v is None or v < best_v:
            best_v, best_g = v, g
    return best_v, best_g

def f_lp(k, N=280, rounds=25, seed=73):
    rng = random.Random(seed)
    ringsum_vec, _ = make_ringsum(k, N)
    ds = THR / N
    # start with structured configs
    configs = []
    for t in range(k):
        s_rest = 1 - t * 0.0
        hi = s_rest / (k - t) if k > t else None
        if hi is not None and hi <= THR:
            configs.append([0.0] * t + [hi] * (k - t))
    configs.append([1.0 / k] * k)
    val_hist = []
    for rd in range(rounds):
        A_ub, b_ub = [], []
        for g in configs:
            A_ub.append(-ringsum_vec(g))
            b_ub.append(-1.0)
        c = np.full(N, ds)
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      bounds=[(0, None)] * N, method="highs")
        if res.status != 0:
            return None, None, None
        fvals = res.x
        intf = float(res.fun)
        floor = 1 - k * (k - 1) * intf
        # adversary
        v, g_bad = adversary_min(fvals, k, N, rng)
        val_hist.append((rd, floor, v))
        if v >= 1 - 1e-6:
            return floor, fvals, val_hist
        configs.append(g_bad)
    return floor, fvals, val_hist

print("=" * 92)
print("THE f-LP (rings included): optimal weight function for the gap-histogram floor")
print("=" * 92)
HONEST = {8: 0.675024, 9: 0.562233, 10: 0.452092, 11: 0.331246, 12: 0.199344, 13: 0.056487}
TENT = {8: 3/4, 9: 31/63, 10: 8/35, 11: 0.0, 12: 0.0, 13: 0.0}
for k in (8, 9, 10, 11):
    floor, fvals, hist = f_lp(k)
    if floor is None:
        print(f"  k={k}: LP failed")
        continue
    conv = "converged" if hist[-1][2] >= 1 - 1e-6 else f"NOT converged (adv found {hist[-1][2]:.4f})"
    print(f"  k={k}: f-LP floor = {floor:.5f}  (tent {TENT[k]:.4f}, honest bar {HONEST[k]:.4f}) "
          f"{'DISCHARGES' if floor > HONEST[k] else 'short'}  [{conv}, {len(hist)} rounds]")
    # describe optimal f coarsely
    N = len(fvals)
    nz = [i for i in range(N) if fvals[i] > 1e-8]
    if nz:
        lo, hi = nz[0] / N * THR, (nz[-1] + 1) / N * THR
        print(f"        supp(f*) ~ [{lo:.4f}, {hi:.4f}] = [{lo*28:.2f}, {hi*28:.2f}]/28; "
              f"max f* = {max(fvals):.4f}")
