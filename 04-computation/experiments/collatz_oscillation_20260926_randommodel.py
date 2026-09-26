#!/usr/bin/env python3
"""collatz_oscillation_20260926_randommodel.py -- shell-revisit multiplicity in the random-address model of a
divergent orbit (session collatz-oscillation-20260926, opus, 2026-09-26).

Toy: heights in bits h_0 = 0, odd step h -> h + log_2(3/2) followed by v - 1 further halvings (h -> h - 1 each),
with v i.i.d.: P(v = 1) = p, P(v = j) = (1 - p) 2^-(j-1) for j >= 2. E[v] = p + 3(1 - p) = 3 - 2p; the walk diverges
(positive drift per odd step, log_2 3 - E[v] > 0) iff p > 0.7075. This is the address of no integer; it is the
'random seed' against which HYP-9161's shell form is calibrated: for X = 2^L and depth D, dippers, landing
points and multiplicities are computed exactly as for orbits (an index i with h_i <= L dips if some h_(i+s) <
h_i - D, 1 <= s <= L; its landing is the least such i + s), on a trajectory of M T-steps, and the mean and maximum
multiplicity are reported for L in {32, 64, 128, 256}, D = ceil(1.05 log_2 L), several p, and 20 trajectories each.
The random model predicts mean multiplicity of order D/(drift) at most; a growth like L would indicate the
hypothesis fails even for random addresses.
Usage: python3 collatz_oscillation_20260926_randommodel.py [M=200000]
"""
import math, random, sys

LOG32 = math.log2(1.5)


def trajectory(p, M, rng):
    h = [0.0]
    while len(h) < M:
        # odd step
        h.append(h[-1] + LOG32)
        # valuation
        if rng.random() < p:
            v = 1
        else:
            v = 2
            while rng.random() < 0.5:
                v += 1
        for _ in range(v - 1):
            h.append(h[-1] - 1.0)
    return h[:M]


def stats(h, L, D):
    mult = {}; tot = 0; nd = 0
    # only indices with h <= L and a full window
    n = len(h)
    for i in range(n - L):
        if h[i] > L:
            continue
        tot += 1
        thr = h[i] - D
        land = None
        for s in range(1, L + 1):
            if h[i + s] < thr:
                land = i + s; break
        if land is None:
            nd += 1
        else:
            mult[land] = mult.get(land, 0) + 1
    ms = list(mult.values())
    return tot, nd, (sum(ms) / len(ms) if ms else 0.0), (max(ms) if ms else 0), len(ms)


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    rng = random.Random(20260926)
    print("random-address model: p = P(v=1); drift per odd step = log_2 3 - (3 - 2p) bits; 20 trajectories of %d T-steps each" % M)
    print(" p     drift   L    D   mean tot<=L  mean ND  mean multiplicity  max multiplicity (over trajectories)  mean landing pts")
    for p in (0.72, 0.75, 0.80, 0.90):
        drift = math.log2(3) - (3 - 2 * p)
        for L in (32, 64, 128, 256):
            D = math.ceil(1.05 * math.log2(L))
            tots, nds, means, maxs, lps = [], [], [], [], []
            for _ in range(20):
                h = trajectory(p, M, rng)
                tot, nd, mean, mx, lp = stats(h, L, D)
                tots.append(tot); nds.append(nd); means.append(mean); maxs.append(mx); lps.append(lp)
            print("%.2f  %+.3f  %3d  %2d   %9.1f  %8.1f   %6.2f             %3d                          %8.1f" % (
                p, drift, L, D, sum(tots) / 20, sum(nds) / 20, sum(means) / 20, max(maxs), sum(lps) / 20))


if __name__ == '__main__':
    main()
