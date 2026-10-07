#!/usr/bin/env python3
"""HYP-9218 conditional test (mac-mini-2026-10-07-oaimath3): are successive FRESH alignment depths of one pair independent
(conditional Haar-ness), away from merges?  Uses run_pair from twoorbit_exponents_20261007.py.
Reports, for consecutive fresh alignments (both with L_s >= 8) within a pair: the joint law of (delta_j, delta_(j+1)) vs the
product of Haar marginals (chi-square), P(delta_(j+1) >= d | delta_j >= 3) / 2^(1-d), and the correlation of 2^(1-delta).
Usage: python3 twoorbit_conditional_haar_20261007.py [NPAIRS BITS SEED NPROC]  (defaults 2000 20000 2031 10)"""
import sys, math
import numpy as np
from multiprocessing import Pool
import twoorbit_exponents_20261007 as T

def fresh_seq(idx):
    r = T.run_pair(idx)
    return [(x[3], x[6]) for x in r['aligns'] if x[9] == 0]       # (delta, L_s) for fresh alignments, in time order

if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:5]]
    NP, BITS, SEED, NPROC = (args + [2000, 20000, 2031, 10][len(args):])[:4]
    T.BITS, T.SEED = BITS, SEED
    with Pool(NPROC) as pool:
        seqs = pool.map(fresh_seq, range(NP))
    pairs = []
    for sq in seqs:
        for (d1, l1), (d2, l2) in zip(sq, sq[1:]):
            if l1 >= 8 and l2 >= 8:
                pairs.append((min(d1, 60), min(d2, 60)))
    P = np.array(pairs)
    n = len(P)
    cap = lambda d: np.minimum(d, 6)
    obs = np.zeros((6, 6))
    for i, j in zip(cap(P[:, 0]), cap(P[:, 1])):
        obs[i - 1, j - 1] += 1
    pm = np.array([2.0 ** -d for d in range(1, 6)] + [2.0 ** -5])
    exp_ = n * np.outer(pm, pm)
    chi2 = float(np.sum((obs - exp_) ** 2 / exp_))
    print(f"consecutive fresh alignments with L >= 8: n = {n}; joint (delta_j, delta_(j+1)) vs Haar x Haar: chi^2 = {chi2:.1f} on 35 dof")
    m = P[:, 0] >= 3
    print("   P(delta' >= d | delta >= 3) / 2^(1-d), d = 1..6: " + " ".join(f"{float((P[m, 1] >= d).mean()) / 2 ** (1 - d):.3f}" for d in range(1, 7))
          + f"   (n = {int(m.sum())})")
    m = P[:, 0] == 1
    print("   P(delta' >= d | delta = 1) / 2^(1-d), d = 1..6: " + " ".join(f"{float((P[m, 1] >= d).mean()) / 2 ** (1 - d):.3f}" for d in range(1, 7))
          + f"   (n = {int(m.sum())})")
    u = 2.0 ** (1 - P[:, 0]); w = 2.0 ** (1 - P[:, 1])
    print(f"   corr(2^(1-delta_j), 2^(1-delta_(j+1))) = {float(np.corrcoef(u, w)[0, 1]):+.4f}  (s.e. ~ {1 / math.sqrt(n):.4f})")
