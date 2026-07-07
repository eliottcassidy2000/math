#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S1 (part 3) -- shift-sum mod-14 + Chung-Erdos as a candidate
PROOF MECHANISM for the k=13 tail floor mu_{1/7}(E) >= c > m_P.

IDEA (HYP-4787 backlog target 4).  Discretize with the 1/14-grid: cells C_j=[j/14,(j+1)/14).
A_j = {x : all 13 points {e_i x} avoid C_j u C_{j+1}} = an ALIGNED empty 1/7-arc.
  mu_{1/7}(E) >= P(U_j A_j)   (a.e.; an aligned empty 2-cell arc forces maxgap >= 1/7)
Chung-Erdos:  P(U A_j) >= S1^2 / (S1 + 2*sum_{j<j'} P(A_j /\ A_{j'})),  S1 = sum_j p_j.

WHY THIS MIGHT BE PROVABLE where HYP-4767's unsigned lattice sums diverged: the SHIFT SUM
S1 = sum_j p_j has Fourier support only on total frequencies == 0 mod 14
(sum_j e(m j/14) = 14*[14|m]), pruning the resonance lattice.  Independence value:
S1_indep = 14*(6/7)^13 ~ 1.871.  If S1 is uniformly ~1.8 and pair sums are tame, CE gives a
uniform mu floor >> m_P = 0.0565 -- turning the k=13 leg into two boundable correlation sums.

This script: numeric test across the strongest known adversaries + an S1-minimizing descent.
Decides whether the analytic program (bound S1 below, pair sum above, via the pruned lattice)
is worth fleet effort.

Tournament Analysis declaration:
  vertices: candidate floor mechanisms (union bound, reverse-Markov, condRM, CE-shiftsum);
  pairwise observable: certified floor value on the adversarial family bank;
  switch/gauge: orient toward the mechanism with the larger worst-case certified floor;
  tie Hamiltonian path: exact ledger -> CE bank -> S1-descent -> verdict.
"""
from fractions import Fraction as F
from math import gcd
import random

import numpy as np

M_P = 14249 / 252252
IND13 = (6 / 7) ** 13


def ce_stats(E, res=160000):
    """Numeric: p_j, S1, pair sum, Chung-Erdos bound, P(union), true mu_{1/7}."""
    xs = (np.arange(res) + 0.5) / res
    ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0)   # res x 13
    # presence[r, j] = some point in C_j u C_{j+1} (arc [j/14, j/14+1/7))
    A = np.empty((14, res), dtype=bool)
    for j in range(14):
        rel = np.mod(ph - j / 14.0, 1.0)
        A[j] = ~(rel < 1.0 / 7.0).any(axis=1)      # A_j holds: arc empty
    p = A.mean(axis=1)
    S1 = p.sum()
    pair = 0.0
    for j in range(14):
        for jp in range(j + 1, 14):
            pair += (A[j] & A[jp]).mean()
    union = A.any(axis=0).mean()
    ce = S1 * S1 / (S1 + 2 * pair) if S1 > 0 else 0.0
    # true mu (any-position gap > 1/7)
    ph.sort(axis=1)
    gaps = np.empty_like(ph)
    gaps[:, :-1] = np.diff(ph, axis=1)
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mu = (gaps.max(axis=1) > 1.0 / 7.0).mean()
    return p, S1, pair, ce, union, mu


def saturated(v): return all(any(x % q == 0 for x in v) for q in range(2, 15))
def primitive(v):
    g = 0
    for x in v: g = gcd(g, x)
    return g == 1
def single_scale(v): return max(v) <= 13 * min(v)
def in_crux(v): return saturated(v) and primitive(v) and single_scale(v)


BANK = {
    "AP {1..13}": list(range(1, 14)),
    "2*{1..12} u {13} (ds record)": [2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],
    "2*{1..11} u {11,13} (my record)": [2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
    "S57 adversarial": [2, 6, 8, 10, 11, 12, 14, 16, 18, 20, 22, 26, 42],
    "consecutive 14..26 (sat block)": list(range(14, 27)),
    "GW {1..11,13,24}": list(range(1, 12)) + [13, 24],
    "random big (min>=60)": [61, 67, 74, 83, 89, 97, 104, 113, 122, 131, 140, 151, 163],
    "random huge (min>=500)": [503, 541, 577, 613, 659, 701, 743, 787, 823, 863, 907, 947, 983],
    "mult-14 heavy": [14, 28, 42, 56, 70, 84, 98, 112, 126, 140, 11, 13, 9],
}

if __name__ == "__main__":
    print("shift-sum mod-14 + Chung-Erdos test;  m_P=%.5f, independence S1=14*(6/7)^13=%.4f"
          % (M_P, 14 * IND13))
    print()
    print(f"{'family':>34} {'S1':>7} {'pairSum':>8} {'CE>=':>7} {'P(univ)':>8} "
          f"{'mu_1/7':>7} {'CE/mu':>6}")
    for nm, E in BANK.items():
        p, S1, pair, ce, un, mu = ce_stats(E)
        print(f"{nm:>34} {S1:>7.4f} {pair:>8.4f} {ce:>7.4f} {un:>8.4f} {mu:>7.4f} "
              f"{ce / mu if mu > 0 else 0:>6.2f}"
              + ("   *** CE < m_P ***" if ce < M_P else ""))

    print()
    print("S1-minimizing adversarial descent (crux class; the mechanism's weak point is S1):")
    random.seed(1414)
    best_S1, best_v, best_ce = 99.0, None, None
    for trial in range(6):
        v = random.choice([list(BANK["2*{1..12} u {13} (ds record)"]),
                           list(BANK["consecutive 14..26 (sat block)"]),
                           list(BANK["AP {1..13}"])])
        _, S1, pair, ce, un, mu = ce_stats(v, 40000)
        for step in range(250):
            i = random.randrange(13)
            cand = random.randrange(max(1, min(v) // 2), min(3 * max(v), 600))
            w = sorted(set(v[:i] + v[i + 1:] + [cand]))
            if len(w) != 13 or not in_crux(w):
                continue
            _, S1b, pairb, ceb, unb, mub = ce_stats(w, 40000)
            if S1b < S1 - 1e-4:
                v, S1, pair, ce, un, mu = w, S1b, pairb, ceb, unb, mub
        if S1 < best_S1:
            best_S1, best_v = S1, v
    p, S1, pair, ce, un, mu = ce_stats(best_v)
    print(f"  min-S1 family: {best_v}")
    print(f"  S1={S1:.4f}  pair={pair:.4f}  CE>={ce:.4f}  P(u)={un:.4f}  mu={mu:.4f}"
          + ("   *** CE < m_P ***" if ce < M_P else "   (CE >= m_P holds)"))
    print(f"  per-j p_j: {[round(float(x), 4) for x in p]}")
    print()
    print("VERDICT inputs: if worst CE across bank+descent stays >= m_P=%.4f with margin,"
          % M_P)
    print("the k=13 leg reduces to TWO boundable correlation sums (S1 below, pairs above)")
    print("on the mod-14-pruned resonance lattice -- an analytic program, not a census.")
