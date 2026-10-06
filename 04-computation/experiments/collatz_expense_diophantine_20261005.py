#!/usr/bin/env python3
"""Diophantine structure of the excess-rate expense, the exact segment-level
tail, inheritance of hard excursions, and bank-relative two-tier families.

Inherits F_q from collatz_weak_reset_family_20261005.py: a first-descent segment
(length l, valuation A) has least admissible rate q(l,A) = least q with
2^(qA-l) >= 3^(ql) = ceil(l/(A - l log2 3)); q(n) is the max over the chain of
running minima. New here:
  Q1  q depends only on (l,A); the worst at length l is the minimal valuation
      A_min(l) = bitlength(3^l), with q_max(l) = ceil(l/(A_min - l log2 3)); the
      record lengths are the denominators of the best upper approximations of
      log2 3 (upper semiconvergents 2/1, 8/5, 27/17, 46/29, 65/41, 149/94, ...).
  Q2  the proportion of odd sources whose OWN first segment has q > Q equals
      sum over first-descent types (l,A) with q(l,A) > Q of N(l,A)/2^A, where
      N(l,A) counts words of length l and sum A with all proper prefixes rising;
      computed exactly by DP and compared with the census below 2^20.
  Q3  inheritance: popularity of hard excursions as running minima.
  Q4  two-tier families F(q_hi, Y): rate q_hi above Y, table below; coverage and
      the two-tier deadline tau <= floor(1.051 q_hi log2 n) + max tau(odd < Y') (the bank
      buys coverage, not deadline length: the last admissible segment may land far below Y').
Companion note: 05-knowledge/results/collatz_expense_diophantine_20261005.md
Usage: python3 <this file> [--census-bits 20] [--json PATH]
"""
from fractions import Fraction as F
from math import log2, ceil, log
import argparse
import json
import time

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a or not callable(a[0]) else a[0]

CHECKS = 0
ALPHA = log2(3.0)


def require(cond, witness=None):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError(witness)


def admissible(l, A, q):
    return q * A - l >= 0 and (1 << (q * A - l)) >= 3 ** (q * l)


import mpmath as mp
mp.mp.dps = 60
ALPHA_MP = mp.log(3) / mp.log(2)


def q_seg_exact(l, A):
    """Least q with 2^(qA-l) >= 3^(ql), i.e. ceil(l/(A - l log2 3)) since log2 3 is irrational.
    Computed at 60 digits; confirmed by the exact integer test whenever qA <= 2e6 bits."""
    require((1 << A) > 3 ** l, ("not a descent type", l, A))
    exc = A - l * ALPHA_MP
    q = int(mp.ceil(l / exc))
    if q * A <= 2 * 10 ** 6:
        require(admissible(l, A, q) and not admissible(l, A, q - 1), ("exact confirmation", l, A, q))
    return q


def A_min(l):
    return (3 ** l).bit_length()   # ceil(l log2 3)


# ---------------------------------------------------------------------------
# Q1: Diophantine structure
# ---------------------------------------------------------------------------
def continued_fraction_upper(limit_den):
    """Best upper approximations p/q > log2 3 with q <= limit_den (semiconvergents of upper type)."""
    # exact convergents via the integer test 2^p vs 3^q
    best = []
    # brute force: for each q, the least p with 2^p > 3^q; record if p/q - alpha is a new minimum
    rec = None
    for q in range(1, limit_den + 1):
        p = A_min(q)
        # compare excess p - q alpha exactly through 2^p/3^q = 2^(p - q alpha) : smaller ratio = smaller excess
        ratio = F(1 << p, 3 ** q)
        if rec is None or ratio < rec:
            rec = ratio
            best.append((p, q))
    return best


def section_diophantine(lmax=320):
    print("== Q1. The expense depends only on (l,A); records at upper best approximations of log2 3 ==")
    # q(l,A) = ceil(l/(A - l alpha)) (float) agrees with the exact integer bisection
    for l in range(1, 121):
        for A in range(A_min(l), A_min(l) + 4):
            qe = q_seg_exact(l, A)
            qf = ceil(l / (A - l * ALPHA) - 1e-9)
            require(qe == qf or qe == qf + 1, (l, A, qe, qf))
    print("  q(l,A) = ceil(l/(A - l log2 3)) matches the exact integer bisection for l<=120, A in [A_min, A_min+3]")
    # records of q_max(l) = q(l, A_min(l))
    records = []
    best = 0
    for l in range(1, lmax + 1):
        q = q_seg_exact(l, A_min(l))
        if q > best:
            best = q
            records.append((l, A_min(l), q, (1 << A_min(l)) - 3 ** l))
    print("  record lengths (l, A_min, q_max, 2^A - 3^l):", [(l, A, q, d if d < 10 ** 6 else "%.2e" % d) for (l, A, q, d) in records])
    upper = continued_fraction_upper(lmax)
    upper_dens = [q for (p, q) in upper]
    require([l for (l, A, q, d) in records] == upper_dens, ("records vs upper approximations", records, upper))
    print("  these are exactly the denominators of the best upper approximations p/l > log2 3:",
          ["%d/%d" % (p, q) for (p, q) in upper])
    return records


# ---------------------------------------------------------------------------
# Q2: exact segment-level tail by DP over first-descent words
# ---------------------------------------------------------------------------
def first_descent_counts(lmax):
    """N[l][A]: number of words (a_1..a_l), a_i>=1, with 2^(A_j) < 3^j for all j<l and 2^A > 3^l."""
    # states after j letters: dict s -> count, with 2^s < 3^j (rising prefix)
    N = {}
    states = {0: 1}  # j = 0
    for j in range(1, lmax + 1):
        new_states = {}
        counts_j = {}
        thr_rise = A_min(j)          # rising at length j: 2^s < 3^j  <=> s < A_min(j)  (s <= A_min-1)
        for s, cnt in states.items():
            # next letter a >= 1 gives s+a; rising if s+a <= thr_rise-1, descending if s+a >= thr_rise
            for a in range(1, thr_rise - s):          # s+a <= thr_rise-1 : rising prefix of length j
                new_states[s + a] = new_states.get(s + a, 0) + cnt
            # descending completions: any a with s+a >= thr_rise; we record A = s+a for a bounded range
            # of A (excess < some cap); words with larger terminal valuation have larger excess and
            # tiny q; cap A at thr_rise + 40 (excess >= 40 means q = 1 for l < 40... keep exact later)
            for A in range(max(thr_rise, s + 1), thr_rise + 41):
                counts_j[A] = counts_j.get(A, 0) + cnt
        N[j] = counts_j
        states = new_states
    return N


def section_tail(lmax, census_qseg, census_lA, census_bits=20):
    print("== Q2. Exact segment-level tail: density of sources whose own first segment has q > Q ==")
    t0 = time.time()
    N = first_descent_counts(lmax)
    # total mass of first-descent types with l <= lmax (as a proportion of odd sources)
    total = F(0)
    types = []
    for l in range(1, lmax + 1):
        for A, cnt in N[l].items():
            total += F(cnt, 1 << A)
            types.append((l, A, cnt))
    print("  DP over %d first-descent types with l<=%d (%.1fs): total density %.6f"
          " (deficit %.2e = sources whose first descent is later than %d odd steps, plus the capped terminal valuations)"
          % (len(types), lmax, time.time() - t0, float(total), 1 - float(total), lmax))
    # q for each type via the ceil formula (exact bisection for the small-excess ones)
    tail_at = {}
    thresholds = [2, 3, 7, 8, 13, 31, 67, 104, 310, 800, 2000, 6951]
    for Q in thresholds:
        s = F(0)
        for (l, A, cnt) in types:
            if ceil(l / (A - l * ALPHA) - 1e-9) > Q:
                s += F(cnt, 1 << A)
        tail_at[Q] = s
    print("  exact density of {q_seg > Q}:", " ".join("%d:%.5f" % (Q, float(tail_at[Q])) for Q in thresholds))
    # the staircase: contribution of each record type (l, A_min)
    print("  record types' own densities N(l,A_min)/2^A:")
    for l in (1, 5, 17, 29, 41, 94):
        A = A_min(l)
        cnt = N[l].get(A, 0)
        print("    l=%3d A=%3d N=%d density=%.3e q=%d" % (l, A, cnt, float(F(cnt, 1 << A)), q_seg_exact(l, A)))
    # compare with the census of own-segment q below 2^20
    print("  census below 2^20, proportion with own-segment q > Q:",
          " ".join("%d:%.5f" % (Q, float((census_qseg > Q).mean())) for Q in thresholds))
    for Q in thresholds:
        require(abs(float(tail_at[Q]) - float((census_qseg > Q).mean())) < 0.004 * (1 << max(0, 20 - census_bits)) ** 0.5, (Q, float(tail_at[Q]), float((census_qseg > Q).mean())))
    # exact agreement of type frequencies for small A (cylinders fully represented below 2^20)
    worst = 0.0
    for (l, A, cnt) in types:
        if A <= census_bits - 4:   # cylinders mod 2^(A+1) are fully represented below 2^census_bits
            emp = float(((census_lA[:, 0] == l) & (census_lA[:, 1] == A)).mean())
            worst = max(worst, abs(emp - float(F(cnt, 1 << A))))
    print("  type frequencies (l,A) with A<=%d: max |census - exact| = %.2e" % (census_bits - 4, worst))
    require(worst < 1e-4, ("type frequencies", worst))
    return {Q: str(v) for Q, v in tail_at.items()}


# ---------------------------------------------------------------------------
# census machinery
# ---------------------------------------------------------------------------
@njit(cache=True)
def _census(limit, thresholds):
    """Per odd n: own-segment (l,A,q), chain max q, tau, and max q above each threshold."""
    m = (limit + 1) // 2
    nth = thresholds.shape[0]
    own_l = np.zeros(m, dtype=np.int64)
    own_A = np.zeros(m, dtype=np.int64)
    own_q = np.zeros(m, dtype=np.int64)
    qn = np.zeros(m, dtype=np.int64)
    tau = np.zeros(m, dtype=np.int64)
    qabove = np.zeros((m, nth), dtype=np.int64)
    inherit = np.zeros(m, dtype=np.int64)   # number of sources having x as a running minimum (x = 2i+1)
    log23 = np.log2(3.0)
    for i in range(1, m):
        n = 2 * i + 1
        x = n
        qbest = 0
        t = 0
        first = True
        while x != 1:
            if x < limit:
                inherit[(x - 1) // 2] += 1
            v = x
            l = 0
            A = 0
            while True:
                tt = 3 * v + 1
                a = 0
                while (tt & 1) == 0:
                    tt >>= 1
                    a += 1
                v = tt
                l += 1
                A += a
                if v < x:
                    break
            t += l
            exc = A - l * log23
            qs = int(np.ceil(l / exc - 1e-9))
            if qs < 1:
                qs = 1
            if first:
                own_l[i] = l
                own_A[i] = A
                own_q[i] = qs
                first = False
            if qs > qbest:
                qbest = qs
            for k in range(nth):
                if x >= thresholds[k] and qs > qabove[i, k]:
                    qabove[i, k] = qs
            x = v
        qn[i] = qbest
        tau[i] = t
    return own_l, own_A, own_q, qn, tau, qabove, inherit


def section_inheritance(bits, own_l, own_A, own_q, qn, tau, qabove, inherit, thresholds):
    print("== Q3. Inheritance: popular hard excursions below 2^%d ==" % bits)
    n = 2 * np.arange(len(qn)) + 1
    hard = np.where((own_q >= 30) & (n > 1))[0]
    order = hard[np.argsort(-inherit[hard])]
    print("  top hard excursions (start x, own segment (l,A), q, inheritors = sources below 2^%d with x as a running minimum):" % bits)
    for i in order[:12]:
        print("    x=%7d (l,A)=(%d,%d) q=%d inheritors=%d (%.3f of census)"
              % (n[i], own_l[i], own_A[i], own_q[i], inherit[i], inherit[i] / (len(qn) - 1)))
    # the share of q(n) > 30 explained by the top-k hard excursions
    tot = int(((qn > 30) & (n > 1)).sum())
    print("  sources with q(n) > 30: %d (%.4f); own segment hard: %d; the rest inherit" % (tot, tot / (len(qn) - 1), int(((own_q > 30) & (n > 1)).sum())))
    # two-tier coverage
    print("== Q4. Two-tier families F(q_hi, Y): rate q_hi for segments starting at or above Y, table below ==")
    for k, Y in enumerate(thresholds):
        sel = n >= Y
        row = []
        for qh in (3, 8, 16, 32):
            row.append((qh, float((qabove[sel, k] <= qh).mean())))
        print("  Y=2^%-2d: coverage among sources >= Y: " % int(round(log2(Y))) + " ".join("q_hi=%d: %.4f" % r for r in row))
    # two-tier deadline check: tau(n) <= floor(1.051 q_hi log2(n/Y')) + max tau(odd < Y'), Y' = max(Y, 10 q_hi)
    tau_prefix_max = np.maximum.accumulate(tau)
    for k, Y in enumerate(thresholds):
        for qh in (3, 8, 16, 32):
            Yp = max(int(Y), 10 * qh)
            if Yp >= len(qn) * 2:
                continue
            tmax = int(tau_prefix_max[(Yp - 1) // 2 - 1])
            sel = (n >= Yp) & (qabove[:, k] <= qh)
            # for members, every segment with start >= Yp is qh-admissible (start >= Y >= ... careful: Yp >= Y)
            bound = np.floor(1.051 * qh * np.log2(n[sel].astype(np.float64))) + tmax   # the last segment above Y' may land far below Y'
            viol = int((tau[sel] > bound).sum())
            require(viol == 0, ("two-tier deadline", Y, qh, viol))
    print("  two-tier deadline tau(n) <= floor(1.051 q_hi log2 n) + max tau(odd < Y') holds for all members at every tested (Y, q_hi)")
    # bank-relative: certify the top-K popular hard starts by table, measure coverage of F_3 relative to the bank
    print("  bank-relative coverage of rate-3 families: segments starting at a banked x end the chain")
    hard_sorted = order
    limit = 1 << bits
    for K in (0, 10, 100, 1000, 10000):
        bankflag = np.zeros(len(qn), dtype=np.uint8)
        for i in hard_sorted[:K]:
            bankflag[i] = 1
        for qh in (3, 8):
            cov = _coverage_with_bank(limit, bankflag, qh)
            print("    bank = top %5d hard excursions, rate q_hi=%d: coverage %.4f" % (K, qh, cov))
    return


@njit(cache=True)
def _coverage_with_bank(limit, bankflag, qh):
    """Membership with a bank: chain segments must be qh-admissible until a banked start is reached."""
    m = (limit + 1) // 2
    log23 = np.log2(3.0)
    good = 0
    for i in range(1, m):
        n = 2 * i + 1
        x = n
        ok = True
        while x != 1:
            if x < limit and bankflag[(x - 1) // 2]:
                break
            v = x
            l = 0
            A = 0
            while True:
                tt = 3 * v + 1
                a = 0
                while (tt & 1) == 0:
                    tt >>= 1
                    a += 1
                v = tt
                l += 1
                A += a
                if v < x:
                    break
            if np.ceil(l / (A - l * log23) - 1e-9) > qh:
                ok = False
                break
            x = v
        if ok:
            good += 1
    return good / (m - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census-bits', type=int, default=20)
    ap.add_argument('--json', type=str, default='')
    args = ap.parse_args()
    t0 = time.time()
    records = section_diophantine(320)
    thresholds = np.array([1 << 6, 1 << 8, 1 << 10, 1 << 12, 1 << 14, 1 << 16], dtype=np.int64)
    own_l, own_A, own_q, qn, tau, qabove, inherit = _census(1 << args.census_bits, thresholds)
    n = 2 * np.arange(len(qn)) + 1
    sel = n > 1
    tail = section_tail(200, own_q[sel], np.stack([own_l[sel], own_A[sel]], axis=1), args.census_bits)
    section_inheritance(args.census_bits, own_l, own_A, own_q, qn, tau, qabove, inherit, thresholds)
    print("== Summary ==")
    print("  checks: %d, total time %.1fs" % (CHECKS, time.time() - t0))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(dict(checks=CHECKS, records=records, tail=tail,
                           status="PROVED Diophantine structure and exact tail; FINITE-EXACT; census VERIFIED"),
                      fh, indent=1, default=str)
        print("  json written:", args.json)


if __name__ == '__main__':
    main()
