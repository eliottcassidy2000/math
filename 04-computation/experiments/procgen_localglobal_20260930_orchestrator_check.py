#!/usr/bin/env python3
"""Orchestrator audit of lane `localglobal` (LRC / Collatz / primes), written from the note's statements; the
lane's code was not read.

  1. Tight LRC sets: for k = 3, 4 (speeds <= 21) and k = 5 (speeds <= 13), the primitive speed sets with
     M(V) = max_t min_i ||t v_i|| = 1/(k+1) are exactly {1..k}, {1,3,4,7} (k=4), {1,3,4,5,9} (k=5); covering
     sets (every q <= k+1 divides a speed) have M >= 2/(2k+1).
  2. kappa(2,3) = 1/5: over t = m/q (q <= 400, gcd(q,6) = 1) the best min over the box {2^s 3^r : s, r <= 14}
     is 1/5 (attained at 1/5).
  3. Gate local densities: for gates (p, a) with 12 <= p <= 22 and a near p*log_3 2, the number of words
     with c_w = 0 mod q, for primes q | D, is close to C(p,a)/q; and N_5 = 0 for every shape a = p-2 with
     5 | D (p <= 22).
  4. The deep well {1..12, 182} is lonely at t = 2/27 at level 1/14.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb, log
from functools import reduce


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def dist(t, v):
    x = (t * v) % 1
    return min(x, 1 - x)


def M(V):
    cands = set()
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            for d in (V[i] + V[j], abs(V[i] - V[j])):
                for m in range(1, d):
                    cands.add(F(m, d))
        for m in range(V[i]):
            cands.add(F(2 * m + 1, 2 * V[i]))
    best = F(0)
    for t in cands:
        val = min(dist(t, v) for v in V)
        if val > best:
            best = val
    return best


def tight_sets(k, vmax):
    out = []
    for V in combinations(range(1, vmax + 1), k):
        if reduce(gcd, V) != 1:
            continue
        if M(V) == F(1, k + 1):
            out.append(V)
    return out


t3 = tight_sets(3, 21)
t4 = tight_sets(4, 21)
t5 = tight_sets(5, 13)
check(t3 == [(1, 2, 3)] and t4 == [(1, 2, 3, 4), (1, 3, 4, 7)] and t5 == [(1, 2, 3, 4, 5), (1, 3, 4, 5, 9)],
      f"tight primitive sets: k=3 {t3}; k=4 {t4}; k=5 (speeds <= 13) {t5}")
viol = 0
for k, vmax in ((3, 21), (4, 21)):
    for V in combinations(range(1, vmax + 1), k):
        if reduce(gcd, V) != 1:
            continue
        if all(any(v % q == 0 for v in V) for q in range(2, k + 2)):
            if M(V) < F(2, 2 * k + 1):
                viol += 1
check(viol == 0, "covering sets (every 2 <= q <= k+1 divides a speed) have M >= 2/(2k+1), k = 3, 4, speeds <= 21")

box = [2 ** s * 3 ** r for s in range(15) for r in range(15)]
best, arg = F(0), None
for q in range(5, 401):
    if gcd(q, 6) != 1:
        continue
    for m in range(1, q):
        if gcd(m, q) != 1:
            continue
        t = F(m, q)
        val = min(dist(t, b) for b in box)
        if val > best:
            best, arg = val, t
check(best == F(1, 5), f"kappa(2,3): best min over the box {{2^s 3^r}} for t = m/q (q <= 400, gcd(q,6)=1) is {best} at t = {arg}")


def factor(n):
    fs, d = [], 2
    while d * d <= n:
        while n % d == 0:
            fs.append(d)
            n //= d
        d += 1
    if n > 1:
        fs.append(n)
    return fs


def carry_counts(p, a, mods):
    """number of words (length p, a ones) with c_w = 0 mod m, for each m in mods (DP over positions)."""
    res = {}
    for m in mods:
        # c_w = sum_i 3^(a-1-i) 2^(s_i): process positions s = 0..p-1, state (ones used i, residue)
        dp = {(0, 0): 1}
        for s in range(p):
            nd = {}
            for (i, r), c in dp.items():
                nd[(i, r)] = nd.get((i, r), 0) + c                 # a zero at position s
                if i < a:
                    r2 = (r + pow(3, a - 1 - i, m) * pow(2, s, m)) % m
                    nd[(i + 1, r2)] = nd.get((i + 1, r2), 0) + c
            dp = nd
        res[m] = dp.get((a, 0), 0)
    return res


ratios = []
for p in range(12, 23):
    for a in range(max(1, int(p * log(2) / log(3)) - 1), int(p * log(2) / log(3)) + 2):
        D = abs(2 ** p - 3 ** a)
        if D < 2:
            continue
        qs = sorted(set(q for q in factor(D) if q > 3))
        if not qs:
            continue
        C = comb(p, a)
        cnt = carry_counts(p, a, qs)
        for q in qs:
            if q <= 2000:
                ratios.append(cnt[q] * q / C)
mean = sum(ratios) / len(ratios)
check(0.9 < mean < 1.1, f"gate local densities at primes q | D (12 <= p <= 22, a near the critical ratio, q <= 2000): mean of N_q q / C = {mean:.3f} over {len(ratios)} (gate, prime) pairs")
zero5 = []
for p in range(4, 23):
    a = p - 2
    D = abs(2 ** p - 3 ** a)
    if D % 5 == 0:
        zero5.append((p, carry_counts(p, a, [5])[5]))
check(all(n == 0 for _, n in zero5) and len(zero5) > 0, f"N_5 = 0 for every shape a = p - 2 with 5 | D (p <= 22): {zero5}")

V = list(range(1, 13)) + [182]
t = F(2, 27)
check(min(dist(t, v) for v in V) == F(2, 27) and F(2, 27) >= F(1, 14), "the deep well {1..12, 182} is lonely at t = 2/27 at level 1/14 (min distance 2/27)")
