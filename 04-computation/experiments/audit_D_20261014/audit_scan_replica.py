#!/usr/bin/env python3
"""Audit D: replicate the A-centric scan of collatz_sources_check.py (c) with an independent implementation.

Safety by the THM-4072 follower graph (integers), agreement by bitmask XOR + trailing-zero count, constants both
in the original float form and exactly (mpmath).  Replicates: number of A <= AMAX with safe depth >= 20, the deepest
ones, the number of (A, u/v) pairs with nonempty agreement (v odd <= 45, |u| <= 12 v, gcd = 1), violations, the
tightest ratio h/bound, and additionally the minimal ADDITIVE slack bound - h.
"""
import math
import sys
import time
import mpmath

mpmath.mp.dps = 50
AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
VMAX = 45
OUT = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_D_20261014/audit_scan_replica.out'
lines = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def greedy_d(n):
    a, d = 1, []
    for k in range(1, n + 1):
        t = 3 * a
        dk = 1 if t >= 2 ** k else 0
        a = t - dk * 2 ** k
        d.append(dk)
    return d


D1 = [None] + greedy_d(500)


def fstep(q, s):
    nxt = D1[q + 1]
    if s == nxt:
        return q + 1
    if s == 0 and nxt == 1:
        return 0
    return None


LOG3F = math.log2(3)
CF = 1 - LOG3F / 3
LOG3 = mpmath.log(3, 2)
CC = 1 - LOG3 / 3

t0 = time.time()
found = []
for A in range(1, AMAX + 1):
    a, q, mask, n = A, 0, 0, 0
    while True:
        b = a & 1
        q2 = fstep(q, b)
        if q2 is None:
            break
        q = q2
        mask |= b << n
        n += 1
        a = (3 * a + b) >> 1
        if n >= 120:
            break
    if n >= 20:
        found.append((A, n, mask))
found.sort(key=lambda t: -t[1])
say(f"A <= {AMAX} with safe depth >= 20: {len(found)}; deepest: " +
    ", ".join(f"A={A} (depth {D})" for A, D, _ in found[:6]) + f"  [{time.time() - t0:.1f}s]")

# rationals
xs = []
for v in range(1, VMAX + 1, 2):
    lv_f = math.floor(math.log2(6.3 * v) / CF) + 1
    lv_x = int(mpmath.floor(mpmath.log(mpmath.mpf('6.3') * v, 2) / CC)) + 1
    for u in range(-12 * v, 12 * v + 1):
        if math.gcd(u, v) != 1:
            continue
        ax = abs(u) / v
        n1_f = 0 if ax <= 1 / 6 else math.ceil(math.log2(6 * ax) / CF)
        n1_x = 0 if 6 * abs(u) <= v else int(mpmath.ceil(mpmath.log(mpmath.mpf(6 * abs(u)) / v, 2) / CC))
        j, mask = u, 0
        for n in range(40):
            b = j & 1
            mask |= b << n
            j = (3 * j + v) // 2 if b else j // 2
        xs.append((u, v, mask, n1_f, n1_x, lv_f, lv_x))
say(f"rationals: {len(xs)}; float/exact constant mismatches: n1 {sum(1 for x in xs if x[3] != x[4])}, "
    f"ell {sum(1 for x in xs if x[5] != x[6])}")


def nfac(code_bits, lo, hi, l):
    s = set()
    for i in range(max(lo, 0), hi - l + 1):
        s.add((code_bits >> i) & ((1 << l) - 1))
    return len(s)


def nfac_full(code_bits, hi, l):     # theorem's definition: factors of e[0, h)
    return nfac(code_bits, 0, hi, l)


pairs = viol = viol_thm = 0
best_ratio = []
min_slack = (1e9, None)
for A, D, E in found:
    LA = math.log2(A + 1)
    full = (1 << D) - 1
    for (u, v, X, n1f, n1x, lvf, lvx) in xs:
        z = (E ^ X) & full
        h = D if z == 0 else (z & -z).bit_length() - 1
        if h == 0:
            continue
        pairs += 1
        st = (21 * v) // 10 + 1
        bstate = LOG3F * (n1f + st) + LA
        lower = min(bstate, LOG3F * (n1f + lvf) + LA)
        if h / lower > 0.5 or lower - h < 9:
            N = nfac(E, n1f, h, lvf)                       # original script's convention: factors in [n1, h)
            bound = min(bstate, LOG3F * (n1f + lvf + N) + LA)
            Nt = nfac_full(E, h, lvx)                        # theorem's convention: factors of e[0, h), exact consts
            bound_thm = float(LOG3 * (n1x + min(st, lvx + Nt))) + LA
            if h >= bound:
                viol += 1
            if h >= bound_thm:
                viol_thm += 1
            best_ratio.append((h / bound, A, u, v, h, round(bound, 2)))
            if bound_thm - h < min_slack[0]:
                min_slack = (bound_thm - h, (A, u, v, h, round(bound_thm, 3)))
best_ratio.sort(reverse=True)
say(f"pairs with nonempty agreement: {pairs}; violations (original convention) {viol}; violations (theorem "
    f"convention, exact constants) {viol_thm}  [{time.time() - t0:.1f}s]")
say(f"tightest h/bound: {[(round(t[0], 3),) + t[1:] for t in best_ratio[:5]]}")
say(f"minimal additive slack bound - h (theorem convention): {min_slack[0]:.4f} at (A, u, v, h, bound) = "
    f"{min_slack[1]}")
with open(OUT, 'w') as fh:
    fh.write("\n".join(lines) + "\n")
