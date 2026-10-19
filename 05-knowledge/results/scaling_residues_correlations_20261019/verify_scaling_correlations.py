#!/usr/bin/env python3
"""Finite exact controls for SCALING_RESIDUES_AND_CORRELATIONS.md.

No floats, network, Git writes, or external packages. All checks are explicit
guards, so python -O preserves them. Writes one JSON receipt beside this file.
The inverse-square controls verify only the transfer algebra, not an analytic
claim about actual escape thresholds.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import prod, isqrt, gcd
from pathlib import Path
import hashlib
import json


COUNTS = Counter()


def check(category, condition, *context):
    COUNTS[category] += 1
    if not condition:
        raise RuntimeError(f"{category}: {context}")


def primes(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def discrepancy(length, modulus):
    t = length % modulus
    return F(t * (modulus - t), modulus * length)


def local_count(p, multiplier, digit):
    u = pow(6, -1, p)
    ell, d = multiplier % p, digit % p
    if ell == 0:
        return (p - 2) * int(d not in {u, (-u) % p})
    if ell in {1, p - 1}:
        return p - (2 if d == 0 else 3 if d in {(2*u) % p, (-2*u) % p} else 4)
    exceptional = {((ell - 1)*u) % p, (-(ell - 1)*u) % p,
                   ((ell + 1)*u) % p, (-(ell + 1)*u) % p}
    check("four_distinct_affine_collisions", len(exceptional) == 4, p, ell)
    return p - (3 if d in exceptional else 4)


# Quantization: a midpoint is compatible with adjacent distinct residues.
for Q in range(3, 60):
    for k in range(1, Q - 1):
        midpoint = F(2*k + 1, 2*Q)
        check("quantization_midpoint", abs(midpoint - F(k, Q)) == F(1, 2*Q)
              and abs(midpoint - F(k + 1, Q)) == F(1, 2*Q), Q, k)
        for W in range(2, 10):
            check("adjacent_residues_differ", k % W != (k + 1) % W, Q, k, W)

# Same macroscopic density, different microscopic pair correlations.
for Q in range(4, 401, 4):
    A = {k for k in range(Q) if k % 2 == 0}
    B = {k for k in range(Q) if k % 4 in {0, 1}}
    check("microscopic_hostile", len(A) == len(B) == Q//2
          and sum((k+1) % Q in A for k in A) == 0
          and sum((k+1) % Q in B for k in B) == Q//4, Q)

# Every row, not merely the initial row, and moduli sharing factors with b.
for b in range(2, 9):
    for W in range(2, 25):
        columns = Counter((b*r+d) % W for r in range(W) for d in range(b))
        check("stationary_columns", all(columns[s] == b for s in range(W)), b, W)
        for r in range(W):
            for s in range(W):
                # bPD=DU checked on every delta-function basis vector.
                left = sum(int((b*r+d) % W == s)-int((b*r+d+1) % W == s)
                           for d in range(b))
                right = int(b*r % W == s)-int((b*(r+1)) % W == s)
                check("difference_intertwiner", left == right, b, W, r, s)
        for h in range(1, 4):
            L = b**h
            for r in range(W):
                actual = Counter((L*r+d) % W for d in range(L))
                tv = sum(abs(F(actual[s], L) - F(1, W)) for s in range(W)) / 2
                check("exact_transition_tv", tv == discrepancy(L, W), b, W, h, r)
                check("sharp_tv_envelope", tv <= F(W, 4*L), b, W, h, r)
                check("uniform_iff_divides", (tv == 0) == (L % W == 0), b, W, h)

# Joint prefix correlation control: direct finite digit enumeration.
for b in range(2, 7):
    for W in range(2, 16):
        Fs = [{0}, {r for r in range(W) if r % 2 == 0}, set(range(W))]
        Gs = [{W-1}, {r for r in range(W) if (r*r+r) % 3 == 0}]
        for n in range(0, 4):
            B = b**n
            for h in range(1, 4):
                L = b**h
                joint = Counter((a % W, (L*a+d) % W)
                                for a in range(B) for d in range(L))
                for f in Fs:
                    for g in Gs:
                        actual = sum(joint[r, s] for r in f for s in g) / F(B*L)
                        mean_f, mean_g = F(len(f), W), F(len(g), W)
                        bound = discrepancy(B, W) + mean_f*discrepancy(L, W)
                        check("joint_correlation_bound", abs(actual - mean_f*mean_g) <= bound,
                              b, W, n, h, sorted(f), sorted(g))

# All affine local types, including primes dividing the multiplier.
for p in primes(71):
    if p < 5:
        continue
    allowed = {r for r in range(p) if (36*r*r - 1) % p}
    for ell in range(p):
        for d in range(p):
            actual = sum((ell*r+d) % p in allowed for r in allowed)
            check("affine_prime_factor", actual == local_count(p, ell, d), p, ell, d)
            if ell in {1, p-1}:
                rho = 2 if d == 0 else 3 if (9*d*d - 1) % p == 0 else 4
                check("inherited_sum_difference_factor", actual == p-rho, p, ell, d)

# Composite CRT products, with both invertible and noninvertible multipliers.
for ps in [(5, 7), (5, 11), (5, 7, 11)]:
    W = prod(ps)
    allowed = {r for r in range(W) if all((36*r*r-1) % p for p in ps)}
    for ell in [0, 1, 2, 5, 6, 7, 10, 11, W-1]:
        for d in range(W):
            actual = sum((ell*r+d) % W in allowed for r in allowed)
            expected = prod(local_count(p, ell, d) for p in ps)
            check("affine_crt_factor", actual == expected, W, ell, d)

# Prefix residues really are the final digit blocks, exactly, including 0s.
digits = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3]
prefix = 0
for n, digit in enumerate(digits, 1):
    prefix = 10*prefix + digit
    for j in range(1, n+1):
        tail = 0
        for x in digits[n-j:n]:
            tail = 10*tail+x
        check("prefix_digit_block", prefix % (10**j) == tail, n, j)

# Exact inverse-square level weights (pi^2 cancels), all subset events of
# a short window. These are not asserted to be actual quadratic thresholds.
for M in [2, 5, 20, 100]:
    for L in range(1, 10):
        weights = [F(1, (n-1)**2) - F(1, n*n) for n in range(M+1, M+L+1)]
        total = sum(weights)
        check("inverse_square_window_telescope", total == F(1, M*M)-F(1, (M+L)**2), M, L)
        ratio = max(weights)/min(weights)
        for mask in range(1<<L):
            count = mask.bit_count()
            measure = sum(w for i, w in enumerate(weights) if mask >> i & 1) / total
            density = F(count, L)
            check("relative_weight_transfer", density/ratio <= measure <= density*ratio,
                  M, L, mask)
        for W in range(2, 10):
            for m in range(W):
                allowed = [int(gcd((36*a*a-1)*(36*(m-a)**2-1), W) == 1)
                           for a in range(M+1, M+L+1)]
                # The transfer bound applies even to rare four-form events.
                density = F(sum(allowed), L)
                measure = sum(w*v for w, v in zip(weights, allowed))/total
                check("relative_transfer_rare_event", density/ratio <= measure <= density*ratio,
                      M, L, W, m)

# Gauge theorem: factorized gauges preserve every Boolean diamond. Independent
# decimal mantissas of actual antichain counts supply the hostile control.
for size in range(1, 8):
    ts = [F(i+2, i+1) for i in range(size)]
    c = {S: F(3, 7)*prod(ts[i] for i in range(size) if S >> i & 1)
         for S in range(1<<size)}
    H = {S: 1+(17*S*S+5*S+size) % 101 for S in range(1<<size)}
    for S in range(1<<size):
        available = [i for i in range(size) if not S >> i & 1]
        for a, b in combinations(available, 2):
            Sa, Sb, Sab = S | 1<<a, S | 1<<b, S | 1<<a | 1<<b
            old = F(H[S]*H[Sab], H[Sa]*H[Sb])
            new = c[S]*H[S]*c[Sab]*H[Sab]/(c[Sa]*H[Sa]*c[Sb]*H[Sb])
            check("factorized_gauge", old == new, size, S, a, b)
    reconstructed = {S: c[0]*prod(c[1<<i]/c[0] for i in range(size) if S >> i & 1)
                     for S in range(1<<size)}
    check("gauge_reconstruction", c == reconstructed, size)

H = (24, 6, 6, 2)
e = (2, 1, 1, 1)
h = tuple(F(v, 10**k) for v, k in zip(H, e))
old, new = F(H[0]*H[3], H[1]*H[2]), h[0]*h[3]/(h[1]*h[2])
check("antichain_mantissa_hostile", old == F(4, 3) and new == F(2, 15)
      and new*10**(e[0]+e[3]-e[1]-e[2]) == old)

receipt = {
    "status": "FINITE-EXACT; all active guards passed; no floating-point claims",
    "counts": dict(sorted(COUNTS.items())),
    "total_guards": sum(COUNTS.values()),
    "antichain_hostile": {"integer_counts": H, "decimal_exponents": e,
                           "integer_ratio": str(old), "mantissa_ratio": str(new)},
    "scope": {
        "transition_bases": "2..8",
        "transition_moduli": "2..24",
        "transition_digit_block_lengths": "1..3",
        "affine_prime_moduli": "all primes 5..71; every multiplier and shift",
        "crt_wheels": [35, 55, 385],
        "cusp_controls": "inverse-square model weights only; phase theorem not numerically certified here"
    },
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
target = Path(__file__).with_name("scaling_correlations_receipt.json")
target.write_text(json.dumps(receipt, indent=2)+"\n")
print(json.dumps(receipt, indent=2))
