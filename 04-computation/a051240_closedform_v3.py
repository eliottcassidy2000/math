#!/usr/bin/env python3
"""
Closed-form c_4 computation for 4-uniform hypergraphs — v3 with fix.

Fix: enumerate all distinct permutations of composition values.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations
from time import perf_counter
import sys


def lcm(a, b):
    return a * b // gcd(a, b)


def lcm_list(lst):
    return reduce(lcm, lst) if lst else 1


def euler_phi(n):
    result = n
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def mobius(n):
    if n == 1:
        return 1
    factors = 0
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            m //= p
            factors += 1
            if m % p == 0:
                return 0
        p += 1
    if m > 1:
        factors += 1
    return (-1) ** factors


def divisors(n):
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)


def binomial(n, k):
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


def compute_ck_burnside(compressed_partition, k):
    """Reference: full Burnside iteration."""
    if not compressed_partition:
        return 0
    L = lcm_list([r for r, m in compressed_partition])
    total = 0
    for t in range(L):
        poly = [0] * (k + 1)
        poly[0] = 1
        for r, m in compressed_partition:
            g = gcd(t, r) if t > 0 else r
            length = r // g
            cnt = m * g
            factor = [0] * (k + 1)
            for j in range(cnt + 1):
                deg = j * length
                if deg > k:
                    break
                factor[deg] = binomial(cnt, j)
            new_poly = [0] * (k + 1)
            for a in range(k + 1):
                if poly[a] == 0:
                    continue
                for b in range(k + 1 - a):
                    if factor[b] == 0:
                        continue
                    new_poly[a + b] += poly[a] * factor[b]
            poly = new_poly
        total += poly[k]
    assert total % L == 0
    return total // L


def count_with_gcd_pattern(active_rs, active_ds, L):
    """Count t in [0,L) with gcd(t, r_a) = d_a for each active part.

    Uses Mobius inversion.
    """
    p = len(active_rs)
    if p == 0:
        return L

    q_values = [active_rs[i] // active_ds[i] for i in range(p)]
    divisor_lists = [divisors(q) for q in q_values]

    total = 0

    def recurse(idx, mu_prod, lcm_so_far):
        nonlocal total
        if idx == p:
            if L % lcm_so_far != 0:
                return
            total += mu_prod * (L // lcm_so_far)
            return
        for e in divisor_lists[idx]:
            mu_e = mobius(e)
            if mu_e == 0:
                continue
            new_lcm = lcm(lcm_so_far, active_ds[idx] * e)
            if new_lcm > L:
                continue
            recurse(idx + 1, mu_prod * mu_e, new_lcm)

    recurse(0, 1, 1)
    return total


def distinct_permutations(seq):
    """Generate all distinct permutations of a sequence."""
    if len(seq) <= 1:
        yield tuple(seq)
        return
    seen = set()
    for i in range(len(seq)):
        if seq[i] in seen:
            continue
        seen.add(seq[i])
        for rest in distinct_permutations(seq[:i] + seq[i+1:]):
            yield (seq[i],) + rest


def compute_c4_closedform(compressed_partition):
    """Compute c_4 using divisor-signature approach."""
    k = 4
    d = len(compressed_partition)
    if d == 0:
        return 0

    L = lcm_list([r for r, m in compressed_partition])

    # For each part, valid (divisor d_val, cycle length l) where l ∈ {1,2,3,4} and l|r
    part_options = []
    for r, m in compressed_partition:
        options = []
        for l in [1, 2, 3, 4]:
            if r % l == 0:
                d_val = r // l
                options.append((d_val, l))
        part_options.append(options)

    total = 0

    # Generate all compositions of 4 into exactly p positive parts, p=1..min(d,4)
    # (order matters — will enumerate permutations separately)
    compositions_by_len = {}
    def gen_comp(remaining, max_parts, min_val, prefix):
        if remaining == 0:
            key = len(prefix)
            if key not in compositions_by_len:
                compositions_by_len[key] = []
            compositions_by_len[key].append(tuple(prefix))
            return
        if max_parts == 0:
            return
        for v in range(min(remaining, 4), min_val - 1, -1):
            prefix.append(v)
            gen_comp(remaining - v, max_parts - 1, v, prefix)
            prefix.pop()

    gen_comp(4, min(d, 4), 1, [])
    # compositions_by_len[p] = sorted compositions of 4 into p parts

    for p in compositions_by_len:
        if p > d:
            continue

        for comp in compositions_by_len[p]:
            # For each distinct permutation of comp
            for perm in distinct_permutations(list(comp)):
                # perm[i] = j value for the i-th active part

                # Choose which p parts are active
                for part_indices in combinations(range(d), p):
                    # For each active part, find valid (d, l) options where l | j
                    valid = True
                    options_per_active = []
                    for idx_in_perm, part_idx in enumerate(part_indices):
                        j = perm[idx_in_perm]
                        r, m = compressed_partition[part_idx]
                        active_options = []
                        for d_val, l in part_options[part_idx]:
                            if j % l == 0:
                                num_cycles = m * d_val
                                num_chosen = j // l
                                if num_chosen <= num_cycles:
                                    bv = binomial(num_cycles, num_chosen)
                                    active_options.append((d_val, bv))
                        if not active_options:
                            valid = False
                            break
                        options_per_active.append(active_options)

                    if not valid:
                        continue

                    # Enumerate all (d_val) choices
                    def enum_opts(idx, rs, ds, bp):
                        nonlocal total
                        if idx == p:
                            cnt = count_with_gcd_pattern(rs, ds, L)
                            total += cnt * bp
                            return
                        pi = part_indices[idx]
                        r = compressed_partition[pi][0]
                        for d_val, bv in options_per_active[idx]:
                            enum_opts(idx + 1, rs + [r], ds + [d_val], bp * bv)

                    enum_opts(0, [], [], 1)

    assert total % L == 0, f"c4 not integer: {total}/{L}"
    return total // L


# Verify
print("=== Verifying c_4 closed-form v3 ===")
errors = 0
tested = 0
t0 = perf_counter()
for n in range(1, 16):
    def gen_parts(rem, maxp, parts):
        if rem == 0:
            yield list(parts)
            return
        for p in range(min(rem, maxp), 0, -1):
            for m in range(1, rem // p + 1):
                parts.append((p, m))
                yield from gen_parts(rem - m * p, p - 1, parts)
                parts.pop()

    for part in gen_parts(n, n, []):
        ref = compute_ck_burnside(part, 4)
        cf = compute_c4_closedform(part)
        tested += 1
        if ref != cf:
            if errors < 5:
                print(f"  MISMATCH n={n} {part}: ref={ref} cf={cf}")
            errors += 1

dt = perf_counter() - t0
if errors == 0:
    print(f"Tested {tested} partitions for n=1..15: ALL OK [{dt:.2f}s]")
else:
    print(f"Tested {tested} partitions: {errors} ERRORS [{dt:.2f}s]")

# Benchmark comparison at n=20
if errors == 0:
    print("\n=== Benchmark: Burnside vs closed-form for n=20 ===")
    parts_20 = list(gen_parts(20, 20, []))

    t0 = perf_counter()
    for part in parts_20:
        compute_ck_burnside(part, 4)
    dt_burn = perf_counter() - t0

    t0 = perf_counter()
    for part in parts_20:
        compute_c4_closedform(part)
    dt_cf = perf_counter() - t0

    print(f"  Burnside:    {dt_burn:.3f}s for {len(parts_20)} partitions")
    print(f"  Closed-form: {dt_cf:.3f}s for {len(parts_20)} partitions")
    print(f"  Speedup: {dt_burn/dt_cf:.1f}x")

print("\nDone.")
