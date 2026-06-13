#!/usr/bin/env python3
"""
Closed-form c_4 computation for 4-uniform hypergraphs (A051240).

Key insight: [x^4] prod_a (1+x^{l_a})^{c_a} is nonzero only when l_a ≤ 4
for the "active" parts. So instead of iterating t=0..L-1, we:

1. Select up to 4 "active" parts (those contributing j_a > 0)
2. For each, choose a cycle length l_a ∈ {1,2,3,4} dividing r_a
3. Find compositions j_1+...+j_p = 4 with l_a | j_a
4. Count t values producing these gcd patterns via Mobius inversion
5. Multiply by the binomial coefficients

This reduces compute_ck from O(L) to O(d^4 * small) where d = #distinct parts
and L = lcm(parts) can be astronomically large.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from functools import reduce
from itertools import combinations, combinations_with_replacement
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
                return 0  # p^2 | n
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


def count_with_gcd_pattern(active_parts, active_divs, L):
    """Count t in [0,L) with gcd(t, r_a) = d_a for each active part a.

    Uses Mobius inversion:
    #{t: gcd(t,r)=d for all (r,d)} = sum_{e_i | (r_i/d_i)} prod mu(e_i) * L/lcm(d_i*e_i)
    """
    p = len(active_parts)
    if p == 0:
        return L

    # Get the q_i = r_i / d_i values
    q_values = [active_parts[i] // active_divs[i] for i in range(p)]

    # Iterate over all divisor tuples (e_1, ..., e_p) with e_i | q_i
    divisor_lists = [divisors(q) for q in q_values]

    total = 0
    # Recursive or iterative enumeration of divisor tuples
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
            new_lcm = lcm(lcm_so_far, active_divs[idx] * e)
            if new_lcm > L:
                continue
            recurse(idx + 1, mu_prod * mu_e, new_lcm)

    recurse(0, 1, 1)
    return total


def compute_c4_closedform(compressed_partition):
    """Compute c_4 using the divisor-signature approach.

    Instead of iterating t=0..L-1, enumerate all valid "active part" selections
    and compute the count via Mobius inversion.
    """
    k = 4
    d = len(compressed_partition)
    if d == 0:
        return 0

    L = lcm_list([r for r, m in compressed_partition])
    n_total = sum(r * m for r, m in compressed_partition)

    # For each part, precompute the valid (divisor, cycle_length) pairs
    # where cycle_length ∈ {1, 2, 3, 4} and divides r
    part_options = []  # part_options[i] = list of (div, length) pairs
    for r, m in compressed_partition:
        options = []
        for l in [1, 2, 3, 4]:
            if r % l == 0:
                d_val = r // l
                options.append((d_val, l))
        part_options.append(options)

    total = 0

    # Enumerate all ways to select up to 4 "active" parts and assign j values.
    # An "active" selection is: for each active part i, choose (d_i, l_i) from options,
    # and j_i (a positive multiple of l_i with j_i ≤ 4).
    # Sum of j_i = 4.

    # We enumerate by "patterns": how many active parts and what (part_idx, d, l, j) each has.
    # Parts can repeat (e.g., two cycles from the same part type).
    # But in our formulation, each part type appears once with multiplicity m_a,
    # and we pick j_a elements from that type's cycles. So each part can be active at most once.

    # Actually, each "part" in the compressed partition is a group of m_a identical cycles.
    # The j_a represents the total contribution from all cycles in that group.
    # So each compressed part appears at most once.

    # Enumerate subsets of {0,...,d-1} with |subset| ≤ 4, and for each, assignments.
    # For efficiency, enumerate by composition of 4.

    # Compositions of 4 into at most d parts (each at most 4):
    # (4), (3,1), (2,2), (2,1,1), (1,1,1,1)
    # For each composition, assign to distinct parts.

    compositions = []
    # Generate all compositions of 4 into at most min(d, 4) positive parts
    def gen_comp(remaining, max_parts, min_val, prefix):
        if remaining == 0:
            compositions.append(tuple(prefix))
            return
        if max_parts == 0:
            return
        for v in range(min(remaining, 4), min_val - 1, -1):
            prefix.append(v)
            gen_comp(remaining - v, max_parts - 1, v, prefix)
            prefix.pop()

    gen_comp(4, min(d, 4), 1, [])

    for comp in compositions:
        p = len(comp)  # number of active parts

        # Choose which parts are active: p distinct parts from d
        if p > d:
            continue

        for part_indices in combinations(range(d), p):
            # For each active part, find valid (d, l) options where l | j
            valid = True
            options_per_active = []
            for idx_in_comp, part_idx in enumerate(part_indices):
                j = comp[idx_in_comp]
                r, m = compressed_partition[part_idx]
                active_options = []
                for d_val, l in part_options[part_idx]:
                    if j % l == 0:
                        num_cycles = m * d_val  # c_a = m * gcd(t, r) = m * d
                        num_chosen = j // l
                        if num_chosen <= num_cycles:
                            binom_val = binomial(num_cycles, num_chosen)
                            active_options.append((d_val, l, binom_val))
                if not active_options:
                    valid = False
                    break
                options_per_active.append(active_options)

            if not valid:
                continue

            # Enumerate all combinations of (d, l) choices for active parts
            def enum_options(idx, active_rs, active_ds, binom_prod):
                nonlocal total
                if idx == p:
                    # Count t values with gcd(t, r_a) = d_a for active parts
                    cnt = count_with_gcd_pattern(active_rs, active_ds, L)
                    total += cnt * binom_prod
                    return
                part_idx = part_indices[idx]
                r = compressed_partition[part_idx][0]
                for d_val, l, bv in options_per_active[idx]:
                    enum_options(idx + 1, active_rs + [r], active_ds + [d_val],
                                binom_prod * bv)

            enum_options(0, [], [], 1)

    # Account for the "no active parts" case: j_a = 0 for all a.
    # This contributes [x^0] product = 1 for each of L values of t.
    # But we need [x^4], not [x^0]. So this contributes 0 (since 4 ≠ 0).

    assert total % L == 0, f"c4 not integer: {total}/{L}"
    return total // L


def a_4uniform(n):
    """Compute a(n) for 4-uniform hypergraphs using closed-form c_4."""
    if n < 4:
        return 1

    total = Fraction(0)
    fact_n = factorial(n)
    count = 0

    def enumerate_partitions(remaining, max_part, parts_so_far):
        nonlocal total, count
        if remaining == 0:
            count += 1
            ck = compute_c4_closedform(parts_so_far)
            z = 1
            for r, m in parts_so_far:
                z *= r ** m * factorial(m)
            pc = fact_n // z
            total += Fraction(pc * (1 << ck), fact_n)
            return

        for p in range(min(remaining, max_part), 0, -1):
            max_m = remaining // p
            for m in range(1, max_m + 1):
                parts_so_far.append((p, m))
                enumerate_partitions(remaining - m * p, p - 1, parts_so_far)
                parts_so_far.pop()

    enumerate_partitions(n, n, [])
    result = int(total)
    assert total == result
    return result, count


# Known A051240 values
known = {0: 1, 1: 1, 2: 1, 3: 1, 4: 2, 5: 6, 6: 156, 7: 7013320,
         8: 29281354514767168}


# Verification
print("=== Verifying c_4 closed-form against Burnside reference ===")
errors = 0
tested = 0
t0 = perf_counter()
for n in range(1, 20):
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
            print(f"  MISMATCH n={n} {part}: ref={ref} cf={cf}")
            errors += 1

dt = perf_counter() - t0
print(f"Tested {tested} partitions for n=1..19: {'ALL OK' if errors == 0 else f'{errors} ERRORS'} [{dt:.2f}s]")

# Now verify full A051240 values
print("\n=== Computing A051240 values ===")
for n in range(0, 15):
    t0 = perf_counter()
    result, nparts = a_4uniform(n)
    dt = perf_counter() - t0
    expected = known.get(n, "?")
    if expected != "?":
        match = "OK" if result == expected else f"FAIL (expected {expected})"
    else:
        match = ""
    print(f"  a({n:2d}) = {result}  [{dt:.4f}s, {nparts} parts]  {match}")
    if dt > 60:
        print("  Stopping: took too long")
        break

print("\nDone.")
