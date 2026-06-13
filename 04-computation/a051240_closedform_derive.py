#!/usr/bin/env python3
"""
Derive closed-form edges formula for 4-uniform hypergraphs (A051240).

For a permutation with cycle type lambda (compressed as [(r_a, m_a)]),
the number of orbits on 4-subsets is:

  c_4(lambda) = (1/L) * sum_{t=0}^{L-1} [x^4] prod_a (1 + x^{r_a/gcd(t,r_a)})^{m_a*gcd(t,r_a)}

A 4-subset fixed by sigma^t must be a union of complete cycles of sigma^t.
The cycles of sigma^t from a cycle of length r_a are: m_a*gcd(t,r_a) cycles each of length r_a/gcd(t,r_a).

Ways to pick 4 elements from these cycles to form a union:
  - 4 from one cycle type (choosing cycles of size 1,2,4 or combinations that sum to 4)

The [x^4] extraction for each partition gives us specific combinatorial terms.

Let's think about what 4-subsets look like:
  From cycles of sigma^t with lengths l_a and counts c_a:
  [x^4] prod_a (1+x^{l_a})^{c_a} = sum over compositions of 4 into parts from {l_a}

For each divisor d of L, we group all t with gcd(t,L)=d (there are phi(L/d) such t).
But different t with same gcd(t,L) can give different gcd(t,r_a) for individual parts.

Wait - that's the issue that made the divisor optimization wrong for general partitions.
Let me instead find N_4(r) = orbits on 4-subsets of a single cycle of length r.

For a single cycle of length r:
  N_4(r) = (1/r) sum_{d|r} phi(d) * C(r/d, 4/d) where the sum is over d|gcd(r,4)
  Wait, that's not right either. Let me think more carefully.

For a single cycle (r), sigma has one cycle of length r.
sigma^t has gcd(t,r) cycles of length r/gcd(t,r).
f_4(t) = [x^4] (1+x^{r/gcd(t,r)})^{gcd(t,r)}

This equals C(gcd(t,r), 4/(r/gcd(t,r))) if r/gcd(t,r) divides 4, else 0 for the single-cycle terms.

Actually, let me just compute it numerically and compare. The goal is to find formulas
for each "interaction type" and verify they match the Burnside sum.

Author: opus-2026-03-08-S48
"""

from math import gcd, factorial
from fractions import Fraction
from functools import reduce
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
    """Compute c_k via full Burnside iteration (reference)."""
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


def N4_single(r):
    """Number of orbits on 4-subsets within a single cycle of length r.

    These are necklaces of r beads with exactly 4 marked, under cyclic rotation.
    = (1/r) * sum_{d|gcd(r,4)} phi(d) * C(r/d, 4/d)
    """
    # Burnside on cyclic group Z_r acting on 4-subsets of [r]
    # sigma^t fixes a 4-subset iff it's a union of orbits of sigma^t
    # sigma^t has gcd(t,r) orbits of size r/gcd(t,r)
    # Need to pick orbits summing to 4

    # More precisely: N4(r) = (1/r) sum_{t=0}^{r-1} f(t)
    # where f(t) = [x^4] (1 + x^{r/gcd(t,r)})^{gcd(t,r)}
    total = 0
    for t in range(r):
        g = gcd(t, r) if t > 0 else r
        length = r // g
        cnt = g
        # [x^4] (1+x^length)^cnt
        if length > 4:
            continue
        # Enumerate ways: choose j copies of x^length such that j*length = 4
        if 4 % length == 0:
            j = 4 // length
            if j <= cnt:
                total += binomial(cnt, j)
    assert total % r == 0, f"N4({r}) not integer: {total}/{r}"
    return total // r


def N4_pair(r, s):
    """Number of orbits on 4-subsets spanning two cycle types r and s.

    Elements from a cycle of length r and a cycle of length s.
    The joint action is Z_{lcm(r,s)} acting on pairs of subsets.

    Actually, for two cycles of lengths r and s, the permutation sigma acts
    as a product of two independent cyclic permutations. sigma^t acts on the
    r-cycle with gcd(t,r) orbits and on the s-cycle with gcd(t,s) orbits.

    We need 4-subsets that use elements from both cycles.
    Specifically, pick j elements from the r-cycle and 4-j from the s-cycle (j=1,2,3).

    For a fixed t, f(t) = sum_{j=1}^{3} [x^j](1+x^{r/g_r})^{g_r} * [x^{4-j}](1+x^{s/g_s})^{g_s}
    where g_r = gcd(t,r), g_s = gcd(t,s).

    Total orbits = (1/lcm(r,s)) * sum_{t=0}^{lcm(r,s)-1} f(t)
    """
    L = lcm(r, s)
    total = 0
    for t in range(L):
        gr = gcd(t, r) if t > 0 else r
        gs = gcd(t, s) if t > 0 else s
        lr = r // gr
        ls = s // gs

        # Compute [x^j](1+x^lr)^gr for j=0..4
        poly_r = [0] * 5
        for j in range(gr + 1):
            deg = j * lr
            if deg > 4:
                break
            poly_r[deg] = binomial(gr, j)

        poly_s = [0] * 5
        for j in range(gs + 1):
            deg = j * ls
            if deg > 4:
                break
            poly_s[deg] = binomial(gs, j)

        # Cross terms: j from r, 4-j from s, with j >= 1 and 4-j >= 1
        for j in range(1, 4):
            total += poly_r[j] * poly_s[4 - j]

    assert total % L == 0, f"N4_pair({r},{s}) not integer: {total}/{L}"
    return total // L


# For the full formula, I need interaction terms for:
# 1. All 4 from one cycle type a: involves N4_single(r_a) and multi-cycle terms
# 2. Split 3+1 across two cycle types
# 3. Split 2+2 across two cycle types
# 4. Split 2+1+1 across three cycle types
# 5. Split 1+1+1+1 across four cycle types

# Actually, the proper decomposition for compressed partitions is:
# Each "part" (r_a, m_a) contributes m_a cycles of length r_a.
# We pick 4 elements total from these cycles.
# The 4 elements come from at most 4 cycles.

# Let me think about this differently. The key formula is:
# c_4(lambda) = (1/L) sum_{t=0}^{L-1} [x^4] prod_a (1+x^{l_a(t)})^{c_a(t)}
# where l_a(t) = r_a/gcd(t,r_a), c_a(t) = m_a*gcd(t,r_a)

# To get a closed form, I can expand [x^4] of the product.
# Let's denote for each part a: P_a(x,t) = (1+x^{l_a})^{c_a}
# Then [x^4] prod = sum of products of [x^{j_a}] P_a where sum j_a = 4.

# The compositions of 4 into non-negative parts are:
# With at most d parts (where d = number of distinct part values):
# (4,0,...), (3,1,0,...), (2,2,0,...), (2,1,1,0,...), (1,1,1,1,0,...), etc.

# For each composition (j_1, ..., j_d) with sum = 4 and each j_a >= 0:
# contribution = prod_a [x^{j_a}] (1+x^{l_a})^{c_a}

# And [x^j] (1+x^l)^c = C(c, j/l) if l|j, else 0.

# So the sum over t becomes:
# sum_t sum_{compositions} prod_a C(c_a(t), j_a/l_a(t)) * [l_a(t)|j_a]

# This is tractable for small k!

# Let me verify: for a single part (r, m), what are the terms?
# j = 4, all from this part type.
# For each t, l = r/gcd(t,r), c = m*gcd(t,r).
# [x^4](1+x^l)^c = C(c, 4/l) if l|4.
# So l must be 1, 2, or 4.

# l=1 means gcd(t,r)=r, i.e., r|t. Count: L/r such values of t (0, r, 2r, ..., L-r).
#   Contribution per t: C(m*r, 4)
# l=2 means gcd(t,r)=r/2, i.e., r/2 exists and gcd(t,r)=r/2.
#   Only if r is even. The count of t with gcd(t,r)=r/2 is phi(2)=1... no wait.
#   gcd(t,r) = d iff d|t and gcd(t/d, r/d)=1.
#   So count of t in [0,L) with gcd(t,r) = d is: (L/r) * phi(r/d) if d|r.
#   Wait, more carefully: for t in [0,L), gcd(t,r) = d iff d|t and gcd(t/d, r/d) = 1.
#   The number of t in [0,L) divisible by d is L/d.
#   Among those, t/d ranges over [0, L/d), and we need gcd(t/d, r/d) = 1.
#   Count = sum_{t'=0}^{L/d-1} [gcd(t', r/d) = 1] = (L/d) * phi(r/d) / (r/d) = (L*phi(r/d))/r

# So for gcd(t,r) = d where d|r:
#   count(d) = L * phi(r/d) / r

# This means we can replace the sum over t with a sum over divisors d of r!
# But ONLY for a single part. For multiple parts, we'd need joint divisor conditions.

# Actually wait — that's the key issue. For a single part (r,m):
# c_4(r,m) = (1/L) sum_{d|r} count(d) * C(m*d, 4/(r/d)) * [r/d | 4]
#          = (1/L) sum_{d|r, (r/d)|4} L*phi(r/d)/r * C(m*d, 4*d/r)
#          = (1/r) sum_{d|r, (r/d)|4} phi(r/d) * C(m*d, 4*d/r)

# This is a single-part formula! Let me verify it.
print("=== Verifying single-part formula for c_4 ===")
for r in range(1, 13):
    for m in range(1, min(8, 20 // r + 1)):
        # Burnside reference
        ref = compute_ck_burnside([(r, m)], 4)

        # Closed form
        cf = 0
        for d in range(1, r + 1):
            if r % d != 0:
                continue
            l = r // d  # cycle length in sigma^t
            if 4 % l != 0:
                continue
            j = 4 // l  # number of cycles to pick
            c = m * d   # number of available cycles
            cf += euler_phi(l) * binomial(c, j)
        assert cf % r == 0, f"r={r}, m={m}: {cf}/{r}"
        cf //= r

        if ref != cf:
            print(f"  MISMATCH r={r} m={m}: ref={ref} cf={cf}")
        # else:
        #     print(f"  OK r={r} m={m}: c4={ref}")

print("Single-part formula verified!")

# Now for two parts (r, mr) and (s, ms), the cross terms.
# We need [x^4] of the product, but only the cross terms (j_r >= 1 and j_s >= 1).
#
# For each t: contributions from (j_r, j_s) = (1,3), (2,2), (3,1)
# where [x^{j_r}](1+x^{l_r})^{c_r} * [x^{j_s}](1+x^{l_s})^{c_s}
#
# For this we need: l_r | j_r and l_s | j_s.
# l_r = r/gcd(t,r), l_s = s/gcd(t,s).
#
# Let d_r = gcd(t,r), d_s = gcd(t,s). We need to count t with given (d_r, d_s).
# t must satisfy: d_r | t and gcd(t/d_r, r/d_r) = 1
#                 d_s | t and gcd(t/d_s, s/d_s) = 1
# So lcm(d_r, d_s) | t, and further coprimality conditions.
#
# The count of t in [0, L) with gcd(t, r) = d_r and gcd(t, s) = d_s is:
#   sum_{t=0}^{L-1} [gcd(t,r)=d_r] * [gcd(t,s)=d_s]
#
# Using Mobius: [gcd(t,r)=d_r] = sum_{e|gcd(t/d_r, r/d_r)} mu(e) (if d_r|t, else 0)
#
# This gets complicated. Let me try a different approach: just compute the pair
# interaction numerically for all (r, s, mr, ms) we need and see if a pattern emerges.

print("\n=== Pair interaction terms for c_4 ===")
print("Testing c_4([(r,mr),(s,ms)]) = single(r,mr) + single(s,ms) + cross(r,s,mr,ms)")

for r in range(1, 8):
    for s in range(r + 1, 8):
        for mr in range(1, 4):
            for ms in range(1, 4):
                full = compute_ck_burnside([(r, mr), (s, ms)], 4)
                single_r = compute_ck_burnside([(r, mr)], 4)
                single_s = compute_ck_burnside([(s, ms)], 4)
                cross = full - single_r - single_s

                # Try to find a formula for cross
                g = gcd(r, s)
                L = lcm(r, s)

                # The cross term should depend on (r, s, mr, ms, g)
                # For k=3, cross was: mr*ms * pair_val(r,s)
                # For k=4, cross involves (1,3), (2,2), (3,1) splits

                if cross != 0 and (mr == 1 and ms == 1):
                    print(f"  r={r} s={s} mr={mr} ms={ms}: cross={cross}, g={g}, L={L}")

# Let me focus on mr=ms=1 first (one cycle of each length)
print("\n=== Cross term for single cycles: c_4_cross(r, s) ===")
for r in range(1, 10):
    for s in range(r + 1, 10):
        full = compute_ck_burnside([(r, 1), (s, 1)], 4)
        sr = compute_ck_burnside([(r, 1)], 4)
        ss = compute_ck_burnside([(s, 1)], 4)
        cross = full - sr - ss
        g = gcd(r, s)
        L = lcm(r, s)

        # Compute via sum over t
        # For (j_r, j_s) with j_r + j_s = 4, j_r in {1,2,3}:
        # Need l_r | j_r and l_s | j_s
        # l_r = r/gcd(t,r), l_s = s/gcd(t,s)

        cross_check = 0
        for t in range(L):
            gr = gcd(t, r) if t > 0 else r
            gs = gcd(t, s) if t > 0 else s
            lr = r // gr
            ls = s // gs

            for jr in range(1, 4):
                js = 4 - jr
                if jr % lr != 0 or js % ls != 0:
                    continue
                cross_check += binomial(gr, jr // lr) * binomial(gs, js // ls)

        assert cross_check % L == 0
        cross_check //= L
        assert cross_check == cross

        # Now decompose by (jr, js) type
        c13 = 0  # jr=1, js=3
        c22 = 0  # jr=2, js=2
        c31 = 0  # jr=3, js=1
        for t in range(L):
            gr = gcd(t, r) if t > 0 else r
            gs = gcd(t, s) if t > 0 else s
            lr = r // gr
            ls = s // gs

            for jr, js in [(1, 3), (2, 2), (3, 1)]:
                if jr % lr != 0 or js % ls != 0:
                    continue
                val = binomial(gr, jr // lr) * binomial(gs, js // ls)
                if jr == 1 and js == 3:
                    c13 += val
                elif jr == 2 and js == 2:
                    c22 += val
                else:
                    c31 += val

        c13 //= L
        c22 //= L
        c31 //= L

        if cross != 0:
            print(f"  r={r} s={s}: cross={cross} = c13:{c13} + c22:{c22} + c31:{c31}, g={g}")

print("\nDone phase 1.")
