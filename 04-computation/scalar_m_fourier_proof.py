#!/usr/bin/env python3
"""
FOURIER PROOF ATTEMPT: M[0,k] = 0 for circulant tournaments at odd n.

KEY INSIGHT from rep_theory.py:
  For Paley T_5, E_a(S+a) summed over all S of size s is CONSTANT across a.
  This means: sum_{S subset [n]\{a}, |S|=s} E_a(S+{a}) = e_s (independent of a).

  This is a MUCH STRONGER statement than just M being circulant.
  If e_s is the same for all a, then the Fourier transform of the function
  a -> sum_{|S|=s} E_a(S+a) is concentrated at mode 0.

  QUESTION: Does this hold for ALL circulant tournaments? All VT?

  If so, then M[0,k] = sum_s sum_{S,R: |S|=s, R=U\S} (-1)^s E_0(S+0) B_k(R+k)
  has a special structure: the E and B components are "spectrally flat."

  But the PRODUCT E_0 * B_k is not simply the product of their marginals.
  The correlation between E_0 and B_k on complementary subsets matters.

APPROACH: Check if the "uniform E/B" property holds for non-Paley circulants,
          then attempt a combinatorial proof.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np
from collections import defaultdict

def make_circulant_tournament(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def E_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def B_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# TEST: Uniform E property for all circulant tournaments
# ============================================================
print("=" * 70)
print("TEST: sum_{|S|=s} E_a(S+a) is constant in a?")
print("=" * 70)

for n in [5, 7]:
    print(f"\n  n={n}:")
    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S_gen in gen_sets[:4]:  # Limit at n=7
        T = make_circulant_tournament(n, S_gen)
        H = count_H(T, n)

        is_uniform = True
        e_values = {}
        for s in range(n):
            e_by_a = []
            for a in range(n):
                others = [v for v in range(n) if v != a]
                e_sum = 0
                for S_tuple2 in combinations(others, s):
                    S_set = sorted(list(S_tuple2) + [a])
                    e_sum += E_v(T, S_set, a)
                e_by_a.append(e_sum)

            if len(set(e_by_a)) > 1:
                is_uniform = False
                print(f"    S={sorted(S_gen)}, |S|={s}: NOT UNIFORM — {e_by_a}")
            else:
                e_values[s] = e_by_a[0]

        if is_uniform:
            print(f"    S={sorted(S_gen)}: H={H}, UNIFORM! e_s = {[e_values[s] for s in range(n)]}")
        else:
            print(f"    S={sorted(S_gen)}: H={H}, NOT uniform")


# ============================================================
# TEST: Uniform E for NON-circulant VT at n=5
# ============================================================
print("\n" + "=" * 70)
print("TEST: Uniform E for non-circulant H=15 tournaments at n=5")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

tested = 0
uniform_count = 0
nonuniform_count = 0

for bits in range(1 << len(pairs)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    H = count_H(T, n)
    if H != 15:
        continue

    tested += 1
    is_uniform = True
    for s in range(n):
        e_by_a = []
        for a in range(n):
            others = [v for v in range(n) if v != a]
            e_sum = 0
            for S_tuple in combinations(others, s):
                S_set = sorted(list(S_tuple) + [a])
                e_sum += E_v(T, S_set, a)
            e_by_a.append(e_sum)
        if len(set(e_by_a)) > 1:
            is_uniform = False
            break

    if is_uniform:
        uniform_count += 1
    else:
        nonuniform_count += 1
        if nonuniform_count <= 3:
            scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
            print(f"  bits={bits}: scores={scores}, NOT uniform")

print(f"\n  H=15 tournaments: {uniform_count} uniform, {nonuniform_count} non-uniform out of {tested}")


# ============================================================
# DEEPER: What property ensures uniform E?
# ============================================================
print("\n" + "=" * 70)
print("What property ensures uniform E_a sums?")
print("=" * 70)

# Check for ALL H values at n=5
n = 5
h_to_uniform = defaultdict(lambda: [0, 0])  # [uniform_count, non_uniform_count]

for bits in range(1 << len(pairs)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    H = count_H(T, n)

    is_uniform = True
    for s in range(n):
        e_by_a = []
        for a in range(n):
            others = [v for v in range(n) if v != a]
            e_sum = 0
            for S_tuple in combinations(others, s):
                S_set = sorted(list(S_tuple) + [a])
                e_sum += E_v(T, S_set, a)
            e_by_a.append(e_sum)
        if len(set(e_by_a)) > 1:
            is_uniform = False
            break

    if is_uniform:
        h_to_uniform[H][0] += 1
    else:
        h_to_uniform[H][1] += 1

print(f"\n  n=5: Uniform E property by H value:")
for H in sorted(h_to_uniform.keys()):
    u, nu = h_to_uniform[H]
    print(f"    H={H:>3}: {u:>3} uniform, {nu:>3} non-uniform")


# ============================================================
# Is "uniform E at all levels" <==> M scalar?
# ============================================================
print("\n" + "=" * 70)
print("uniform E <==> M scalar?")
print("=" * 70)

# Check a few non-H=15 tournaments
def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            val += sign * E_v(T, S_set, a) * B_v(T, R_set, b)
        return val

n = 5
scalar_but_not_uniform = 0
uniform_but_not_scalar = 0
both = 0
neither = 0

for bits in range(1 << len(pairs)):
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    # Check uniform E
    is_uniform = True
    for s in range(n):
        e_by_a = []
        for a in range(n):
            others = [v for v in range(n) if v != a]
            e_sum = 0
            for S_tuple in combinations(others, s):
                S_set = sorted(list(S_tuple) + [a])
                e_sum += E_v(T, S_set, a)
            e_by_a.append(e_sum)
        if len(set(e_by_a)) > 1:
            is_uniform = False
            break

    # Check scalar M
    is_scalar = True
    for a in range(n):
        for b in range(a+1, n):
            if compute_M_entry(T, n, a, b) != 0:
                is_scalar = False
                break
        if not is_scalar:
            break

    if is_scalar and is_uniform:
        both += 1
    elif is_scalar and not is_uniform:
        scalar_but_not_uniform += 1
    elif not is_scalar and is_uniform:
        uniform_but_not_scalar += 1
    else:
        neither += 1

print(f"\n  n=5 results:")
print(f"    Both scalar M and uniform E: {both}")
print(f"    Scalar M but NOT uniform E:  {scalar_but_not_uniform}")
print(f"    Uniform E but NOT scalar M:  {uniform_but_not_scalar}")
print(f"    Neither:                      {neither}")

if scalar_but_not_uniform == 0 and uniform_but_not_scalar == 0:
    print(f"\n  PERFECT EQUIVALENCE: M scalar <==> uniform E at n=5!")
elif uniform_but_not_scalar == 0:
    print(f"\n  uniform E => M scalar (but not converse)")
elif scalar_but_not_uniform == 0:
    print(f"\n  M scalar => uniform E (but not converse)")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
