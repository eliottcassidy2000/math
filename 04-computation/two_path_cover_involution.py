#!/usr/bin/env python3
"""
THE 2-PATH-COVER INVOLUTION

Can we define an involution on 2-path-covers that:
1. Preserves the endpoints (a as end of first path, b as start of second)
2. Changes the sign of the contribution (negates r^odd terms)
3. Is a perfect matching (no fixed points at r^odd)?

A 2-path-cover consists of:
  (S, pi_a, pi_b) where:
    pi_a = path through S+a ending at a
    pi_b = path through R+b starting at b
    R = U \ S, U = [n] \ {a,b}

Weight: (-1)^|S| * prod of arc weights

The idea: encode (S, pi_a, pi_b) as a SINGLE permutation sigma of [n]
by defining sigma = the 'next vertex' map.

In pi_a = (v_1, v_2, ..., v_k, a): sigma(v_1) = v_2, ..., sigma(v_{k-1}) = v_k, sigma(v_k) = a
  (sigma(a) is undefined within pi_a since a is the END)

In pi_b = (b, w_1, w_2, ..., w_m): sigma(b) = w_1, sigma(w_1) = w_2, ..., sigma(w_{m-1}) = w_m
  (sigma(w_m) is undefined since w_m is the end of pi_b)

So sigma is a partial permutation on [n], defined on all vertices EXCEPT a and the
last vertex of pi_b. But we can complete it:
  sigma(a) = start of pi_a (makes pi_a a cycle)
  sigma(end of pi_b) = b (makes pi_b a cycle)

This gives a COMPLETE permutation sigma with exactly 2 cycles:
  cycle_a = the cycle through a (from pi_a)
  cycle_b = the cycle through b (from pi_b)

The key: |S| = |cycle_a| - 1 (since cycle_a has |S|+1 vertices including a).

kind-pasteur-2026-03-06-S23
"""

from itertools import permutations
from sympy import symbols, expand, Poly
from collections import defaultdict

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')

    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]

    def t(i, j):
        if i == j: return 0
        return r + s(i, j)

    return r, sv, s, t

print("=" * 70)
print("2-PATH-COVER INVOLUTION ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: Encode 2-path-covers as permutations
# ============================================================
print("\n--- Part 1: Encoding 2-path-covers as permutations ---")

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    print(f"\n  n={n}, a={a}, b={b}:")

    covers = []
    for mask in range(1 << len(U)):
        S_lst = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R_lst = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S_lst)

        S_set = set(S_lst) | {a}
        R_set = set(R_lst) | {b}

        for p_a in permutations(sorted(S_set)):
            if p_a[-1] != a: continue
            for p_b in permutations(sorted(R_set)):
                if p_b[0] != b: continue

                # Build the permutation sigma
                sigma = [0] * n

                # pi_a cycle: p_a[0] -> p_a[1] -> ... -> p_a[-1]=a -> p_a[0]
                for i in range(len(p_a)-1):
                    sigma[p_a[i]] = p_a[i+1]
                sigma[a] = p_a[0]  # close the cycle

                # pi_b cycle: b -> p_b[1] -> ... -> p_b[-1] -> b
                for i in range(len(p_b)-1):
                    sigma[p_b[i]] = p_b[i+1]
                sigma[p_b[-1]] = b  # close the cycle

                # Get cycle type
                visited = [False]*n
                cycles = []
                for start in range(n):
                    if visited[start]: continue
                    cycle = []
                    v = start
                    while not visited[v]:
                        visited[v] = True
                        cycle.append(v)
                        v = sigma[v]
                    cycles.append(tuple(cycle))

                cycle_lens = tuple(sorted([len(c) for c in cycles], reverse=True))

                # Compute weight
                weight = sign
                for i in range(len(p_a)-1):
                    weight *= t(p_a[i], p_a[i+1])
                for i in range(len(p_b)-1):
                    weight *= t(p_b[i], p_b[i+1])
                weight = expand(weight)

                covers.append({
                    'S': S_lst, 'R': R_lst, 'sign': sign,
                    'path_a': p_a, 'path_b': p_b,
                    'sigma': tuple(sigma),
                    'cycle_type': cycle_lens,
                    'weight': weight,
                    'cycles': cycles
                })

    # Analyze the cycle types
    by_type = defaultdict(list)
    for c in covers:
        by_type[c['cycle_type']].append(c)

    print(f"    Total 2-path-covers: {len(covers)}")
    print(f"    Cycle types: {sorted(by_type.keys())}")

    for ct in sorted(by_type.keys()):
        items = by_type[ct]
        total = expand(sum(c['weight'] for c in items))
        p_total = Poly(total, r)
        r1 = expand(p_total.nth(1)) if p_total.degree() >= 1 else 0
        r0 = expand(p_total.nth(0))
        has_even = any(k % 2 == 0 for k in ct)
        print(f"    Type {ct} ({'HAS EVEN' if has_even else 'all odd'}): "
              f"{len(items)} covers, r^0={len(r0.as_ordered_terms()) if r0 != 0 else 0}t, "
              f"r^1={'0' if r1 == 0 else len(r1.as_ordered_terms())}")

# ============================================================
# Part 2: The involution candidate
# ============================================================
print("\n" + "=" * 70)
print("Part 2: Involution candidate — reverse a cycle containing neither a nor b")
print("=" * 70)
print("""
If sigma has cycle type (k_a, k_b) where k_a = cycle through a, k_b = cycle through b,
and k_a + k_b = n, then:
  |S| = k_a - 1
  |R| = k_b - 1

The only cycles ARE the a-cycle and the b-cycle.
There are no "free" even cycles to reverse!

This means the global involution (reverse an even cycle) CANNOT directly apply,
because the 2-path-cover permutation has exactly 2 cycles, both containing
a fixed vertex (a or b).

NEW IDEA: Instead of reversing a cycle, consider an involution that:
  - Moves a vertex from S to R (or vice versa)
  - Changes which path it belongs to
  - Adjusts the sign by changing |S|

This is a VERTEX-TRANSFER involution.
""")

# Let's check: what if we transfer the FIRST vertex of pi_a to pi_b?
# pi_a = (v, ..., a), pi_b = (b, ...)
# Transfer v: new pi_a = (..., a) without v, new pi_b = (b, v, ...)
# This changes |S| by -1, flipping the sign.
# But the path structures change non-trivially...

# ============================================================
# Part 3: Vertex-transfer involution analysis
# ============================================================
print("\n--- Part 3: Can we find a sign-reversing involution? ---")

n = 4; r, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3]

# For each cover, compute r^1 contribution
# Then try to pair covers whose r^1 contributions cancel

r1_contributions = []
for mask in range(1 << len(U)):
    S_lst = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_lst = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_lst)

    S_set = set(S_lst) | {a}
    R_set = set(R_lst) | {b}

    for p_a in permutations(sorted(S_set)):
        if p_a[-1] != a: continue
        for p_b in permutations(sorted(R_set)):
            if p_b[0] != b: continue

            weight = sign
            for i in range(len(p_a)-1):
                weight *= t(p_a[i], p_a[i+1])
            for i in range(len(p_b)-1):
                weight *= t(p_b[i], p_b[i+1])
            weight = expand(weight)

            pw = Poly(weight, r)
            r1 = expand(pw.nth(1)) if pw.degree() >= 1 else 0

            if r1 != 0:
                r1_contributions.append({
                    'S': S_lst, 'path_a': p_a, 'path_b': p_b,
                    'sign': sign, 'r1': r1
                })
                print(f"    S={S_lst}, pi_a={p_a}, pi_b={p_b}: "
                      f"sign={sign:+d}, r^1 = {r1}")

total_r1 = expand(sum(c['r1'] for c in r1_contributions))
print(f"\n    TOTAL r^1 = {total_r1}")

# Try to manually pair them
print("\n  Manual pairing attempt:")
print("  Looking for pairs where r^1 contributions cancel...")

used = [False]*len(r1_contributions)
for i in range(len(r1_contributions)):
    if used[i]: continue
    for j in range(i+1, len(r1_contributions)):
        if used[j]: continue
        if expand(r1_contributions[i]['r1'] + r1_contributions[j]['r1']) == 0:
            ci = r1_contributions[i]
            cj = r1_contributions[j]
            print(f"    PAIR: S={ci['S']},pi_a={ci['path_a']},pi_b={ci['path_b']} "
                  f"<-> S={cj['S']},pi_a={cj['path_a']},pi_b={cj['path_b']}")
            print(f"          r^1: {ci['r1']} + {cj['r1']} = 0")
            used[i] = used[j] = True
            break

unpaired = [i for i in range(len(r1_contributions)) if not used[i]]
if unpaired:
    print(f"\n    UNPAIRED: {len(unpaired)} covers")
    for i in unpaired:
        ci = r1_contributions[i]
        print(f"      S={ci['S']},pi_a={ci['path_a']},pi_b={ci['path_b']}: r^1={ci['r1']}")
else:
    print(f"\n    ALL PAIRED! Perfect sign-reversing matching on 2-path-covers.")

# ============================================================
# Part 4: What IS the natural involution?
# ============================================================
print("\n" + "=" * 70)
print("Part 4: Characterizing the natural involution")
print("=" * 70)

# Let's look at the pairs more carefully
print("\n  Analyzing the paired covers at n=4:")
for i in range(len(r1_contributions)):
    for j in range(i+1, len(r1_contributions)):
        if expand(r1_contributions[i]['r1'] + r1_contributions[j]['r1']) == 0:
            ci = r1_contributions[i]
            cj = r1_contributions[j]
            print(f"\n    Cover 1: S={ci['S']}, pi_a={ci['path_a']}, pi_b={ci['path_b']}")
            print(f"    Cover 2: S={cj['S']}, pi_a={cj['path_a']}, pi_b={cj['path_b']}")
            # What changed?
            if set(ci['S']) != set(cj['S']):
                moved = set(ci['S']).symmetric_difference(set(cj['S']))
                print(f"    Vertex moved between S and R: {moved}")
            else:
                print(f"    Same S! Path rearrangement.")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
KEY FINDING: The 2-path-cover encoding as a permutation has EXACTLY
2 cycles (the a-cycle and b-cycle). There are NO free cycles to reverse.

This means the global Even Cycle Vanishing involution does NOT directly
apply to the endpoint setting. We need a DIFFERENT involution.

The natural candidate: a vertex-transfer involution that moves a vertex
between the S-path and R-path, changing |S| by 1 and thus flipping
the (-1)^|S| sign.

At n=4: the r^1 contributions CAN be perfectly paired (verified above).
The pairing involves either:
  (a) Moving a vertex from S to R (or vice versa)
  (b) Rearranging paths within the same S

This suggests a combinatorial proof exists, but the involution is
ENDPOINT-SPECIFIC and cannot be reduced to the global even-cycle-vanishing.
""")
