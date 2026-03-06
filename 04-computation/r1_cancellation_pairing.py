#!/usr/bin/env python3
"""
Investigate the PAIRING STRUCTURE in r^1 cancellation.

At n=4, each s_ij term in the r^1 coefficient cancels between
two subsets S that differ by adding/removing one vertex.

Can this pairing be made into a general involution that proves
r^1 = 0 for all n? And does it extend to all odd r-powers?

opus-2026-03-06-S21
"""

from itertools import permutations, combinations
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
print("r^1 CANCELLATION PAIRING ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: Enumerate ALL contributions to r^1 at n=4
# ============================================================
print("\n--- Part 1: n=4 detailed r^1 contributions ---")

n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3]
m = n - 2  # = 2

# Each cover: path P_a through S+a ending at a, path P_b through R+b starting at b
# Product of arc weights: (r+s_{e1})(r+s_{e2}) for m=2 arcs
# r^1 coefficient: s_{e1} + s_{e2} (= e_1 of skew arc values)

contributions = []  # (sign, S, P_a_arcs, P_b_arcs, r1_coeff)
for mask in range(1 << len(U)):
    S_l = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_l = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_l)
    S_set = set(S_l) | {a}
    R_set = set(R_l) | {b}

    for p1 in permutations(sorted(S_set)):
        if p1[-1] != a: continue
        for p2 in permutations(sorted(R_set)):
            if p2[0] != b: continue
            arcs_a = [(p1[i], p1[i+1]) for i in range(len(p1)-1)]
            arcs_b = [(p2[i], p2[i+1]) for i in range(len(p2)-1)]
            all_arcs = arcs_a + arcs_b
            s_vals = [s(e[0], e[1]) for e in all_arcs]

            # r^1 = sum of s values
            r1 = sum(s_vals)
            r1 = expand(r1)

            # For each ARC, its individual contribution to r^1
            # r^1 = sum_j prod_{i!=j} s_i evaluated as sum of s values when m=2
            # Actually r^1 coeff of (r+s1)(r+s2) = s1+s2 (for m=2)
            # For general m: r^1 coeff = e_{m-1}(s1,...,sm)

            contributions.append((sign, tuple(S_l), arcs_a, arcs_b, r1))
            print(f"  S={S_l}, P_a={'->'.join(str(x) for x in p1)}, P_b={'->'.join(str(x) for x in p2)}")
            print(f"    sign={sign:+d}, arcs={all_arcs}, r^1_contrib = {sign} * ({r1}) = {expand(sign * r1)}")

total = sum(expand(c[0] * c[4]) for c in contributions)
print(f"\n  TOTAL r^1 = {total}")

# ============================================================
# Part 2: Decompose by individual s_ij
# ============================================================
print("\n--- Part 2: Which covers contribute which s_ij? ---")

s_contributions = defaultdict(list)  # s_ij -> list of (sign, S, cover)
for sign, S_l, arcs_a, arcs_b, r1 in contributions:
    all_arcs = list(arcs_a) + list(arcs_b)
    for arc in all_arcs:
        s_val = s(arc[0], arc[1])
        # The r^1 contribution from this arc is:
        # sign * (product of OTHER arcs' s-values)
        # For m=2: the other arc's s-value
        other_arcs = [e for e in all_arcs if e != arc]
        other_product = 1
        for e in other_arcs:
            other_product *= s(e[0], e[1])
        other_product = expand(other_product)

        key = str(s_val)
        s_contributions[key].append({
            'sign': sign,
            'S': S_l,
            'arc': arc,
            'other_product': other_product,
            'total_contrib': expand(sign * other_product)
        })

print(f"\n  Contributions by arc s-value:")
for key in sorted(s_contributions.keys()):
    entries = s_contributions[key]
    total = sum(e['total_contrib'] for e in entries)
    print(f"    {key}:")
    for e in entries:
        print(f"      S={e['S']}, arc={e['arc']}, sign={e['sign']:+d}, "
              f"other={e['other_product']}, contrib={e['total_contrib']}")
    print(f"      NET = {expand(total)}")

# ============================================================
# Part 3: n=5 — same analysis
# ============================================================
print("\n--- Part 3: n=5 r^1 pairing ---")
n = 5
r, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3, 4]
m = 3  # arcs per cover

# For m=3, r^1 coeff of (r+s1)(r+s2)(r+s3) = e_2(s1,s2,s3) = s1*s2 + s1*s3 + s2*s3

total_r1 = 0
per_arc_contribs = defaultdict(lambda: 0)
per_S_r1 = {}

for mask in range(1 << len(U)):
    S_l = tuple(U[i] for i in range(len(U)) if mask & (1 << i))
    R_l = tuple(U[i] for i in range(len(U)) if not (mask & (1 << i)))
    sign = (-1)**len(S_l)
    S_set = set(S_l) | {a}
    R_set = set(R_l) | {b}

    S_r1 = 0
    for p1 in permutations(sorted(S_set)):
        if p1[-1] != a: continue
        for p2 in permutations(sorted(R_set)):
            if p2[0] != b: continue
            arcs_a = [(p1[i], p1[i+1]) for i in range(len(p1)-1)]
            arcs_b = [(p2[i], p2[i+1]) for i in range(len(p2)-1)]
            all_arcs = arcs_a + arcs_b
            s_vals = [s(e[0], e[1]) for e in all_arcs]

            # r^1 coeff of product = e_{m-1}(s_vals)
            # For m=3: e_2 = s1*s2 + s1*s3 + s2*s3
            r1_coeff = 0
            for i in range(m):
                for j in range(i+1, m):
                    r1_coeff += s_vals[i] * s_vals[j]
            r1_coeff = expand(r1_coeff)
            S_r1 += expand(sign * r1_coeff)

    per_S_r1[S_l] = expand(S_r1)
    total_r1 += S_r1

total_r1 = expand(total_r1)
print(f"  r^1 contributions by S:")
for S_l in sorted(per_S_r1.keys(), key=len):
    print(f"    S={list(S_l)}: r^1 = {per_S_r1[S_l]}")
print(f"  TOTAL r^1 = {total_r1}")

# ============================================================
# Part 4: Complement pairing for r^1
# ============================================================
print("\n--- Part 4: Complement pairing S <-> U\\S for r^1 ---")
for S_l in sorted(per_S_r1.keys(), key=len):
    R_l = tuple(u for u in U if u not in S_l)
    if S_l <= R_l:  # avoid double-counting
        pair_sum = expand(per_S_r1[S_l] + per_S_r1[R_l])
        print(f"    S={list(S_l)} + S={list(R_l)}: pair sum r^1 = {pair_sum}")

# ============================================================
# Part 5: Single-element toggle pairing for r^1
# ============================================================
print("\n--- Part 5: Toggle pairing S <-> S xor {u} for r^1 ---")
for u in U:
    print(f"\n  Toggling vertex {u}:")
    for S_l in sorted(per_S_r1.keys(), key=len):
        if u in S_l:
            S_without = tuple(v for v in S_l if v != u)
        else:
            S_without = S_l
            S_l_toggled = tuple(sorted(list(S_l) + [u]))
            pair = expand(per_S_r1[S_l] + per_S_r1.get(S_l_toggled, 0))
            print(f"    S={list(S_l)} + S={list(S_l_toggled)}: {pair}")

# ============================================================
# Part 6: r^3 at n=5 (another odd power)
# ============================================================
print("\n--- Part 6: r^3 cancellation at n=5 ---")

per_S_r3 = {}
total_r3 = 0
for mask in range(1 << len(U)):
    S_l = tuple(U[i] for i in range(len(U)) if mask & (1 << i))
    R_l = tuple(U[i] for i in range(len(U)) if not (mask & (1 << i)))
    sign = (-1)**len(S_l)
    S_set = set(S_l) | {a}
    R_set = set(R_l) | {b}

    S_r3 = 0
    for p1 in permutations(sorted(S_set)):
        if p1[-1] != a: continue
        for p2 in permutations(sorted(R_set)):
            if p2[0] != b: continue
            arcs_a = [(p1[i], p1[i+1]) for i in range(len(p1)-1)]
            arcs_b = [(p2[i], p2[i+1]) for i in range(len(p2)-1)]
            all_arcs = arcs_a + arcs_b
            s_vals = [s(e[0], e[1]) for e in all_arcs]

            # r^3 coeff for m=3: e_0 = 1 (product of zero s-values)
            # Wait: product (r+s_i) = r^3 + r^2*e_1 + r*e_2 + e_3
            # So r^3 coeff = 1 (coefficient of r^3 is 1 for each cover)
            r3_coeff = 1
            S_r3 += sign * r3_coeff

    per_S_r3[S_l] = S_r3
    total_r3 += S_r3

print(f"  r^3 contributions by S (each is sign * count_of_covers):")
for S_l in sorted(per_S_r3.keys(), key=len):
    print(f"    S={list(S_l)}: {per_S_r3[S_l]}")
print(f"  TOTAL r^3 = {total_r3}")

# r^3 coefficient = sum_S (-1)^|S| * (number of covers for S)
# = sum_S (-1)^|S| * |S|! * |R|!

print(f"\n  Check: sum (-1)^|S| |S|! |R|! = ?")
check = 0
for mask in range(1 << len(U)):
    S_l = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_l = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    k = len(S_l)
    import math
    term = (-1)**k * math.factorial(k) * math.factorial(len(U) - k)
    check += term
    print(f"    |S|={k}: (-1)^{k} * {k}! * {len(U)-k}! = {term}")
print(f"  Sum = {check}")
print(f"  = (n-2)! * sum_k (-1)^k = {math.factorial(len(U))} * "
      f"{'0' if len(U) % 2 == 1 else '1'}")

# ============================================================
# Part 7: The KEY QUESTION — does the toggle pairing work in general?
# ============================================================
print("\n" + "=" * 70)
print("Part 7: GENERAL TOGGLE INVOLUTION")
print("=" * 70)
print("""
IDEA: For a fixed vertex u in U, define the involution:
  I_u: (S, cover C) -> (S xor {u}, cover C')

where C' is obtained by "transferring vertex u" between the two paths.

When u is in S (in P_a's territory), move it to R (to P_b's territory).
This changes (-1)^|S| to (-1)^{|S|-1} = -(-1)^|S|.

For this to cancel the r^1 contribution, we need:
  r^1(C') = r^1(C)   (same r^1 magnitude, opposite sign from (-1)^|S|)

But C and C' have DIFFERENT arc sets! So the r^1 contributions are
different polynomials in s. For the cancellation, we need the SUM
over all covers with u "toggled" to have the same r^1 magnitude.

This is NOT a simple involution on individual covers, but rather
a BULK cancellation across all covers.
""")

# Let me verify: for each vertex u, does sum over S containing u of
# (-1)^|S| * (r^1 for S) = -sum over S not containing u of (-1)^|S| * (r^1 for S)?

for u in U:
    with_u = 0
    without_u = 0
    for S_l, val in per_S_r1.items():
        if u in S_l:
            with_u += val
        else:
            without_u += val
    with_u = expand(with_u)
    without_u = expand(without_u)
    print(f"  u={u}: sum(S contains u) = {with_u}, sum(S not contains u) = {without_u}")
    print(f"    Equal in magnitude? {expand(with_u + without_u) == 0}")

# ============================================================
# Part 8: The DERIVATIVE approach
# ============================================================
print("\n" + "=" * 70)
print("Part 8: DERIVATIVE APPROACH")
print("=" * 70)
print("""
r^1 coeff of M = dM/dr at r=0.

M(r) = sum_S (-1)^|S| E_a(S+a; r) * B_b(R+b; r)
M'(r) = sum_S (-1)^|S| [E_a'(r) B_b(r) + E_a(r) B_b'(r)]
M'(0) = sum_S (-1)^|S| [E_a'(0) B_b(0) + E_a(0) B_b'(0)]

Now use path reversal at c=0:
  B_b(R+b; 0) = (-1)^|R| E_b(R+b; 0)

  B_b'(r) = d/dr B_b(R+b; r)
  At r=0: what is B_b'(0)?

  B_b(R+b; r) = sum_paths prod(r + s_e)
  B_b(R+b; r) at r=0 uses path reversal: B_b(0) = (-1)^|R| E_b(0)
  d/dr B_b(r) = sum_paths sum_j prod_{i!=j} (r + s_{e_i})
  At r=0: B_b'(0) = sum_paths sum_j prod_{i!=j} s_{e_i}

  Hmm, can we relate B_b'(0) to E_b'(0)?

  Actually: B_b(r) = E_b(r at -s) * (-1)^|R| ... no.
  B_b(R+b; r,s) = E_b(R+b; r,-s) ... wait, the path reversal is:
  B_v(S+v; c,s) = E_v(S+v; c,-s)

  So d/dr B_v(r,s) = d/dr E_v(r,-s).
  At r=0: B_v'(0,s) = E_v'(0,-s).

  Therefore:
  M'(0) = sum_S (-1)^|S| [E_a'(0,s) * (-1)^|R| E_b(0,s) +
                            E_a(0,s) * E_b'(0,-s) * ... ]

  Hmm, this is getting messy. But the key: E_v'(0,-s) is E_v' with s negated.
  And E_v(0) is a homogeneous polynomial of degree |S| in s.
  E_v'(0) is a sum of degree-(|S|-1) terms.

  For M'(0) to vanish, we need cancellation across the S-sum.
""")

# Compute E_a'(0) and E_a(0) for each S at n=5
n = 5
r, sv, s, t = setup(n)
a, b = 0, 1
U = [2, 3, 4]

def E_val(v, S_set, r_val):
    """E_v(S+v) at given r."""
    vset = S_set | {v}
    total = 0
    for p in permutations(sorted(vset)):
        if p[-1] != v: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= (r_val + s(p[i], p[i+1]))
        total += expand(prod)
    return expand(total)

subs_neg = {sv[k]: -sv[k] for k in sv}

print(f"\n  n=5: E_a'(0), B_b(0), E_a(0), B_b'(0) for each S:")
M_prime_0 = 0
for mask in range(1 << len(U)):
    S_l = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_l = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_l)

    Ea_r = E_val(a, set(S_l), r)
    Bb_r = E_val(b, set(R_l), r)  # B_b(R+b;r) as a function

    # Actually B_b = paths starting at b, not ending at b.
    # Let me compute B directly.
    vset_R = set(R_l) | {b}
    Bb_total = 0
    for p in permutations(sorted(vset_R)):
        if p[0] != b: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= (r + s(p[i], p[i+1]))
        Bb_total += expand(prod)

    Ea_r_poly = Poly(Ea_r, r)
    Bb_r_poly = Poly(expand(Bb_total), r)

    Ea_0 = expand(Ea_r_poly.nth(0)) if Ea_r_poly.degree() >= 0 else 0
    Ea_1 = expand(Ea_r_poly.nth(1)) if Ea_r_poly.degree() >= 1 else 0  # coefficient of r
    Bb_0 = expand(Bb_r_poly.nth(0)) if Bb_r_poly.degree() >= 0 else 0
    Bb_1 = expand(Bb_r_poly.nth(1)) if Bb_r_poly.degree() >= 1 else 0

    contrib = expand(sign * (Ea_1 * Bb_0 + Ea_0 * Bb_1))
    M_prime_0 += contrib
    if len(S_l) <= 1 or len(S_l) >= len(U) - 1:  # only print small/large S
        print(f"    S={S_l}: sign*[E_a'(0)*B_b(0) + E_a(0)*B_b'(0)] = {contrib}")

M_prime_0 = expand(M_prime_0)
print(f"\n  M'(0) = {M_prime_0}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
The r^1 cancellation at n=4 works by a TERM-BY-TERM pairing where
each s_ij appears exactly twice with opposite signs in the r^1 sum.

At n=5, the per-S r^1 contributions DON'T pair individually by
toggling a single vertex. The cancellation is GLOBAL (across all
subsets simultaneously).

The toggle pairing (S contains u vs S doesn't contain u) splits
the sum into two halves that cancel. This works for EACH choice of u.
This is related to the fact that sum (-1)^|S| * f(S) = 0 when
f is "balanced" with respect to each element.

KEY INSIGHT: M'(0) = 0 is a consequence of the MORE GENERAL fact
that M(r) is even in r. But proving M'(0) = 0 directly (without
knowing M is even) might be possible through the specific structure
of E and B at c=0.
""")
