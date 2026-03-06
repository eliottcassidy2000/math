#!/usr/bin/env python3
"""
STRUCTURE OF M_n - M_{n-1}: WHY DOES IT HAVE EVEN r-POWERS?

We decompose M_n[a,b] - M_{n-1}[a,b] (deleting vertex v = n-1) into
terms where v contributes 1 or 2 arcs, and analyze each.

KEY DECOMPOSITION:
  M_n = sum_{S subset U} (-1)^|S| E_a(S+a) B_b(R+b)

Split into v in S and v not in S:
  Part_R = sum_{S' subset U'} (-1)^|S'| E_a(S'+a) B_b(R'+v+b)
  Part_S = sum_{S' subset U'} -(-1)^|S'| E_a(S'+v+a) B_b(R'+b)

  M_n = Part_R + Part_S
  M_{n-1} = sum_{S' subset U'} (-1)^|S'| E_a(S'+a) B_b(R'+b)

  M_n - M_{n-1} = sum_{S'} (-1)^|S'| E_a(S'+a) [B_b(R'+v+b) - B_b(R'+b)]
                 - sum_{S'} (-1)^|S'| E_a(S'+v+a) B_b(R'+b)

The first sum: v is added to the b-path (enlarging it by v).
The second sum: v joins the a-path (with sign flip).

Now: B_b(R'+v+b) - B_b(R'+b) is the "marginal value" of adding v to the b-path.
Can we express this as t(?,v)*B_b(something) + B_b(something)*t(v,?) ?

Actually, B_b(W) for the b-path starting at b counts Hamiltonian paths
b -> ... through W. Adding v to W:
B_b(W+v) = sum over positions of v in the path * contribution.
Specifically: v can be inserted at any position in the path.

Similarly, E_a(V+v) where v is added to the a-path ending at a.

QUESTION: Is there a "transfer matrix" formulation where adding v
produces a product that's manifestly even in r?

kind-pasteur-2026-03-06-S23b
"""
from itertools import permutations
from sympy import symbols, expand, Poly, factor
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

def hp_sum(t_fn, vset, start=None, end=None):
    """Sum of Hamiltonian path weights through vset."""
    vl = sorted(vset)
    if len(vl) == 0:
        return 0
    if len(vl) == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for p in permutations(vl):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= t_fn(p[i], p[i+1])
        total += prod
    return expand(total)

print("=" * 70)
print("STRUCTURE OF M_n - M_{n-1}")
print("=" * 70)

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    v = n - 1
    U = [u for u in range(n) if u != a and u != b]
    Uprime = [u for u in U if u != v]

    print(f"\nn={n}, v={v}, U={U}, U'={Uprime}")

    # Compute the two parts of the difference
    # Part 1: v added to b-path
    # sum_{S' subset U'} (-1)^|S'| E_a(S'+a) [B_b(R'+v+b) - B_b(R'+b)]
    part1 = 0
    # Part 2: v joins a-path
    # -sum_{S' subset U'} (-1)^|S'| E_a(S'+v+a) B_b(R'+b)
    part2 = 0

    for mask in range(1 << len(Uprime)):
        Sp = [Uprime[i] for i in range(len(Uprime)) if mask & (1 << i)]
        Rp = [Uprime[i] for i in range(len(Uprime)) if not (mask & (1 << i))]
        sign = (-1)**len(Sp)

        Ea = hp_sum(t, set(Sp) | {a}, end=a)
        Bb_with_v = hp_sum(t, set(Rp) | {v, b}, start=b)
        Bb_without_v = hp_sum(t, set(Rp) | {b}, start=b)
        Ea_with_v = hp_sum(t, set(Sp) | {v, a}, end=a)
        Bb_basic = hp_sum(t, set(Rp) | {b}, start=b)

        part1 += sign * Ea * (Bb_with_v - Bb_without_v)
        part2 += -sign * Ea_with_v * Bb_basic

    part1 = expand(part1)
    part2 = expand(part2)
    diff = expand(part1 + part2)

    print(f"\n  Part 1 (v added to b-path): {part1}")
    print(f"  Part 2 (v joins a-path):     {part2}")
    print(f"  Total difference:             {diff}")

    # Check r-parity of each part
    for label, expr in [("Part1", part1), ("Part2", part2), ("Diff", diff)]:
        if expr != 0:
            p = Poly(expr, r)
            odd_coeffs = {}
            even_coeffs = {}
            for k in range(p.degree() + 1):
                c = expand(p.nth(k))
                if c != 0:
                    if k % 2 == 1:
                        odd_coeffs[k] = c
                    else:
                        even_coeffs[k] = c
            print(f"\n  {label}: degree {p.degree()}")
            for k in sorted(even_coeffs):
                print(f"    r^{k} [even]: {even_coeffs[k]}")
            for k in sorted(odd_coeffs):
                print(f"    r^{k} [odd]:  {odd_coeffs[k]}")

    # KEY CHECK: does Part1 at -r equal -Part2 at r (or some relation)?
    # M(-r) = M(b,a), so the substitution r -> -r transforms Part1 and Part2.
    print(f"\n  --- Part1(-r) vs Part2(r) relationship ---")
    part1_neg = expand(part1.subs(r, -r))
    print(f"  Part1(-r) = {part1_neg}")
    print(f"  Part2(r)  = {part2}")
    print(f"  Part1(-r) + Part2(r) = {expand(part1_neg + part2)}")
    print(f"  Part1(-r) - Part2(r) = {expand(part1_neg - part2)}")

    # What about Part1(-r) and Part2's "swapped" version?
    # Under r -> -r, the FULL M[a,b] becomes M[b,a].
    # Part1 (v in b-path for M[a,b]) under r -> -r should relate to
    # v in a-path for M[b,a] (which is Part2 for M[b,a]).

    # Let's compute Part1 and Part2 for M[b,a] and see the relationship
    part1_ba = 0
    part2_ba = 0

    for mask in range(1 << len(Uprime)):
        Sp = [Uprime[i] for i in range(len(Uprime)) if mask & (1 << i)]
        Rp = [Uprime[i] for i in range(len(Uprime)) if not (mask & (1 << i))]
        sign = (-1)**len(Sp)

        Eb = hp_sum(t, set(Sp) | {b}, end=b)
        Ba_with_v = hp_sum(t, set(Rp) | {v, a}, start=a)
        Ba_without_v = hp_sum(t, set(Rp) | {a}, start=a)
        Eb_with_v = hp_sum(t, set(Sp) | {v, b}, end=b)
        Ba_basic = hp_sum(t, set(Rp) | {a}, start=a)

        part1_ba += sign * Eb * (Ba_with_v - Ba_without_v)
        part2_ba += -sign * Eb_with_v * Ba_basic

    part1_ba = expand(part1_ba)
    part2_ba = expand(part2_ba)

    print(f"\n  --- M[b,a] parts ---")
    print(f"  Part1_ba (v in a-path of M[b,a]): {part1_ba}")
    print(f"  Part2_ba (v in b-path of M[b,a]): {part2_ba}")

    # Check: Part1_ab(-r) = Part2_ba(r) ?
    print(f"\n  Part1_ab(-r) = Part2_ba(r)? {expand(part1_neg - part2_ba) == 0}")
    part2_neg = expand(part2.subs(r, -r))
    print(f"  Part2_ab(-r) = Part1_ba(r)? {expand(part2_neg - part1_ba) == 0}")

    # If Part1_ab(-r) = Part2_ba(r), then the "v in b-path" at -r maps to
    # "v in a-path" of M[b,a] at r. And vice versa.
    # Then Diff_ab(-r) = Part1_ab(-r) + Part2_ab(-r) = Part2_ba(r) + Part1_ba(r) = Diff_ba(r).
    # So Diff_ab(-r) = Diff_ba(r), same as the full M identity.
    # Diff is even iff Diff_ab = Diff_ba.
    diff_ba = expand(part1_ba + part2_ba)
    print(f"\n  Diff_ab = Diff_ba? {expand(diff - diff_ba) == 0}")

print("\n" + "=" * 70)
print("Part 2: FACTORING the difference into even-r-powers")
print("=" * 70)

# At n=4: M_4 - M_3 should factor nicely
n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
v = 3

# M_3 on {0,1,2}
M3 = hp_sum(t, {0, 2}, end=0) * hp_sum(t, {1}, start=1) - hp_sum(t, {0,2}, end=0) * 1 + 1 * hp_sum(t, {1,2}, start=1)  # Wrong, let me recompute

# Actually: M_3[0,1] on vertices {0,1,2}
# U = {2}, S subset {2}
# S={}: E_0({0}) * B_1({1,2}) = 1 * t(1,2)
# S={2}: -E_0({0,2}) * B_1({1}) = -t(2,0) * 1
M3 = expand(1 * t(1,2) - t(2,0) * 1)
print(f"\nn=4: M_3[0,1] = {M3}")

# M_4 on {0,1,2,3}
Uprime = [2]
M4 = 0
for mask in range(1 << 2):  # U = {2, 3}
    S_list = [u for i, u in enumerate([2,3]) if mask & (1 << i)]
    R_list = [u for i, u in enumerate([2,3]) if not (mask & (1 << i))]
    sign = (-1)**len(S_list)
    ea = hp_sum(t, set(S_list) | {0}, end=0)
    bb = hp_sum(t, set(R_list) | {1}, start=1)
    M4 += sign * ea * bb
M4 = expand(M4)
print(f"M_4[0,1] = {M4}")

diff = expand(M4 - M3)
print(f"Diff = {diff}")

# Try to factor
print(f"Factor(Diff) = {factor(diff)}")

# Check: is diff = f(r^2) * g(s) + h(r^2) * k(s)?
p_diff = Poly(diff, r)
print(f"  r^0: {expand(p_diff.nth(0))}")
print(f"  r^2: {expand(p_diff.nth(2))}")

# The r^0 part of the difference
r0 = expand(p_diff.nth(0))
# Factor r0
print(f"  Factor r^0: {factor(r0)}")

# Can we write diff = r^2 * C + D where C, D are functions of s only?
# diff = 2*r^2 + (r^0 stuff)
# The r^0 stuff is: s02*s13 + s02*s23 - s02 + s03*s12 - s03*s23 + s12*s23 - s12 - s13*s23

# Actually: let me check if r^0 of diff factors as (s02+s12)(something) - (something)
# r^0 = s02*s13 + s02*s23 - s02 + s03*s12 - s03*s23 + s12*s23 - s12 - s13*s23
# = s02*(s13 + s23 - 1) + s12*(s03 + s23 - 1) - s23*(s03 + s13)
# = (s02 + s12)(s23 - 1) + s02*s13 + s03*s12 - s23*(s03 + s13) + s12*s23 + s02*s23
# Hmm, that doesn't simplify nicely.

# Let me try: r^0 with s03 = -s30, etc.
# Remember M_3 = s12 + s02. So M_3 - (something) might relate.
# r^0 of diff = r^0(M_4) - M_3 = (s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23) - (s12 + s02)
# = s02*s13 + s02*s23 - s02 + s03*s12 - s03*s23 + s12*s23 - s12 - s13*s23
# = s02*(s13 + s23 - 1) + s12*(s03 + s23 - 1) - s23*(s03 + s13)
# = (s02 + s12)*(s23 - 1) + s02*s13 + s03*s12 - s23*s03 - s23*s13 + s02*s23 + s12*s23
# Hmm... = (s02+s12)(s23-1) + s13(s02-s23) + s03(s12-s23) + s23(s02+s12)
# = (s02+s12)(s23-1+s23) + s13(s02-s23) + s03(s12-s23)
# = (s02+s12)(2*s23-1) + s13(s02-s23) + s03(s12-s23)
# Not obviously factorable.

# KEY INSIGHT CHECK: is r^0 of diff = M_3 * (s23 - 1) + something?
# M_3 * (s23 - 1) = (s02+s12)*(s23-1) = s02*s23 - s02 + s12*s23 - s12
# r^0 diff - M_3*(s23-1) = s02*s13 + s03*s12 - s03*s23 - s13*s23
# = s13*(s02-s23) + s03*(s12-s23)
# = s13*(-s(0,2)+s(0,2)) ... wait, s02 is just s02 and s23 is s23.
# Hmm, no factoring here.

print("\n" + "=" * 70)
print("Part 3: Alternative decomposition using v's CONTRIBUTION")
print("=" * 70)

# For n=5: the b-path includes v at some position.
# B_b(W+v) - B_b(W) = sum over paths through W+v starting at b that USE v
# minus the paths that DON'T use v. But all paths through W+v must use v.
# So B_b(W+v) - B_b(W) doesn't make sense for subsets: ALL paths through W+v use v.
#
# Wait: B_b(W) counts paths through ALL vertices in W. B_b(W+v) counts
# paths through ALL vertices in W+v. So B_b(W+v) - B_b(W) is NOT the
# "v-contribution" — it's the difference of two sums over different vertex sets.
#
# The right decomposition: in B_b(W+v), condition on v's position.
# v can be at position 2, 3, ..., |W+v| (since b is always at position 1).
# If v is at position k (1-indexed), then:
# - The first k-1 arcs form a path b -> ... -> prev(v) through a (k-1)-subset
# - Arc prev(v) -> v
# - The last |W+v|-k arcs form a path v -> ... through the remaining vertices
#
# So B_b(W+v) = sum_{u in W} t(u,v) * [paths b->...->u through subset] * [paths v->... through rest]
#             + sum_{w in W} [paths b->...->v->w->... with v->w arc]
# This is complicated.
#
# SIMPLER: condition on v's NEIGHBORS (predecessor and successor).
# v has predecessor u and successor w in the path (or v is at an end).
#
# For the b-path (starting at b):
# - v at the end: ...-> u -> v. Predecessor = u, no successor.
#   Contribution: B_b(W, end=u) * t(u,v)
# - v in the middle: ...-> u -> v -> w ->...
#   Contribution: sum over (u,w) of B_b(W1, start=b, end=u) * t(u,v) * t(v,w) * B_b(W2, start=w)
#   where W1 + W2 = W (partition of vertices other than v in the b-path)

# For n=4: U={2,3}, v=3.
n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
v = 3

print(f"\nn=4: Decomposing B_b(R+v+b) by v's neighbors")
# R' subset U' = {2}. Cases:
# R' = {}: b-path through {v, b} = {3, 1}. Path: 1->3.
#   v is at end. Contribution: t(1,3).
# R' = {2}: b-path through {v, 2, b} = {3, 2, 1}. Paths: 1->2->3 or 1->3->2.
#   For 1->2->3: v=3 at end, pred=2. Contribution: t(1,2)*t(2,3).
#   For 1->3->2: v=3 in middle, pred=1, succ=2. Contribution: t(1,3)*t(3,2).

# Let me compute t(u,v)*t(v,w) for specific u,w:
for u in range(n):
    for w in range(n):
        if u == v or w == v or u == w:
            continue
        prod = expand(t(u,v) * t(v,w))
        # t(u,v)*t(v,w) = (r+s_{uv})(r+s_{vw}) = r^2 + r(s_{uv}+s_{vw}) + s_{uv}*s_{vw}
        # Under r -> -r: (-r+s_{uv})(-r+s_{vw}) = r^2 - r(s_{uv}+s_{vw}) + s_{uv}*s_{vw}
        # = t(v,u)*t(w,v) at r  [since t(v,u) = r+s(v,u) = r-s(u,v)]
        # So: [t(u,v)*t(v,w)](-r) = t(v,u)*t(w,v)(r)
        # The ODD part is: r*(s_{uv}+s_{vw}) -> under r->-r, this flips sign.
        # The EVEN part is: r^2 + s_{uv}*s_{vw} -> unchanged.
        #
        # For the sum over all covers to have even r-powers, the r*(s_{uv}+s_{vw})
        # terms from v-in-middle must cancel the r*s_{*} terms from v-at-end.
        print(f"  t({u},{v})*t({v},{w}) = {prod}")

# For a single arc (endpoint): t(u,v) at -r = -(r-s(u,v)) = -t(v,u).
# Sign flip on the arc level means v-at-end contributions get an extra -.

print("\n  KEY: Under r -> -r, the arc products transform as:")
print("  t(u,v)*t(v,w)|_{-r} = t(v,u)*t(w,v)|_r")
print("  t(u,v)|_{-r} = -t(v,u)|_r")
print("  So pairs of arcs through v are EVEN in r (no sign),")
print("  while single arcs involving v are ODD (sign flip).")
print("  The cancellation of odd powers requires matching single-arc")
print("  and double-arc contributions involving v.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
