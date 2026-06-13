#!/usr/bin/env python3
"""
VERTEX TRANSFER PROOF ATTEMPT for even r-powers in M[a,b].

KEY INSIGHT (from wronskian_even_proof.py at n=4):
At [r^1], each s_{ij} appears exactly twice with opposite signs in the
Wronskian sum. The mechanism: moving vertex u between S and R flips (-1)^|S|
while preserving the s-product.

THIS SCRIPT: Formalizes and tests the vertex transfer at [r^1] and [r^3].

THE PROOF STRATEGY:
1. Express [r^k] G as a sum indexed by (S, path pair, edge subset T)
   where T is the set of edges contributing r (|T| = k) and the rest
   contribute s.
2. Show that for odd k, there exists a sign-reversing involution on
   these index tuples.

The involution: for a tuple (S, Pa, Pb, T), identify a "transfer vertex"
w that can be moved between S and R. This changes:
  - (-1)^|S| -> -(-1)^|S| (sign flip)
  - The path structure Pa, Pb changes (w moves to other path)
  - The r-contributing edge set T changes
  - But the s-PRODUCT remains the same

If for every tuple with odd |T|, we can find such a transfer, we're done.

kind-pasteur-2026-03-06-S24
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Symbol, Poly
from collections import defaultdict

def make_tournament(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    def t(i, j):
        if i == j: return 0
        return r + s[(i,j)]
    return r, s, t

def compute_G_terms(n, a, b, r, s):
    """
    Compute G(r) = sum_S E_a(S+a; r) * E_b(R+b; -r) term by term.

    Returns list of (S, path_a, path_b, weight) tuples.
    """
    U = [v for v in range(n) if v != a and v != b]
    terms = []

    for mask in range(2**len(U)):
        S = tuple(sorted([U[i] for i in range(len(U)) if (mask >> i) & 1]))
        R = tuple(sorted([u for u in U if u not in S]))
        Sa = set(S) | {a}
        Rb = set(R) | {b}

        # Paths ending at a through Sa
        paths_a = []
        if len(Sa) == 1:
            paths_a = [()]
        else:
            for p in permutations(sorted(Sa)):
                if p[-1] == a:
                    paths_a.append(tuple((p[k], p[k+1]) for k in range(len(p)-1)))

        # Paths ending at b through Rb
        paths_b = []
        if len(Rb) == 1:
            paths_b = [()]
        else:
            for p in permutations(sorted(Rb)):
                if p[-1] == b:
                    paths_b.append(tuple((p[k], p[k+1]) for k in range(len(p)-1)))

        for pa in paths_a:
            for pb in paths_b:
                # Weight = prod(r+s_e for e in pa) * prod(s_e - r for e in pb)
                # = prod(r+s_e for e in pa) * prod(-r+s_e for e in pb)
                weight = 1
                for u, v in pa:
                    weight *= (r + s[(u,v)])
                for u, v in pb:
                    weight *= (-r + s[(u,v)])

                terms.append({
                    'S': S, 'R': R,
                    'pa': pa, 'pb': pb,
                    'pa_edges': pa, 'pb_edges': pb,
                    'weight': expand(weight)
                })

    return terms

def extract_r_monomials(weight, r):
    """Extract (r_degree, s_monomial, coefficient) triples from a term."""
    p = Poly(weight, r)
    monomials = []
    for (k,), coeff in p.as_dict().items():
        monomials.append((k, expand(coeff)))
    return monomials

print("=" * 70)
print("VERTEX TRANSFER ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: n=5, analyze [r^1] and [r^3] cancellation
# ============================================================
n = 5
r, s, t = make_tournament(n)
a, b = 0, 1

terms = compute_G_terms(n, a, b, r, s)
print(f"\nn={n}: {len(terms)} total (S, path_a, path_b) triples")

# Group terms by r-degree
r_groups = defaultdict(list)
for term in terms:
    monos = extract_r_monomials(term['weight'], r)
    for k, coeff in monos:
        r_groups[k].append({
            'S': term['S'], 'R': term['R'],
            'pa': term['pa'], 'pb': term['pb'],
            's_coeff': coeff
        })

print(f"\n  Terms by r-degree:")
for k in sorted(r_groups.keys()):
    entries = r_groups[k]
    total = expand(sum(e['s_coeff'] for e in entries))
    print(f"  [r^{k}]: {len(entries)} terms, total = {total}")

# For [r^1]: analyze the s-monomials
print(f"\n  [r^1] detailed cancellation at n=5:")
r1_entries = r_groups.get(1, [])

# Group by s-monomial
from collections import Counter
s_monomial_groups = defaultdict(list)
for e in r1_entries:
    s_monomial_groups[str(e['s_coeff'])].append(e)

# Show groups that cancel
nonzero_groups = 0
for mono_str, entries in sorted(s_monomial_groups.items()):
    vals = [e['s_coeff'] for e in entries]
    total = expand(sum(vals))
    if total != 0:
        nonzero_groups += 1

print(f"  {len(s_monomial_groups)} distinct s-monomials, {nonzero_groups} with nonzero sum")

# Better: group by ABSOLUTE VALUE of s-monomial
abs_groups = defaultdict(list)
for e in r1_entries:
    coeff = e['s_coeff']
    # Use the string representation after removing sign
    abs_groups[str(abs(coeff)) if str(coeff).startswith('-') else str(coeff)].append(e)

# Actually, just collect all and check total
total_r1 = expand(sum(e['s_coeff'] for e in r1_entries))
print(f"  Total [r^1] = {total_r1}")

# For [r^3]: analyze at n=5 (if present)
if 3 in r_groups:
    r3_entries = r_groups[3]
    total_r3 = expand(sum(e['s_coeff'] for e in r3_entries))
    print(f"\n  [r^3] at n=5: {len(r3_entries)} terms, total = {total_r3}")

# ============================================================
# Part 2: The VERTEX TRANSFER involution
# ============================================================
print("\n" + "=" * 70)
print("Part 2: VERTEX TRANSFER INVOLUTION at [r^1]")
print("=" * 70)

# At [r^1], each term has ONE edge contributing r and the rest contributing s.
# The term is: (product of s-values of non-promoted edges) with some sign.
#
# The "promoted" edge is the one contributing r.
# It's in Pa if it contributes (r+s_e) -> r part.
# It's in Pb if it contributes (-r+s_e) -> -r part.
#
# So from Pa: the r^1 contribution is +r * (prod of s of other Pa edges) * (prod of Pb s-edges)
# From Pb: the r^1 contribution is -r * (prod of Pa s-edges) * (prod of s of other Pb edges)
# The sign from (-r) in Pb edges: if promoted edge is in Pb, it contributes -r.
# Actually, EACH Pb edge contributes (s_e - r), so the r part is -r.

# Let's enumerate: for each term, which edge is promoted?
print(f"\nn=5, [r^1]: Edge-by-edge decomposition")

r1_by_edge = defaultdict(list)

for term in terms:
    pa_edges = list(term['pa'])
    pb_edges = list(term['pb'])
    all_edges = pa_edges + pb_edges
    m = len(all_edges)
    if m == 0:
        continue

    # For each promoted edge:
    for j in range(m):
        is_pb = j >= len(pa_edges)
        promoted = all_edges[j]
        other_edges = all_edges[:j] + all_edges[j+1:]

        # s-product of other edges
        s_prod = 1
        for u, v in other_edges:
            s_prod *= s[(u,v)]

        # Sign from promoted edge: +r from Pa edge, -r from Pb edge
        sign = -1 if is_pb else 1

        r1_contrib = expand(sign * s_prod)

        r1_by_edge[(term['S'], promoted, is_pb)].append({
            'other_edges': other_edges,
            'contrib': r1_contrib,
            'full_pa': term['pa'],
            'full_pb': term['pb']
        })

# Now: for each PROMOTED edge, list all contributions
promoted_edge_totals = defaultdict(list)
for (S, promoted, is_pb), entries in r1_by_edge.items():
    for e in entries:
        promoted_edge_totals[(promoted, S)].append(e['contrib'])

# Actually, let's just look at which edges get promoted and from which S
print(f"\n  Total [r^1] contributions by promoted edge and subset S:")

# Collect contributions more carefully
edge_S_contrib = defaultdict(lambda: 0)
for term in terms:
    pa_edges = list(term['pa'])
    pb_edges = list(term['pb'])
    all_edges = pa_edges + pb_edges
    m = len(all_edges)
    if m == 0: continue

    for j in range(m):
        is_pb = j >= len(pa_edges)
        promoted = all_edges[j]
        other_edges = all_edges[:j] + all_edges[j+1:]

        s_prod = 1
        for u, v in other_edges:
            s_prod *= s[(u,v)]

        sign = -1 if is_pb else 1
        r1_contrib = expand(sign * s_prod)

        key = (term['S'], promoted, tuple(other_edges))
        edge_S_contrib[key] = expand(edge_S_contrib[key] + r1_contrib)

# For EACH s-monomial in [r^1], find which (S, promoted edge) pairs produce it
# and verify they cancel
print(f"\n  Collecting cancellation partners:")

# Group by s-monomial (the product of non-promoted edge weights)
s_mono_to_sources = defaultdict(list)
for (S, promoted, others), contrib in edge_S_contrib.items():
    if contrib == 0: continue
    s_mono_to_sources[str(expand(contrib))].append({
        'S': S, 'promoted': promoted, 'others': others,
        'value': contrib
    })

cancellation_count = 0
shown = 0
for mono_str in sorted(s_mono_to_sources.keys()):
    sources = s_mono_to_sources[mono_str]
    total = expand(sum(s['value'] for s in sources))

    # Check if there's a partner with negative value
    neg_str = str(expand(-sources[0]['value']))
    if neg_str in s_mono_to_sources:
        partner = s_mono_to_sources[neg_str]
        combined = expand(sum(s['value'] for s in sources) + sum(s['value'] for s in partner))
        if combined == 0:
            cancellation_count += 1
            if shown < 5:
                print(f"\n  Canceling pair:")
                for s_entry in sources[:2]:
                    print(f"    S={s_entry['S']}, promoted={s_entry['promoted']}, value={s_entry['value']}")
                for s_entry in partner[:2]:
                    print(f"    S={s_entry['S']}, promoted={s_entry['promoted']}, value={s_entry['value']}")
                shown += 1

print(f"\n  Total canceling pairs found: {cancellation_count}")

# ============================================================
# Part 3: Direct proof structure — vertex u transfer
# ============================================================
print("\n" + "=" * 70)
print("Part 3: TRANSFER VERTEX identification")
print("=" * 70)

# For each [r^1] contributing tuple, can we identify a "transfer vertex"?
# The transfer operation: take vertex u from S and move it to R (or vice versa),
# producing a new (S', Pa', Pb') with the same s-product but opposite sign.

n = 4
r, s, t = make_tournament(n)
a, b = 0, 1
U = [2, 3]

# At n=4, [r^1] G:
# The Wronskian sum gave us:
# s12+s13 (from S={}) + (s02-s13) (from S={2}) + (s03-s12) (from S={3}) + (-s02-s03) (from S={2,3})
#
# The pairings:
# s02: +1 from S={2}, -1 from S={2,3} -> transfer vertex 3 (move 3 into S)
# s03: +1 from S={3}, -1 from S={2,3} -> transfer vertex 2 (move 2 into S)
# s12: +1 from S={}, -1 from S={3} -> transfer vertex 3 (move 3 into S)
# s13: +1 from S={}, -1 from S={2} -> transfer vertex 2 (move 2 into S)

print("""
n=4 [r^1] pairings via vertex transfer:
  s02: S={2} (+) <-->[add 3]-> S={2,3} (-)
  s03: S={3} (+) <-->[add 2]-> S={2,3} (-)
  s12: S={}  (+) <-->[add 3]-> S={3}   (-)
  s13: S={}  (+) <-->[add 2]-> S={2}   (-)

Pattern: each s_{xy} with x in {a,b} and y in U is paired
by ADDING vertex y (or another U-vertex) to S.

The transfer vertex is NOT necessarily the vertex appearing in s_{xy}.
For s_{12}: the transfer vertex is 3 (not 2).
For s_{13}: the transfer vertex is 2 (not 3).
""")

# At n=5: check if similar vertex transfers explain the [r^1] cancellation
n = 5
r, s, t = make_tournament(n)
a, b = 0, 1
U = [2, 3, 4]

# Compute the Wronskian at [r^1] for each S
# From the previous analysis:
# S={}: Wronskian r^1 part involves s12, s13, s14 (and higher-order s-products)
# S={u}: involves s0u and others
# etc.

# Let me compute the FULL r^1 coefficient of G for each S separately
from sympy import sympify

def E_v(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1: return sympify(1)
    result = sympify(0)
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = sympify(1)
        for k in range(len(path)-1):
            u, v = path[k], path[k+1]
            w *= (r + s[(u,v)])
        result += w
    return expand(result)

print("n=5: [r^1] G contribution from each S (Wronskian)")
total_r1 = sympify(0)

for mask in range(2**len(U)):
    S = tuple(sorted([U[i] for i in range(len(U)) if (mask >> i) & 1]))
    R = tuple(sorted([u for u in U if u not in S]))

    f_S = E_v(set(S) | {a}, a, r, s)
    g_S = E_v(set(R) | {b}, b, r, s)

    # G contribution: f_S(r) * g_S(-r)
    g_neg = expand(g_S.subs(r, -r))
    product = expand(f_S * g_neg)

    # Extract [r^1]
    p = Poly(product, r)
    r1 = expand(p.nth(1)) if p.degree() >= 1 else sympify(0)

    total_r1 = expand(total_r1 + r1)

    if r1 != 0:
        # Decompose r1 into individual s-monomials
        print(f"\n  S={set(S)}: [r^1] = {r1}")

        # List individual terms
        terms_list = str(r1).replace(' - ', ' + -').split(' + ')
        for term_str in terms_list:
            term_str = term_str.strip()
            if term_str and term_str != '0':
                print(f"    {term_str}")

print(f"\n  Total [r^1] = {total_r1}")
print(f"  Is zero: {total_r1 == 0}")

# Now: for EACH s-monomial in the total, trace which S-subsets contribute
print(f"\n  Tracing individual s-monomials:")
from sympy import Add, Mul

# Compute [r^1] for each S and collect
S_contributions = {}
for mask in range(2**len(U)):
    S = tuple(sorted([U[i] for i in range(len(U)) if (mask >> i) & 1]))
    R = tuple(sorted([u for u in U if u not in S]))

    f_S = E_v(set(S) | {a}, a, r, s)
    g_S = E_v(set(R) | {b}, b, r, s)
    g_neg = expand(g_S.subs(r, -r))
    product = expand(f_S * g_neg)

    p = Poly(product, r)
    r1 = expand(p.nth(1)) if p.degree() >= 1 else sympify(0)
    S_contributions[S] = r1

# For each s-variable pair (i,j), find which S contribute it
for i in range(n):
    for j in range(i+1, n):
        sij = s[(i,j)]
        sij_name = f"s{i}{j}"

        contributors = []
        for S, r1 in S_contributions.items():
            if r1 == 0: continue
            # Check if sij appears as a term (at degree 1)
            coeff_pos = expand(r1.coeff(sij))
            if coeff_pos != 0:
                contributors.append((S, coeff_pos))

        if contributors:
            total_coeff = expand(sum(c for _, c in contributors))
            print(f"\n  {sij_name}: total coeff = {total_coeff}")
            for S, c in contributors:
                print(f"    S={set(S)}: coeff = {c}")

# ============================================================
# Part 4: The [r^3] coefficient at n=6 (numerical)
# ============================================================
print("\n" + "=" * 70)
print("Part 4: [r^3] cancellation structure (numerical, n=6)")
print("=" * 70)

import random
random.seed(123)

n_test = 6
s_num = {}
for i in range(n_test):
    for j in range(n_test):
        if i < j:
            val = random.uniform(-1, 1)
            s_num[(i,j)] = val
            s_num[(j,i)] = -val

def E_v_num(vertex_set, target, r_val, s_num):
    vs = list(vertex_set)
    if len(vs) == 1: return 1.0
    total = 0.0
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = 1.0
        for k in range(len(path)-1):
            u, v = path[k], path[k+1]
            w *= (r_val + s_num[(u,v)])
        total += w
    return total

a, b = 0, 1
U = [v for v in range(n_test) if v != a and v != b]

# Compute G(r) for several r values and fit polynomial
import numpy as np

r_values = np.linspace(-2, 2, 20)
G_values = []

for r_val in r_values:
    G = 0.0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        f = E_v_num(set(S + [a]), a, r_val, s_num)
        g = E_v_num(set(R + [b]), b, -r_val, s_num)
        G += f * g
    G_values.append(G)

G_values = np.array(G_values)

# Fit polynomial and check coefficients
from numpy.polynomial import polynomial as P
coeffs = np.polyfit(r_values, G_values, n_test - 2)
coeffs = coeffs[::-1]  # low to high degree

print(f"\nn=6: Polynomial fit of G(r):")
for k, c in enumerate(coeffs):
    parity = "EVEN" if k % 2 == 0 else "ODD"
    print(f"  [r^{k}] = {c:+.10f}  ({parity})")
    if k % 2 == 1:
        print(f"    ^ Should be ~0. Is it? |c| = {abs(c):.2e}")

# Double-check: G(r) vs G(-r) at specific points
print(f"\n  Verification G(r) = G(-r):")
for r_val in [0.5, 1.0, 1.5]:
    G_pos = 0.0
    G_neg = 0.0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        G_pos += E_v_num(set(S+[a]), a, r_val, s_num) * E_v_num(set(R+[b]), b, -r_val, s_num)
        G_neg += E_v_num(set(S+[a]), a, -r_val, s_num) * E_v_num(set(R+[b]), b, r_val, s_num)
    print(f"  r={r_val}: |G(r)-G(-r)| = {abs(G_pos - G_neg):.2e}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
