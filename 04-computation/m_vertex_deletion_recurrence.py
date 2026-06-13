#!/usr/bin/env python3
"""
VERTEX DELETION RECURRENCE for M[a,b]

Goal: Express M_n[a,b] in terms of M_{n-1}[a,b] (with vertex v deleted)
plus correction terms involving t(v, *) and t(*, v).

If the recurrence preserves even-r-powers, this gives an inductive proof.

Approach: Delete vertex v from U. Each 2-path-cover either has v in the
a-path or the b-path. Condition on v's neighbors in the path to get
a recurrence.

Also: test whether M_{n-1}[a,b] (with v deleted) is a "leading term" of M_n.

kind-pasteur-2026-03-06-S23b
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

def transfer_M(t_fn, n_verts, vertex_set, a, b):
    """Compute M[a,b] for a tournament on vertex_set."""
    U = [v for v in vertex_set if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        S_set = set(S) | {a}
        R_set = set(R) | {b}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1

        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("VERTEX DELETION RECURRENCE FOR M[a,b]")
print("=" * 70)

# ============================================================
# Part 1: Compare M_n[0,1] with M_{n-1}[0,1] (deleting vertex n-1)
# ============================================================

for n in [4, 5, 6]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    v = n - 1  # vertex to delete

    # Full M_n[0,1]
    full_verts = list(range(n))
    M_n = transfer_M(t, n, full_verts, a, b)

    # Reduced M_{n-1}[0,1] (delete vertex v)
    reduced_verts = [i for i in range(n) if i != v]
    M_n1 = transfer_M(t, n, reduced_verts, a, b)

    # Difference
    diff = expand(M_n - M_n1)

    print(f"\nn={n}, deleting vertex {v}:")
    print(f"  M_n[0,1] = {M_n}")
    print(f"  M_{{n-1}}[0,1] = {M_n1}")
    print(f"  Difference = {diff}")

    # Check: is the difference even in r?
    if diff != 0:
        p_diff = Poly(diff, r)
        all_even = True
        for k in range(p_diff.degree() + 1):
            coeff = expand(p_diff.nth(k))
            if k % 2 == 1 and coeff != 0:
                all_even = False
                print(f"    diff r^{k} = {coeff}")
        if all_even:
            print(f"    Difference has ONLY EVEN r-powers!")
        else:
            print(f"    Difference has ODD r-powers.")

    # Check: is M_{n-1} even in r?
    if M_n1 != 0:
        p_m = Poly(M_n1, r)
        all_even_m = True
        for k in range(p_m.degree() + 1):
            if k % 2 == 1 and expand(p_m.nth(k)) != 0:
                all_even_m = False
        print(f"    M_{{n-1}} even in r: {all_even_m}")

# ============================================================
# Part 2: Decompose the difference by WHERE v appears
# ============================================================
print("\n" + "=" * 70)
print("Part 2: Decomposition by v's role")
print("=" * 70)

for n in [4, 5]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    v = n - 1
    U = [u for u in range(n) if u != a and u != b]
    U_no_v = [u for u in U if u != v]

    # Terms where v is in the a-path (v in S)
    terms_v_in_S = 0
    # Terms where v is in the b-path (v in R)
    terms_v_in_R = 0

    for mask in range(1 << len(U)):
        S_list = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R_list = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S_list)
        S_set = set(S_list) | {a}
        R_set = set(R_list) | {b}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1

        contrib = expand(sign * ea * bb)

        if v in S_list:
            terms_v_in_S += contrib
        else:
            terms_v_in_R += contrib

    terms_v_in_S = expand(terms_v_in_S)
    terms_v_in_R = expand(terms_v_in_R)

    print(f"\nn={n}, v={v}:")
    print(f"  Terms with v in a-path: {terms_v_in_S}")
    print(f"  Terms with v in b-path: {terms_v_in_R}")
    print(f"  Sum = {expand(terms_v_in_S + terms_v_in_R)}")

    # Check r-parity of each part
    for label, terms in [("v in a-path", terms_v_in_S), ("v in b-path", terms_v_in_R)]:
        if terms != 0:
            p = Poly(terms, r)
            odd_terms = []
            for k in range(p.degree() + 1):
                if k % 2 == 1 and expand(p.nth(k)) != 0:
                    odd_terms.append(k)
            print(f"    {label}: odd r-powers at {odd_terms if odd_terms else 'NONE'}")

    # KEY: Does terms_v_in_S + terms_v_in_R cancel odd powers?
    # Check: terms_v_in_S at -r vs terms_v_in_R at r (or vice versa)
    # By the M(-r) = M(b,a) identity, the "v in S" part at -r should relate to
    # "v in R" part of M[b,a].

# ============================================================
# Part 3: Alternative decomposition - condition on v's neighbors
# ============================================================
print("\n" + "=" * 70)
print("Part 3: Condition on v's neighbors in path")
print("=" * 70)

n = 5
r, sv, s, t = setup(n)
a, b = 0, 1
v = 4  # vertex to focus on
U = [u for u in range(n) if u != a and u != b]

# For each cover, record:
# - Whether v is in a-path or b-path
# - v's predecessor and successor in the path
# - The contribution

neighbor_groups = defaultdict(lambda: 0)

for mask in range(1 << len(U)):
    S_list = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R_list = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1)**len(S_list)
    S_set = set(S_list) | {a}
    R_set = set(R_list) | {b}

    for p1 in permutations(sorted(S_set)):
        if p1[-1] != a: continue
        for p2 in permutations(sorted(R_set)):
            if p2[0] != b: continue

            arcs = []
            for i in range(len(p1)-1):
                arcs.append((p1[i], p1[i+1]))
            for i in range(len(p2)-1):
                arcs.append((p2[i], p2[i+1]))

            prod_val = sign
            for (u, w) in arcs:
                prod_val *= t(u, w)
            prod_val = expand(prod_val)

            # Find v in p1 or p2
            if v in S_set:
                path = p1
                path_type = 'a'
            else:
                path = p2
                path_type = 'b'

            idx = list(path).index(v)
            pred = path[idx-1] if idx > 0 else None
            succ = path[idx+1] if idx < len(path)-1 else None

            key = (path_type, pred, succ)
            neighbor_groups[key] = expand(neighbor_groups[key] + prod_val)

print(f"\nn=5, v=4: Contributions grouped by v's neighbors")
for key in sorted(neighbor_groups.keys(), key=str):
    val = neighbor_groups[key]
    if val != 0:
        path_type, pred, succ = key
        print(f"  path={path_type}, pred={pred}, succ={succ}: {val}")

        # Check r-parity
        p = Poly(val, r)
        has_odd = any(expand(p.nth(k)) != 0 for k in range(p.degree()+1) if k % 2 == 1)
        if has_odd:
            for k in range(p.degree()+1):
                if k % 2 == 1:
                    c = expand(p.nth(k))
                    if c != 0:
                        print(f"       r^{k}: {c}")

# Check: do neighbor groups pair up to cancel odd powers?
print(f"\n  Pairing analysis:")
keys = sorted(neighbor_groups.keys(), key=str)
used = set()
for key in keys:
    if key in used or neighbor_groups[key] == 0:
        continue
    path_type, pred, succ = key
    # Natural partner: swap path type, swap pred/succ roles
    partner_candidates = []

    # If v is in a-path between pred->v->succ, the partner might be
    # v in b-path between succ'->v->pred' (reversed role)
    other_type = 'b' if path_type == 'a' else 'a'
    # Try: partner = (other_type, succ, pred)
    partner = (other_type, succ, pred)
    if partner in neighbor_groups and partner not in used:
        combined = expand(neighbor_groups[key] + neighbor_groups[partner])
        p_combined = Poly(combined, r)
        has_odd = any(expand(p_combined.nth(k)) != 0 for k in range(p_combined.degree()+1) if k % 2 == 1)
        print(f"  {key} + {partner}: odd r-powers = {'YES' if has_odd else 'NO'}, total = {combined}")
        if not has_odd:
            used.add(key)
            used.add(partner)

for key in keys:
    if key not in used and neighbor_groups[key] != 0:
        print(f"  UNPAIRED: {key}: {neighbor_groups[key]}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
