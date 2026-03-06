#!/usr/bin/env python3
"""
KEY IDENTITY: odd(B_b(V\{a})) = r * sum_v M^{(a)}[v,b]

This identity is the LAST STEP needed for the inductive proof that M[a,b]
has only even powers of r.

WHAT IT SAYS:
  B_b(V\{a}) = sum of Hamiltonian paths starting at b through V\{a}
  M^{(a)}[v,b] = transfer matrix entry on V\{a} (vertex a deleted)

  The odd-in-r part of B_b equals r times the column sum of M^{(a)}[*,b].

WHY IT MIGHT BE TRUE:
  B_b(V\{a}) counts all Hamiltonian paths b -> ... through V\{a}.
  M^{(a)}[v,b] counts signed 2-path-covers of V\{a} with endpoints v,b.

  The column sum sum_v M^{(a)}[v,b] is a "boundary" quantity.

APPROACH 1: Expand both sides and match.
APPROACH 2: Find a combinatorial/algebraic proof.
APPROACH 3: Relate to known identities (e.g., cofactor expansion).

Instance: opus-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Poly, Symbol, factor

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

def ham_paths_from(vertex_set, start, t_fn):
    """Sum of all Hamiltonian paths starting at 'start' through vertex_set."""
    vs = sorted(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    for perm in permutations(vs):
        if perm[0] != start:
            continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, vertex_set, a, b):
    """Compute M[a,b] for tournament on vertex_set."""
    V = sorted(vertex_set)
    U = [v for v in V if v != a and v != b]
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

def odd_part(expr, r):
    """Extract the odd-in-r part of expr."""
    if expr == 0:
        return 0
    p = Poly(expr, r)
    d = p.as_dict()
    return expand(sum(c * r**deg[0] for deg, c in d.items() if deg[0] % 2 == 1))

def even_part(expr, r):
    """Extract the even-in-r part of expr."""
    if expr == 0:
        return 0
    p = Poly(expr, r)
    d = p.as_dict()
    return expand(sum(c * r**deg[0] for deg, c in d.items() if deg[0] % 2 == 0))


# ============================================================
# TEST 1: Verify the key identity at n=4,5,6
# ============================================================
print("=" * 70)
print("KEY IDENTITY TEST: odd(B_b(V\\{a})) = r * sum_v M^{(a)}[v,b]")
print("=" * 70)

for n in [3, 4, 5, 6]:
    r, sv, s, t = setup(n)
    a, b = 0, 1
    V = set(range(n))
    V_minus_a = V - {a}
    U = sorted(V - {a, b})

    # LHS: odd part of B_b(V\{a})
    Bb = ham_paths_from(V_minus_a, b, t)
    lhs = odd_part(Bb, r)

    # RHS: r * sum_v M^{(a)}[v,b] for v in U
    col_sum = 0
    for v in U:
        M_sub = transfer_M(t, V_minus_a, v, b)
        col_sum = expand(col_sum + M_sub)
    rhs = expand(r * col_sum)

    match = expand(lhs - rhs) == 0
    print(f"\nn={n}: KEY IDENTITY holds = {match}")

    if n <= 5:
        print(f"  B_b(V\\{{a}}) = {Bb}")
        print(f"  odd(B_b) = {lhs}")
        print(f"  col_sum = {col_sum}")
        print(f"  r * col_sum = {rhs}")

    if not match:
        print(f"  DIFFERENCE = {expand(lhs - rhs)}")

    # ALSO: Check the EVEN part identity
    # even(B_b) - even(s_part) should give even(M_n[a,b])
    # Actually let's check: what IS even(B_b)?
    even_Bb = even_part(Bb, r)
    if n <= 5:
        print(f"  even(B_b) = {even_Bb}")


# ============================================================
# TEST 2: Can we PROVE the key identity?
# ============================================================
print("\n" + "=" * 70)
print("STRUCTURAL ANALYSIS OF THE KEY IDENTITY")
print("=" * 70)

# The key identity says:
# odd(B_b(W)) = r * sum_{v in W\{b}} M_W[v,b]
# where W = V\{a}.
#
# This is a statement PURELY about W (no reference to a anymore).
# Let's rename: for a set W with distinguished element b,
#   odd(B_b(W)) = r * sum_{v in W\{b}} M_W[v,b]
#
# This means: the odd-r-part of the total Hamiltonian path weight from b
# equals r times the column-b sum of the transfer matrix on W.

# Let's test this for ALL choices of b in W (not just b=1):
print("\nTest for all choices of b:")
for n_W in [3, 4, 5]:
    r, sv, s, t = setup(n_W + 1)  # need enough symbols
    W = set(range(1, n_W + 1))  # W = {1, ..., n_W}

    for b in sorted(W):
        Bb = ham_paths_from(W, b, t)
        lhs = odd_part(Bb, r)

        col_sum = 0
        for v in sorted(W - {b}):
            M_sub = transfer_M(t, W, v, b)
            col_sum = expand(col_sum + M_sub)
        rhs = expand(r * col_sum)

        match = expand(lhs - rhs) == 0
        print(f"  |W|={n_W}, b={b}: {match}")
        if not match:
            print(f"    DIFF = {expand(lhs - rhs)}")


# ============================================================
# TEST 3: What is the column sum of M?
# ============================================================
print("\n" + "=" * 70)
print("COLUMN SUMS OF M")
print("=" * 70)

# sum_v M[v,b] where v ranges over W\{b}
# = sum_v sum_S (-1)^|S| E_v(S+{v}) B_b(R+{b})
#
# Interchanging sums: for each partition S,R of W\{v,b},
# sum over v not in S,R... no, v is always an endpoint.
#
# Actually M[v,b] uses 2-path-covers of W where v is the "a" endpoint.
# So sum_v M[v,b] sums over all v != b.
#
# Let's compute and see if there's a pattern.

for n_W in [3, 4, 5]:
    r, sv, s, t = setup(n_W + 1)
    W = set(range(1, n_W + 1))
    b = 1

    print(f"\n|W|={n_W}, b={b}:")
    col_sum = 0
    for v in sorted(W - {b}):
        M_vb = transfer_M(t, W, v, b)
        print(f"  M[{v},{b}] = {M_vb}")
        col_sum = expand(col_sum + M_vb)

    print(f"  Column sum = {col_sum}")

    # Is the column sum related to B_b / r?
    Bb = ham_paths_from(W, b, t)
    # odd(Bb) / r should equal col_sum
    odd_Bb = odd_part(Bb, r)
    if odd_Bb != 0:
        ratio_check = expand(odd_Bb - r * col_sum)
        print(f"  odd(B_b)/r = col_sum? {ratio_check == 0}")

    # What about: col_sum = d/dr B_b evaluated somehow?
    # B_b has degree n_W - 1 in r (n_W - 1 edges in path through n_W vertices)
    # d/dr B_b has degree n_W - 2
    # col_sum has degree n_W - 3 (from M which has degree n_W - 2, and sum doesn't increase)
    from sympy import diff as sym_diff
    dBb = expand(sym_diff(Bb, r))
    print(f"  d/dr B_b = {dBb}")
    print(f"  col_sum = {col_sum}")

    # At r=0:
    col_sum_0 = expand(col_sum.subs(r, 0))
    dBb_0 = expand(dBb.subs(r, 0))
    print(f"  col_sum(0) = {col_sum_0}")
    print(f"  d/dr B_b(0) = {dBb_0}")
    print(f"  Match at r=0? {expand(col_sum_0 - dBb_0) == 0}")


# ============================================================
# TEST 4: Reformulate using B_b(W) decomposition
# ============================================================
print("\n" + "=" * 70)
print("REFORMULATION: B_b as column expansion")
print("=" * 70)

# B_b(W) = sum over all Hamiltonian paths b -> v1 -> v2 -> ... -> v_{m}
# where m = |W| - 1 and {v1,...,v_m} = W\{b}.
#
# Group by the LAST vertex v_m:
# B_b(W) = sum_{v in W\{b}} (sum of paths b -> ... -> v) * 1
#         = sum_v E_v(W) ... wait, E_v counts paths ENDING at v.
#
# Actually B_b counts paths BEGINNING at b.
# If we group by the last vertex v:
# B_b(W) = sum_v [paths b -> ... -> v through W]
#
# A path b -> ... -> v = B_b starting at b, ending at v, through W
# = sum of products t(b,w1)*t(w1,w2)*...*t(w_{m-1},v)
#
# This is NOT directly M[v,b] (which involves 2-path-covers, not single paths).

# But the KEY IDENTITY relates B_b to M. Can we find an intermediate step?

# IDEA: Write B_b(W) = sum_v H(b,v,W) where H(b,v,W) = Hamiltonian path weight b->v through W
# Then: odd(B_b) = sum_v odd(H(b,v,W))
# And: r * sum_v M[v,b] = r * sum_v M[v,b]
# So: sum_v odd(H(b,v,W)) = r * sum_v M[v,b]
# Does this hold TERMWISE? i.e., odd(H(b,v,W)) = r * M[v,b]?

print("\nTest: odd(H(b,v,W)) = r * M[v,b] termwise?")
for n_W in [3, 4, 5]:
    r, sv, s, t = setup(n_W + 1)
    W = set(range(1, n_W + 1))
    b = 1

    for v in sorted(W - {b}):
        # H(b,v,W) = Hamiltonian paths from b to v through W
        H_bv = 0
        for perm in permutations(sorted(W)):
            if perm[0] != b or perm[-1] != v:
                continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            H_bv += prod
        H_bv = expand(H_bv)

        odd_H = odd_part(H_bv, r)
        M_vb = transfer_M(t, W, v, b)
        rM = expand(r * M_vb)

        match = expand(odd_H - rM) == 0
        print(f"  |W|={n_W}, b={b}, v={v}: termwise match = {match}")
        if not match:
            print(f"    odd(H) = {odd_H}")
            print(f"    r*M[v,b] = {rM}")
            print(f"    diff = {expand(odd_H - rM)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
