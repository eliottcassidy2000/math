#!/usr/bin/env python3
"""
Attacking the LAST GAP: why do odd r-powers vanish in M[a,b]?

M[a,b] = sum_{S subset U} (-1)^|S| E_a(S+a) * B_b(R+b)

Each arc weight is r + s_ij where s_ij = -s_ji.
M is polynomial of degree n-2 in r.

CLAIM: odd powers of r vanish.

APPROACH: Expand each E_a * B_b as a product of n-2 arc weights.
Each product contributes to r^k via elementary symmetric functions
of the s-values along the arcs.

The r^k coefficient of a single 2-path-cover (path ending at a through S+a,
path starting at b through R+b) with arcs having skew values s_1,...,s_{n-2}
is: C(n-2,k) * e_{n-2-k}(s_1,...,s_{n-2}) ... no, it's actually the
elementary symmetric polynomial evaluated at s_1,...,s_{n-2}.

Wait: product_{i=1}^{m} (r + s_i) = sum_{k=0}^{m} r^k * e_{m-k}(s_1,...,s_m)

So the r^k coeff of a cover with arcs s_1,...,s_m (where m=n-2) is e_{m-k}(s_1,...,s_m).

For odd k, we need: sum over all 2-path-covers, weighted by (-1)^|S|,
of e_{m-k}(arc skew values) = 0.

Let's investigate this at small n.

opus-2026-03-06-S21
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, Rational, Symbol
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

    return r, sv, s

print("=" * 70)
print("EVEN r-POWERS: STRUCTURAL ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: Enumerate all 2-path-covers and their arc sets
# ============================================================
print("\n--- Part 1: 2-path-cover enumeration ---")

for n in [4, 5]:
    r, sv, s = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = n - 2  # number of arcs in each 2-path-cover

    print(f"\n  n={n}, m={m} arcs per cover:")

    # A "2-path-cover" consists of:
    # - A subset S of U
    # - A Hamiltonian path P1 through S+a ending at a
    # - A Hamiltonian path P2 through R+b starting at b
    # Together P1 and P2 use exactly n-2 arcs covering all n vertices.

    all_covers = []
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        # Enumerate paths ending at a through S+a
        S_set = set(S) | {a}
        for p1 in permutations(sorted(S_set)):
            if p1[-1] != a: continue
            # Enumerate paths starting at b through R+b
            R_set = set(R) | {b}
            for p2 in permutations(sorted(R_set)):
                if p2[0] != b: continue
                # Collect arcs
                arcs = []
                for i in range(len(p1)-1):
                    arcs.append((p1[i], p1[i+1]))
                for i in range(len(p2)-1):
                    arcs.append((p2[i], p2[i+1]))
                all_covers.append((sign, arcs, S))

    print(f"  Total signed covers: {len(all_covers)}")

    # For each cover, compute arc skew values and the r^k contributions
    r1_total = 0
    for sign, arcs, S in all_covers:
        arc_s_vals = [s(i, j) for (i, j) in arcs]
        # r^1 coefficient = e_{m-1}(s_1,...,s_m) = sum of products of m-1 of the s_i
        # = sum_i product_{j != i} s_j
        # Actually: product (r + s_i) = r^m + r^{m-1} e_1 + r^{m-2} e_2 + ... + e_m
        # So r^k coeff = e_{m-k}
        # r^1 coeff = e_{m-1} = sum_i (product of all s_j except s_i)

        # For m=2 (n=4): e_1(s1,s2) = s1 + s2. r^1 coeff = e_1.
        # For m=3 (n=5): e_2(s1,s2,s3) = s1*s2 + s1*s3 + s2*s3. r^1 coeff = e_2.

        from sympy import prod as sprod
        total_prod = 1
        for sv_val in arc_s_vals:
            total_prod *= sv_val
        total_prod = expand(total_prod)

        # e_{m-1} = sum_i product_{j!=i} s_j
        e_m_minus_1 = 0
        for i in range(len(arc_s_vals)):
            term = 1
            for j in range(len(arc_s_vals)):
                if j != i:
                    term *= arc_s_vals[j]
            e_m_minus_1 += term
        e_m_minus_1 = expand(e_m_minus_1)

        r1_total += sign * e_m_minus_1

    r1_total = expand(r1_total)
    print(f"  r^1 coefficient (sum of signed e_{{m-1}}): {r1_total}")

# ============================================================
# Part 2: Group covers by their ARC MULTISETS
# ============================================================
print("\n--- Part 2: Arc multiset analysis ---")

for n in [4]:
    r, sv, s = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    # Group covers by arc tuple (sorted)
    arc_groups = defaultdict(list)
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        S_set = set(S) | {a}
        R_set = set(R) | {b}
        for p1 in permutations(sorted(S_set)):
            if p1[-1] != a: continue
            for p2 in permutations(sorted(R_set)):
                if p2[0] != b: continue
                arcs = []
                for i in range(len(p1)-1):
                    arcs.append((p1[i], p1[i+1]))
                for i in range(len(p2)-1):
                    arcs.append((p2[i], p2[i+1]))
                key = tuple(sorted(arcs))
                arc_groups[key].append((sign, S, list(p1), list(p2)))

    print(f"\n  n={n}: {len(arc_groups)} distinct arc multisets")
    for key, entries in sorted(arc_groups.items()):
        total_sign = sum(e[0] for e in entries)
        arc_vals = [s(i,j) for (i,j) in key]
        print(f"    arcs={key}: {len(entries)} covers, net sign={total_sign:+d}")
        for sign, S, p1, p2 in entries:
            print(f"      S={S}, P1={p1}, P2={p2}, sign={sign:+d}")

# ============================================================
# Part 3: Key insight — what is the r^1 coefficient REALLY?
# ============================================================
print("\n--- Part 3: What is the r^1 coefficient? ---")
print("""
The r^1 coefficient of M[a,b] is:
  sum_{S} (-1)^|S| * sum_{covers C of (S,R)} e_{m-1}(s(C))

where e_{m-1}(s(C)) = sum_i prod_{j!=i} s_j for the arcs of cover C.

This is equivalent to: for each cover C, "remove one arc" and take
the product of the remaining arcs' s-values, summed over all choices
of which arc to remove.

Removing arc (u,v) from cover C and taking the product of the
remaining s-values gives: (product of all arc s-values) / s_{uv}.

So e_{m-1} = (total product) * sum_i (1/s_i).

This is a "logarithmic derivative" structure.
""")

# ============================================================
# Part 4: Involution approach
# ============================================================
print("\n--- Part 4: Involution that cancels odd-r terms ---")
print("""
IDEA: Find an involution on the set of (S, cover) pairs that:
1. Preserves the arc multiset (so same product of arc weights)
2. Changes (-1)^|S| by -1 (so signs cancel)
3. Preserves even-degree r contributions but flips odd-degree ones

OR: Find a sign-reversing involution on the r^{odd} contributions.

The simplest involution: S <-> U\\S (complement).
This changes sign by (-1)^{|U|} = (-1)^{n-2}.
And it swaps the E_a and B_b roles.

For n even: (-1)^{n-2} = +1, so complement preserves sign.
  Can't use complement alone to cancel.
For n odd: (-1)^{n-2} = -1, so complement flips sign.
  Complement pairing would cancel EVERYTHING, but M is nonzero.

So complement alone doesn't work. We need a more refined involution.
""")

# ============================================================
# Part 5: The r -> -r substitution directly
# ============================================================
print("\n--- Part 5: What does r -> -r do to each factor? ---")

for n in [4, 5]:
    r, sv, s = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    print(f"\n  n={n}:")
    # M(r) = sum_S (-1)^|S| E_a(S+a; r,s) * B_b(R+b; r,s)
    # M(-r) = sum_S (-1)^|S| E_a(S+a; -r,s) * B_b(R+b; -r,s)

    # Now E_a(S+a; -r,s) = E_a with arc weights -r+s_ij = -(r-s_ij) = -(r+s_ji)
    # So E_a(S+a; -r,s) = (-1)^|S| * E_a(S+a; r, -s)
    # Wait, each arc weight becomes -r + s_ij. There are |S| arcs in the path through S+a.
    # Product of |S| factors (-r + s_i) = (-1)^|S| * product(r - s_i)
    # = (-1)^|S| * product(r + (-s_i))
    # = (-1)^|S| * (path weight with r, -s)

    # Hmm, that's not quite right. Let me be more careful.
    # Arc weight at (u,v) is r + s(u,v). Under r -> -r: becomes -r + s(u,v).
    # This is -(r - s(u,v)) = -(r + s(v,u)).
    # So it's -1 times the arc weight of the REVERSED arc with the SAME r.

    # A path v0->v1->...->vk has weight prod_{i} (r + s(vi, v_{i+1})).
    # Under r -> -r: prod (-r + s(vi, v_{i+1})) = prod (-(r - s(vi, v_{i+1})))
    #   = (-1)^k * prod (r - s(vi, v_{i+1}))
    #   = (-1)^k * prod (r + s(v_{i+1}, vi))
    #   = (-1)^k * (weight of reversed path with same r)

    # So: E_a(S+a; -r, s) = sum over paths ending at a of (-1)^|S| * (reversed path weight at r)
    # Reversed paths ending at a = paths starting at a
    # So: E_a(S+a; -r, s) = (-1)^|S| * B_a(S+a; r, s)

    # Similarly: B_b(R+b; -r, s) = (-1)^|R| * E_b(R+b; r, s)

    # Therefore:
    # M(-r, s) = sum_S (-1)^|S| * (-1)^|S| B_a(S+a; r, s) * (-1)^|R| E_b(R+b; r, s)
    #          = sum_S (-1)^{2|S|+|R|} * B_a(S+a; r, s) * E_b(R+b; r, s)
    #          = sum_S (-1)^{|R|} * B_a(S+a) * E_b(R+b)
    #          = sum_S (-1)^{n-2-|S|} * B_a(S+a) * E_b(R+b)
    #          = (-1)^{n-2} * sum_S (-1)^|S| * B_a(S+a) * E_b(R+b)

    # But sum_S (-1)^|S| B_a(S+a) E_b(R+b) — what is this?
    # It's M[b,a] with E and B swapped!
    # Wait: M[b,a] = sum_S (-1)^|S| E_b(S+b) B_a(R+a)

    # Our expression is sum_S (-1)^|S| B_a(S+a) E_b(R+b)
    # This is NOT M[b,a]. It has B_a on S+a (not R+a) and E_b on R+b (not S+b).

    # Relabel S <-> R:
    # sum_R (-1)^{|U|-|R|} B_a(R+a) E_b(S+b)
    # = (-1)^{n-2} sum_S (-1)^|S| E_b(S+b) B_a(R+a)
    # = (-1)^{n-2} M[b,a]

    # So: sum_S (-1)^|S| B_a(S+a) E_b(R+b) = (-1)^{n-2} M[b,a]

    # Therefore: M(-r, s) = (-1)^{n-2} * (-1)^{n-2} * M[b,a](r,s)
    #                     = (-1)^{2(n-2)} * M[b,a](r,s)
    #                     = M[b,a](r,s)

    # WAIT! This gives M[a,b](-r,s) = M[b,a](r,s).

    # For M to be even in r: M(r) = M(-r)
    # <=> M[a,b](r,s) = M[b,a](r,s)
    # <=> M is symmetric!

    # So "M is even in r" is DIRECTLY EQUIVALENT to "M[a,b] = M[b,a]".
    # This is just another restatement, not a new proof.

    # BUT WAIT — let me verify this derivation!!!

    def t_fn(rv):
        def f(i, j):
            if i == j: return 0
            return rv + s(i, j)
        return f

    def hp_fn(t, vset, start=None, end=None):
        vl = sorted(vset)
        k = len(vl)
        if k == 0: return 0
        if k == 1:
            if start is not None and vl[0] != start: return 0
            if end is not None and vl[0] != end: return 0
            return 1
        total = 0
        for p in permutations(vl):
            if start is not None and p[0] != start: continue
            if end is not None and p[-1] != end: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            total += prod
        return expand(total)

    def M_fn(t, n_v, a_v, b_v):
        U_v = [v for v in range(n_v) if v != a_v and v != b_v]
        result = 0
        for mask in range(1 << len(U_v)):
            S = [U_v[i] for i in range(len(U_v)) if mask & (1 << i)]
            R = [U_v[i] for i in range(len(U_v)) if not (mask & (1 << i))]
            sign_v = (-1)**len(S)
            ea = hp_fn(t, set(S)|{a_v}, end=a_v)
            bb = hp_fn(t, set(R)|{b_v}, start=b_v)
            result += sign_v * ea * bb
        return expand(result)

    M_ab_r = M_fn(t_fn(r), n, a, b)
    M_ab_neg_r = M_fn(t_fn(-r), n, a, b)
    M_ba_r = M_fn(t_fn(r), n, b, a)

    print(f"    M[a,b](-r,s) = M[b,a](r,s)? {expand(M_ab_neg_r - M_ba_r) == 0}")
    print(f"    This confirms: M even in r <=> M[a,b] = M[b,a]")

# ============================================================
# Part 6: THE KEY DERIVATION (clean)
# ============================================================
print("\n" + "=" * 70)
print("Part 6: CLEAN DERIVATION — M[a,b](-r,s) = M[b,a](r,s)")
print("=" * 70)
print("""
THEOREM: M[a,b](-r, s) = M[b,a](r, s).

PROOF:
  Arc weight t(u,v) = r + s(u,v) where s(u,v) = -s(v,u).
  Under r -> -r: t(u,v) becomes -r + s(u,v) = -(r - s(u,v)) = -(r + s(v,u)).

  For a path P: v0 -> v1 -> ... -> vk (using k arcs):
    w(P; -r) = prod_{i=0}^{k-1} (-r + s(vi, v_{i+1}))
             = (-1)^k * prod (r + s(v_{i+1}, vi))
             = (-1)^k * w(P^rev; r)

  where P^rev: vk -> ... -> v1 -> v0 is the reversed path.

  Therefore:
    E_a(S+a; -r, s) = (-1)^|S| * B_a(S+a; r, s)   [paths of length |S|]
    B_b(R+b; -r, s) = (-1)^|R| * E_b(R+b; r, s)    [paths of length |R|]

  Now:
    M[a,b](-r, s) = sum_S (-1)^|S| * E_a(S+a; -r) * B_b(R+b; -r)
                   = sum_S (-1)^|S| * [(-1)^|S| B_a(S+a; r)] * [(-1)^|R| E_b(R+b; r)]
                   = sum_S (-1)^{|S| + |S| + |R|} * B_a(S+a; r) * E_b(R+b; r)
                   = sum_S (-1)^{2|S| + n-2-|S|} * B_a(S+a; r) * E_b(R+b; r)
                   = (-1)^{n-2} * sum_S (-1)^{|S|} * B_a(S+a; r) * E_b(R+b; r)

  Relabel S <-> R (replace S by U\\S):
    = (-1)^{n-2} * sum_S (-1)^{n-2-|S|} * B_a(R+a; r) * E_b(S+b; r)
    = (-1)^{2(n-2)} * sum_S (-1)^{|S|} * E_b(S+b; r) * B_a(R+a; r)
    = sum_S (-1)^|S| * E_b(S+b; r) * B_a(R+a; r)
    = M[b,a](r, s)

  QED.

COROLLARY: M[a,b] = M[b,a] (i.e., M is symmetric) if and only if
M[a,b](r,s) is even in r. And the PROOF of M[a,b](-r) = M[b,a](r)
is COMPLETE — it uses only path reversal and S<->R relabeling.

*** THIS IS NOT CIRCULAR! ***

The previous "circularity" arose from trying to prove BOTH
  (i) M[b,a] = (-1)^{n-2} M[a,b](c,-s)   [from path reversal on s]
  (ii) M[a,b] = M[b,a]                     [from (i) + parity]

But the CORRECT route is:
  (i) M[a,b](-r,s) = M[b,a](r,s)    [from path reversal on r]

This is a DIRECT EQUIVALENCE, not circular. However, it still doesn't
PROVE symmetry — it just says symmetry <=> even in r.

Wait... let me re-examine. We proved:
  M[a,b](-r,s) = M[b,a](r,s)

This is a THEOREM, not an equivalence. It holds for ALL r,s.
But it doesn't by itself prove M[a,b] = M[b,a] OR that M is even in r.
It RELATES the two.

So we're still stuck: we need an independent argument for one of:
  - M[a,b](r) = M[a,b](-r)  [even in r]
  - M[a,b] = M[b,a]          [symmetry]

These are equivalent by the theorem above, but neither is proved.
""")

# Verify the theorem
print("Verification:")
for n in [4, 5]:
    r, sv, s = setup(n)
    a, b = 0, 1

    def t_fn_r(rv):
        def f(i, j):
            if i == j: return 0
            return rv + s(i, j)
        return f

    def hp_f(t, vset, start=None, end=None):
        vl = sorted(vset)
        k = len(vl)
        if k == 0: return 0
        if k == 1:
            if start is not None and vl[0] != start: return 0
            if end is not None and vl[0] != end: return 0
            return 1
        total = 0
        for p in permutations(vl):
            if start is not None and p[0] != start: continue
            if end is not None and p[-1] != end: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            total += prod
        return expand(total)

    def M_f(t, n_v, a_v, b_v):
        U_v = [v for v in range(n_v) if v != a_v and v != b_v]
        result = 0
        for mask in range(1 << len(U_v)):
            S_l = [U_v[i] for i in range(len(U_v)) if mask & (1 << i)]
            R_l = [U_v[i] for i in range(len(U_v)) if not (mask & (1 << i))]
            sign_v = (-1)**len(S_l)
            ea = hp_f(t, set(S_l)|{a_v}, end=a_v)
            bb = hp_f(t, set(R_l)|{b_v}, start=b_v)
            result += sign_v * ea * bb
        return expand(result)

    M_ab_r = M_f(t_fn_r(r), n, a, b)
    M_ab_neg = M_f(t_fn_r(-r), n, a, b)
    M_ba_r = M_f(t_fn_r(r), n, b, a)

    check = expand(M_ab_neg - M_ba_r)
    print(f"  n={n}: M[a,b](-r,s) = M[b,a](r,s)? {check == 0}")

# ============================================================
# Part 7: NEW IDEA — interpolation proof
# ============================================================
print("\n" + "=" * 70)
print("Part 7: INTERPOLATION APPROACH")
print("=" * 70)
print("""
We know:
  (A) M[a,b](-r,s) = M[b,a](r,s)  [proved]
  (B) M[a,b](0,s) = M[b,a](0,s)   [proved: c=0 symmetry]

From (A): M[a,b](-r) = M[b,a](r).
From (B): at r=0, M[a,b](0) = M[b,a](0).

Define f(r) = M[a,b](r) - M[b,a](r) = M[a,b](r) - M[a,b](-r) [by (A)].
So f(r) is the ODD part of M[a,b](r) in r.
We know f(0) = 0 from (B).

BUT f(r) is a polynomial of degree n-2. If f is odd: f(r) = c_1 r + c_3 r^3 + ...
And f(0) = 0 is automatic for any odd polynomial!

So (B) gives no new information beyond what (A) already implies.
We need MORE.

Alternative: what other special values of r give information?
""")

# Check M at r = s_ij for specific arc values
print("What happens at r = some specific s-value?")
for n in [4]:
    r, sv, s = setup(n)
    a, b = 0, 1

    # At r = s02, the arc (0,2) has weight r + s02 = 2*s02
    # and the arc (2,0) has weight r + s20 = r - s02 = s02 - s02 = 0
    # So arc (2,0) vanishes!

    # This means: any path using arc (2,0) contributes 0.
    # The tournament becomes "partially directed" at this r value.

    print(f"\n  n={n}: at r = s02, arc (2,0) has weight 0")
    print("  This doesn't obviously help with parity...")

# ============================================================
# Part 8: Determinantal / matrix identity approach
# ============================================================
print("\n" + "=" * 70)
print("Part 8: Can M be expressed as a determinant?")
print("=" * 70)

for n in [4]:
    r, sv, s = setup(n)
    a, b = 0, 1

    def t_fn(i, j):
        if i == j: return 0
        return r + s(i, j)

    from sympy import Matrix, eye, ones, det

    # Build the n x n arc matrix A
    A = Matrix(n, n, lambda i, j: t_fn(i, j) if i != j else 0)
    I_n = eye(n)

    # The transfer matrix relates to det(I - zA) somehow
    # At z=1: det(I - A) and its adjugate
    IminA = I_n - A

    # Check: is M[a,b] related to a minor/cofactor of (I-A)?
    M_ab = M_f(t_fn, n, a, b) if 'M_f' in dir() else None

    def hp_f2(vset, start=None, end=None):
        vl = sorted(vset)
        k = len(vl)
        if k == 0: return 0
        if k == 1:
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

    def M_direct(n_v, a_v, b_v):
        U_v = [v for v in range(n_v) if v != a_v and v != b_v]
        result = 0
        for mask in range(1 << len(U_v)):
            S_l = [U_v[i] for i in range(len(U_v)) if mask & (1 << i)]
            R_l = [U_v[i] for i in range(len(U_v)) if not (mask & (1 << i))]
            sign_v = (-1)**len(S_l)
            ea = hp_f2(set(S_l)|{a_v}, end=a_v)
            bb = hp_f2(set(R_l)|{b_v}, start=b_v)
            result += sign_v * ea * bb
        return expand(result)

    M_ab = M_direct(n, a, b)

    # Try: M[a,b] = cofactor(I-A, b, a)?
    # cofactor(b,a) = (-1)^{a+b} det(minor obtained by deleting row b, col a)
    minor_ba = IminA.minor_submatrix(b, a)
    cof_ba = (-1)**(a+b) * det(minor_ba)
    cof_ba = expand(cof_ba)

    print(f"\n  n={n}:")
    print(f"    M[0,1] = {M_ab}")
    print(f"    cofactor(I-A, 1,0) = {cof_ba}")
    print(f"    Equal? {expand(M_ab - cof_ba) == 0}")

    # Try adj(I-A)[a,b]
    adj_ab = expand(IminA.adjugate()[a,b])
    print(f"    adj(I-A)[0,1] = {adj_ab}")
    print(f"    Equal? {expand(M_ab - adj_ab) == 0}")

    # What IS M[a,b]?
    # Hmm, let's compare with the PERMANENT of submatrices
    # or other algebraic objects

    # Actually, for the transfer matrix of a tournament,
    # M[a,b] counts 2-path-covers with inclusion-exclusion.
    # This is the (a,b) entry of the "path matrix" of the tournament.

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
STATUS: The even-r-powers property remains unproved.

We have:
  M[a,b](-r,s) = M[b,a](r,s)  [PROVED]
  M[a,b](0,s) = M[b,a](0,s)   [PROVED, c=0 case]

These are equivalent restatements. We need a NEW ingredient.

REMAINING APPROACHES TO TRY:
1. Direct algebraic proof that M(r) is even in r
2. Find a matrix identity expressing M as a manifestly even function
3. Induction on n
4. Use the specific structure of the inclusion-exclusion sum
5. Relate to a known result in algebraic combinatorics
""")
