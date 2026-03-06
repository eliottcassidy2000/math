#!/usr/bin/env python3
"""
WHY does odd(B_b(W)) = r * sum_v M_W[v,b]?

REFORMULATION: Let W be a set of m vertices with distinguished element b.
Let t(i,j) = r + s(i,j), s antisymmetric.

B_b(W) = sum of Hamiltonian paths starting at b through W
       = sum_{sigma in S_{W\{b}}} prod_{i=0}^{m-2} t(sigma_i, sigma_{i+1})
         where sigma_0 = b

M_W[v,b] = transfer matrix entry (2-path-cover sum with inclusion-exclusion)

KEY IDENTITY: odd(B_b(W)) = r * sum_{v in W\{b}} M_W[v,b]

APPROACH 1: Use the symmetry M[v,b] = M[b,v] (ASSUMING even-r, which is
            what we're trying to prove -- CIRCULAR!)

APPROACH 2: Direct algebraic manipulation.
  B_b(W) = sum_sigma prod (r + s_{sigma_i, sigma_{i+1}})
  Expand each product: B_b = sum_{k=0}^{m-1} r^k * C_k
  where C_k = sum over paths, sum over choosing k edges to contribute r.

  odd(B_b) = sum_{k odd} r^k C_k

  Meanwhile, M_W[v,b] = sum_S (-1)^|S| E_v(S+{v}) B_b(R+{b})

  sum_v M[v,b] = sum_v sum_S (-1)^|S| E_v(S+{v}) B_b(R+{b})

APPROACH 3: Relate to the matrix (I - T) where T_{ij} = t(i,j).
  From lgv_two_path_analysis.py, we showed M[a,b] = cofactor(I-T, a, b)?
  Let me re-check this. If M = cofactor, then sum_v M[v,b] = sum_v cof(I-T,v,b)
  = the column-b sum of the cofactor matrix = related to adj(I-T) * e_b.
  And adj(I-T) * (I-T) = det(I-T) * I, so sum_v (I-T)_{v,w} * cof(I-T,v,b) = det(I-T) * delta_{w,b}.
  This is the cofactor expansion along column b.

  But is M[v,b] really the cofactor of (I-T)?

APPROACH 4: Path reversal identity.
  B_b(W; r, s) = E_b(W; r, -s) (reversing all paths)
  = E_b(W; -r, s) ... NO. Let's be careful.

  A path b -> v1 -> ... -> v_{m-1} has weight
  t(b,v1) * t(v1,v2) * ... * t(v_{m-2}, v_{m-1})
  = prod (r + s(v_i, v_{i+1}))

  Reversing: path v_{m-1} -> ... -> v1 -> b has weight
  t(v_{m-1}, v_{m-2}) * ... * t(v1, b)
  = prod (r + s(v_{i+1}, v_i))
  = prod (r - s(v_i, v_{i+1}))

  So B_b(W; r, s) = sum over sigma: prod(r + s) = E_b(W; r, -s)
  where E_b counts paths ENDING at b.

  Now E_b(W; r, -s) at -r: E_b(W; -r, -s) = paths ending at b with weight
  prod(-r - s) = (-1)^{m-1} prod(r + s) = (-1)^{m-1} E_b(W; r, s)

  So B_b(W; -r, s) = E_b(W; -r, -s) = (-1)^{m-1} E_b(W; r, s)

Actually, let me focus on the matrix approach. Let me verify the
cofactor connection first.

opus-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Poly, Matrix, det, eye

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
    vs = sorted(vertex_set)
    if len(vs) == 1: return 1
    total = 0
    for perm in permutations(vs):
        if perm[0] != start: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, vertex_set, a, b):
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
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1
        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1
        result += sign * ea * bb
    return expand(result)

def odd_part(expr, r):
    if expr == 0: return 0
    p = Poly(expr, r)
    d = p.as_dict()
    return expand(sum(c * r**deg[0] for deg, c in d.items() if deg[0] % 2 == 1))


# ============================================================
# STEP 1: Is M[v,b] the cofactor of (I-T)?
# ============================================================
print("=" * 70)
print("STEP 1: Cofactor relationship M[v,b] vs cof(I-T, v, b)")
print("=" * 70)

for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = list(range(m))
    b = 0  # Use 0 as distinguished

    # Build T matrix on W
    T = Matrix(m, m, lambda i, j: t(W[i], W[j]) if W[i] != W[j] else 0)
    I_T = eye(m) - T

    for v_idx in range(1, m):
        v = W[v_idx]
        M_vb = transfer_M(t, set(W), v, b)

        # Cofactor of (I-T) at position (v_idx, b_idx=0)
        # cof = (-1)^(v_idx + 0) * det(minor)
        minor = I_T.minor_submatrix(v_idx, 0)
        cof = expand((-1)**(v_idx) * det(Matrix(minor)))

        match = expand(M_vb - cof) == 0
        neg_match = expand(M_vb + cof) == 0
        print(f"  m={m}, v={v}, b={b}: M = cof? {match}, M = -cof? {neg_match}")
        if not match and not neg_match:
            print(f"    M = {M_vb}")
            print(f"    cof = {cof}")
            # Try other signs
            for sign in [1, -1]:
                for power in [0, 1]:
                    test = expand(M_vb - sign * (-1)**power * cof)
                    if test == 0:
                        print(f"    FOUND: M = {sign}*(-1)^{power}*cof")


# ============================================================
# STEP 2: Try cofactors of (T - I) instead
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Try various matrix formulations")
print("=" * 70)

for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = list(range(m))
    b = 0

    T = Matrix(m, m, lambda i, j: t(W[i], W[j]) if W[i] != W[j] else 0)

    # Try: M[v,b] = (-1)^{m-1} * cofactor of (T-I)
    T_I = T - eye(m)

    for v_idx in range(1, m):
        v = W[v_idx]
        M_vb = transfer_M(t, set(W), v, b)

        # Try all sign combinations
        minor = T_I.minor_submatrix(v_idx, 0)
        base_det = expand(det(Matrix(minor)))

        for overall_sign in [1, -1]:
            for cof_sign_power in range(4):
                test_val = expand(overall_sign * (-1)**cof_sign_power * base_det)
                if expand(M_vb - test_val) == 0:
                    print(f"  m={m}, v={v}: M = {overall_sign}*(-1)^{cof_sign_power} * det(minor(T-I,{v_idx},0))")
                    break

    # Try: cofactors of T itself
    for v_idx in range(1, min(3, m)):
        v = W[v_idx]
        M_vb = transfer_M(t, set(W), v, b)
        minor_T = T.minor_submatrix(v_idx, 0)
        det_minor_T = expand(det(Matrix(minor_T)))
        for s_val in [1, -1]:
            if expand(M_vb - s_val * det_minor_T) == 0:
                print(f"  m={m}: M[{v},{b}] = {s_val} * det(minor(T,{v_idx},0))")


# ============================================================
# STEP 3: The COLUMN SUM as a cofactor sum
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Column sum interpretation")
print("=" * 70)

# If M[v,b] is related to cofactors, then sum_v M[v,b] relates to
# a row/column of the adjugate.

for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = list(range(m))
    b = 0

    # Compute column sum of M
    col_sum = 0
    for v in W:
        if v == b: continue
        M_vb = transfer_M(t, set(W), v, b)
        col_sum = expand(col_sum + M_vb)

    # Compute B_b
    Bb = ham_paths_from(set(W), b, t)
    odd_Bb = odd_part(Bb, r)

    print(f"\nm={m}:")
    print(f"  col_sum = {col_sum}")
    print(f"  B_b = {Bb}")

    # IDENTITY: odd(B_b) = r * col_sum
    # Equivalently: (B_b(r) - B_b(-r)) / 2 = r * col_sum
    # i.e., B_b(r) - B_b(-r) = 2r * col_sum

    Bb_neg = expand(Bb.subs(r, -r))
    lhs = expand(Bb - Bb_neg)
    rhs = expand(2 * r * col_sum)
    print(f"  B_b(r) - B_b(-r) = {lhs}")
    print(f"  2r * col_sum = {rhs}")
    print(f"  Match: {expand(lhs - rhs) == 0}")

    # KEY INSIGHT: B_b(-r) = ?
    # B_b(r,s) = sum_sigma prod(r + s_{sigma})
    # B_b(-r,s) = sum_sigma prod(-r + s_{sigma})
    # = sum_sigma (-1)^{m-1} prod(r - s_{sigma})
    # = (-1)^{m-1} B_b(r, -s)
    # = (-1)^{m-1} E_b(r, s)  [by path reversal: B with -s = E]

    # So: B_b(r) - B_b(-r) = B_b(r) - (-1)^{m-1} E_b(r)
    # And the identity becomes:
    # B_b(r) - (-1)^{m-1} E_b(r) = 2r * col_sum

    # Let's verify E_b:
    def ham_paths_to(vertex_set, target, t_fn):
        vs = sorted(vertex_set)
        if len(vs) == 1: return 1
        total = 0
        for perm in permutations(vs):
            if perm[-1] != target: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t_fn(perm[i], perm[i+1])
            total += prod
        return expand(total)

    Eb = ham_paths_to(set(W), b, t)
    alt_rhs = expand(Bb - (-1)**(m-1) * Eb)
    print(f"  E_b = {Eb}")
    print(f"  B_b - (-1)^{{m-1}} E_b = {alt_rhs}")
    print(f"  = 2r * col_sum? {expand(alt_rhs - rhs) == 0}")

    # So the EQUIVALENT identity is:
    # B_b(W) + (-1)^m E_b(W) = 2r * sum_v M_W[v,b]
    # i.e., B_b + (-1)^m E_b = 2r * col_sum_b

    print(f"\n  EQUIVALENT: B_b + (-1)^m * E_b = 2r * col_sum")
    equiv = expand(Bb + (-1)**m * Eb - 2 * r * col_sum)
    print(f"  Check: {equiv == 0}")

    # For EVEN m: B_b + E_b = 2r * col_sum
    # For ODD m:  B_b - E_b = 2r * col_sum


# ============================================================
# STEP 4: What is B_b + (-1)^m E_b combinatorially?
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: B_b + (-1)^m E_b = ?")
print("=" * 70)

# B_b(W) = sum_{sigma: b->...} prod t(sigma_i, sigma_{i+1})
# E_b(W) = sum_{sigma: ...->b} prod t(sigma_i, sigma_{i+1})
#
# For a path sigma = (v_0, v_1, ..., v_{m-1}):
# B_b: v_0 = b, weight = prod_{i=0}^{m-2} t(v_i, v_{i+1})
# E_b: v_{m-1} = b, weight = prod_{i=0}^{m-2} t(v_i, v_{i+1})
#
# B_b + (-1)^m E_b:
# = sum_{sigma: b starts} prod t + (-1)^m sum_{sigma: b ends} prod t
#
# For each path sigma starting at b:
#   sigma = (b, v_1, ..., v_{m-1})
#   its reverse sigma^R = (v_{m-1}, ..., v_1, b) ends at b
#   weight of sigma^R = prod t(v_{i+1}, v_i) = prod (r - s(v_i, v_{i+1}))
#
# So E_b = sum_{sigma: b starts} prod(r - s_{edges})
# (Summing over the REVERSES of paths ending at b, which ARE paths starting at b)
#
# Wait, more carefully:
# E_b = sum over Ham paths ending at b = sum over sigma with sigma_{m-1} = b
# = sum over sigma^R with sigma^R_0 = b, weight = prod t(sigma^R_{m-1-i}, sigma^R_{m-i})
# ... this reversal sends E_b(r,s) to B_b(r,-s)

# So: E_b(W; r, s) = B_b(W; r, -s)
# (Reverse every path: a path ending at b becomes a path starting at b,
#  with all arc directions flipped, i.e., s -> -s)

# Therefore:
# B_b(r,s) + (-1)^m E_b(r,s) = B_b(r,s) + (-1)^m B_b(r,-s)
#
# Write B_b = sum_k r^k C_k(s) where C_k has parity (-1)^{m-1-k} in s
# (C_k involves choosing k edges for r, the other m-1-k edges contribute s terms)
# Actually C_k has mixed s-parity in general.
#
# B_b(r,-s) = sum_k r^k C_k(-s)
# Each C_k(-s): for a path of m-1 edges, choosing k to contribute r,
# the remaining m-1-k contribute s. So C_k(-s) = (-1)^{m-1-k} C_k(s)
# IF each s monomial in C_k has degree exactly m-1-k. YES! Each term
# in C_k is a product of exactly m-1-k s-values (one per non-r edge).
#
# So C_k(-s) = (-1)^{m-1-k} C_k(s)
# Therefore:
# B_b(r,-s) = sum_k r^k (-1)^{m-1-k} C_k(s) = (-1)^{m-1} sum_k (-r)^k C_k(s) ... wait
#
# No: B_b(r,-s) = sum_k r^k * (-1)^{m-1-k} C_k(s)
#
# B_b(r,s) + (-1)^m B_b(r,-s)
# = sum_k r^k C_k + (-1)^m sum_k r^k (-1)^{m-1-k} C_k
# = sum_k r^k C_k [1 + (-1)^m (-1)^{m-1-k}]
# = sum_k r^k C_k [1 + (-1)^{2m-1-k}]
# = sum_k r^k C_k [1 + (-1)^{-1-k}]  (since (-1)^{2m} = 1)
# = sum_k r^k C_k [1 - (-1)^k]
#
# For k even: 1 - 1 = 0
# For k odd: 1 - (-1) = 2
#
# So B_b + (-1)^m E_b = 2 * sum_{k odd} r^k C_k = 2 * odd(B_b)
#
# And the identity says: 2 * odd(B_b) = 2r * col_sum
# i.e., odd(B_b) = r * col_sum.  ✓ (This is just a consistency check.)

# BUT WAIT — this derivation shows:
# B_b(r,s) + (-1)^m B_b(r,-s) = 2 * odd_r(B_b)
# Which is TRIVIALLY TRUE from the definition of odd part!
# So the identity B_b + (-1)^m E_b = 2r * col_sum is the SAME as
# odd(B_b) = r * col_sum, just written differently.

# The SUBSTANTIVE content is: odd(B_b) = r * sum_v M[v,b]
# We need a DIFFERENT angle to prove this.

print("The reformulation B_b + (-1)^m E_b = 2*odd(B_b) is tautological.")
print("Need a different approach to prove odd(B_b) = r * col_sum.")

# ============================================================
# STEP 5: Row-sum / column-sum of M
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Row and column sums of M — algebraic structure")
print("=" * 70)

# What is sum_v M[v,b]?
# sum_{v != b} M_W[v,b]
# = sum_{v != b} sum_{S subset W\{v,b}} (-1)^|S| E_v(S+{v}) B_b(R+{b})
#
# Exchange sums: for each partition of W\{b} into {v} + S + R:
# sum_v sum_S (-1)^|S| E_v(S+{v}) B_b(R+{b})
#
# For a fixed set S' = S + {v} and R, where v ranges over S':
# = sum_{S' subset W\{b}} sum_{v in S'} (-1)^{|S'|-1} E_v(S') B_b(R+{b})
#   where R = (W\{b}) \ S'
#
# = sum_{S' subset W\{b}} (-1)^{|S'|-1} B_b(R+{b}) * sum_{v in S'} E_v(S')

# KEY: sum_{v in S'} E_v(S') = sum of all Ham paths through S' (any endpoint)
# For each Ham path through S', it ends at SOME v in S', so this is the total
# Hamiltonian path sum through S' = T(S') (total Ham path weight, sum over ALL orderings)

# Wait: E_v(S') = sum of Ham paths through S' ending at v.
# So sum_v E_v(S') = sum of ALL Hamiltonian paths through S' (regardless of endpoint)
# = total Hamiltonian path weight of S'.
#
# For a SINGLE vertex S' = {v}: E_v({v}) = 1, so total = 1.
# For S' = {u,v}: E_u({u,v}) = t(v,u), E_v({u,v}) = t(u,v),
# total = t(v,u) + t(u,v) = 2r.
# For |S'| = 3: total = sum over all 3! = 6 permutations of products of 2 edges,
# but grouped by endpoint.

for m in [3, 4, 5]:
    r_sym, sv, s, t_fn = setup(m)
    W = list(range(m))

    for size in range(1, m+1):
        from itertools import combinations
        for S_prime in [list(range(size))]:  # just test one
            total_ham = 0
            for perm in permutations(S_prime):
                prod = 1
                for i in range(len(perm)-1):
                    prod *= t_fn(perm[i], perm[i+1])
                total_ham += prod
            total_ham = expand(total_ham)
            print(f"  m={m}, |S'|={size}, S'={S_prime}: total Ham paths = {total_ham}")

print("\n  Pattern: total Ham paths through S' of size k = k! * r^{k-1}?")
print("  Check: |S'|=1: 1 = 1!*r^0 ✓")
print("  Check: |S'|=2: 2r = 2!*r^1 ✓")
print("  But |S'|=3: total ≠ 6*r^2 in general (has s-terms)")

# So sum_v E_v(S') ≠ |S'|! * r^{|S'|-1} in general.
# BUT: what about the LEADING r-term?
# The coefficient of r^{k-1} in total Ham paths through k vertices:
# Each path has k-1 edges, each contributing r. So r^{k-1} coefficient
# = number of Ham paths = k!.
# So total Ham paths = k! * r^{k-1} + lower order terms.

# Actually wait. Each permutation of k elements gives a Hamiltonian path
# through those k vertices. And the r^{k-1} coefficient of each path is 1
# (all edges contribute r). So the r^{k-1} coeff of total is k!.

# Now: col_sum = sum_{S'} (-1)^{|S'|-1} B_b(R+{b}) * T(S')
# where T(S') = total Ham path weight through S'.

# This is interesting. Can we relate T(S') to something simpler?

# Actually, T(S') = sum_v E_v(S').
# And we know E_v(S') has a specific r-parity structure.

# Let me check: is T(S') always a multiple of r?
# |S'| = 1: T = 1 (NOT a multiple of r)
# |S'| = 2: T = 2r (multiple of r)
# |S'| >= 2: T = k! r^{k-1} + ...

# Hmm, for |S'| = 1, T = 1. So the formula
# col_sum = sum_{S'} (-1)^{|S'|-1} B_b(R+{b}) T(S')
# has both T=1 terms (when |S'|=1) and T with r-factors.

# Let me compute this explicitly for small cases.

print("\n" + "=" * 70)
print("STEP 6: Explicit col_sum decomposition")
print("=" * 70)

for m in [3, 4]:
    r_sym, sv, s, t_fn = setup(m)
    W = set(range(m))
    b = 0

    print(f"\nm={m}, b={b}:")

    # col_sum = sum_v M[v,b] = sum_{S' subset W\{b}} (-1)^{|S'|-1} T(S') B_b(R+{b})
    # where R = (W\{b}) \ S', and S' is non-empty

    W_minus_b = sorted(W - {b})

    total = 0
    for mask in range(1, 1 << len(W_minus_b)):  # non-empty S'
        S_prime = [W_minus_b[i] for i in range(len(W_minus_b)) if mask & (1 << i)]
        R = [w for w in W_minus_b if w not in S_prime]
        R_b = set(R) | {b}
        sign = (-1)**(len(S_prime) - 1)

        # T(S') = total Ham path weight
        T_S = 0
        for perm in permutations(S_prime):
            prod = 1
            for i in range(len(perm)-1):
                prod *= t_fn(perm[i], perm[i+1])
            T_S += prod
        T_S = expand(T_S)

        # B_b(R+{b})
        Bb_R = ham_paths_from(R_b, b, t_fn)

        contrib = expand(sign * T_S * Bb_R)
        print(f"  S'={S_prime}, R={R}: sign={sign}, T(S')={T_S}, B_b(R+b)={Bb_R}, contrib={contrib}")
        total = expand(total + contrib)

    print(f"  Total = {total}")

    # Compare with direct col_sum
    direct_col = 0
    for v in W_minus_b:
        direct_col = expand(direct_col + transfer_M(t_fn, W, v, b))
    print(f"  Direct col_sum = {direct_col}")
    print(f"  Match: {expand(total - direct_col) == 0}")


# ============================================================
# STEP 7: Focus on the identity at r=0
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: The identity at r=0: d/dr B_b |_{r=0} = col_sum(0)")
print("=" * 70)

# At r=0: t(i,j) = s(i,j)
# B_b(W; 0, s) = sum_sigma prod s(sigma_i, sigma_{i+1})
# This is the "pure skew" permanent-like quantity.
#
# d/dr B_b |_{r=0}: differentiate each factor (r + s) and set r=0:
# = sum_sigma sum_{edge e} [prod_{e' != e} s(e')] * 1
# = sum_sigma (m-1) * [1/product of one edge] * [product of all edges]
# ... no, more carefully:
#
# B_b = sum_sigma prod_{i=0}^{m-2} (r + s_i)
# d/dr B_b = sum_sigma sum_{j=0}^{m-2} prod_{i != j} (r + s_i)
# At r=0: = sum_sigma sum_j prod_{i != j} s_i
#          = sum_sigma [sum_j prod_{i != j} s_i]
#          = sum_sigma [prod_all s_i] * [sum_j 1/s_j]  (when no s_j = 0)

# This is related to the "edge deletion" operation: for each edge in the path,
# remove it and take the product of s-weights of the remaining edges.

# Meanwhile, col_sum(0) = sum_v M_W[v,b](0)
# = sum_v M[v,b] evaluated with t -> s (pure skew tournament).

# The identity says these are equal. This is a purely antisymmetric identity!

print("At r=0, the identity becomes:")
print("  sum_sigma sum_j prod_{i!=j} s_{sigma_i,sigma_{i+1}} = sum_v M[v,b]|_{r=0}")
print("  where sigma ranges over Ham paths from b.")
print("  This is a statement about skew-symmetric matrices.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
