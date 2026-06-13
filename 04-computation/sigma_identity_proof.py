#!/usr/bin/env python3
"""
PROVING THE SIGMA IDENTITY: the last piece for even r-powers.

NEEDED:
  |W| odd:  Sigma(W) = 0
  |W| even: r * Sigma(W) = T(W)

where Sigma = sum_{a!=w} M[a,w], T = total Hamiltonian path weight.

APPROACH: Use the TRANSFER MATRIX recurrence to express Sigma(W)
in terms of Sigma(W') for W' = W\{v} (deleting a vertex).

The column recurrence gives:
  M_W[a,b] = (-1)^{m-2} E_a(W\{b}) + sum_w t(b,w) M_{W\{b}}[a,w]

For Sigma, sum over ALL (a,b) pairs:
  Sigma = sum_{a!=b} M[a,b]
  = sum_{a!=b} [(-1)^{m-2} E_a(W\{b}) + sum_w t(b,w) M_{W\{b}}[a,w]]

This gives a recurrence for Sigma.

opus-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Poly

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

def total_ham(vertex_set, t_fn):
    vs = sorted(vertex_set)
    if len(vs) <= 1: return 1
    total = 0
    for perm in permutations(vs):
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)


# ============================================================
# APPROACH 1: Express Sigma using column recurrence, fixing b
# ============================================================
print("=" * 70)
print("APPROACH 1: Sigma via column recurrence (fixing b)")
print("=" * 70)

# Sigma(W) = sum_{b in W} col_sum(b)
#           = sum_{b in W} sum_{a != b} M[a,b]
#
# Using column recurrence for each M[a,b]:
# M[a,b] = (-1)^{m-2} E_a(W\{b}) + sum_{w != a,b} t(b,w) M_{W\{b}}[a,w]
#
# Sum over a != b (fixed b):
# col_sum(b) = (-1)^{m-2} sum_{a != b} E_a(W\{b}) + sum_{a != b} sum_w t(b,w) M_{W\{b}}[a,w]
# = (-1)^{m-2} T(W\{b}) + sum_w t(b,w) sum_{a != b, a != w} M_{W\{b}}[a,w]
# = (-1)^{m-2} T(W\{b}) + sum_w t(b,w) cs_{W\{b}}(w)
#
# where cs_{W\{b}}(w) = sum_{a in W\{b}\{w}} M_{W\{b}}[a,w] = column-w sum on W\{b}.

# Now sum over b:
# Sigma(W) = sum_b [(-1)^{m-2} T(W\{b}) + sum_w t(b,w) cs_{W\{b}}(w)]
# = (-1)^{m-2} sum_b T(W\{b}) + sum_b sum_w (r + s_{bw}) cs_{W\{b}}(w)
# = (-1)^{m-2} sum_b T(W\{b}) + r sum_b sum_w cs_{W\{b}}(w) + sum_b sum_w s_{bw} cs_{W\{b}}(w)
#
# Now: sum_w cs_{W\{b}}(w) = Sigma(W\{b})  [total off-diag sum on W\{b}]
# So: r sum_b Sigma(W\{b}) appears.
#
# And sum_b T(W\{b}) = sum_b (total Ham path weight through W\{b}).
# For each b, T(W\{b}) is the total on m-1 vertices.
# sum_b T(W\{b}): for each (b, ordered path through W\{b}), sum weight.
# Equivalently: for each path through m-1 vertices, it appears once for the missing vertex.
# If we let T_k(W) = sum over subsets of size k of total Ham path, then
# sum_b T(W\{b}) = sum_{S of size m-1} T(S) = "total paths on m-1 vertices" summed.

# Let me compute numerically.
for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))

    Sigma = 0
    for a in sorted(W):
        for w in sorted(W):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W, a, w))

    sum_T = 0
    for b in sorted(W):
        Wb = W - {b}
        sum_T = expand(sum_T + total_ham(Wb, t))

    sum_Sigma_sub = 0
    for b in sorted(W):
        Wb = W - {b}
        S_sub = 0
        for a in sorted(Wb):
            for w in sorted(Wb):
                if a == w: continue
                S_sub = expand(S_sub + transfer_M(t, Wb, a, w))
        sum_Sigma_sub = expand(sum_Sigma_sub + S_sub)

    # sum_b sum_w s_{bw} cs_{W\{b}}(w)
    s_weighted_total = 0
    for b in sorted(W):
        Wb = W - {b}
        for w in sorted(Wb):
            cs_w = 0
            for a in sorted(Wb):
                if a == w: continue
                cs_w = expand(cs_w + transfer_M(t, Wb, a, w))
            s_weighted_total = expand(s_weighted_total + s(b, w) * cs_w)

    # Check: Sigma = (-1)^{m-2} sum_T + r * sum_Sigma_sub + s_weighted
    rhs = expand((-1)**(m-2) * sum_T + r * sum_Sigma_sub + s_weighted_total)

    print(f"\nm={m}:")
    print(f"  Sigma(W) = {Sigma}")
    print(f"  (-1)^(m-2) sum_b T(W\\b) = {expand((-1)**(m-2) * sum_T)}")
    print(f"  r * sum_b Sigma(W\\b) = {expand(r * sum_Sigma_sub)}")
    print(f"  sum_b sum_w s_bw cs_w = {s_weighted_total}")
    print(f"  RHS = {rhs}")
    print(f"  Match: {expand(Sigma - rhs) == 0}")

    # What is s_weighted_total?
    # sum_b sum_w s_{bw} cs_{W\{b}}(w)
    # = sum_{b!=w} s_{bw} cs_{W\{b}}(w)
    # cs_{W\{b}}(w) = sum_{a in W\{b,w}} M_{W\{b}}[a,w]

    T_W = total_ham(W, t)
    print(f"  T(W) = {T_W}")

    # Does sum_b T(W\{b}) have a pattern?
    print(f"  sum_b T(W\\b) = {sum_T}")


# ============================================================
# APPROACH 2: M[a,b] + M[b,a] = known even quantity
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 2: Analyze M[a,b]+M[b,a] directly")
print("=" * 70)

# M[a,b] + M[b,a] has even r-powers. What IS it?
# M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# M[b,a] = sum_S' (-1)^|S'| E_b(S'+b) B_a(R'+a)
# where S' ranges over subsets of W\{a,b}, R' = W\{a,b}\S'.
# Setting S' = R (reverse the partition!):
# M[b,a] = sum_R (-1)^|R| E_b(R+b) B_a(S+a)
#         = sum_S (-1)^{m-2-|S|} E_b((W\{a,b}\S)+b) B_a(S+a)
#
# Wait, if S' = W\{a,b}\S = R, then |S'| = m-2-|S|. And R' = S.
# So M[b,a] = sum_S (-1)^{m-2-|S|} E_b(R+b) B_a(S+a)
# = (-1)^{m-2} sum_S (-1)^{-|S|} (-1)^{-(-|S|)}... let me be more careful.
# (-1)^{m-2-|S|} = (-1)^{m-2} * (-1)^{-|S|} = (-1)^{m-2} * (-1)^{|S|}.
# So M[b,a] = (-1)^{m-2} sum_S (-1)^{|S|} E_b(R+b) B_a(S+a)
#
# Compare with M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b).
# So: M[a,b] + M[b,a] = sum_S (-1)^|S| [E_a(S+a) B_b(R+b) + (-1)^{m-2} E_b(R+b) B_a(S+a)]
#
# For m even: M[a,b] + M[b,a] = sum_S (-1)^|S| [E_a B_b + E_b B_a] (each product).
# This is SYMMETRIC in the partition! E_a through S, B_b through R, AND E_b through R, B_a through S.

# For m even, each term is:
# (-1)^|S| [E_a(S+a) B_b(R+b) + E_b(R+b) B_a(S+a)]
# This is the "symmetrized" 2-path-cover.

# For m odd: M[a,b] + M[b,a] = sum_S (-1)^|S| [E_a B_b - E_b B_a]
# (antisymmetrized!)

# Let me verify:
for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))
    a, b = 0, 1

    Mab = transfer_M(t, W, a, b)
    Mba = transfer_M(t, W, b, a)

    print(f"\nm={m}:")
    print(f"  M[{a},{b}] = {Mab}")
    print(f"  M[{b},{a}] = {Mba}")
    print(f"  Sum = {expand(Mab + Mba)}")
    print(f"  Diff = {expand(Mab - Mba)}")

    # Check the formula for M+M
    U = sorted(W - {a, b})
    formula_sum = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R_list = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        # E_a(S+a)
        Sa = set(S) | {a}
        ea = 0
        for p in permutations(sorted(Sa)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1): prod *= t(p[i], p[i+1])
            ea += prod
        if len(Sa) == 1: ea = 1

        # B_b(R+b)
        Rb = set(R_list) | {b}
        bb_val = 0
        for p in permutations(sorted(Rb)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1): prod *= t(p[i], p[i+1])
            bb_val += prod
        if len(Rb) == 1: bb_val = 1

        # E_b(R+b) and B_a(S+a)
        eb = 0
        for p in permutations(sorted(Rb)):
            if p[-1] != b: continue
            prod = 1
            for i in range(len(p)-1): prod *= t(p[i], p[i+1])
            eb += prod
        if len(Rb) == 1: eb = 1

        ba_val = 0
        for p in permutations(sorted(Sa)):
            if p[0] != a: continue
            prod = 1
            for i in range(len(p)-1): prod *= t(p[i], p[i+1])
            ba_val += prod
        if len(Sa) == 1: ba_val = 1

        formula_sum += sign * (ea * bb_val + (-1)**(m-2) * eb * ba_val)

    formula_sum = expand(formula_sum)
    print(f"  Formula = {formula_sum}")
    print(f"  Match sum: {expand(formula_sum - (Mab + Mba)) == 0}")


# ============================================================
# APPROACH 3: Sigma from 2-path-cover summed over all endpoints
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 3: What is Sigma combinatorially?")
print("=" * 70)

# Sigma = sum_{a!=w} M[a,w]
# = sum_{a!=w} sum_S (-1)^|S| E_a(S+a) B_w(R+w)
# where S ⊆ W\{a,w}, R = W\{a,w}\S.
#
# Rewrite: for each PARTITION of W into S'+a' (with endpoint a') and R'+w' (with start w'):
# where |S'|+|R'|+2 = m (a and w are separate from S and R in W\{a,w}).
#
# Actually, S ⊆ W\{a,w}, so S and R partition W\{a,w}.
# S+a has |S|+1 vertices, R+w has m-2-|S|+1 = m-1-|S| vertices.
#
# Total vertices covered: (|S|+1) + (m-1-|S|) = m. Good.
#
# Sigma = sum over all 2-path-covers of W where:
#   - P1 ends at some a, goes through S+{a}
#   - P2 starts at some w, goes through R+{w}
#   - sign = (-1)^{|P1|-1} = (-1)^|S|
#
# This is the TOTAL SIGNED 2-PATH-COVER sum, summed over ALL endpoint pairs.

# Now: what is this in terms of T(W)?
# T(W) = sum of all Hamiltonian paths through W (any order) = sum_sigma prod t(sigma_i, sigma_{i+1}).
# T(W) involves single paths (1-path-covers).
# Sigma involves 2-path-covers with inclusion-exclusion.

# IDEA: Is there a matrix identity connecting Sigma to T?
# If M were a matrix with M_ij = cofactor(A, i, j) / det(A), then sum of entries = ?
# For the adjugate: adj(A) has entries = cofactors, and adj(A) * A = det(A) * I.
# sum of all entries of adj(A) = 1^T adj(A) 1 = (1^T adj(A)) * 1.
# From A * adj(A) = det(A) I: adj(A) = det(A) A^{-1}.
# So sum entries = det(A) * 1^T A^{-1} 1.
# Not obviously helpful since M is NOT adj(A).

# INSTEAD: Let me check if Sigma / T is a CONSTANT for small cases.
print("Sigma / T ratio:")
for m in [2, 4, 6]:
    r, sv, s, t = setup(m)
    W = set(range(m))

    Sigma = 0
    for a in sorted(W):
        for w in sorted(W):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W, a, w))

    T_W = total_ham(W, t)
    # r * Sigma should equal T for |W| even
    print(f"  |W|={m}: r*Sigma = T? {expand(r * Sigma - T_W) == 0}")
    if m == 2:
        print(f"    Sigma = {Sigma}, T = {T_W}")


# ============================================================
# APPROACH 4: Check r*Sigma = T as a STANDALONE identity
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 4: r*Sigma(W) = T(W) for even |W| (standalone)")
print("=" * 70)

# For even |W|:
# r * sum_{a!=w} M[a,w] = T(W)
#
# LHS = r * sum_{a!=w} sum_S (-1)^|S| E_a(S+a) B_w(R+w)
# RHS = sum_sigma prod t(sigma_i, sigma_{i+1})
#
# Can we write T(W) as r * (some inclusion-exclusion)?
# T(W) = sum_sigma prod_{i=0}^{m-2} (r + s_i)
# = sum_{k=0}^{m-1} r^k * T_k
# where T_k is the k-th coefficient (sum over choosing k edges for r).
#
# Since T has parity (-1)^{m-1}: for m even, parity is -1 (odd r-powers only).
# So T = r * (even function) when m is even.
# Specifically, T = r * [T_1 + T_3 r^2 + T_5 r^4 + ...]
# = r * sum_j T_{2j+1} r^{2j}
#
# So T/r is well-defined and even in r for even m.
# We need: Sigma = T/r.

# Let me compute T/r and compare with Sigma for m=4:
m = 4
r, sv, s, t = setup(m)
W = set(range(m))

T_W = total_ham(W, t)
p_T = Poly(T_W, r)
# T/r: divide each coefficient by r
T_over_r = expand(sum(c * r**(d[0]-1) for d, c in p_T.as_dict().items()))

Sigma = 0
for a in sorted(W):
    for w in sorted(W):
        if a == w: continue
        Sigma = expand(Sigma + transfer_M(t, W, a, w))

print(f"m=4: T/r = {T_over_r}")
print(f"m=4: Sigma = {Sigma}")
print(f"m=4: T/r = Sigma? {expand(T_over_r - Sigma) == 0}")


# ============================================================
# APPROACH 5: Direct check of r*Sigma = T at DEFINITION level
# ============================================================
print("\n" + "=" * 70)
print("APPROACH 5: What IS T/r for even m?")
print("=" * 70)

# T(W) = sum over all orderings of W, product of edge weights.
# For |W| = m (even), T has only odd r-powers.
# T/r = sum over all orderings, sum_j [product of m-2-j s-values] * r^{2j}
# ... where j = number of r-edges minus 1 (choosing from m-1 edges).
# Actually: T = sum_k r^k T_k, where k ranges over odd values only.
# T_k = coefficient of r^k = sum_sigma C(sigma, k)
# where C(sigma, k) = sum over choosing k edges out of m-1 to contribute r.

# For k=1 (r^1 coefficient):
# T_1 = sum_sigma (m-1) * prod_{all edges except one} s_e
# ... no. More carefully:
# coefficient of r^k = sum_sigma sum_{|E|=k subsets of edges} prod_{e not in E} s_e
# = sum_sigma e_{m-1-k}(s_{e_1}, ..., s_{e_{m-1}})  [elementary symmetric poly in s-values]

# HMMMM. For m=2: T = 2r. T/r = 2. Sigma = 2 for |W|=2.
# For m=4: T/r = Sigma (verified above).

# What is the DEGREE 0 term of T/r?
# Coefficient of r^0 in T/r = coefficient of r^1 in T.
# Coefficient of r^1 in T = sum_sigma sum_j prod_{i!=j} s_{e_i}
# = sum_sigma e_{m-2}(s_{e_1}, ..., s_{e_{m-1}})
# = sum_sigma [prod_all s_e] * [sum_j 1/s_{e_j}]  (when s_e != 0)

# And Sigma(0) = sum_{a!=w} M[a,w]|_{r=0}.
# M[a,w] at r=0 uses t(i,j) = s(i,j), which is a skew-symmetric "tournament."

# OK, I think the key insight needed is:
# Sigma = T/r for even m is a KNOWN identity in the theory of
# transfer matrices / 2-path-covers.

# Actually, let me check: is Sigma = (m-1) * Sigma_smaller or something?
# For m=2: Sigma = 2
# For m=4: Sigma has degree 2 in r (leading term 24r^2)
#   T has degree 3 in r (leading term 24r^3). T/r leading = 24r^2. ✓

# The leading coefficient of T is m! r^{m-1} (all edges contribute r).
# So T/r leading = m! r^{m-2}.
# Sigma leading: the leading r-power of Sigma should be m! r^{m-2} as well.
# This is the highest power. For m=4: 4! r^2 = 24 r^2. ✓

print("For even m: Sigma(W) = T(W)/r")
print("This is the identity we need to prove.")
print()
print("T(W) = total signed Hamiltonian path weight through W")
print("Sigma(W) = total off-diagonal sum of transfer matrix M")
print()
print("Both are computable from the tournament weights.")
print()
print("The identity T = r*Sigma can potentially be proved by:")
print("1. Expanding both in terms of path products")
print("2. Using the path reversal symmetry")
print("3. Direct algebraic argument")

print("\n" + "=" * 70)
print("APPROACH 6: Check if T(W) = r * Sigma(W) relates to d/dr")
print("=" * 70)

# We showed earlier: d/dr B_b|_{r=0} = col_sum_b(0).
# What about d/dr T|_{r=0}?

for m in [2, 3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))

    T_W = total_ham(W, t)
    from sympy import diff
    dT = expand(diff(T_W, r))
    dT_0 = expand(dT.subs(r, 0))

    Sigma = 0
    for a in sorted(W):
        for w in sorted(W):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W, a, w))
    Sigma_0 = expand(Sigma.subs(r, 0))

    # (m-1) * Sigma at r=0?
    print(f"\nm={m}: dT/dr|_0 = {dT_0}, Sigma(0) = {Sigma_0}")
    print(f"  dT/dr|_0 = (m-1)*Sigma(0)? {expand(dT_0 - (m-1)*Sigma_0) == 0}")
    print(f"  dT/dr|_0 = Sigma(0)? {expand(dT_0 - Sigma_0) == 0}")

    # Actually: T = sum_sigma prod (r + s_e)
    # dT/dr = sum_sigma sum_j prod_{i!=j} (r + s_{e_i})
    # This is a sum over (sigma, j) pairs: delete one edge from path sigma.
    # At r=0: = sum_sigma sum_j prod_{i!=j} s_{e_i}
    #
    # If T = r * Sigma (for even m):
    # dT/dr = Sigma + r * dSigma/dr
    # At r=0: dT/dr|_0 = Sigma(0).

    # For odd m: T is even, so T = Sigma_0 + Sigma_2 r^2 + ...
    # dT/dr = 2 Sigma_2 r + ... At r=0: 0.
    # Meanwhile Sigma = 0 for odd m.
    # So dT/dr|_0 = 0 and Sigma(0) = 0. Both zero. ✓


print("\n" + "=" * 70)
print("SUMMARY OF PROOF STATUS")
print("=" * 70)

print("""
PROVEN (computationally verified n=2,...,7):
  1. KEY IDENTITY: odd(B_b(W)) = r * col_sum(b)
  2. COLUMN RECURRENCE: M[a,b] = (-1)^{m-2} E_a(W') + sum t(b,w) M'[a,w]
  3. COL_SUM RECURRENCE: cs_W(b) = alpha*T(W') + sum s_bw cs'(w)
  4. SIGMA IDENTITIES: Sigma=0 for odd, r*Sigma=T for even

PROOF CHAIN (each step verified computationally):
  Base: M even for |W| <= 3.
  Inductive step (|W| = m >= 4):
    a) By induction, M on size < m has even r-powers (hence symmetric).
    b) Column recurrence (ALGEBRAIC IDENTITY from definition).
    c) Sum over a -> col_sum decomposition.
    d) KEY IDENTITY at size m-1 (inductive hypothesis).
    e) Sigma identity at size m-1.
    f) => KEY IDENTITY at size m => even r-powers at size m. QED.

REMAINING GAP:
  Algebraic proof of either:
    - The KEY IDENTITY directly, OR
    - The SIGMA IDENTITY (r*Sigma = T for even m, Sigma = 0 for odd m)
  These are equivalent given the column recurrence.

NOTE: Sigma having even r-powers follows from the PROVEN identity
  M(-r) = M^T (re-indexing). This is not circular — it uses T^op,
  not even r-powers. But going from "even r-powers" to "Sigma = 0
  for odd m" requires MORE than just parity.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
