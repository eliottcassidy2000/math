#!/usr/bin/env python3
"""
INDUCTIVE PROOF OF KEY IDENTITY via col_sum recurrence.

DERIVATION: We want to prove B_b(W) - (-1)^{m-1} E_b(W) = 2r * col_sum_W(b)
where m = |W|, col_sum_W(b) = sum_{v != b} M_W[v,b].

Decompose:
  B_b(W) = sum_v t(b,v) B_v(W\{b})  [first step of path from b]
  E_b(W) = sum_v t(v,b) E_v(W\{b})  [last step of path to b]

  LHS = sum_v [(r+s_bv) B_v - (-1)^{m-1}(r-s_bv) E_v]

Using -(-1)^{m-1} = (-1)^m = (-1)^{m-2}:
  = r * sum_v [B_v + (-1)^{m-2} E_v] + sum_v s_bv [B_v - (-1)^{m-2} E_v]

By INDUCTION (key identity at size m-1):
  B_v(W') - (-1)^{m-2} E_v(W') = 2r * cs_v

where cs_v = col_sum_{W'}(v), W' = W\{b}.

Also: sum_v E_v(W') = sum_v B_v(W') = T(W') [total Ham path weight].

So sum_v [B_v + (-1)^{m-2} E_v] = [1 + (-1)^{m-2}] T(W')
  = {2*T if m even, 0 if m odd}

Therefore:
  LHS = r * [1+(-1)^{m-2}] T(W') + 2r * sum_v s_bv * cs_v

Setting LHS = 2r * col_sum_W(b):
  col_sum_W(b) = alpha * T(W') + sum_v s_bv * cs_v

where alpha = [1+(-1)^{m-2}]/2 = {1 if m even, 0 if m odd}.

THIS IS THE RECURRENCE. If it holds, the induction goes through.

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
    """T(W) = total Hamiltonian path weight through W (all orderings)."""
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
# TEST THE RECURRENCE
# ============================================================
print("=" * 70)
print("COL_SUM RECURRENCE TEST")
print("col_sum_W(b) = alpha*T(W') + sum_v s_bv * col_sum_{W'}(v)")
print("where alpha = 1 if |W| even, 0 if |W| odd, W' = W\\{b}")
print("=" * 70)

for m in range(2, 7):
    # Need enough symbols
    n_syms = max(m, 6)
    r, sv, s, t = setup(n_syms)
    W = set(range(m))
    b = 0
    W_prime = W - {b}

    # LHS: col_sum_W(b) = sum_{v != b} M_W[v,b]
    col_sum_W = 0
    for v in sorted(W - {b}):
        col_sum_W = expand(col_sum_W + transfer_M(t, W, v, b))

    # RHS: alpha * T(W') + sum_v s_bv * cs_v
    alpha = 1 if m % 2 == 0 else 0
    T_Wp = total_ham(W_prime, t)

    s_weighted = 0
    for v in sorted(W_prime):
        cs_v = 0
        for u in sorted(W_prime - {v}):
            cs_v = expand(cs_v + transfer_M(t, W_prime, u, v))
        s_weighted = expand(s_weighted + s(b, v) * cs_v)

    rhs = expand(alpha * T_Wp + s_weighted)

    match = expand(col_sum_W - rhs) == 0
    print(f"\nm={m}: alpha={alpha}")
    if m <= 5:
        print(f"  col_sum_W(b) = {col_sum_W}")
        print(f"  alpha*T(W') = {expand(alpha * T_Wp)}")
        print(f"  sum s_bv cs_v = {s_weighted}")
        print(f"  RHS = {rhs}")
    print(f"  RECURRENCE HOLDS: {match}")

    if not match:
        print(f"  DIFF = {expand(col_sum_W - rhs)}")


# ============================================================
# FULL INDUCTIVE PROOF CHAIN
# ============================================================
print("\n" + "=" * 70)
print("FULL INDUCTIVE PROOF STRUCTURE")
print("=" * 70)

print("""
THEOREM: M_W[a,b] has only even powers of r, for all sets W, all a,b in W.

PROOF BY STRONG INDUCTION ON |W|.

Base case |W| = 2: M[a,b] = 1. Even. ✓
Base case |W| = 3: M[a,b] = s_{ac} + s_{bc} where c is the third vertex.
  This has r-degree 0. Even. ✓

Inductive step: Assume all M_{W'} with |W'| < |W| have even r-powers.
Need: M_W[a,b] has even r-powers.

Step 1: Row recurrence.
  M_W[a,b] = B_b(W\\{a}) - sum_v t(v,a) M^{(a)}[v,b]
  = B_b(W\\{a}) - r * col_sum^{(a)} - sum_v s_{va} M^{(a)}[v,b]

  where col_sum^{(a)} = sum_v M^{(a)}[v,b] and W^{(a)} = W\\{a}.

Step 2: By induction, M^{(a)} has even r-powers.
  So sum_v s_{va} M^{(a)}[v,b] is even in r.
  And r * col_sum^{(a)} is odd in r.

Step 3: odd(M_W[a,b]) = odd(B_b(W\\{a})) - r * col_sum^{(a)}.

Step 4: KEY IDENTITY: odd(B_b(W')) = r * col_sum_{W'}(b).
  This says: odd(B_b(W\\{a})) = r * col_sum^{(a)}.
  Therefore: odd(M_W[a,b]) = 0. QED.

PROOF OF KEY IDENTITY by strong induction on |W'| = m:

Base: m=1: B_b({b}) = 1 (even), col_sum = 0. odd(1) = 0 = r*0. ✓
Base: m=2: B_b({b,v}) = r + s_{bv}. odd = r. col_sum = M[v,b] = 1. r*1 = r. ✓

Inductive step: Assume key identity for all sets of size < m.

Reformulation: B_b(W) - (-1)^{m-1} E_b(W) = 2r * col_sum_W(b).

Decompose using first/last step of paths:
  B_b(W) = sum_v t(b,v) B_v(W')    [W' = W\\{b}]
  E_b(W) = sum_v t(v,b) E_v(W')

  LHS = r * [1+(-1)^{m-2}] T(W') + 2r * sum_v s_{bv} cs_v

where cs_v = col_sum_{W'}(v) (by inductive hypothesis at size m-1).

So: col_sum_W(b) = alpha * T(W') + sum_v s_{bv} * cs_v
where alpha = [1+(-1)^{m-2}]/2.

THIS RECURRENCE IS WHAT WE NEED TO VERIFY INDEPENDENTLY.
If it holds from the DEFINITION of M, the proof is complete.
""")


# ============================================================
# ATTEMPT TO PROVE THE RECURRENCE FROM DEFINITION
# ============================================================
print("=" * 70)
print("PROVING THE RECURRENCE FROM DEFINITION OF M")
print("=" * 70)

# col_sum_W(b) = sum_{v in W\\{b}} M_W[v,b]
# Using the row recurrence for each M_W[v,b]:
# M_W[v,b] = B_b(W\\{v}) - sum_{u in W\\{v,b}} t(u,v) M_{W\\{v}}[u,b]
#
# Sum over v:
# col_sum_W(b) = sum_v B_b(W\\{v}) - sum_v sum_{u!=v,b} t(u,v) M_{W\\{v}}[u,b]
#
# First term: sum_v B_b(W\\{v}) = sum_{v in W\\{b}} B_b(W\\{v})
# This counts pairs (v, Ham path from b through W\\{v}).
# Each such path has m-2 edges through m-1 vertices (W minus v).
#
# Second term: sum_v sum_{u!=v,b} (r + s_{uv}) M^{(v)}[u,b]
# = r * sum_v sum_u M^{(v)}[u,b] + sum_v sum_u s_{uv} M^{(v)}[u,b]
# = r * sum_v col_sum_{W\\{v}}(b) + sum_v sum_u s_{uv} M^{(v)}[u,b]
#
# Wait, col_sum_{W\\{v}}(b) = sum_{u in (W\\{v})\\{b}} M_{W\\{v}}[u,b]
# = sum_{u in W\\{v,b}} M_{W\\{v}}[u,b]
# Yes!
#
# So: col_sum_W(b) = sum_v B_b(W\\{v}) - r * sum_v cs_{W\\{v}}(b) - sum_v sum_u s_{uv} M^{(v)}[u,b]
#
# The recurrence we want to prove is:
# col_sum_W(b) = alpha * T(W') + sum_v s_{bv} * cs_{W'}(v)
# where W' = W\\{b}.
#
# These are DIFFERENT: the first expansion deletes v (for each v != b),
# the recurrence deletes b.

print("Need different approach. Let me try the COLUMN recurrence instead.")
print()

# COLUMN recurrence: expand M[v,b] by the B_b path.
# M_W[v,b] = sum_S (-1)^|S| E_v(S+{v}) B_b(R+{b})
# Expand B_b(R+{b}) by its first step:
# B_b(R+{b}) = sum_{w in R} t(b,w) B_w((R\\{w})+{w})... hmm, more complex.
#
# Actually, let's try: M_W[v,b] is defined by 2-path-covers of W.
# For the 2-path cover: path1 ending at v through S+{v}, path2 from b through R+{b}.
# If we fix b's first step to go to some w, then...
# This gets complicated.

# Let me try directly from the DEFINITION and see if the recurrence holds algebraically.
# For small m, we verified it computationally. The question is: WHY?

# KEY INSIGHT: The recurrence is:
# col_sum_W(b) = alpha * T(W') + sum_v s_{bv} * cs_{W'}(v)
# where cs_{W'}(v) = sum_{u != v} M_{W'}[u,v]
#
# Can we prove this using the relationship between M_W and M_{W'}?
#
# M_W[v,b] for v in W' involves 2-path-covers of W.
# M_{W'}[u,v] involves 2-path-covers of W' = W\\{b}.
#
# The difference: M_W has vertex b available, M_{W'} doesn't.
# In M_W[v,b]: b is the START of the second path.
# In M_{W'}[u,v]: v is the END of the first path (and u is also an endpoint).

# SUBSTITUTE the row recurrence FOR M_W:
# M_W[v,b] = B_b(W\\{v}) - sum_{u != v,b} t(u,v) M^{(v)}[u,b]
#
# Now B_b(W\\{v}) = B_b(W'\\{v} + {b})... wait, W\\{v} = (W'\\{v}) + {b}
# since W = W' + {b}.
# So B_b(W\\{v}) = Ham paths from b through W\\{v}.
# And W\\{v} = W'\\{v} + {b}.
# So B_b(W\\{v}) = sum_{w in W'\\{v}} t(b,w) B_w(W'\\{v}\\{w} + {w})
# = sum_{w in W'\\{v}} t(b,w) B_w(W'\\{v,w} + {w})
# Hmm, that's wrong too. Let me be careful.
# W\\{v} contains b and all vertices of W' except v.
# B_b(W\\{v}) = sum over permutations of W\\{v} starting at b.
# = sum_{w in W\\{v,b}} t(b,w) B_w(W\\{v,b})
#             ... wait, W\\{v,b} = W'\\{v}
# So B_b(W\\{v}) = sum_{w in W'\\{v}} t(b,w) B_w(W'\\{v})

# Now: sum_v B_b(W\\{v}) = sum_{v in W'} sum_{w in W'\\{v}} t(b,w) B_w(W'\\{v})
# = sum_{v in W'} sum_{w in W'\\{v}} (r + s_{bw}) B_w(W'\\{v})
# = r * sum_{v,w, v!=w} B_w(W'\\{v}) + sum_{v,w,v!=w} s_{bw} B_w(W'\\{v})

# Hmm, this is getting complicated but might lead somewhere.
# B_w(W'\\{v}) counts Ham paths from w through W'\\{v} (m-2 vertices, m-3 edges).

# For each (v,w) pair with v,w in W', v != w:
# B_w(W'\\{v}) is a path from w through all of W' except v.

# sum_{v != w} B_w(W'\\{v}) = sum over Ham paths from w through W'\\{v}
# = sum_v (paths from w through W' that DON'T visit v)
# ... this is related to inclusion-exclusion on visits.

# Actually for each w, sum_{v != w} B_w(W'\\{v}) counts (v, path from w through W'\\{v})
# pairs. This is like "delete one vertex and count paths from w".
# It's (m-2) * average + corrections...

print("The direct algebraic proof via row recurrence expansion")
print("leads to sums over B_w(W'\\{v}) that need further analysis.")
print()
print("Let me verify one more time that the recurrence holds at m=6,7.")

for m in [6, 7]:
    n_syms = max(m, 7)
    r, sv, s, t = setup(n_syms)
    W = set(range(m))
    b = 0
    W_prime = W - {b}

    col_sum_W = 0
    for v in sorted(W - {b}):
        col_sum_W = expand(col_sum_W + transfer_M(t, W, v, b))

    alpha = 1 if m % 2 == 0 else 0
    T_Wp = total_ham(W_prime, t)

    s_weighted = 0
    for v in sorted(W_prime):
        cs_v = 0
        for u in sorted(W_prime - {v}):
            cs_v = expand(cs_v + transfer_M(t, W_prime, u, v))
        s_weighted = expand(s_weighted + s(b, v) * cs_v)

    rhs = expand(alpha * T_Wp + s_weighted)
    match = expand(col_sum_W - rhs) == 0
    print(f"m={m}: RECURRENCE HOLDS: {match}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
