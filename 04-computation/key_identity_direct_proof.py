#!/usr/bin/env python3
"""
DIRECT PROOF ATTEMPT for KEY IDENTITY: odd(B_b(W)) = r * sum_v M_W[v,b]

STRATEGY: Expand both sides and show they are equal.

RHS = r * sum_v M_W[v,b]
    = r * sum_{v != b} sum_{S subset W\{v,b}} (-1)^|S| E_v(S+{v}) B_b(R+{b})

Exchange sums (v and S):
    = r * sum_{S' non-empty subset W\{b}} (-1)^{|S'|-1} T(S') B_b(R+{b})
where T(S') = sum_{v in S'} E_v(S') = total Ham path weight through S',
and R = (W\{b}) \ S'.

LHS = odd(B_b(W)) = [B_b(W) - B_b(W;-r)] / 2 = [B_b - (-1)^{m-1} E_b] / 2

So we need: B_b(W) - (-1)^{m-1} E_b(W) = 2r * sum_{S'} (-1)^{|S'|-1} T(S') B_b(R+{b})

NEW IDEA: Maybe we can prove B_b(W) = sum_{S'} f(S') * T(S') * B_b(R+{b})
for some coefficients f(S'), and similarly for E_b(W).

Actually, let's think about B_b(W) as a Hamiltonian path from b through W.
For each path b -> v1 -> ... -> v_{m-1}, partition the vertices by
"prefix up to some point" and "suffix from that point."

If we split the path at position k (after visiting k+1 vertices including b):
prefix = (b, v1, ..., v_k) and suffix = (v_k, v_{k+1}, ..., v_{m-1})
But the prefix doesn't end at a fixed vertex.

BETTER IDEA: Use the DELETION formula.
B_b(W) = sum_{v in W\{b}} H(b,v; W)
where H(b,v; W) = Hamiltonian paths b -> v through W.

Now, H(b,v; W) relates to M[v,b] somehow. We showed:
odd(H(b,v;W)) != r * M[v,b]  (not termwise)
But sum_v odd(H) = r * sum_v M[v,b] (at sum level).

The difference: odd(H(b,v)) - r*M[v,b] sums to zero over v.
What IS odd(H(b,v)) - r*M[v,b]?

Let's analyze this "error term" e(v) = odd(H(b,v)) - r*M[v,b].
We need: sum_v e(v) = 0.

opus-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Poly, factor, collect

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

def ham_paths_between(vertex_set, start, end, t_fn):
    vs = sorted(vertex_set)
    if len(vs) == 1:
        return 1 if start == end else 0
    total = 0
    for perm in permutations(vs):
        if perm[0] != start or perm[-1] != end:
            continue
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
# PART 1: Analyze the "error term" e(v) = odd(H(b,v)) - r*M[v,b]
# ============================================================
print("=" * 70)
print("ERROR TERM ANALYSIS: e(v) = odd(H(b,v;W)) - r*M[v,b]")
print("=" * 70)

for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))
    b = 0

    print(f"\nm={m}, b={b}:")
    errors = {}
    for v in sorted(W - {b}):
        H_bv = ham_paths_between(W, b, v, t)
        M_vb = transfer_M(t, W, v, b)
        e_v = expand(odd_part(H_bv, r) - r * M_vb)
        errors[v] = e_v
        print(f"  e({v}) = {e_v}")

    total_error = expand(sum(errors.values()))
    print(f"  sum e(v) = {total_error}")

    # Can we express e(v) in terms of M entries?
    # e(v) involves odd parts of H and r*M.
    # Note: H(b,v) = B_b restricted to paths ending at v
    # M[v,b] = inclusion-exclusion over 2-path-covers with endpoints v,b

    # KEY QUESTION: Is e(v) antisymmetric in some way?
    # e.g., e(v) = -e(w) for some pairing v <-> w?

    if m == 4:
        print(f"\n  Pairwise sums:")
        vs = sorted(W - {b})
        for i in range(len(vs)):
            for j in range(i+1, len(vs)):
                print(f"    e({vs[i]}) + e({vs[j]}) = {expand(errors[vs[i]] + errors[vs[j]])}")

    # Factor e(v):
    for v in sorted(W - {b}):
        if errors[v] != 0:
            # Check if e(v) is divisible by r
            p = Poly(errors[v], r)
            min_deg = min(d[0] for d in p.as_dict().keys())
            print(f"  e({v}): min r-degree = {min_deg}, factored = {factor(errors[v])}")


# ============================================================
# PART 2: Alternative — prove via row recurrence APPLIED TWICE
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Row recurrence applied at different levels")
print("=" * 70)

# The row recurrence for B_b(W) itself:
# B_b(W) = sum_v t(b,v) * B_v(W\{b})  ... no, this isn't right.
# B_b(W) = sum_{v != b} t(b,v) * H(v, ..., W\{b})
# where H(v, ..., W\{b}) counts Ham paths from v through W\{b} (ending anywhere).

# Actually: B_b(W) = sum_{v in W\{b}} t(b,v) * B_v(W\{b})... NO.
# t(b,v) * B_v(W\{b}) would count paths b -> v -> (path from v through W\{b}).
# But B_v(W\{b}) counts paths from v through W\{b}, which is correct!
# So: B_b(W) = sum_{v in W\{b}} t(b,v) * B_v(W\{b})  ... WAIT.
# B_v(W\{b}) = Ham paths from v through W\{b} (all m-1 vertices in W\{b}).
# Then t(b,v)*B_v(W\{b}) = paths b -> v -> ... through all of W.
# But these end at ANY vertex in W\{b,v}... or at v? No, v is visited first.
# So: b -> v -> w1 -> w2 -> ... where {w1,...,w_{m-2}} = W\{b,v}.
# This is a Ham path through W starting at b with b's successor = v.
# Summing over v: we get ALL Ham paths from b through W. YES!

# So: B_b(W) = sum_{v != b} t(b,v) * B_v(W\{b})

# Verify:
for m in [3, 4]:
    r, sv, s, t = setup(m)
    W = set(range(m))
    b = 0

    # Direct B_b
    Bb_direct = 0
    for perm in permutations(sorted(W)):
        if perm[0] != b: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t(perm[i], perm[i+1])
        Bb_direct += prod
    Bb_direct = expand(Bb_direct)

    # Via recurrence
    Bb_rec = 0
    for v in sorted(W - {b}):
        Bv = 0
        W_minus_b = W - {b}
        for perm in permutations(sorted(W_minus_b)):
            if perm[0] != v: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            Bv += prod
        Bv = expand(Bv)
        Bb_rec = expand(Bb_rec + t(b, v) * Bv)

    print(f"  m={m}: B_b recurrence check = {expand(Bb_direct - Bb_rec) == 0}")

# So: B_b(W) = sum_v t(b,v) * B_v(W\{b})
#            = sum_v (r + s_{bv}) * B_v(W\{b})
#            = r * sum_v B_v(W\{b}) + sum_v s_{bv} * B_v(W\{b})
#
# Now: odd(B_b(W)) = r * even(sum_v B_v(W\{b})) + odd(sum_v s_{bv} * B_v(W\{b}))
#                   + r * odd(...) + even(...) ... actually this isn't clean.
#
# Let's be precise. Let B_v = B_v(W\{b}).
# B_b(W) = r * sum_v B_v + sum_v s_{bv} B_v
#
# odd(B_b) = r * even(sum_v B_v) + odd(sum_v s_{bv} B_v)
#
# And we need this to equal r * col_sum where col_sum = sum_v M_W[v,b].
#
# So: r * [even(sum_v B_v) - col_sum] + odd(sum_v s_{bv} B_v) = 0
# Both terms must vanish separately (different r-parities in the bracket):
# (A) even(sum_v B_v) = col_sum (matching even parts)
# (B) odd(sum_v s_{bv} B_v) = 0

print("\n--- Testing conditions (A) and (B) ---")

for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))
    b = 0

    # sum_v B_v(W\{b})
    sum_Bv = 0
    for v in sorted(W - {b}):
        Bv = 0
        for perm in permutations(sorted(W - {b})):
            if perm[0] != v: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            Bv += prod
        sum_Bv = expand(sum_Bv + Bv)

    # This is T(W\{b}), total Ham path weight through W\{b}!
    # We showed T has definite parity (-1)^{|W\{b}|-1} = (-1)^{m-2}
    print(f"\nm={m}: sum_v B_v(W\\{{b}}) = T(W\\{{b}}) = {sum_Bv}")

    from sympy import Poly as P
    if sum_Bv != 0:
        p = P(sum_Bv, r)
        parities = set(d[0] % 2 for d, c in p.as_dict().items() if c != 0)
        print(f"  T(W\\{{b}}) r-parities: {parities}")

    # even part of T
    if sum_Bv != 0:
        p_dict = P(sum_Bv, r).as_dict()
        even_T = expand(sum(c * r**d[0] for d, c in p_dict.items() if d[0] % 2 == 0))
    else:
        even_T = 0

    # col_sum = sum_v M[v,b]
    col_sum = 0
    for v in sorted(W - {b}):
        col_sum = expand(col_sum + transfer_M(t, W, v, b))

    print(f"  even(T) = {even_T}")
    print(f"  col_sum = {col_sum}")
    print(f"  (A) even(T) = col_sum? {expand(even_T - col_sum) == 0}")

    # sum_v s_{bv} B_v(W\{b})
    s_weighted = 0
    for v in sorted(W - {b}):
        Bv = 0
        for perm in permutations(sorted(W - {b})):
            if perm[0] != v: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            Bv += prod
        Bv = expand(Bv)
        s_weighted = expand(s_weighted + s(b, v) * Bv)

    odd_s_weighted = odd_part(s_weighted, r)
    print(f"  sum_v s_bv B_v = {s_weighted}")
    print(f"  (B) odd(sum_v s_bv B_v) = {odd_s_weighted}")
    print(f"  (B) = 0? {odd_s_weighted == 0}")


# ============================================================
# PART 3: Recursive structure — if (A) and (B) hold, what do they mean?
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Understanding conditions (A) and (B)")
print("=" * 70)

# T(W') has definite parity (-1)^{|W'|-1}.
# For W' = W\{b} with |W'| = m-1:
#   If m is odd: |W'| = m-1 is even, so T has parity (-1)^{m-2} = (-1)^{m}.
#     For m odd: T has parity (-1)^{m} = -1, so T is odd.
#     Then even(T) = 0.
#     Condition (A): col_sum = 0.
#
#   If m is even: |W'| = m-1 is odd, so T has parity (-1)^{m-2} = 1.
#     T is entirely even.
#     Then even(T) = T itself.
#     Condition (A): col_sum = T(W\{b}).

# Condition (B): odd(sum_v s_{bv} B_v) = 0.
# s_{bv} doesn't change r-degree. B_v has parity (-1)^{m-2} (same as T).
# So s_{bv} * B_v has the same r-parity as B_v.
# Wait, B_v(W\{b}) is a SINGLE-STARTING-POINT path sum, not total.
# Does B_v have definite parity? Not necessarily!

# Actually, T(W') = sum_v B_v(W') has definite parity,
# but individual B_v(W') terms need not.
# For m=3, W\{b} = {1,2}: B_1 = t(1,2) = r+s12, B_2 = t(2,1) = r-s12.
# Both have r^1 and r^0, mixed parity!
# But T = B_1 + B_2 = 2r, purely odd.

# So condition (B) says the odd part of sum_v s_{bv} B_v vanishes.
# Since B_v has both parities, s_{bv} B_v also has both, but the sum
# over v somehow kills the odd part.

# For m=3: sum_v s_{0v} B_v = s01*(r+s12) + s02*(r-s12)
#         = r*(s01+s02) + s01*s12 - s02*s12
# odd part = r*(s01+s02). For this to be 0, we need s01+s02 = 0? NO!
# But we computed odd(sum s_bv B_v) = 0 above... let me recheck.

print("\nRecheck m=3:")
r, sv, s, t = setup(3)
W = {0, 1, 2}
b = 0
for v in [1, 2]:
    Bv = ham_paths_between(W - {b}, v, None, t)  # Wait, need B_v not H
    # B_v(W\{b}) = Ham paths starting at v through W\{b}
print("  Need to recompute B_v properly...")

# B_v(W\{b}) for W\{b} = {1,2}:
# B_1({1,2}) = path starting at 1 through {1,2} = t(1,2) = r + s12
# B_2({1,2}) = path starting at 2 through {1,2} = t(2,1) = r - s12
Bv1 = expand(t(1, 2))
Bv2 = expand(t(2, 1))
print(f"  B_1({{1,2}}) = {Bv1}")
print(f"  B_2({{1,2}}) = {Bv2}")
print(f"  s01*B_1 + s02*B_2 = {expand(s(0,1)*Bv1 + s(0,2)*Bv2)}")
sw = expand(s(0,1)*Bv1 + s(0,2)*Bv2)
print(f"  odd part = {odd_part(sw, r)}")
# s01*(r+s12) + s02*(r-s12) = r(s01+s02) + s01*s12 - s02*s12
# odd part = r(s01+s02)
# This is NOT zero unless s01+s02=0!

# Hmm, but above we computed (B) = 0 for m=3. Let me recheck...
# Wait, the computation above was:
# B_b(W) = r * T(W\{b}) + sum_v s_{bv} B_v(W\{b})
# For m=3: B_0({0,1,2}) = t(0,1)*t(1,2) + t(0,2)*t(2,1)
B0_direct = expand(t(0,1)*t(1,2) + t(0,2)*t(2,1))
print(f"\n  B_0({{0,1,2}}) = {B0_direct}")
# = (r+s01)(r+s12) + (r+s02)(r-s12)
# = r^2 + r*s01 + r*s12 + s01*s12 + r^2 + r*s02 - r*s12 - s02*s12
# = 2r^2 + r(s01+s02) + s01*s12 - s02*s12

# r*T = r * 2r = 2r^2
# sum s_bv B_v = s01*(r+s12) + s02*(r-s12) = r(s01+s02) + s01*s12 - s02*s12
# Total: 2r^2 + r(s01+s02) + s01*s12 - s02*s12 ✓

# odd(B_0) = r(s01+s02)
# r * col_sum: col_sum = M[1,0] + M[2,0]
M10 = transfer_M(t, {0,1,2}, 1, 0)
M20 = transfer_M(t, {0,1,2}, 2, 0)
print(f"  M[1,0] = {M10}")
print(f"  M[2,0] = {M20}")
print(f"  col_sum = {expand(M10+M20)}")
print(f"  r*col_sum = {expand(r*(M10+M20))}")
print(f"  odd(B_0) = {odd_part(B0_direct, r)}")

# So odd(B_0) = r*(s01+s02) and col_sum = s01+s02.
# KEY IDENTITY: r*(s01+s02) = r*(s01+s02). ✓

# And the decomposition:
# odd(B_b) = r*even(T) + odd(sum s_bv B_v) = r*0 + r*(s01+s02)
# Wait: T = 2r which is ODD, so even(T) = 0.
# odd(sum s_bv B_v) = r*(s01+s02) which is ODD.
# So odd(B_b) = 0 + r*(s01+s02) = r*(s01+s02).
# And r*col_sum = r*(s01+s02). ✓

# BUT condition (B) said odd(sum s_bv B_v) = 0, which is FALSE here!
# The error is in the decomposition above.

# Let me redo: B_b = r * sum_v B_v + sum_v s_bv B_v
# odd(B_b) = r * even(sum_v B_v) + odd(sum_v s_bv B_v)
# = r * 0 + r*(s01+s02) = r*(s01+s02)
# r * col_sum = r*(s01+s02)
# So the identity holds, but NOT because (B)=0.
# It holds because: r*even(T) + odd(s_part) = r*col_sum
# => odd(s_part) = r*(col_sum - even(T))
# => odd(s_part) = r*col_sum when T is purely odd (even(T) = 0)

# So the identity becomes:
# r * even(T(W\{b})) + odd(sum s_{bv} B_v(W\{b})) = r * col_sum
# r * [col_sum - even(T)] = odd(sum s_{bv} B_v)
# BOTH sides should be equal!

print("\n--- Corrected analysis ---")
for m in [3, 4, 5]:
    r, sv, s, t = setup(m)
    W = set(range(m))
    b = 0

    # T(W\{b})
    T_total = 0
    Bv_dict = {}
    for v in sorted(W - {b}):
        Bv = 0
        for perm in permutations(sorted(W - {b})):
            if perm[0] != v: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            Bv += prod
        Bv = expand(Bv)
        Bv_dict[v] = Bv
        T_total = expand(T_total + Bv)

    from sympy import Poly as P
    if T_total != 0:
        p_dict = P(T_total, r).as_dict()
        even_T = expand(sum(c * r**d[0] for d, c in p_dict.items() if d[0] % 2 == 0))
    else:
        even_T = 0

    # s-weighted sum
    s_part = 0
    for v in sorted(W - {b}):
        s_part = expand(s_part + s(b, v) * Bv_dict[v])
    odd_s = odd_part(s_part, r)

    # col_sum
    col_sum = 0
    for v in sorted(W - {b}):
        col_sum = expand(col_sum + transfer_M(t, W, v, b))

    # Check: r*(col_sum - even(T)) = odd(s_part)?
    lhs = expand(r * (col_sum - even_T))
    rhs = odd_s
    print(f"\nm={m}:")
    print(f"  even(T) = {even_T}")
    print(f"  col_sum = {col_sum}")
    print(f"  r*(col_sum - even(T)) = {lhs}")
    print(f"  odd(s_part) = {rhs}")
    print(f"  Match: {expand(lhs - rhs) == 0}")

    # Also check the full identity
    full_lhs = expand(r * even_T + odd_s)
    full_rhs = expand(r * col_sum)
    print(f"  Full identity: r*even(T) + odd(s) = r*col_sum? {expand(full_lhs - full_rhs) == 0}")


# ============================================================
# PART 4: Recursive formulation
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Can we prove this recursively?")
print("=" * 70)

# We have: odd(B_b(W)) = r * col_sum_W(b)
# Using: B_b(W) = r * T(W') + s_part, where W' = W\{b}
# This gives: odd(B_b) = r * even(T(W')) + odd(s_part)
# And we need this = r * col_sum.

# Now, T(W') itself satisfies the key identity for SMALLER sets!
# T(W') = sum_v B_v(W') where each B_v is a path sum.
# And T(W') has definite parity.

# Can we relate col_sum to T and s_part?
# col_sum = sum_v M_W[v,b]
# Using row recurrence on M_W[v,b] (with "a" = v):
# M_W[v,b] = B_b(W\{v}) - sum_{u != v,b} t(u,v) M^{(v)}[u,b]

# So: col_sum = sum_v B_b(W\{v}) - sum_v sum_u t(u,v) M^{(v)}[u,b]
# This introduces M^{(v)} which is on W\{v}... getting complicated.

# INSTEAD: let's check if there's a DIRECT relation between col_sum
# and the boundary terms.

# At r=0: col_sum(0) = d/dr B_b(W)|_{r=0}
# This was verified. What about higher order?
# col_sum = sum_k r^k D_k, B_b = sum_k r^k C_k
# Identity: C_{2j+1} = D_{2j} for all j.
# So D_0 = C_1, D_2 = C_3, D_4 = C_5, ...
# Each even coefficient of col_sum equals the next odd coefficient of B_b.
#
# C_1 = [r^1] of B_b = sum_sigma sum_j prod_{i!=j} s_{edges}
# D_0 = [r^0] of col_sum = col_sum at r=0 = sum_v M_W[v,b]|_{r=0}

# We already verified C_1 = D_0 (the derivative identity at r=0).
# C_3 = D_2 is a NEW identity at higher degree.

# Let me verify C_3 = D_2 at m=5:
m = 5
r, sv, s, t = setup(m)
W = set(range(m))
b = 0

Bb = 0
for perm in permutations(sorted(W)):
    if perm[0] != b: continue
    prod = 1
    for i in range(len(perm)-1):
        prod *= t(perm[i], perm[i+1])
    Bb += prod
Bb = expand(Bb)

col_sum = 0
for v in sorted(W - {b}):
    col_sum = expand(col_sum + transfer_M(t, W, v, b))

Bb_poly = Poly(Bb, r)
cs_poly = Poly(col_sum, r)

C1 = expand(Bb_poly.nth(1))
C3 = expand(Bb_poly.nth(3))
D0 = expand(cs_poly.nth(0))
D2 = expand(cs_poly.nth(2))

print(f"\nm=5: C_1 = [r^1]B_b = {C1}")
print(f"     D_0 = [r^0]col = {D0}")
print(f"     C_1 = D_0? {expand(C1 - D0) == 0}")

print(f"\n     C_3 = [r^3]B_b = {C3}")
print(f"     D_2 = [r^2]col = {D2}")
print(f"     C_3 = D_2? {expand(C3 - D2) == 0}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
