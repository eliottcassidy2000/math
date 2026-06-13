#!/usr/bin/env python3
"""
COLUMN RECURRENCE for M and algebraic proof of the key identity.

DERIVATION: By expanding B_b in the 2-path-cover formula by its first step:
  B_b(R+{b}) = sum_{w in R} t(b,w) B_w(R)  for |R| >= 1

This gives the COLUMN RECURRENCE:
  M_W[a,b] = (-1)^{m-2} E_a(W') + sum_{w in U} t(b,w) M_{W'}[a,w]

where W' = W\{b}, U = W\{a,b}, m = |W|.

CONSEQUENCE: Summing over a (to get col_sum):
  col_sum_W(b) = (-1)^{m-2} T(W') + r * Sigma + sum_w s_bw * cs_{W'}(w)

where Sigma = total off-diagonal sum of M on W'.

Comparing with the RECURRENCE we need:
  col_sum_W(b) = alpha * T(W') + sum_w s_bw * cs_{W'}(w)
  where alpha = [1+(-1)^{m-2}]/2

This requires:
  r * Sigma = [1-(-1)^{m-2}]/2 * T(W')

  For |W'| odd (m even): Sigma = 0
  For |W'| even (m odd): r * Sigma = T(W')

These are NEW IDENTITIES about the transfer matrix!

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


# ============================================================
# STEP 1: Verify the COLUMN RECURRENCE
# ============================================================
print("=" * 70)
print("STEP 1: Column recurrence verification")
print("M_W[a,b] = (-1)^{m-2} E_a(W') + sum_w t(b,w) M_{W'}[a,w]")
print("=" * 70)

for m in [3, 4, 5, 6]:
    r, sv, s, t = setup(m)
    W = set(range(m))

    for a, b in [(0, 1)]:
        W_prime = W - {b}
        U = sorted(W - {a, b})

        M_direct = transfer_M(t, W, a, b)

        # Column recurrence
        Ea = ham_paths_to(W_prime, a, t)
        rec_sum = 0
        for w in U:
            M_aw = transfer_M(t, W_prime, a, w)
            rec_sum = expand(rec_sum + t(b, w) * M_aw)

        M_rec = expand((-1)**(m-2) * Ea + rec_sum)
        match = expand(M_direct - M_rec) == 0
        print(f"  m={m}, a={a}, b={b}: column recurrence holds = {match}")


# ============================================================
# STEP 2: Total off-diagonal sum of M (Sigma)
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Sigma = total off-diagonal sum of M")
print("=" * 70)

for m_prime in range(2, 7):
    r, sv, s, t = setup(m_prime)
    W_prime = set(range(m_prime))

    Sigma = 0
    for a in sorted(W_prime):
        for w in sorted(W_prime):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W_prime, a, w))

    T_Wp = total_ham(W_prime, t)

    print(f"\n|W'|={m_prime}:")
    print(f"  Sigma = {Sigma}")
    print(f"  T(W') = {T_Wp}")

    if m_prime % 2 == 1:  # |W'| odd
        print(f"  |W'| odd => need Sigma = 0: {Sigma == 0}")
    else:  # |W'| even
        check = expand(r * Sigma - T_Wp)
        print(f"  |W'| even => need r*Sigma = T: {check == 0}")
        if check != 0:
            print(f"    r*Sigma - T = {check}")


# ============================================================
# STEP 3: Understanding Sigma
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Why does Sigma behave this way?")
print("=" * 70)

# Sigma = sum_{a != w} M[a,w]
# = sum_w (sum_a M[a,w]) = sum_w cs_w (sum of column sums)
# = sum_a (sum_w M[a,w]) = sum_a rs_a (sum of row sums)
#
# Row sum: sum_w M[a,w] = ?
# Each M[a,w] involves 2-path-covers of W' with endpoints a,w.
# sum_w M[a,w] sums over all second endpoints w.
#
# Let's compute row sums separately:

for m_prime in [3, 4, 5]:
    r, sv, s, t = setup(m_prime)
    W_prime = set(range(m_prime))

    print(f"\n|W'|={m_prime}:")
    for a in sorted(W_prime):
        row_sum = 0
        for w in sorted(W_prime):
            if w == a: continue
            row_sum = expand(row_sum + transfer_M(t, W_prime, a, w))
        print(f"  row_sum(a={a}) = {row_sum}")

    # Is row_sum = B_a(W') (total Ham path weight from a)?
    for a in sorted(W_prime):
        Ba = 0
        for perm in permutations(sorted(W_prime)):
            if perm[0] != a: continue
            prod = 1
            for i in range(len(perm)-1):
                prod *= t(perm[i], perm[i+1])
            Ba += prod
        Ba = expand(Ba)

        row_sum = 0
        for w in sorted(W_prime):
            if w == a: continue
            row_sum = expand(row_sum + transfer_M(t, W_prime, a, w))

        print(f"  B_{a}(W') = {Ba}, row_sum(a={a}) = {row_sum}")
        print(f"    row_sum = B_a? {expand(Ba - row_sum) == 0}")


# ============================================================
# STEP 4: Does row_sum(a) relate to B_a?
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Relationship between row sums and Ham path sums")
print("=" * 70)

# From the data above, it seems row_sum != B_a in general.
# Let me check: is row_sum always even in r (by induction)?
# If so, and Sigma = sum of row sums, then Sigma is even in r.
# For |W'| odd: T is even in r, so r*Sigma = T requires Sigma odd.
# But if Sigma is even, r*Sigma is odd, and T is even => r*Sigma != T unless both = 0.
# Hmm, wait.

# For |W'| odd: T has parity (-1)^{|W'|-1} = (-1)^{even} = 1 (even r-powers only).
# We need Sigma = 0.

# For |W'| even: T has parity (-1)^{|W'|-1} = -1 (odd r-powers only).
# We need r * Sigma = T. Since T is odd, T = r * T' where T' is even.
# So Sigma = T' = T/r.

# Let's check the r-parity of Sigma:
for m_prime in range(2, 7):
    r, sv, s, t = setup(m_prime)
    W_prime = set(range(m_prime))

    Sigma = 0
    for a in sorted(W_prime):
        for w in sorted(W_prime):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W_prime, a, w))

    if Sigma == 0:
        print(f"|W'|={m_prime}: Sigma = 0")
    else:
        p = Poly(Sigma, r)
        parities = set(d[0] % 2 for d, c in p.as_dict().items() if c != 0)
        print(f"|W'|={m_prime}: Sigma r-parities = {parities}, Sigma = {Sigma}")


# ============================================================
# STEP 5: Can we prove Sigma identities WITHOUT assuming even r?
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Direct proof of Sigma identities")
print("=" * 70)

# Sigma = sum_{a != w} M[a,w]
# = sum_{a != w} sum_{S ⊆ W'\{a,w}} (-1)^|S| E_a(S+{a}) B_w(R+{w})
# = sum_{ordered pairs (a,w)} sum_{S ⊆ W'\{a,w}} (-1)^|S| E_a(S+{a}) B_w(R+{w})
#
# For each 2-path-cover: S+{a} is the "ending at a" path,
# R+{w} = (W'\{a,w})\S + {w} is the "starting at w" path.
# Together S+{a} and R+{w} partition W' (since S ∪ R = W'\{a,w}, plus a and w).
#
# So: each contribution is a PAIR of ordered Hamiltonian paths
# (P1 ending at a, P2 starting at w) that together cover all of W',
# with sign (-1)^{|P1|-1}.
#
# Sigma sums this over ALL choices of a and w.
#
# For a FIXED pair of paths (P1, P2):
# P1 = (v1, ..., v_k, a) through S+{a}
# P2 = (w, u1, ..., u_l) through R+{w}
# |P1| = k+1 = |S|+1, |P2| = l+1 = |R|+1.
# Sign = (-1)^|S| = (-1)^{k} = (-1)^{|P1|-1}.
#
# In Sigma, this term appears for the specific pair (a, w).
# But different (a,w) pairs give different 2-path-covers.
#
# So Sigma = sum over ALL 2-path-covers (P1 ending at some a, P2 starting at some w)
# of (-1)^{|P1|-1} * weight(P1) * weight(P2).

# Now consider the REVERSAL of P1: P1^R = (a, v_k, ..., v1) starts at a.
# The reversed cover has P1^R starting at a and P2 starting at w.
# weight(P1^R) = prod t(v_{i+1}, v_i) = prod (r - s_{v_i, v_{i+1}}) = E_a(weight; r, -s)...
# This doesn't obviously simplify.

# KEY INSIGHT: In Sigma, we sum M[a,w] over ALL (a,w).
# By the re-indexing identity: M[a,w](T) = M[a,w](r,s).
# And M[w,a](r,s) = M[a,w](r,s)??? (THIS IS WHAT WE'RE TRYING TO PROVE!)
#
# If M were symmetric (M[a,w] = M[w,a]), then:
# Sigma = 2 * sum_{a < w} M[a,w] and the identity would follow from symmetry.
# But we can't assume symmetry.
#
# HOWEVER: M[a,w] + M[w,a] has even r-powers
# (since M[a,w](-r) = M[w,a](r), so M[a,w]+M[w,a] is even in r).
# Similarly, M[a,w] - M[w,a] has ODD r-powers.
#
# Sigma = sum_{a!=w} M[a,w] = sum_{a<w} (M[a,w] + M[w,a]).
# So Sigma has EVEN r-powers (each pair M[a,w]+M[w,a] is even).

# WAIT! This is PROVABLE without assuming M is symmetric!
# M[a,w](-r) = M[w,a](r) is ALREADY PROVEN (the re-indexing identity).
# So M[a,w](r) + M[w,a](r) = M[a,w](r) + M[a,w](-r) [by re-indexing]
# which has even r-powers.
# Therefore Sigma has even r-powers.

print("KEY OBSERVATION: M[a,w] + M[w,a] has even r-powers (PROVED).")
print("So Sigma = sum_{a<w} (M[a,w]+M[w,a]) has even r-powers.")
print()

# For |W'| ODD: T(W') has EVEN r-powers (parity (-1)^{|W'|-1} = 1).
# Need: Sigma = 0.
# We know Sigma is even in r. But is it zero?
# NOT OBVIOUS from parity alone!

# For |W'| EVEN: T(W') has ODD r-powers (parity -1).
# Need: r * Sigma = T.
# Sigma is even, so r*Sigma is odd. T is odd. Compatible parities.
# But need equality.

# Let me check: for |W'| = 3, Sigma = 0. Why?
print("Checking Sigma for small |W'|:")
for m_prime in [2, 3, 4, 5]:
    r, sv, s, t = setup(m_prime)
    W_prime = set(range(m_prime))
    Sigma = 0
    for a in sorted(W_prime):
        for w in sorted(W_prime):
            if a == w: continue
            Sigma = expand(Sigma + transfer_M(t, W_prime, a, w))
    print(f"  |W'|={m_prime}: Sigma = {Sigma}")


# ============================================================
# STEP 6: Prove Sigma = 0 for odd |W'| and r*Sigma = T for even |W'|
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Direct computation of Sigma")
print("=" * 70)

# For |W'| = 2: W' = {0,1}. M[0,1] = 1, M[1,0] = 1.
# Sigma = 2. T = 2r. r*Sigma = 2r = T. ✓

# For |W'| = 3: W' = {0,1,2}.
# Sigma = sum of all 6 M-entries.
# M[0,1] = s02+s12, M[1,0] = s02-s12
# M[0,2] = s01-s12, M[2,0] = s01+s12
# M[1,2] = s01+s02, M[2,1] = -(s01+s02)??? Let me compute.

print("\n|W'|=3 detailed:")
r, sv, s, t = setup(3)
W = {0, 1, 2}
for a in sorted(W):
    for w in sorted(W):
        if a == w: continue
        M_aw = transfer_M(t, W, a, w)
        print(f"  M[{a},{w}] = {M_aw}")

# Check: M[a,w] + M[w,a]
print("  Symmetric parts:")
for a in sorted(W):
    for w in sorted(W):
        if a >= w: continue
        M_aw = transfer_M(t, W, a, w)
        M_wa = transfer_M(t, W, w, a)
        print(f"  M[{a},{w}]+M[{w},{a}] = {expand(M_aw + M_wa)}")
        print(f"  M[{a},{w}]-M[{w},{a}] = {expand(M_aw - M_wa)}")


# ============================================================
# STEP 7: Row sum = column sum?
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Row sums vs column sums of M")
print("=" * 70)

for m_prime in [3, 4, 5]:
    r, sv, s, t = setup(m_prime)
    W = set(range(m_prime))
    print(f"\n|W|={m_prime}:")
    for v in sorted(W):
        rs = 0  # row sum: sum_w M[v,w]
        cs = 0  # col sum: sum_a M[a,v]
        for w in sorted(W):
            if w == v: continue
            rs = expand(rs + transfer_M(t, W, v, w))
            cs = expand(cs + transfer_M(t, W, w, v))
        print(f"  v={v}: row_sum = {rs}, col_sum = {cs}")
        print(f"       row=col? {expand(rs-cs)==0}, sum={expand(rs+cs)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
