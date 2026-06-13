"""
sigma_a2_connection.py -- kind-pasteur-2026-03-13-S61

Deep dive into Sigma = (n-2)*J_off - sym(A^2)_off

Since Sigma encodes the complement of A^2 symmetrization, and A^2 entries
are literally delta (forward paths) and lambda (reverse paths), we can
relate tr(Sigma^k) to traces of A^k via polynomial expansion.

Key identity: Sigma = c*J_off - M, where c = n-2, M = sym(A^2)_off
M[u,v] = A^2[u,v] + A^2[v,u] for u != v, M[u,u] = 0

Since M encodes lambda + delta = n-2-sigma per pair, and
M[u,v] = lambda(u,v) + delta(u,v), we have:
Sigma = c*J_off - M => Sigma + M = c*J_off

This script:
1. Expresses tr(Sigma^k) in terms of tr(M^k) and tr(J^k)
2. Relates M to A^2 explicitly
3. Finds what 48 corresponds to in the A^2 language
4. Tests whether M determines c7 directly
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"SIGMA-A^2 CONNECTION AT n={n}")
print("=" * 60)

# The matrix M = sym(A^2)_off = A^2 + (A^2)^T with zeroed diagonal
# M[u,v] = lambda(u,v) + delta(u,v) = n-2-sigma(u,v) for u != v
# This is also M[u,v] = #{2-paths from u to v} + #{2-paths from v to u}

# Now: Sigma = (n-2)*J_off - M
# => Sigma^2 = (n-2)^2 * J_off^2 - 2*(n-2)*J_off*M + M^2
# J_off^2[u,v] = sum_w J_off[u,w]*J_off[w,v] = (n-1) if u=v (sum of n-1 ones excluding u,w=u)
#              = (n-2) if u!=v (sum of n-2 ones, excluding w=u and w=v from J_off)
# Wait: J_off[u,w] = 1 if u!=w, 0 if u=w
# J_off^2[u,v] = sum_w 1[u!=w]*1[w!=v]
#              = n if u=v (sum over w!=u of 1[w!=u] = wait...
# Let me be more careful: J_off = J - I where J is all-ones.
# J_off^2 = (J-I)^2 = J^2 - 2J + I = nJ - 2J + I = (n-2)J + I
# So tr(J_off^2) = tr((n-2)J + I) = (n-2)*n + n = n*(n-1) = 42

np.random.seed(42)

# Verify the identity algebraically
test_bits = np.random.randint(0, 1 << total_bits)
A = bits_to_adj(test_bits, n)
A2 = A @ A
M = A2 + A2.T
np.fill_diagonal(M, 0)

J_off = np.ones((n, n), dtype=int) - np.eye(n, dtype=int)
Sigma = (n-2) * J_off - M

# Check traces
print(f"\ntr(J_off^k) values:")
for k in range(1, 6):
    Jk = np.linalg.matrix_power(J_off, k)
    print(f"  k={k}: tr(J_off^{k}) = {int(np.trace(Jk))}")

# 1. Express tr(Sigma^k) in terms of M
# Sigma = c*J_off - M, c = n-2 = 5
# tr(Sigma^2) = c^2*tr(J_off^2) - 2c*tr(J_off*M) + tr(M^2)
# tr(Sigma^3) = c^3*tr(J_off^3) - 3c^2*tr(J_off^2*M) + 3c*tr(J_off*M^2) - tr(M^3)

print(f"\n--- tr(Sigma^k) vs tr(M^k) ---")

data = []
np.random.seed(42)
for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    A2 = A @ A

    M = A2 + A2.T
    np.fill_diagonal(M, 0)
    Sigma = (n-2) * J_off - M

    c7 = count_directed_k_cycles(A, n, 7)
    c3 = count_directed_k_cycles(A, n, 3)

    tr_s2 = int(np.trace(Sigma @ Sigma))
    tr_s3 = int(np.trace(Sigma @ Sigma @ Sigma))
    tr_m2 = int(np.trace(M @ M))
    tr_m3 = int(np.trace(M @ M @ M))

    # Also: tr(A^k) for the TOURNAMENT adjacency matrix
    tr_a2 = int(np.trace(A2))  # = 0 always (tournament has no 2-cycles)
    tr_a3 = int(np.trace(A2 @ A))  # = 3*c3
    tr_a4 = int(np.trace(np.linalg.matrix_power(A, 4)))
    tr_a5 = int(np.trace(np.linalg.matrix_power(A, 5)))
    tr_a6 = int(np.trace(np.linalg.matrix_power(A, 6)))
    tr_a7 = int(np.trace(np.linalg.matrix_power(A, 7)))  # = 7*c7

    data.append({
        'c7': c7,
        'c3': c3,
        'tr_s2': tr_s2,
        'tr_s3': tr_s3,
        'tr_m2': tr_m2,
        'tr_m3': tr_m3,
        'tr_a3': tr_a3,
        'tr_a4': tr_a4,
        'tr_a5': tr_a5,
        'tr_a6': tr_a6,
        'tr_a7': tr_a7,
    })

# 2. Check: does tr(A^7) = 7*c7 trivially?
print(f"  tr(A^7) = 7*c7: {all(d['tr_a7'] == 7*d['c7'] for d in data)}")
print(f"  tr(A^3) = 3*c3: {all(d['tr_a3'] == 3*d['c3'] for d in data)}")

# 3. Now the key: express tr(A^7) in terms of tr(Sigma^k)
# tr(A^7) = 7*c7, so we need c7 from sigma traces
# But Sigma = f(A^2), so tr(Sigma^k) depends on A through A^2

# Let's check: does tr(M^k) determine c7?
print(f"\n--- Does tr(M^k) determine c7? ---")
mk_groups = defaultdict(set)
for d in data:
    mk_groups[(d['tr_m2'], d['tr_m3'])].add(d['c7'])
ambig = sum(1 for c7s in mk_groups.values() if len(c7s) > 1)
print(f"  (tr_m2, tr_m3) -> c7: {len(mk_groups)} groups, {ambig} ambiguous")

# 4. Cross-check: tr(Sigma^k) and tr(M^k) must be related by binomial expansion
# Sigma = c*J - I*c - M + I*c ? No: Sigma = c*J_off - M = c*(J-I) - M
# Since J_off = J - I:
# Sigma = c*J - c*I - M
# Hmm but Sigma has 0 diagonal, and c*J has c on diagonal, c*I has c on diagonal
# c*J[u,u] = c, c*I[u,u] = c, M[u,u] = 0
# So Sigma[u,u] = c - c - 0 = 0. Good.
# For u!=v: Sigma[u,v] = c*1 - 0 - M[u,v] = c - M[u,v]. Good.

# So Sigma = c*J - c*I - M
# = c*(J - I) - M
# = c*J_off - M

# Sigma^3 = (c*J_off - M)^3
# = c^3*J_off^3 - 3c^2*J_off^2*M + 3c*J_off*M^2 - M^3
# (actually need to be careful about non-commutative!)
# J_off*M != M*J_off in general!

# But J_off = J - I, and J commutes with everything symmetric
# J*M = M*J (since JX = X*J for any symmetric X? No! JX[u,v] = sum_w X[w,v])
# Actually J commutes with EVERY matrix: J*X[u,v] = sum_w X[w,v] = column_sum(X, v)
# X*J[u,v] = sum_w X[u,w] = row_sum(X, u)
# These are NOT equal in general!

# However, M is SYMMETRIC, so row_sum = column_sum transpose.
# JM[u,v] = sum_w M[w,v] = column_sum_v(M) = row_sum_v(M) (M symmetric)
# MJ[u,v] = sum_w M[u,w] = row_sum_u(M)
# These are NOT equal unless all row sums are equal.

# When are all row sums of M equal? M[u] = sum_v (A^2[u,v] + A^2[v,u]) = total 2-paths through u
# = sum_v (#{u->w->v} + #{v->w->u})
# = out-deg_2(u) + in-deg_2(u)
# = sum_v #{u->w->v} + sum_v #{v->w->u}
# First sum = s(u)*(n-1-s(u)) where s(u) = score
# Wait: sum_v A^2[u,v] = sum_v sum_w A[u,w]*A[w,v] = sum_w A[u,w] * sum_v A[w,v]
# = sum_w A[u,w] * s(w) = out-weighted-degree of u
# This is NOT s(u)*(n-1-s(u)) in general.

# Let me just compute
print(f"\n--- M row sums ---")
test_A = bits_to_adj(np.random.randint(0, 1 << total_bits), n)
test_M = test_A @ test_A + (test_A @ test_A).T
np.fill_diagonal(test_M, 0)
row_sums = [int(test_M[i].sum()) for i in range(n)]
print(f"  Row sums of M: {row_sums}")
print(f"  Are they all equal? {len(set(row_sums)) == 1}")
# M row sum for u = sum_v M[u,v] = sum_v (A^2[u,v] + A^2[v,u])
# = sum_v A^2[u,v] + sum_v A^2[v,u]
# = (A^2 * ones)[u] + ((A^2)^T * ones)[u]
# = (A^2 + (A^T)^2) * ones)[u]
# = ((A + A^T)^2 - A*A^T - A^T*A) * ones)[u]... getting complicated

# Row sum of M for u = sum_{v!=u} (n-2) - sigma(u,v) = (n-1)*(n-2) - sigma_deg(u)
# = 30 - sigma_deg(u) at n=7
test_Sigma = (n-2) * J_off - test_M
sigma_degs = [int(test_Sigma[i].sum()) for i in range(n)]
print(f"  Sigma degree seq: {sorted(sigma_degs)}")
print(f"  M row sums = 30 - sigma_deg: {all(row_sums[i] == 30 - sigma_degs[i] for i in range(n))}")

# 5. The BIG check: are sigma degree sequences equivalent to score sequences?
print(f"\n--- Sigma Degree vs Score Sequence ---")
sigma_score_map = defaultdict(set)
np.random.seed(42)
for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))

    # sigma_deg(u) = sum_v sigma(u,v) = sum_v (n-2 - A^2[u,v] - A^2[v,u])
    A2 = A @ A
    sig_degs = []
    for u in range(n):
        sd = sum((n-2) - int(A2[u][v]) - int(A2[v][u]) for v in range(n) if v != u)
        sig_degs.append(sd)
    sig_deg_sorted = tuple(sorted(sig_degs))

    sigma_score_map[scores].add(sig_deg_sorted)

n_score_multi = sum(1 for v in sigma_score_map.values() if len(v) > 1)
print(f"  Score sequences: {len(sigma_score_map)}")
print(f"  With multiple sigma deg seqs: {n_score_multi}")

if n_score_multi == 0:
    print(f"  => Sigma degree sequence DETERMINED by score sequence!")
else:
    for scores, sig_degs in sorted(sigma_score_map.items()):
        if len(sig_degs) > 1:
            print(f"  scores={scores}: {len(sig_degs)} different sigma deg seqs")
            for sd in sorted(sig_degs):
                print(f"    {sd}")

# 6. Sigma degree formula
# sigma_deg(u) = sum_{v!=u} sigma(u,v)
# = sum_{v!=u} (#{w: u->w, v->w} + #{w: w->u, w->v})
# = sum_{v!=u} sum_{w!=u,v} (A[u,w]*A[v,w] + A[w,u]*A[w,v])
# = sum_w sum_{v!=u,w} A[u,w]*A[v,w] + sum_w sum_{v!=u,w} A[w,u]*A[w,v]
# For first term: fix w. A[u,w] is 0 or 1.
#   If A[u,w]=1 (u->w): sum_{v!=u,w} A[v,w] = (n-1-s(w)) - A[u,w] = (n-1-s(w))-1
#     Wait: sum_v A[v,w] = in-degree of w = n-1-s(w). But we exclude v=u and v=w:
#     sum_{v!=u,w} A[v,w] = (n-1-s(w)) - A[u,w] = (n-1-s(w)) - 1
#   If A[u,w]=0 (w->u): contribution is 0
# First term = sum_{w: u->w} ((n-1-s(w)) - 1) = sum_{w: u->w} (n-2-s(w))
# Similarly: second term = sum_{w: w->u} ((s(w)) - 1 + ...) -- getting complicated

# Actually: sigma_deg(u) = sum_v sigma(u,v)
# total_sigma = (1/2) * sum_u sigma_deg(u) = sum(s^2) - n(n-1)/2
# So sum_u sigma_deg(u) = 2*(sum(s^2) - n(n-1)/2) = 2*sum(s^2) - n(n-1)
# For regular n=7: 2*7*9 - 42 = 126 - 42 = 84. sigma_deg(u) = 84/7 = 12 each.

print(f"\n--- Sigma Degree Formula ---")
print(f"  sum sigma_deg = 2*(sum(s^2) - n(n-1)/2)")
print(f"  For regular n=7: sum = 2*(63 - 21) = 84, each vertex = 12")

# More detailed: sigma_deg(u) = C(s(u), 2) + C(n-1-s(u), 2) for tournament?
# sigma_deg(u) = sum_{v!=u} (common succ + common pred between u and v)
# = sum_{v!=u} #{w: u->w, v->w} + sum_{v!=u} #{w: w->u, w->v}
# = sum_w #{v: u->w, v->w, v!=u,w} + sum_w #{v: w->u, w->v, v!=u,w}
# = sum_{w: u->w} (n-2-s(w)) + sum_{w: w->u} (s(w) - ... hmm)

# OK let me just verify: sigma_deg(u) = sum over v of sigma(u,v)
# = sum_v (common succ(u,v) + common pred(u,v))
# Common successors of u and v = #{w!=u,v: u->w AND v->w}
# Summing over all v: sum_v #{w: u->w, v->w, v!=u,w}
# = sum_w sum_{v!=u,w: v->w} A[u,w]
# = sum_{w: u->w} in_deg(w, excluding u and w) = sum_{w: u->w} ((n-1-s(w)) - A[u,w])
# Wait: in_deg of w excluding u means sum_{v!=u} A[v,w]
# Actually I also need v!=w, but A[w,w]=0 so:
# = sum_{w: u->w} (in_deg(w) - A[u,w]) = sum_{w: u->w} ((n-1-s(w)) - 1)
# Hmm, A[u,w] here means does u beat w as a "v" --- no, I'm computing from v's perspective.
# Let me restart more carefully.

# sigma_deg(u) = sum_{v!=u} common_succ(u,v) + sum_{v!=u} common_pred(u,v)

# Term 1: sum_{v!=u} common_succ(u,v)
# = sum_{v!=u} sum_{w!=u,v} A[u,w]*A[v,w]
# = sum_{w!=u} A[u,w] * sum_{v!=u,w} A[v,w]
# A[v,w] summed over v!=u,w = in_deg(w) - A[u,w] (subtracting w->w=0 and v=u case)
# Wait: in_deg(w) = sum_{v!=w} A[v,w]. If we also exclude v=u:
# sum_{v!=u,w} A[v,w] = in_deg(w) - A[u,w]
# = (n-1-s(w)) - A[u,w]
# But A[u,w] = 1 (since we're in the case w: u->w, so A[u,w]=1)
# = n-2-s(w)

# So Term 1 = sum_{w: u->w} (n-2-s(w)) = s(u)*(n-2) - sum_{w: u->w} s(w)
# Let me define: S_out(u) = sum of scores of vertices u beats = sum_{w: u->w} s(w)

# Term 2: sum_{v!=u} common_pred(u,v) = sum_{w: w->u} (s(w) - (1 - A[u,w]))
# Hmm, by symmetry: common_pred(u,v) = #{w: w->u AND w->v}
# = sum_{w!=u,v} A[w,u]*A[w,v]
# sum_{v!=u} of this = sum_{w!=u} A[w,u] * sum_{v!=u,w} A[w,v]
# = sum_{w: w->u} (s(w) - A[w,u]) = sum_{w: w->u} (s(w) - 1)
# Wait: A[w,u] = 0 when w->u (since A[u,w] would be... hmm)
# For w->u: that means A[w,u]=1. sum_{v!=u,w} A[w,v] = out_deg(w) - A[w,u] = s(w) - 1
# But we need A[w,u]=1 from "w beats u". If w->u, then A[w][u] = 1.
# sum_{v!=u,w} A[w,v] = s(w) - A[w,u] ... no, sum_{v!=w} A[w,v] = s(w). Exclude v=u:
# = s(w) - A[w,u] = s(w) - 1 (since w->u means A[w,u]=1)

# So Term 2 = sum_{w: w->u} (s(w) - 1) = S_in(u) - (n-1-s(u))
# where S_in(u) = sum of scores of vertices that beat u

# Total: sigma_deg(u) = [s(u)*(n-2) - S_out(u)] + [S_in(u) - (n-1-s(u))]

# Let me verify this formula
np.random.seed(42)
match = 0
for trial in range(500):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    A2 = A @ A
    scores = [int(sum(A[i])) for i in range(n)]

    for u in range(n):
        # Compute sigma_deg directly
        sd = sum((n-2) - int(A2[u][v]) - int(A2[v][u]) for v in range(n) if v != u)

        # Formula
        s_u = scores[u]
        S_out_u = sum(scores[w] for w in range(n) if w != u and A[u][w])
        S_in_u = sum(scores[w] for w in range(n) if w != u and A[w][u])

        formula = s_u*(n-2) - S_out_u + S_in_u - (n-1-s_u)
        if sd == formula:
            match += 1
        else:
            print(f"  MISMATCH: u={u}, sd={sd}, formula={formula}")
            break

print(f"\n  sigma_deg(u) = s(u)*(n-2) - S_out(u) + S_in(u) - (n-1-s(u))")
print(f"  Verified: {match}/{500*n}")

# Simplify: sigma_deg(u) = s(u)*(n-2) - S_out(u) + S_in(u) - n + 1 + s(u)
#                         = s(u)*(n-1) - n + 1 - S_out(u) + S_in(u)
# Note: S_out(u) + S_in(u) = sum_{w!=u} s(w) = n*(n-1)/2 - s(u) [total score minus own]
# Also: S_in(u) = n*(n-1)/2 - s(u) - S_out(u)
# So: S_in - S_out = n*(n-1)/2 - s(u) - 2*S_out(u)

# sigma_deg(u) = s(u)*(n-1) - n + 1 + n*(n-1)/2 - s(u) - 2*S_out(u)
#              = s(u)*(n-2) + 1 - n + n*(n-1)/2 - 2*S_out(u)
#              = s(u)*(n-2) + n(n-3)/2 + 1 - 2*S_out(u)

# At n=7: sigma_deg(u) = 5*s(u) + 14 + 1 - 2*S_out(u) = 5*s(u) + 15 - 2*S_out(u)

# For REGULAR tournament (all s(u)=3):
# S_out(u) = sum of scores of 3 beaten opponents
# sigma_deg(u) = 15 + 15 - 2*S_out(u) = 30 - 2*S_out(u)
# If S_out ranges... for regular T_7 (Paley): each vertex beats 3 with scores 3 each
# S_out = 9. sigma_deg = 30 - 18 = 12. Correct!

print(f"\n  Simplified: sigma_deg(u) = s(u)*(n-2) + n(n-3)/2 + 1 - 2*S_out(u)")
print(f"  At n=7: sigma_deg(u) = 5*s(u) + 15 - 2*S_out(u)")
print(f"  For regular (s=3, S_out=9): sigma_deg = 15 + 15 - 18 = 12. Correct!")

print(f"\n  Key insight: sigma_deg depends on score AND S_out (out-score sum).")
print(f"  S_out is NOT determined by score alone — it encodes 'who you beat'.")
print(f"  This is the Level 1.5 information: beyond scores, below full adjacency.")

print("\nDone.")
