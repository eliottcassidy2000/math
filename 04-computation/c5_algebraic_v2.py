"""
c5_algebraic_v2.py -- kind-pasteur-2026-03-13-S61

Correct algebraic analysis of c5 and lambda.

KEY IDENTITIES (for tournament, A[u][v]+A[v][u]=1, u!=v):

For each w in V\{u,v}, exactly ONE of four categories:
  d+(w): u->w->v   [contributes to (A^2)[u][v]]
  d-(w): v->w->u   [contributes to (A^2)[v][u]]
  out(w): u->w, v->w  [common successor]
  in(w):  w->u, w->v  [common predecessor]

So: d+ + d- + out + in = n-2
And: lambda(u,v) = d- when A[u][v]=1, or d+ when A[u][v]=0

Define:
  mu(u,v) = #{common successors} = out(u,v) = sum_w A[u][w]*A[v][w]
  (NOT symmetric: mu(u,v) counts w where both u->w and v->w)

Wait, actually mu IS symmetric since it's about common successors, which
doesn't care about direction between u and v.

So: (A^2)[u][v] + (A^2)[v][u] + mu + mu' = n-2
where mu = common successors, mu' = common predecessors.

This is getting complex. Let me just directly check whether c5 can be
expressed as a symmetric polynomial in the lambda entries.
"""

import numpy as np
from itertools import combinations, permutations
from numpy.linalg import lstsq, matrix_rank
from collections import defaultdict

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

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

n = 5
total_bits = n * (n-1) // 2
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]

# Correct identity check
print("=" * 60)
print("CORRECT IDENTITIES: (A^2)[u][v] decomposition")
print("=" * 60)

for bits in [0, 1, 42, 100, 1023]:
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    print(f"\nbits={bits}, scores={sorted(int(sum(A[i])) for i in range(n))}")
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            # Count the 4 categories
            dp = sum(1 for w in range(n) if w != u and w != v
                     and A[u][w] and A[w][v])
            dm = sum(1 for w in range(n) if w != u and w != v
                     and A[v][w] and A[w][u])
            out_c = sum(1 for w in range(n) if w != u and w != v
                        and A[u][w] and A[v][w])
            in_c = sum(1 for w in range(n) if w != u and w != v
                       and A[w][u] and A[w][v])

            assert dp + dm + out_c + in_c == n - 2
            assert A2[u][v] == dp

            if A[u][v] == 1:
                assert L[u][v] == dm, f"lambda check failed"
            else:
                assert L[u][v] == dp, f"lambda check failed"

print("All identities verified.")

# KEY INSIGHT: When A[u][v]=0 (v beats u), lambda(u,v) = (A^2)[u][v] = d+.
# When A[u][v]=1 (u beats v), lambda(u,v) = (A^2)[v][u] = d-.
# So: lambda(u,v) = A[v][u]*(A^2)[u][v] + A[u][v]*(A^2)[v][u]

# This means: (A^2)[u][v] = lambda(u,v) when A[v][u]=1 (i.e., A[u][v]=0)
# But when A[u][v]=1: lambda = d-, and d+ = n-2 - lambda - out - in.
# So (A^2)[u][v] = n-2 - lambda(u,v) - out(u,v) - in(u,v).

# The "out" and "in" counts are NOT determined by lambda alone.
# Define: sigma(u,v) = out(u,v) + in(u,v) = #{common successors} + #{common predecessors}
# sigma(u,v) = n-2 - (A^2)[u][v] - (A^2)[v][u]
# Note: sigma IS symmetric!

# Is sigma determined by lambda?
print(f"\n{'='*60}")
print("IS sigma(u,v) DETERMINED BY lambda(u,v)?")
print(f"{'='*60}")

lam_sigma_map = defaultdict(set)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    for u in range(n):
        for v in range(u+1, n):
            lam = L[u][v]
            sigma = n - 2 - A2[u][v] - A2[v][u]
            lam_sigma_map[lam].add(sigma)

print("lambda -> sigma values:")
for lam in sorted(lam_sigma_map.keys()):
    vals = sorted(lam_sigma_map[lam])
    print(f"  lambda={lam}: sigma in {vals}")

# sigma is NOT determined by lambda alone!
# But maybe (lambda, score_pair) or (lambda, vertex_degrees) determines sigma?

# Actually, what matters is whether tr(A^5) can be expressed using lambda.
# Let's try a completely different approach: SYMMETRIC polynomial in lambda entries.

print(f"\n{'='*60}")
print("FINDING c5 AS SYMMETRIC POLYNOMIAL IN LAMBDA ENTRIES")
print(f"{'='*60}")

# The S_5 action on pairs of vertices induces an action on the 10 lambda values.
# A symmetric polynomial must be invariant under this action.
#
# Degree-1 S_5-invariants: S1 = sum lambda[i][j]
# Degree-2 S_5-invariants: sum lambda^2, and products by overlap type
#   Type 0 (disjoint): pairs {i,j}, {k,l} with {i,j} cap {k,l} = empty
#     There are C(5,2)*C(3,2)/2 = 15 such unordered pair-pairs
#     But at n=5, two disjoint pairs use 4 of 5 vertices.
#   Type 1 (sharing): pairs {i,j}, {i,k} sharing exactly 1 vertex
#     There are 5*C(4,2) = 30 such ordered pairs, 30/... actually unordered = 30
#     Wait: for each vertex i, there are C(4,2)=6 pairs of pairs sharing vertex i.
#     Total: 5*6 = 30. But each unordered pair-pair is counted once. So 30.
#
# So the S_5-invariant quadratic ring has 3 generators:
#   sum L[i][j]^2 = S2
#   sum_{disjoint} L[ij]*L[kl] = D_prod
#   sum_{sharing_1} L[ij]*L[ik] = S_prod

# Note: sum_{all pairs of pairs} L[ij]*L[kl] = S1^2 - S2 (when including same-pair)
# Actually: (sum L)^2 = sum L^2 + 2*sum_{i<j} L_i*L_j
# where the second sum is over all unordered pairs of distinct elements.
# So sum_{distinct pairs of pairs} L[ij]*L[kl] = (S1^2 - S2)/2
# And D_prod + S_prod = (S1^2 - S2)/2

# So there are really only 3 independent degree-2 invariants: S1, S2, D_prod (or S_prod)
# Plus the constant 1. Total degree-2 space: 4 dimensions? Wait let's count properly.

# Degree-0 invariants: 1
# Degree-1 invariants: S1
# Degree-2 invariants: S2, D_prod (S_prod = (S1^2-S2)/2 - D_prod)
# So the degree <=2 invariant space has dimension 4.

# But we have 9 distinct groups! So degree <=2 can have at most 4 independent values,
# meaning up to 4^1 groups distinguished. Not enough for 9.

# Degree-3 invariants:
#   S3 = sum L[ij]^3
#   sum L[ij]^2 * L[kl] where {ij} disjoint {kl}
#   sum L[ij]^2 * L[ik] (sharing)
#   sum L[ij]*L[jk]*L[ik] (triangle in pairs = triangle in K5)
#   sum L[ij]*L[kl]*L[mn] where {ij},{kl},{mn} pairwise disjoint
#     -- IMPOSSIBLE at n=5 (would need 6 vertices)
#   sum L[ij]*L[ik]*L[jl] where {ij},{ik} share i; {jl} disjoint from {ik}?
#     -- these various overlap patterns

# Let me just build all symmetric features systematically.

# Actually, since the 9 groups correspond to distinct lambda multisets,
# and the multiset is determined by (n0, n1, n2, n3), and n0+n1+n2+n3=10,
# we can express c5 as a function of (n0, n1, n2, n3).
# And nk = sum_{i<j} [L[i][j] == k] = sum_{i<j} I_k(L[i][j])
# where I_k(x) = 1 if x=k, 0 otherwise.
# I_k is NOT polynomial... but we can write:
#   n0 = #{pairs with L=0}
#   n1 = #{pairs with L=1}
#   n2 = #{pairs with L=2}
#   n3 = #{pairs with L=3}
# And from the power sums: S1 = n1+2n2+3n3, S2 = n1+4n2+9n3, etc.

# The nk can be recovered from S1,...,S3 (and n0+...+n3=10):
# We have 4 unknowns (n0,...,n3) and 4 equations (sum=10, S1, S2, S3).
# The Vandermonde system:
# [1  1  1  1] [n0]   [10]
# [0  1  2  3] [n1] = [S1]
# [0  1  4  9] [n2]   [S2]
# [0  1  8 27] [n3]   [S3]

# This is a Vandermonde matrix with nodes 0,1,2,3 — invertible!
# So (n0,n1,n2,n3) <-> (10, S1, S2, S3) is a bijection.
# And c5 = f(n0,n1,n2,n3) = g(S1, S2, S3).

# But we showed that g is NOT polynomial of degree <= 2 in (S1,S2,S3).
# It IS determined by (S1,S2,S3) since the Vandermonde map is invertible.

# The non-polynomial nature is because indicator functions I_k(x) = [x==k]
# are not polynomial. But they CAN be expressed as Lagrange interpolation
# in the lambda value (which takes integer values 0,1,2,3):

# I_0(x) = (x-1)(x-2)(x-3)/(-6)
# I_1(x) = x(x-2)(x-3)/2
# I_2(x) = x(x-1)(x-3)/(-2)
# I_3(x) = x(x-1)(x-2)/6

# So n_k = sum_{i<j} I_k(L[i][j]) is a polynomial of DEGREE 3 in the L entries.
# And c5 = f(n0,...,n3) where f is a function we need to determine.

# From the 9-group table:
# (n0, n1, n2, n3, c5) =
# (10, 0, 0, 0, 0)
# ( 7, 3, 0, 0, 0)
# ( 5, 4, 1, 0, 0)
# ( 3, 5, 2, 0, 1)
# ( 3, 6, 0, 1, 1)
# ( 1, 6, 3, 0, 2)
# ( 0, 9, 0, 1, 3)
# ( 2, 4, 4, 0, 1)
# ( 0, 5, 5, 0, 2)

# Let's find c5 = f(n0,n1,n2,n3) as a polynomial.
# Using the constraint n0 = 10-n1-n2-n3, we work in (n1,n2,n3).

groups = [
    (10, 0, 0, 0, 0),
    (7, 3, 0, 0, 0),
    (5, 4, 1, 0, 0),
    (3, 5, 2, 0, 1),
    (3, 6, 0, 1, 1),
    (1, 6, 3, 0, 2),
    (0, 9, 0, 1, 3),
    (2, 4, 4, 0, 1),
    (0, 5, 5, 0, 2),
]

# Using Lagrange interpolation on the lambda indicator functions,
# c5 = sum over 10 pairs of f(L[i][j]) where f is at most degree 3.
# BUT c5 is NOT a separable function (sum of single-variable functions).
# The multiset structure means it's a SYMMETRIC function of the 10 values.

# From Newton's identities, any symmetric polynomial in the L values
# can be expressed in terms of the power sums S1, S2, ..., S10.
# But at n=5, L values in {0,1,2,3}, so S_k for k>=4 can be expressed
# via S1, S2, S3 and the nk.

# The key insight: c5 viewed as a function of (n0,n1,n2,n3) must be
# expressible in a way that works at ALL n (since THM-172 holds for all n).

# At n=5: lambda ranges 0..3 (since 3 witnesses).
# At n=7: lambda ranges 0..5.
# At n=9: lambda ranges 0..7.

# So the formula changes with n! Not a universal polynomial.

# Let's just find the exact formula at n=5 via interpolation of the 9 points.
print("\n--- Lagrange interpolation approach ---")
print("The 9 groups define c5 as a function of (n1, n2, n3):")

for n0, n1, n2, n3, c5 in groups:
    print(f"  f({n1}, {n2}, {n3}) = {c5}")

# Check: is there a simple closed form?
# Pattern: c5 ~ n2 for n3=0:
#   n2=0: c5=0 (two cases)
#   n2=1: c5=0
#   n2=2: c5=1
#   n2=3: c5=2
#   n2=4: c5=1
#   n2=5: c5=2
# Not simply n2-dependent.

# c5 = n3 + something?
#   n3=0: c5 in {0,0,0,1,2,1,2} — no
#   n3=1: c5 in {1,3} — varies

# Let me try: at n=5, c5 = number of directed 5-cycles = number of
# Hamiltonian directed cycles. For tournament on 5 vertices.

# Key formula: c5 = (1/5)*tr(A^5)
# = (1/5)*(sum of products along all closed 5-walks that are cycles)
# = (1/5)*(5! * fraction)

# Actually, let me try something more direct.
# Can we express c5 in terms of n3 and (n2 - something)?

# The c5=3 case: (n0,n1,n2,n3)=(0,9,0,1). 9 pairs with lambda=1, one with lambda=3.
# The lambda=3 pair is in ALL 4 three-cycles. The tournament has 3 directed 5-cycles.

# The c5=0 cases all have c3 <= 2 (S1 <= 6).
# c3 = S1/3 = (n1+2n2+3n3)/3

# c5 seems to jump at c3=3 (S1=9). Let me check:
# c3=0: c5=0, c3=1: c5=0, c3=2: c5=0, c3=3: c5=1, c3=4: c5={1,2,3}, c3=5: c5=2

# For c3=4: the discriminator is S3 (or equivalently n3).
# n3=0, n2=3: c5=2 (star arrangement)
# n3=0, n2=4: c5=1 (distributed)
# n3=1, n2=0: c5=3 (hub pair)

# So at c3=4: c5 = 3*n3 + max(0, n2-2) - [n2==4]? No...

# Let me just try a direct cubic fit on the full dataset
print(f"\n{'='*60}")
print("CUBIC FIT ON FULL DATASET")
print(f"{'='*60}")

all_feats = []
all_c5 = []
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5

    lam = [L[i][j] for i, j in pairs]

    # Symmetric features
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)

    # Vertex degree sums
    d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
    D2 = sum(dv**2 for dv in d_v)
    D3 = sum(dv**3 for dv in d_v)

    # Triangle sum = (1/6)*tr(L^3)
    trL3 = int(np.trace(L @ L @ L))

    # Disjoint pair product
    disj = sum(lam[a]*lam[b] for a in range(10) for b in range(a+1, 10)
               if len(set(pairs[a]) & set(pairs[b])) == 0)

    # Sharing pair product
    share = sum(lam[a]*lam[b] for a in range(10) for b in range(a+1, 10)
                if len(set(pairs[a]) & set(pairs[b])) == 1)

    # All cubics that might help:
    # S1*S2, S1*D2, S1*disj, S1*share, S1*trL3, S3
    feats = [1, S1, S2, S3, D2, D3, trL3, disj, share,
             S1*S2, S1*D2, S1**2, S2**2, S1*S3,
             S1*disj, S1*share, S1*trL3//6,
             D2*S1, disj*S1]
    all_feats.append(feats)
    all_c5.append(c5)

X = np.array(all_feats, dtype=float)
y = np.array(all_c5, dtype=float)
names = ['1', 'S1', 'S2', 'S3', 'D2', 'D3', 'trL3', 'disj', 'share',
         'S1*S2', 'S1*D2', 'S1^2', 'S2^2', 'S1*S3',
         'S1*disj', 'S1*share', 'S1*tri', 'D2*S1', 'disj*S1']

rank = matrix_rank(X)
print(f"Features: {len(names)}, rank: {rank}")

coeffs, _, _, _ = lstsq(X, y, rcond=None)
err = np.max(np.abs(X @ coeffs - y))
print(f"Max error: {err:.6f}")

if err < 0.001:
    # Backward elimination
    active = list(range(len(names)))
    while True:
        best_remove = None
        best_err = float('inf')
        for i in active:
            trial = [j for j in active if j != i]
            X_t = X[:, trial]
            c_t, _, _, _ = lstsq(X_t, y, rcond=None)
            e_t = np.max(np.abs(X_t @ c_t - y))
            if e_t < best_err:
                best_err = e_t
                best_remove = i
        if best_err < 0.001:
            active.remove(best_remove)
        else:
            break

    print(f"\nMinimal set ({len(active)} features):")
    min_names = [names[i] for i in active]
    print(f"  {min_names}")
    X_min = X[:, active]
    c_min, _, _, _ = lstsq(X_min, y, rcond=None)
    err_min = np.max(np.abs(X_min @ c_min - y))
    print(f"  Error: {err_min:.6f}")

    # Rational form
    for denom in range(1, 2001):
        int_c = np.round(c_min * denom)
        err_r = np.max(np.abs(X_min @ (int_c / denom) - y))
        if err_r < 0.001:
            print(f"\n  EXACT (denom={denom}):")
            terms = []
            for i, name in enumerate(min_names):
                c = int(int_c[i])
                if c != 0:
                    terms.append(f"({c:+d})*{name}")
            print(f"  5*c5 = ({' '.join(terms)}) / {denom/5:.0f}")
            print(f"  c5 = ({' '.join(terms)}) / {denom}")

            # Verify
            success = sum(1 for idx in range(len(y))
                         if abs(sum(int_c[i]*X_min[idx, i] for i in range(len(active)))/denom - y[idx]) < 0.001)
            print(f"  Verified: {success}/1024")
            break

# THE REALLY BIG QUESTION: does the formula generalize?
# At n=6, do we get c5_total = sum_V f(lambda_V)?
# Where lambda_V is the restricted lambda?

# We showed: lambda_V(u,v) = lambda(u,v) - witness(k,u,v)
# where k is the omitted vertex.

# So c5_total(T, n=6) = sum_{k=0}^{5} c5(T[V_k])
# where V_k = [6]\{k} and c5(T[V_k]) = f(lambda_V_k entries).

# And lambda_{V_k}(u,v) = lambda(u,v) - [k witnesses {u,v}]

# This is a beautiful structure! The total c5 at n=6 is a SUM of
# functions of RESTRICTED lambda, each of which differs from the
# FULL lambda by indicator corrections.

print(f"\n{'='*60}")
print("THE WITNESS STRUCTURE")
print(f"{'='*60}")

# At n=5: each pair {u,v} has n-2=3 potential witnesses.
# lambda(u,v) = number that ARE witnesses.
# For restricted lambda at n=6 (omitting vertex k):
# lambda_V(u,v) = lambda(u,v) - w(k,u,v) where w is the witness indicator.

# w(k,u,v) = 1 iff k is a witness for {u,v}
# = 1 iff (u,v,k) or (v,u,k) forms a directed 3-cycle
# = A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v]

# So lambda(u,v) = sum_{k != u,v} w(k,u,v).
# And w(k,u,v) depends on A, not just lambda.
# But the SUM is lambda. And the pattern of witnesses
# carries more information than lambda alone.

# At n=5: the witness MATRIX W[k][{u,v}] is a 5 x 10 binary matrix.
# Each row sum = score-related quantity.
# Each column sum = lambda(u,v).
# The row space might determine c5?

print("Witness matrix W[k][(u,v)] at n=5:")
for bits in [40, 120, 341]:  # interesting cases
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5

    W = np.zeros((n, 10), dtype=int)
    for pi, (u, v) in enumerate(pairs):
        for k in range(n):
            if k == u or k == v:
                continue
            w = int(A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v])
            W[k][pi] = w

    print(f"\nbits={bits}, c5={c5}, scores={sorted(int(sum(A[i])) for i in range(n))}")
    print(f"  lambda: {[int(L[i][j]) for i,j in pairs]}")
    print(f"  W (5x10 witness matrix):")
    for k in range(n):
        row = [int(W[k][pi]) for pi in range(10)]
        print(f"    k={k}: {row}  sum={sum(row)}")
    print(f"  Col sums = lambda: {[int(sum(W[k][pi] for k in range(n))) for pi in range(10)]}")

    # Check: col sums = lambda
    for pi, (u, v) in enumerate(pairs):
        col_sum = sum(W[k][pi] for k in range(n) if k != u and k != v)
        assert col_sum == L[u][v], f"Col sum check failed"

    # Row sums: for each vertex k, how many pairs does k witness?
    # w(k,u,v) = 1 iff k is in a 3-cycle with {u,v}
    # So row_sum(k) = #{pairs {u,v} that k is a witness for}
    #              = #{3-cycles containing k} * 2? No...
    # Each 3-cycle {k,u,v} has k witnessing {u,v}.
    # So row_sum(k) = #{3-cycles containing k} = delta(k).
    # And sum delta(k) = 3*c3 (each 3-cycle counted 3 times).
    # This is consistent with sum of all W entries = sum lambda = 3*c3 = S1.

    for k in range(n):
        row_sum = sum(int(W[k][pi]) for pi in range(10))
        # delta(k) = #{3-cycles containing k}
        delta_k = sum(1 for u in range(n) for v in range(u+1, n)
                      if u != k and v != k
                      and (A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v]))
        assert row_sum == delta_k, f"Row sum check failed for k={k}"
    print(f"  Row sums = delta(k): verified")

print("\nDone.")
