#!/usr/bin/env python3
"""
PROOF ATTEMPT: Uniform position distribution => M[a,b] = 0 for a != b at odd n.

KNOWN (diagonal):
  M[a,a] = sum_P (-1)^{pos(a,P)}
  If P[a,k] = H/n for all k, then M[a,a] = (H/n) sum_k (-1)^k.
  At odd n: sum_k (-1)^k = 1, so M[a,a] = H/n. ✓

UNKNOWN (off-diagonal):
  M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)  for a != b.
  Why does uniform position force this to vanish?

KEY OBSERVATION: E_a(S+a) counts paths in subtournament T[S+a] ending at a.
B_b(R+b) counts paths in T[R+b] starting at b. The product E*B counts
"split path pairs" — one ending at a, one starting at b, on complementary vertex sets.

The inclusion-exclusion sum adds these with alternating signs.

APPROACH 1: Express M[a,b] in terms of (v,k)-pair position counts.

Define: for a != b, and for each Ham path P, let j = pos(a,P), k = pos(b,P).
The path P naturally splits into:
  - A prefix of length j+1 ending at a (vertices at positions 0,...,j)
  - A suffix of length n-j-1... no, wait, a and b are both on the path.

Actually, the path visits vertices in order perm[0], perm[1], ..., perm[n-1].
Vertex a is at position j, vertex b is at position k. If we remove vertex a from
the path, the prefix (positions 0 to j-1) and suffix (positions j+1 to n-1)
are connected through a. Similarly for b.

But this doesn't directly relate to the E*B splitting, because in M[a,b],
the sets S and R partition [n]\\{a,b} — so BOTH a and b are "removed" and
placed in their respective sets.

APPROACH 2: Try a generating function argument.

Define f(a, x) = sum_S x^|S| E_a(S+a) = sum of x^|S| * (paths ending at a in T[S+a]).
Then M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) where R = U\\S.

This is a TWISTED convolution because S+R = U = [n]\\{a,b}, and the
convolution pairs (S, U\\S).

APPROACH 3: Compute M[a,b] directly for position-uniform tournaments
and show cancellation using position matrix uniformity.

Let's try: M[a,b] = sum over all "split decompositions" where
vertex a ends a path in S+a and vertex b starts a path in R+b.
Each such decomposition has a natural "position" interpretation:
the vertices in S+a precede a, the vertices in R+b follow b.

For a Hamiltonian path P = (v_0, v_1, ..., v_{n-1}):
  If we split at some position j (so v_j = a):
    The first j+1 vertices form a path in T ending at a.
    The remaining n-j-1 vertices form a path in T starting at v_{j+1}.

  But this gives a split (S+a, {v_{j+1},...,v_{n-1}}) where the second part
  starts at v_{j+1}, NOT at b.

The M[a,b] formula doesn't correspond to splitting a single path.
It's an INDEPENDENT pair of paths on complementary vertex sets.

APPROACH 4: Direct computation of the signed sum using position pairing.

For paths P in T[S+a] ending at a and paths Q in T[R+b] starting at b:
  Each term in M[a,b] is (-1)^|S| * (# such P) * (# such Q).

  The signed sum telescopes if we can find a sign-reversing involution
  on the pairs (P, Q, S, R) that maps (S, R) -> (S', R') with |S'| = |S| ± 1.

  Such an involution would move a single vertex u from R to S (or vice versa),
  and show that the new (P', Q') counts match.

  TOGGLE INVOLUTION: Pick vertex u in U. Define tau_u: S -> S xor {u}.
  If u in S: remove u from S, add to R. New: E_a(S\\{u}+a) vs E_a(S+a).
  The path count in S+a changes by removing one vertex (u).
  The path count in R+b changes by adding one vertex (u).

  For this to cancel: E_a(S+a) * B_b(R+b) = E_a(S\\{u}+a) * B_b(R+{u}+b)
  up to sign.

  This works if E_a(S+a) = E_a(S\\{u}+a) and B_b(R+b) = B_b(R+{u}+b),
  which is NEVER true (changing the vertex set changes the path count).

  So toggle involution on a single vertex doesn't directly cancel.

APPROACH 5: Examine the structure at n=3 (smallest odd n).

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np

def E_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def B_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count


# ============================================================
# n=3: Full analysis of M[a,b] for the 3-cycle tournament
# ============================================================
print("=" * 70)
print("n=3: The 3-cycle tournament (position-uniform)")
print("=" * 70)

n = 3
# 3-cycle: 0->1->2->0
T = {(0,1): 1, (1,0): 0, (1,2): 1, (2,1): 0, (2,0): 1, (0,2): 0}

# Ham paths
paths = []
for perm in permutations(range(n)):
    prod = 1
    for k in range(n-1):
        prod *= T.get((perm[k], perm[k+1]), 0)
    if prod > 0:
        paths.append(perm)

print(f"\n  Hamiltonian paths: {paths}")
print(f"  H = {len(paths)}")

# Position matrix
P = np.zeros((n, n), dtype=int)
for perm in paths:
    for k in range(n):
        P[perm[k], k] += 1
print(f"\n  Position matrix:")
for v in range(n):
    print(f"    v={v}: {list(P[v,:])}")

# Check uniform
is_uniform = all(P[v,k] == len(paths) // n for v in range(n) for k in range(n))
print(f"  Uniform: {is_uniform}")

# M[a,b] decomposition for a=0, b=1
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]
print(f"\n  M[{a},{b}] decomposition:")
print(f"  U = {U}")

total = 0
for mask in range(1 << len(U)):
    S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
    sign = (-1)**len(S_list)
    S_set = sorted(set(S_list) | {a})
    R_set = sorted(set(R) | {b})

    ea = E_v(T, S_set, a)
    bb = B_v(T, R_set, b)
    contrib = sign * ea * bb
    total += contrib
    print(f"    S={sorted(S_list)}: S+a={S_set}, R+b={R_set}, "
          f"E_{a}={ea}, B_{b}={bb}, sign={sign}, contrib={contrib}")

print(f"  Total M[{a},{b}] = {total}")

# At n=3, U has only 1 element!
# S = {} or S = {u}
# S={}: E_a({a})=1, B_b({b,u})=B_1({1,2})
#   There's a path 1->2 (T[1,2]=1). So B_1({1,2})=1.
#   contrib = (+1)*1*1 = 1
# S={u}: E_a({a,u})=E_0({0,2})
#   Path ending at 0 in {0,2}: 2->0 (T[2,0]=1). E_0({0,2})=1.
#   B_b({b})=B_1({1})=1.
#   contrib = (-1)*1*1 = -1
# Total = 1 - 1 = 0. ✓


# ============================================================
# Why does this cancel? Structural analysis.
# ============================================================
print("\n" + "=" * 70)
print("n=3: Why does M[0,1] = 0?")
print("=" * 70)

print("""
At n=3 with 3-cycle 0->1->2->0:
  M[0,1] = E_0({0}) * B_1({1,2}) - E_0({0,2}) * B_1({1})
         = 1 * 1 - 1 * 1 = 0

This cancels because:
  E_0({0}) = 1 (trivially)
  B_1({1,2}) = 1 (path 1->2)
  E_0({0,2}) = 1 (path 2->0)
  B_1({1}) = 1 (trivially)

The cancellation is between:
  Term 1: empty split — a alone, {b,u} together. B from b = path b->u = 1 path.
  Term 2: full split — {a,u} together, b alone. E to a = path u->a = 1 path.

The counts are both 1 because the tournament is a 3-cycle:
  - 0->1->2->0, so 1->2 exists (B_1({1,2})=1) and 2->0 exists (E_0({0,2})=1).

For the TRANSITIVE tournament 0->1, 0->2, 1->2:
  M[0,1]: E_0({0}) * B_1({1,2}) - E_0({0,2}) * B_1({1})
        = 1 * B_1({1,2}) - E_0({0,2}) * 1
  B_1({1,2}) = 1 (path 1->2)
  E_0({0,2}) = 1 (path 2->0? No! T[2,0]=0. So E_0({0,2})=0.)
  M[0,1] = 1*1 - 0*1 = 1

  Position matrix for transitive n=3:
  Only path: 0->1->2
  P = [[1,0,0],[0,1,0],[0,0,1]] — NOT uniform!

So transitive has non-uniform position => M not scalar => M[0,1] != 0.
3-cycle has uniform position => M scalar => M[0,1] = 0.

The KEY: For uniform position, there are EXACTLY as many paths with
u preceding a as there are with u following b. This makes the E*B
products match in absolute value, and the alternating sign kills them.
""")


# ============================================================
# n=5: Can we find a direct formula M[a,b] = f(position info)?
# ============================================================
print("=" * 70)
print("n=5: Direct proof attempt")
print("=" * 70)

# For n=5 Paley with uniform positions:
# At each |S|=s, the layer sum sum_{|S|=s} (-1)^s E_0(S+0) B_k(R+k) should vanish.
# But we showed in scalar_m_vt_proof.py that the layers are [1,-2,2,-1] and sum to 0.
# The individual layers are NOT zero, but the alternating sum of layers IS zero.

# For uniform position: P[v,k] = 3 for all v,k at n=5, H=15.
# The diagonal M[a,a] = 3 comes from 3*(1-1+1-1+1) = 3.

# For M[a,b] with a != b:
# Each Hamiltonian path P contributes to M[a,b] via the subset decomposition.
# The contribution of path P to M[a,b] depends on WHERE a and b appear in P.

# Claim: M[a,b] = sum_{j < k} [paths with a at j, b at k] * (-1)^j
#                - sum_{j > k} [paths with a at j, b at k] * (-1)^{j-1}
# ... or some similar formula. Let me derive it.

# From the definition: M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b).
# Each path P in T[S+a] ending at a has |S+a| vertices.
# Each path Q in T[R+b] starting at b has |R+b| vertices.
# |S+a| + |R+b| = n. The pair (P,Q) uses ALL vertices, with a at position |S|
# in P (the last position), and b at position 0 in Q (the first position).

# But (P,Q) is NOT a single Hamiltonian path! It's two separate paths
# on complementary vertex sets. There's no edge connecting P to Q.

# However, if we CONCATENATE P and Q (with an imaginary edge a->b),
# we'd get a sequence v_1,...,v_{|S|}, a, b, w_1,...,w_{|R|}.
# This uses all n vertices but is not necessarily a Hamiltonian path of T.

# The concatenation is a Ham path of T iff T[a,b]=1 AND the internal
# edges of P and Q are edges of T.

# So the "split path pairs" include SOME actual Ham paths (those where a->b)
# and SOME non-paths (where a doesn't beat b, or the concatenation has
# non-edges at the junction).

# This means M[a,b] is NOT simply a sum over Hamiltonian paths.
# It's a fundamentally different object.

# HOWEVER: at n=3, the total sum over split pairs DID correspond to
# a difference of 1-1=0. And for all n, if we SUM M[a,b] over all b,
# we get row sums.

# Let me compute: for each Hamiltonian path P of T, what is its
# "contribution" to M[a,b] for EACH pair (a,b)?

n = 5
from itertools import combinations as comb

def make_circulant_tournament(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

T = make_circulant_tournament(n, {1, 2})

# For each "split pair" (S, path_in_S+a, path_in_R+b),
# count how many correspond to actual Ham paths of T.
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

print(f"\n  Paley T_5, M[{a},{b}] split pairs:")
total_actual_paths = 0
total_non_paths = 0

for mask in range(1 << len(U)):
    S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
    sign = (-1)**len(S_list)
    S_set = sorted(set(S_list) | {a})
    R_set = sorted(set(R) | {b})

    # Find all (P,Q) pairs
    for p in permutations(S_set):
        if p[-1] != a: continue
        p_valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                p_valid = False; break
        if not p_valid: continue

        for q in permutations(R_set):
            if q[0] != b: continue
            q_valid = True
            for k in range(len(q)-1):
                if T.get((q[k], q[k+1]), 0) != 1:
                    q_valid = False; break
            if not q_valid: continue

            # This is a valid split pair
            concatenated = list(p) + list(q)
            # Check if concatenation is a Ham path
            is_ham = True
            for k in range(len(concatenated)-1):
                if T.get((concatenated[k], concatenated[k+1]), 0) != 1:
                    is_ham = False; break

            if is_ham:
                total_actual_paths += 1
                print(f"    S={sorted(S_list)}: P={p}->Q={q} IS Ham path, sign={sign}")
            else:
                total_non_paths += 1

print(f"\n  Split pairs that are Ham paths: {total_actual_paths}")
print(f"  Split pairs that are NOT Ham paths: {total_non_paths}")


# ============================================================
# KEY: What fraction of split pairs are actual Ham paths?
# ============================================================
print("\n  For EACH S-size, how many split pairs form Ham paths?")
for s_size in range(len(U)+1):
    actual = 0
    total = 0
    for S_list in combinations(U, s_size):
        S_list = list(S_list)
        R = [u for u in U if u not in S_list]
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        for p in permutations(S_set):
            if p[-1] != a: continue
            p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
            if not p_valid: continue

            for q in permutations(R_set):
                if q[0] != b: continue
                q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                if not q_valid: continue

                total += 1
                concatenated = list(p) + list(q)
                is_ham = all(T.get((concatenated[k], concatenated[k+1]), 0) == 1 for k in range(len(concatenated)-1))
                if is_ham:
                    actual += 1

    print(f"    |S|={s_size}: {actual} Ham paths out of {total} split pairs, "
          f"sign=(-1)^{s_size}={(-1)**s_size}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
