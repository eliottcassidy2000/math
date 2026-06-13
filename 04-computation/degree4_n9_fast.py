#!/usr/bin/env python3
"""
Fast degree-4 Fourier analysis at n=9.

Strategy: instead of computing over all 9! = 362880 permutations,
sample random tournaments and compute H(T) values, then extract
Fourier coefficients via Walsh-Hadamard on small subsets of edges.

Actually, the key question is: what is the DIMENSION of the degree-4
Fourier space of H(T) at n=9?

At n=7: dim = 2 (spanned by type P and type Q).
At n=9: we expect dim >= 2. The question is whether it's 2, 3, or more.

Approach:
1. Sample many random tournaments
2. For each, compute H(T) and the s_e = A_ij - 1/2 values
3. Extract degree-4 Walsh coefficients by computing H_hat(S) for |S|=4
4. Check rank of the resulting coefficient matrix

Since C(36,4) = 58905 is large, we'll check rank by sampling
a subset of degree-4 monomials.

Actually, the fastest approach: use the fact that degree-4 Fourier
coefficients can be expressed as:
  H_hat(S) = (1/2^m) sum_T (-1)^{sum_{e in S} A_e} H(T)

where the sum is over all 2^m tournaments. This is intractable.

Better: for each tournament T, compute
  w_4(T) = sum_{|S|=4} H_hat(S) * prod_{e in S} s_e

This is the degree-4 component of H(T) in the multilinear expansion.

For a SPECIFIC tournament T, w_4(T) can be computed from the definition:
  H(T) = sum_P product_{i=0}^{n-2} (1/2 + s_{P_i -> P_{i+1}})

Expand the product: H(T) = sum_{k=0}^{n-1} sum_{|S|=k} c_S(T) * prod_{e in S} s_e

where c_S(T) is the coefficient. For degree 4:
  w_4(T) = sum_{|S|=4} c_S * prod_{e in S} s_e

And c_S = number of Hamiltonian paths P such that exactly the edges in S
are "selected" from the path's edges.

Actually c_S = sum over all perms P of: (product of 1/2 for non-S edges on path)
times (product of s_e for S edges on path)... this is getting complicated.

Simpler: compute w_{n-1-2k} directly.
w_k(T) = sum_{|S|=k, S subset of path-edge-sets} product s_e

Wait, let's be precise. For each Hamiltonian path P = (v_0, ..., v_{n-1}),
define edge set E(P) = {(v_i, v_{i+1}) : i = 0,...,n-2}.
The expansion is:
H(T) = sum_P prod_{i=0}^{n-2} (1/2 + s_{E(P)_i})
     = sum_P sum_{S subset E(P)} (1/2)^{n-1-|S|} prod_{e in S} s_e

So the degree-k term is:
w_k(T) = (1/2)^{n-1-k} * sum_P sum_{S subset E(P), |S|=k} prod_{e in S} s_e

For each S of size k, the coefficient is:
c_S = (1/2)^{n-1-k} * #{paths P : S subset E(P)}

So the degree-4 Fourier space is determined by the function
S -> #{paths P : S subset E(P)} for 4-element edge sets S.

The rank of the degree-4 space = rank of the matrix
M[T, S] = #{paths of T containing all edges in S} for |S|=4.

But this depends on T! The FOURIER coefficients of H as a function
on the Boolean hypercube {0,1}^m are:
H_hat(S) = (1/2^m) * sum_T H(T) * chi_S(T)

These are FIXED numbers (not depending on T).

For degree 4, we need rank of H_hat restricted to |S|=4.

Let me compute this by sampling: pick random T, compute H(T), build
up the Walsh transform incrementally.

Hmm, 2^36 is too many. Let me think differently.

For n=9, there are C(36,4) = 58905 degree-4 Walsh coefficients.
Each H_hat(S) for |S|=4 is a real number. The space they span
(as functions of the tournament) has some dimension d.

Actually what we want is: expand H(T) as multilinear polynomial.
The degree-4 part is a homogeneous degree-4 polynomial in the 36
variables s_e. We want the dimension of this space.

But H(T) is a SINGLE fixed polynomial — it's not a space. The question
is: what is the dimension of the space of degree-4 monomials that
appear with nonzero coefficient?

No — the question from INV-050 is: the degree-4 part of the OCF
identity. If we write:
  w_{n-1-2k}(T) for the degree-(n-1-2k) Walsh component,
then w_4 at n=9 means the degree-4 Walsh component of H(T).

This is a multilinear polynomial in 36 variables. We want to decompose
it into "types" analogous to the n=7 case (Type P, Type Q).

The dimension question is: how many independent "types" span this polynomial.

Let me compute the degree-4 coefficients directly for small tournaments.

opus-2026-03-07-S37
"""
from itertools import combinations, permutations
from collections import defaultdict
import time

n = 9
edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
m = len(edges)  # 36
edge_idx = {e: k for k, e in enumerate(edges)}

def tournament_from_bits(bits):
    """Create adjacency from bit vector (0 = i->j forward, 1 = j->i backward)."""
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(edges):
        if bits[k]:
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A

def count_ham_paths(A):
    """Count Hamiltonian paths using DP (Held-Karp)."""
    # dp[mask][v] = number of paths visiting vertices in mask, ending at v
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))

def degree4_coefficients(A):
    """Compute degree-4 Walsh-Fourier coefficients of H(T).

    For each 4-element subset S of edges, compute the coefficient
    by enumerating all Hamiltonian paths.
    """
    # First, enumerate all Hamiltonian paths and their edge sets
    s_vals = [A[i][j] - 0.5 for i, j in edges]

    # For each path, compute its contribution to each degree-4 monomial
    # This is too slow for all C(36,4) = 58905 monomials

    # Instead, compute w_4(T) = sum over paths of sum over 4-subsets of path edges
    # of product of s values

    # For each path P, path has 8 edges. Choose 4 of them.
    # Contribution = (1/2)^4 * product of s values for those 4 edges
    # Actually (1/2)^{8-4} = (1/2)^4 times product of 4 s-values

    # Total: w_4(T) = (1/2)^4 * sum_P sum_{S subset E(P), |S|=4} prod_{e in S} s_e

    # But s_e in {-0.5, 0.5}, so prod of 4 s_e's = +/- (0.5)^4 = +/- 1/16
    # And the outer (1/2)^4 = 1/16
    # So each 4-subset contributes +/- 1/256 to w_4

    # Let's track which monomials (4-edge-subsets) have nonzero coefficients
    coeffs = defaultdict(int)  # monomial (tuple of 4 edge indices) -> signed count

    # Enumerate Hamiltonian paths via backtracking
    path_count = 0

    def backtrack(path, visited):
        nonlocal path_count
        if len(path) == n:
            path_count += 1
            # Get edge indices for this path
            path_edges = []
            for i in range(n-1):
                u, v = path[i], path[i+1]
                if u < v:
                    path_edges.append(edge_idx[(u, v)])
                else:
                    path_edges.append(edge_idx[(v, u)])

            # Sign of each s_e: +1 if forward (A[i][j]=1), -1 if backward
            signs = []
            for i in range(n-1):
                u, v = path[i], path[i+1]
                if u < v:
                    signs.append(1 if A[u][v] else -1)
                else:
                    signs.append(1 if A[v][u] else -1)
            # Wait: s_e = A[i][j] - 0.5 for edge (i,j) with i<j.
            # If path goes i->j (forward for edge), s_e = A[i][j] - 0.5
            #   If A[i][j]=1 (i beats j), s_e = 0.5
            #   If A[i][j]=0 (j beats i), s_e = -0.5
            # But the path goes i->j means i BEATS j (that's how we traverse).
            # So for a Hamiltonian path of T, each step is along an arc of T.
            # Thus A[path[i]][path[i+1]] = 1, always.
            # So s_{canonical_edge} = 0.5 if the arc direction matches canonical,
            #                         -0.5 otherwise.

            # Actually: for undirected edge (u,v) with u<v:
            # s_{u,v} = A[u][v] - 0.5
            # The path traverses u->v or v->u. Either way, the arc exists in T.
            # If u->v: A[u][v]=1, s = 0.5. If v->u: A[u][v]=0, s = -0.5.

            # For the path edge at position i: path[i]->path[i+1] is an arc of T.
            # If path[i] < path[i+1]: canonical edge is (path[i], path[i+1]),
            #   A[path[i]][path[i+1]]=1, s = 0.5. Sign = +1.
            # If path[i] > path[i+1]: canonical edge is (path[i+1], path[i]),
            #   A[path[i+1]][path[i]]=0 (since path[i]->path[i+1] means path[i] beats path[i+1]),
            #   s = -0.5. Sign = -1.

            path_signs = []
            for i in range(n-1):
                if path[i] < path[i+1]:
                    path_signs.append(1)  # forward = +
                else:
                    path_signs.append(-1)  # backward = -

            # Choose 4 edges from the 8 path edges
            for combo in combinations(range(n-1), 4):
                mono = tuple(sorted(path_edges[i] for i in combo))
                sign = 1
                for i in combo:
                    sign *= path_signs[i]
                coeffs[mono] += sign

            return

        v = path[-1]
        for u in range(n):
            if u not in visited and A[v][u]:
                path.append(u)
                visited.add(u)
                backtrack(path, visited)
                path.pop()
                visited.remove(u)

    for start in range(n):
        backtrack([start], {start})

    return dict(coeffs), path_count

# Test with a specific tournament
import random
random.seed(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

print("Computing degree-4 coefficients for random tournament on n=9...")
print(f"This enumerates all {n}! = {362880} orderings via backtracking")
t0 = time.time()
coeffs, H = degree4_coefficients(A)
t1 = time.time()

print(f"H(T) = {H}")
print(f"Time: {t1-t0:.1f}s")
print(f"Number of nonzero degree-4 monomials: {len(coeffs)}")

# Now classify monomials by graph structure
def classify_monomial(mono_edges):
    """Classify a 4-edge set by its graph structure."""
    elist = [edges[i] for i in mono_edges]
    verts = set()
    deg = defaultdict(int)
    for u, v in elist:
        verts.add(u)
        verts.add(v)
        deg[u] += 1
        deg[v] += 1
    nv = len(verts)
    deg_seq = tuple(sorted([deg[v] for v in verts], reverse=True))
    return (nv, deg_seq)

# Classify all nonzero monomials
type_counts = defaultdict(int)
type_coeffs = defaultdict(list)
for mono, coeff in coeffs.items():
    if coeff != 0:
        t = classify_monomial(mono)
        type_counts[t] += 1
        type_coeffs[t].append(coeff)

print(f"\nMonomial types (nonzero only):")
for t in sorted(type_counts.keys()):
    abs_coeffs = sorted(set(abs(c) for c in type_coeffs[t]))
    print(f"  {t}: {type_counts[t]} monomials, |coeffs| in {abs_coeffs[:5]}")
