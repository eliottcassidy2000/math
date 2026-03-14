#!/usr/bin/env python3
"""
cayley_independence.py — opus-2026-03-14-S75

The Hamiltonian paths of a tournament T form an INDEPENDENT SET
in the Cayley graph Γ_n = Cay(S_n, {s₁,...,s_{n-1}}).

This is the "permutohedron graph" — each vertex is a permutation,
adjacent iff they differ by an adjacent transposition.

Key questions:
1. What is the independence number α(Γ_n)?
2. Does max H = α(Γ_n) for tournaments?
3. What is the structure of maximum independent sets?
4. How does this relate to the OCF and the keys 2,3?
"""

from itertools import permutations
from math import factorial

def build_cayley_graph(n):
    """Build Cayley graph of S_n with adjacent transpositions."""
    perms = list(permutations(range(n)))
    idx = {p: i for i, p in enumerate(perms)}
    adj = [[] for _ in range(len(perms))]
    for i, p in enumerate(perms):
        for k in range(n-1):
            q = list(p)
            q[k], q[k+1] = q[k+1], q[k]
            j = idx[tuple(q)]
            adj[i].append(j)
    return perms, idx, adj

def max_ham_paths(n):
    """Find max H among all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    all_perms = list(permutations(range(n)))
    max_h = 0
    for bits in range(2**len(edges)):
        t = [[False]*n for _ in range(n)]
        for idx_e, (i,j) in enumerate(edges):
            if bits & (1 << idx_e):
                t[i][j] = True
            else:
                t[j][i] = True
        h = sum(1 for p in all_perms if all(t[p[k]][p[k+1]] for k in range(n-1)))
        max_h = max(max_h, h)
    return max_h

print("=" * 70)
print("PART 1: CAYLEY GRAPH PROPERTIES")
print("=" * 70)
print()

for n in range(2, 7):
    nf = factorial(n)
    perms, idx, adj = build_cayley_graph(n)
    degree = n - 1

    # Chromatic number of Cayley graph?
    # This is the Cayley graph of S_n with generators s_1,...,s_{n-1}
    # It's the 1-skeleton of the permutohedron

    print(f"  n={n}: |S_n|={nf}, degree={degree}")

    # Independence number by brute force (only feasible for small n)
    if nf <= 120:  # n≤5
        best = 0
        best_sets = []
        # For n≤4: brute force all subsets
        if nf <= 24:
            for mask in range(1 << nf):
                verts = [v for v in range(nf) if mask & (1 << v)]
                indep = True
                for v in verts:
                    for u in adj[v]:
                        if mask & (1 << u) and u > v:
                            indep = False
                            break
                    if not indep:
                        break
                if indep and len(verts) > best:
                    best = len(verts)
                    best_sets = [set(verts)]
                elif indep and len(verts) == best:
                    best_sets.append(set(verts))

            print(f"    α(Γ_{n}) = {best}")
            print(f"    Number of maximum independent sets: {len(best_sets)}")

            # Check which are achievable by tournaments
            max_h = max_ham_paths(n)
            print(f"    max H (tournaments) = {max_h}")
            print(f"    α(Γ_{n}) = max H? {best == max_h}")

        else:
            # n=5: too many subsets (2^120), use greedy
            max_h = max_ham_paths(n)
            print(f"    max H (tournaments) = {max_h}")
            # Lower bound: max_h ≤ α(Γ_n)
            # Upper bound: use fractional relaxation
            # Lovász theta?
            print(f"    α(Γ_{n}) ≥ {max_h} (from tournaments)")

    print()

print("=" * 70)
print("PART 2: STRUCTURE OF THE CAYLEY GRAPH")
print("=" * 70)
print()

# The Cayley graph Γ_n has some known properties:
# - It's vertex-transitive (S_n acts on itself)
# - It's (n-1)-regular
# - Its eigenvalues come from the representation theory of S_n

# The adjacency matrix eigenvalues of Cay(S_n, {s_1,...,s_{n-1}}):
# For the generator set {s_1,...,s_{n-1}} in S_n,
# eigenvalue for irrep λ is: Σ_i χ_λ(s_i) / dim(λ)
# This is: (n-1) times the content polynomial evaluation

# Actually, for Cayley graphs, the eigenvalues are:
# λ_ρ = (1/dim(ρ)) Σ_{g in S} χ_ρ(g)
# where S = generating set

# For S = {s_1,...,s_{n-1}} in S_n:
# χ_λ(s_i) where s_i = (i, i+1) is a transposition

import numpy as np

for n in range(2, 6):
    nf = factorial(n)
    perms, idx, adj_list = build_cayley_graph(n)

    # Build adjacency matrix
    A = np.zeros((nf, nf))
    for i in range(nf):
        for j in adj_list[i]:
            A[i][j] = 1

    # Eigenvalues
    eigvals = sorted(np.linalg.eigvalsh(A))
    print(f"  n={n}: eigenvalues of Cay(S_{n}): {[round(e,4) for e in eigvals]}")

    # Hoffman bound: α ≤ n * (-λ_min) / (λ_max - λ_min)
    # where λ_max = degree = n-1, λ_min = min eigenvalue
    lmax = max(eigvals)
    lmin = min(eigvals)
    hoffman = nf * (-lmin) / (lmax - lmin)
    print(f"    Hoffman bound: α ≤ {hoffman:.2f}")
    print()

print()

print("=" * 70)
print("PART 3: MAX H SEQUENCE AND INDEPENDENCE NUMBERS")
print("=" * 70)
print()

# Known max H values:
# n=1: H=1
# n=2: H=1
# n=3: H=3
# n=4: H=5
# n=5: H=15
# n=6: H=45
# n=7: H=189 (from OEIS A000568 context? or need to check)

max_h_known = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}
print("  Max H values (known):")
for n in sorted(max_h_known):
    print(f"    n={n}: max H = {max_h_known[n]}, n! = {factorial(n)}, ratio = {max_h_known[n]/factorial(n):.6f}")

print()
print("  RATIOS max_H/n!:")
for n in sorted(max_h_known):
    r = max_h_known[n] / factorial(n)
    print(f"    n={n}: {r:.6f}")

print()
print("  PATTERN: max H = (2n-1)!! / 2^{n-1}? Let me check.")
for n in range(1, 8):
    # Double factorial (2n-1)!!
    df = 1
    for k in range(1, 2*n, 2):
        df *= k
    val = df / (2**(n-1))
    print(f"    n={n}: (2n-1)!!={df}, (2n-1)!!/2^(n-1)={val}, max H={max_h_known.get(n,'?')}")

print()

# Actually, the max H for regular tournaments is known:
# For n odd (n=2m+1): max H = (2m-1)!! * ???
# Let me check: n=3: max H=3, n=5: max H=15, n=7: max H=189
# 3, 15, 189: ratios 15/3=5, 189/15=12.6
# Not obvious pattern. Let me factor:
# 3 = 3
# 15 = 3·5
# 189 = 3·63 = 3·7·9 = 3·7·3² = 27·7

print("  Factorizations of max H:")
from sympy import factorint
for n in sorted(max_h_known):
    h = max_h_known[n]
    f = factorint(h)
    fstr = " · ".join(f"{p}^{e}" if e > 1 else str(p) for p,e in sorted(f.items()))
    print(f"    n={n}: max H = {h} = {fstr}")

# Check: n=5 max H. Is it really 15?
# Regular tournament on 5 vertices: score sequence (2,2,2,2,2)
# H = 15 for the doubly-regular tournament? Yes, verified.

print()
print("  Max H for n=3,5,7 (odd):")
for n in [3, 5, 7]:
    h = max_h_known[n]
    nf = factorial(n)
    print(f"    n={n}: max H = {h}, n!/H = {nf//h}, n!/max_H = {nf/h:.2f}")

print()
print("  n!/max_H: 2, 8, 26.67...")
print("  Hmm: 6/3=2, 120/15=8, 5040/189≈26.67")
print("  2, 8, 26.67 → not obvious.")
print()
print("  But: 3 = 1·3, 15 = 3·5, 189 = 15·12.6?")
print("  Ratios: 15/3=5, 189/15=12.6")
print("  Not a clean pattern.")
print()

# Check: is max H at n=7 really 189?
# From exhaustive computation at n=7, or from literature?
# At n=7, the max score tournament has score sequence (3,3,3,3,3,3,3)
# For the quadratic residue tournament QR_7: H = ???
# Let me compute for the Paley tournament on 7 vertices

# QR_7: vertex set {0,...,6}, i→j iff j-i is a QR mod 7
# QRs mod 7: 1²=1, 2²=4, 3²=2. So QRs = {1,2,4}
# i→j iff j-i ∈ {1,2,4} mod 7

print("  Paley tournament QR_7:")
n = 7
qrs = {1, 2, 4}  # quadratic residues mod 7
adj = [[False]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % n in qrs:
            adj[i][j] = True

# Count Ham paths
from itertools import permutations as perms_iter
h = sum(1 for p in perms_iter(range(n)) if all(adj[p[k]][p[k+1]] for k in range(n-1)))
print(f"    H(QR_7) = {h}")
print(f"    Is {h} the maximum? (Expected 189)")

# Also check score sequence
scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
print(f"    Score sequence: {sorted(scores)}")

print()

# What about the doubly-regular tournament at n=7?
# There are multiple regular tournaments at n=7
# Let me try a few random ones

import random
random.seed(42)

max_h_found = h
for trial in range(100):
    # Generate a random regular tournament
    # Actually, let's generate all tournaments and find max H
    # No, 2^21 is too many for brute force in Python
    # Let me just check a few constructions

    # Rotational tournament: i→j iff j-i mod n ∈ S for fixed S
    # with |S| = (n-1)/2 = 3
    pass

# Check some specific regular tournaments
# S = {1,2,3}: i→j iff j-i ∈ {1,2,3} mod 7
for S_set in [{1,2,3}, {1,2,4}, {1,3,5}, {1,2,5}, {1,3,4}]:
    adj2 = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj2[i][j] = True
    scores2 = [sum(adj2[i][j] for j in range(n)) for i in range(n)]
    if all(s == 3 for s in scores2):
        h2 = sum(1 for p in perms_iter(range(n)) if all(adj2[p[k]][p[k+1]] for k in range(n-1)))
        print(f"  Rotational T with S={S_set}: H={h2}, scores={sorted(scores2)}")
        max_h_found = max(max_h_found, h2)

print(f"\n  Max H found at n=7: {max_h_found}")

print()
print("=" * 70)
print("PART 4: THE CAYLEY GRAPH AND THE 2-3 RECURRENCE")
print("=" * 70)
print()
print("  The eigenvalues of Cay(S_n, {s₁,...,s_{n-1}}) are determined by")
print("  the irreducible representations of S_n.")
print()
print("  For the standard representation (Young diagram [n-1,1]):")
print("    eigenvalue = (n-1) - 2·(content sum)")
print("    Hmm, let me think about this more carefully.")
print()
print("  The adjacency matrix of Cay(S_n, C) where C = {s_1,...,s_{n-1}}:")
print("  A = Σ_i L(s_i) where L(g) is the left-regular representation matrix of g")
print()
print("  Eigenvalues: For each irrep ρ of S_n of dimension d_ρ,")
print("  the eigenvalue is λ_ρ = Σ_i χ_ρ(s_i) / d_ρ")
print("  with multiplicity d_ρ².")
print()
print("  For the trivial rep: λ = (n-1)·1/1 = n-1 (the degree)")
print("  For the sign rep: λ = (n-1)·(-1)/1 = -(n-1)")
print()

# Compute character values of transpositions for irreps of S_n
# For S_3: irreps are [3] (trivial), [2,1] (standard), [1,1,1] (sign)
# χ_{[3]}(s_i) = 1, χ_{[2,1]}(s_i) = 0, χ_{[1,1,1]}(s_i) = -1
# So eigenvalues: [3]: 2, [2,1]: 0 (mult 4), [1,1,1]: -2

# For S_4: irreps [4],[3,1],[2,2],[2,1,1],[1,1,1,1]
# dims: 1, 3, 2, 3, 1
# χ for transposition (i,i+1):
# [4]: 1, [3,1]: 1 (content!), [2,2]: 0, [2,1,1]: -1, [1,1,1,1]: -1
# eigenvalue = 3·χ/d
# [4]: 3·1/1 = 3
# [3,1]: 3·1/3 = 1
# [2,2]: 3·0/2 = 0
# [2,1,1]: 3·(-1)/3 = -1
# [1,1,1,1]: 3·(-1)/1 = -3

print("  Eigenvalues of Cay(S_n, adjacent transpositions):")
print()
print("  S_3: irreps [3],[2,1],[1,1,1]")
print("    dims: 1, 2, 1")
print("    χ(s_i): 1, 0, -1")
print("    eigenvalues: 2, 0, -2 (with mults 1, 4, 1)")
print()
print("  S_4: irreps [4],[3,1],[2,2],[2,1,1],[1,1,1,1]")
print("    dims: 1, 3, 2, 3, 1")
print("    χ(s_i): 1, 1, 0, -1, -1")
print("    λ = (n-1)χ/d: 3, 1, 0, -1, -3")
print("    mults: 1, 9, 4, 9, 1")
print()

# Hoffman bound: α ≤ N · (-λ_min) / (λ_max - λ_min)
for n, lmin, lmax, N in [(3, -2, 2, 6), (4, -3, 3, 24), (5, -4, 4, 120)]:
    hoffman = N * (-lmin) / (lmax - lmin)
    print(f"  n={n}: Hoffman ≤ {hoffman:.0f}, max H = {max_h_known[n]}")

print()
print("  Hoffman gives α ≤ n!/2 for all n (since λ_min = -(n-1), λ_max = n-1).")
print("  This is a WEAK bound. The actual max H is much smaller.")
print()
print("  Better eigenvalue-based bound uses Lovász theta:")
print("  θ = N · (-λ_min) / (λ_max - λ_min)")
print("  But the Cayley graph has symmetric eigenvalues, so θ = N/2.")
print()

# The actual max H values grow much slower than n!/2:
print("  max H / (n!/2):")
for n in sorted(max_h_known):
    r = max_h_known[n] / (factorial(n) / 2)
    print(f"    n={n}: {r:.6f}")

print()
print("  The fraction max_H/(n!/2) decreases rapidly.")
print("  By Stirling: max_H ≈ n! / (e · ...) but the exact asymptotics")
print("  are not known for tournament max H.")
print()

print("=" * 70)
print("PART 5: INDEPENDENCE NUMBER = MAX H?")
print("=" * 70)
print()
print("  Question: Does every maximum independent set in Cay(S_n)")
print("  come from a tournament?")
print()
print("  A tournament selects H permutations with the constraint:")
print("  For each pair (i,j), EITHER all selected perms have i before j,")
print("  OR all have j before i.")
print()
print("  Wait, that's NOT what tournament selection does.")
print("  A selected perm σ satisfies σ(k)→σ(k+1) for all k.")
print("  Different selected perms can have different relative orders of i,j!")
print()
print("  The constraint is: no two selected perms differ by an adjacent swap.")
print("  This is the INDEPENDENT SET condition, which is WEAKER than")
print("  'coming from a tournament'.")
print()
print("  So: not every max independent set in Cay(S_n) comes from a tournament.")
print("  Tournament independent sets are a SUBSET of all independent sets.")
print()
print("  Is α(Cay(S_n)) > max_H?")
print()

# At n=3: α = 3 = max_H ✓ (verified by brute force)
# At n=4: let's verify α = 5
# We verified α = 5 for tournaments. Is there a non-tournament independent set of size 6?

n = 4
perms_4, idx_4, adj_4 = build_cayley_graph(4)
nf = factorial(4)

# Check: is there an independent set of size 6?
# Use greedy search
best_found = 0
best_set = None

import random
random.seed(0)

for trial in range(10000):
    # Random greedy
    remaining = list(range(nf))
    random.shuffle(remaining)
    selected = []
    for v in remaining:
        conflict = False
        for u in selected:
            if u in adj_4[v]:
                conflict = True
                break
        if not conflict:
            selected.append(v)
    if len(selected) > best_found:
        best_found = len(selected)
        best_set = selected

print(f"  n=4: greedy max independent set found: {best_found}")
# Verify it's actually independent
is_indep = True
for i in range(len(best_set)):
    for j in range(i+1, len(best_set)):
        if best_set[j] in adj_4[best_set[i]]:
            is_indep = False
print(f"    Is independent? {is_indep}")
print(f"    max H = 5")
print(f"    α(Cay(S_4)) ≥ {best_found}")

# Brute force at n=4 (2^24 subsets... too many)
# But we can use max independent set algorithm
# Actually for n=4, we already computed α=5 in the previous script
# Let's verify more carefully

# Check if size 6 is achievable
# The permutohedron at n=4 is 3-regular on 24 vertices
# Max independent set of a 3-regular graph on 24 vertices:
# by Brook's theorem, χ ≤ 3 (unless complete or odd cycle)
# So α ≥ 24/3 = 8
# But we found only 5 from tournaments!

# Wait: α ≥ ceiling(N/Δ) where Δ = max degree? No.
# α ≥ N/(Δ+1) by greedy. Here: 24/4 = 6.
# So α(Cay(S_4)) ≥ 6!

# But we found 5 from tournaments. So maybe α > max_H?

# Let's check: our greedy found 6? Let me verify.
if best_found >= 6:
    print(f"\n  FOUND independent set of size {best_found} > max_H=5!")
    print(f"    Vertices: {[perms_4[v] for v in best_set[:8]]}")

    # Check: is this achievable by a tournament?
    # For each pair (i,j), check relative orders in selected perms
    # If both orders appear, it's NOT from a tournament
    for v1 in range(n):
        for v2 in range(v1+1, n):
            orders = set()
            for vi in best_set:
                p = perms_4[vi]
                pos1 = p.index(v1) if v1 in p else -1  # always in p
                pos2 = p.index(v2) if v2 in p else -1
                # wait, p is a permutation so find positions
                pos1 = list(p).index(v1)
                pos2 = list(p).index(v2)
                if pos1 < pos2:
                    orders.add('before')
                else:
                    orders.add('after')
            if len(orders) == 2:
                # v1 appears both before and after v2
                print(f"    Pair ({v1},{v2}): both orders appear → NOT from tournament")
                break
        else:
            continue
        break

print()
print("=" * 70)
print("PART 6: MAX H IS NOT α(Cayley)!")
print("=" * 70)
print()
print("  The independence number α(Cay(S_n, adjacent transpositions))")
print("  is LARGER than max H from tournaments for n ≥ 4.")
print()
print("  Tournament Ham paths form a STRUCTURED independent set:")
print("  the constraint that arcs are consistent is ADDITIONAL")
print("  beyond just being non-adjacent in the Cayley graph.")
print()
print("  The GAP between max_H and α(Cayley) measures the 'cost'")
print("  of the tournament consistency constraint.")
print()
print("  α(Cay(S_n)) ≥ n!/(n) = (n-1)! (by greedy, since degree = n-1)")
print("  max H ≤ (2n-1)!! or similar (tournament-specific bound)")
print()
print("  The tournament constraint is like a COLORING constraint:")
print("  not just any independent set, but one that's 'arc-consistent'.")
print()

# Summary table
print("  Summary:")
print(f"  {'n':>3s} {'n!':>6s} {'max_H':>6s} {'α≥':>6s} {'n!/(n)':>8s}")
for n in range(2, 8):
    nf = factorial(n)
    mh = max_h_known.get(n, '?')
    alpha_lb = nf // n  # greedy lower bound
    print(f"  {n:3d} {nf:6d} {str(mh):>6s} {alpha_lb:6d} {nf/n:8.1f}")

print()
print("  The greedy bound α ≥ (n-1)! is much larger than max H for n ≥ 5.")
print("  So tournaments use only a SMALL fraction of the Cayley graph's")
print("  independence capacity.")
