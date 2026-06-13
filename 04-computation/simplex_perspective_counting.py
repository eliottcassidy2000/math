#!/usr/bin/env python3
"""
Explore the "perspective counting" observation:
  - A tournament on n vertices is a SIMPLEX with binary edge labels
  - Each vertex has a "perspective" = its view of the tournament
  - # distinct perspectives per tournament = # vertex orbits under Aut(T)

User's observation:
  n=3: transitive has 3 perspectives, cyclic has 1. 3+1=4.
  n=4: paired classes contribute 2 each, unpaired contribute 4 each. 2+2+4+4=12.

Question: does this sum give a nice formula? Does it connect n to n+1?

Also explore: tournaments as simplices in cubes, packing triangles, tetrahedra.

kind-pasteur-2026-03-06-S25g
"""

from itertools import permutations, combinations
from collections import defaultdict
import sys

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def tournament_canonical(A):
    """Return a canonical form for isomorphism testing."""
    n = len(A)
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def automorphism_group_size(A):
    """Count automorphisms of tournament A."""
    n = len(A)
    count = 0
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            count += 1
    return count

def vertex_orbits(A):
    """Compute number of vertex orbits under Aut(T)."""
    n = len(A)
    # Find all automorphisms
    auts = []
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] != A[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            auts.append(perm)

    # Find orbits using union-find
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for aut in auts:
        for v in range(n):
            union(v, aut[v])

    return len(set(find(v) for v in range(n)))

def vertex_perspective(A, v):
    """The 'perspective' of vertex v: its out-degree sequence to each orbit."""
    n = len(A)
    # Simple perspective: sorted tuple of (out-neighbor scores)
    out_neighbors = tuple(sorted(sum(A[j]) for j in range(n) if A[v][j] == 1))
    in_neighbors = tuple(sorted(sum(A[j]) for j in range(n) if A[j][v] == 1))
    return (sum(A[v]), out_neighbors, in_neighbors)

def is_self_converse(A):
    """Check if T is isomorphic to T^op."""
    n = len(A)
    # Build T^op
    Aop = [[A[j][i] for j in range(n)] for i in range(n)]
    return tournament_canonical(A) == tournament_canonical(Aop)

print("=" * 70)
print("PERSPECTIVE COUNTING: TOURNAMENTS AS SIMPLICES")
print("=" * 70)

for n in range(3, 8):
    m = n * (n-1) // 2  # number of edges
    total_labeled = 2**m

    # Enumerate isomorphism classes
    canon_to_rep = {}
    for bits in range(total_labeled):
        A = tournament_from_bits(n, bits)
        c = tournament_canonical(A)
        if c not in canon_to_rep:
            canon_to_rep[c] = A

    iso_classes = list(canon_to_rep.values())

    total_orbits = 0
    total_labeled_check = 0

    if n <= 6:  # Only feasible for small n
        print(f"\n--- n={n} ---")
        print(f"{'class':>5} {'|Aut|':>5} {'orbits':>6} {'n!/|Aut|':>8} {'SC':>3} {'scores':>20}")

        for idx, A in enumerate(iso_classes):
            aut_size = automorphism_group_size(A)
            n_orbits = vertex_orbits(A)
            n_labelings = len(list(permutations(range(n)))) // aut_size
            sc = is_self_converse(A)
            scores = tuple(sorted(sum(A[i]) for i in range(n)))

            total_orbits += n_orbits
            total_labeled_check += n_labelings

            print(f"{idx+1:5d} {aut_size:5d} {n_orbits:6d} {n_labelings:8d} {'Y' if sc else 'N':>3} {str(scores):>20}")

        print(f"\nTotal vertex orbits across all iso classes: {total_orbits}")
        print(f"Total labeled tournaments (check): {total_labeled_check} = 2^{m} = {total_labeled}")
        print(f"Number of iso classes: {len(iso_classes)}")

        # Key ratios
        print(f"orbits/classes = {total_orbits}/{len(iso_classes)} = {total_orbits/len(iso_classes):.4f}")
        print(f"orbits * something = ???")

        # Check user's claim
        if n == 3:
            print(f"\nUser says: transitive=3, cyclic=1, sum=4. Got: {total_orbits}")
        elif n == 4:
            print(f"\nUser says: 2+2+4+4=12. Got: {total_orbits}")
    else:
        print(f"\n--- n={n}: {len(iso_classes)} iso classes (not printing details) ---")

# Now compute the sequence of total vertex orbits
print(f"\n{'='*70}")
print("SEQUENCE OF TOTAL VERTEX ORBITS")
print(f"{'='*70}")
orbit_sequence = []
for n in range(2, 7):
    m = n*(n-1)//2
    canon_to_rep = {}
    for bits in range(2**m):
        A = tournament_from_bits(n, bits)
        c = tournament_canonical(A)
        if c not in canon_to_rep:
            canon_to_rep[c] = A

    total = sum(vertex_orbits(A) for A in canon_to_rep.values())
    orbit_sequence.append(total)
    print(f"  n={n}: {total} total vertex orbits across {len(canon_to_rep)} classes")

print(f"\nSequence: {orbit_sequence}")
print(f"Ratios: {[orbit_sequence[i+1]/orbit_sequence[i] for i in range(len(orbit_sequence)-1)]}")

# Check: is this sequence in OEIS?
print(f"\nLook for sequence in OEIS: {orbit_sequence}")

# Also: count total DISTINCT perspectives across all classes
print(f"\n{'='*70}")
print("DISTINCT PERSPECTIVES (unique vertex views)")
print(f"{'='*70}")
for n in range(3, 7):
    m = n*(n-1)//2
    canon_to_rep = {}
    for bits in range(2**m):
        A = tournament_from_bits(n, bits)
        c = tournament_canonical(A)
        if c not in canon_to_rep:
            canon_to_rep[c] = A

    # For each iso class, compute the set of distinct vertex perspectives
    all_perspectives = set()
    for A in canon_to_rep.values():
        for v in range(n):
            persp = vertex_perspective(A, v)
            all_perspectives.add((tournament_canonical(A), persp))

    # Total distinct perspectives = total vertex orbits (by definition!)
    total_orbits = sum(vertex_orbits(A) for A in canon_to_rep.values())
    print(f"  n={n}: {total_orbits} vertex orbits = {len(all_perspectives)} distinct (canon, perspective) pairs")

print("\nDONE")
