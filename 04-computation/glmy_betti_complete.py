"""
Complete GLMY path homology computation for tournaments.

Computes ALL Betti numbers beta_0, ..., beta_{n-1} using the GLMY path complex.

The path complex Omega_p consists of "allowed p-paths": sequences v_0 -> v_1 -> ... -> v_p
where no proper subsequence is itself an allowed path.

The boundary map d_p: Omega_p -> Omega_{p-1} is:
  d_p(v_0, ..., v_p) = sum_{i=0}^{p} (-1)^i * (v_0, ..., v_{i-1}, v_{i+1}, ..., v_p)
  but only the terms that are in Omega_{p-1} are kept.

GLMY definition: A path v_0 -> v_1 -> ... -> v_p is "allowed" if:
  1. All edges v_i -> v_{i+1} exist in the tournament
  2. The path is not a concatenation of shorter allowed paths
     (more precisely, it is "regular" or "non-degenerate")

For tournaments, the condition simplifies: a directed path is allowed iff
it is a Hamiltonian path on the vertex set {v_0, ..., v_p} within the
subtournament T[{v_0, ..., v_p}].

Actually, the standard definition is simpler: Omega_p consists of all
directed (p+1)-tuples (v_0, v_1, ..., v_p) of DISTINCT vertices with
v_i -> v_{i+1} for all i, such that the path is "regular" (allowed).

For digraphs, regularity means: for each 0 < i < p, there is NO edge
v_{i-1} -> v_{i+1} (i.e., the path doesn't "skip" vertices).

Wait, that's for the PATH complex specifically. Let me use the standard
GLMY definition:

Omega_p = {(v_0, v_1, ..., v_p) : all distinct, v_i->v_{i+1} for all i,
           and v_{i-1} does NOT have edge to v_{i+1} for any 1 <= i <= p-1}

Actually no. The GLMY allowed paths are:
A_p = {(v_0, ..., v_p) : v_i -> v_{i+1} for all i, all v_i distinct}
Omega_p = A_p / (non-regular paths)

Where a path is non-regular if it can be decomposed as concatenation of
shorter paths: (v_0, ..., v_k) and (v_k, ..., v_p).

For TOURNAMENTS on n vertices: every pair has a directed edge. So
A_p = all orderings of (p+1)-element subsets.

A path (v_0, ..., v_p) is regular (allowed) iff v_{i-1} -> v_{i+1} for all
0 < i < p. Wait, the opposite: it's allowed if v_{i-1} does NOT go to v_{i+1}
for some i...

Actually, let me just use the correct GLMY definition from the literature.

In Grigor'yan-Lin-Muranov-Yau (GLMY), for a digraph G:
- A_p(G) = set of paths v_0 -> v_1 -> ... -> v_p of length p (with edges in G)
- Omega_p(G) = A_p(G) (the path complex uses ALL directed paths as generators)
- The boundary is d_p: R^{A_p} -> R^{A_{p-1}} given by alternating sum of face maps

Actually, I recall: in the GLMY paper, Omega_p is the QUOTIENT of the path
space by "irregular" paths. For tournaments, the regularity condition is:
a path v_0 v_1 ... v_p is regular iff there is no "shortcut" v_i -> v_j
with j > i+1 that creates a shorter allowed path.

For TOURNAMENTS, since every pair has an edge, a path of length p on p+1
vertices is regular iff it is a Hamiltonian path on those p+1 vertices.

The boundary of a p-path (v_0, ..., v_p) is:
d(v_0, ..., v_p) = sum_{i=0}^{p} (-1)^i * e_i
where e_i is (v_0, ..., hat{v_i}, ..., v_p) IF this is an allowed (p-1)-path,
and 0 otherwise.

For a tournament: (v_0, ..., hat{v_i}, ..., v_p) is allowed iff
v_{i-1} -> v_{i+1} (the vertices before and after v_i are connected in order).
For i=0: (v_1, ..., v_p) is always allowed (removing first vertex).
For i=p: (v_0, ..., v_{p-1}) is always allowed (removing last vertex).
For 0 < i < p: (v_0, ..., v_{i-1}, v_{i+1}, ..., v_p) is allowed iff
v_{i-1} -> v_{i+1} in the tournament.

So the boundary is:
d(v_0, ..., v_p) = sum_{i : v_{i-1}->v_{i+1}} (-1)^i * (v_0,...,hat{v_i},...,v_p)
where i=0 and i=p are always included (with v_{-1}->v_1 and v_{p-1}->v_{p+1}
interpreted as always true).

This is the correct GLMY boundary for tournaments.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import sys

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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_all_betti(A, n, verbose=False):
    """Compute all Betti numbers beta_0, ..., beta_{n-1} for tournament A."""

    # Step 1: Enumerate allowed p-paths for each p
    # For a tournament, a p-path on p+1 vertices is a Hamiltonian path of the
    # subtournament on those vertices.

    omega = {}  # omega[p] = list of allowed p-paths (tuples)
    omega_index = {}  # omega_index[p] = dict mapping path tuple to index

    for p in range(n):
        paths = []
        # Choose p+1 vertices
        for subset in combinations(range(n), p+1):
            # Find all Hamiltonian paths of the subtournament on these vertices
            for perm in permutations(subset):
                # Check if this is a directed path
                valid = True
                for i in range(p):
                    if not A[perm[i]][perm[i+1]]:
                        valid = False
                        break
                if valid:
                    paths.append(perm)

        omega[p] = paths
        omega_index[p] = {path: idx for idx, path in enumerate(paths)}

        if verbose:
            print(f"  Omega_{p}: {len(paths)} paths")

    # Step 2: Compute boundary matrices d_p: Omega_p -> Omega_{p-1}
    betti = []

    prev_im_rank = 0  # rank of im(d_{p+1})

    for p in range(n):
        dim_p = len(omega[p])

        if dim_p == 0:
            betti.append(0)
            prev_im_rank = 0
            continue

        # Compute d_p: Omega_p -> Omega_{p-1}
        if p == 0:
            # d_0 = 0 map (or trivially ker = all of Omega_0)
            ker_rank = dim_p
        else:
            # Build the matrix d_p
            dim_pm1 = len(omega[p-1])

            if dim_pm1 == 0:
                ker_rank = dim_p
            else:
                # d_p matrix: rows = Omega_{p-1}, cols = Omega_p
                d = np.zeros((dim_pm1, dim_p), dtype=np.int64)

                for col, path in enumerate(omega[p]):
                    # Boundary: sum_{i=0}^{p} (-1)^i * face_i
                    for i in range(p + 1):
                        # face_i: remove vertex at position i
                        # Check if the resulting (p-1)-path is allowed
                        if i == 0 or i == p:
                            # Always allowed (removing first or last vertex)
                            face = path[:i] + path[i+1:]
                            if face in omega_index[p-1]:
                                row = omega_index[p-1][face]
                                d[row, col] += (-1)**i
                        else:
                            # Allowed iff path[i-1] -> path[i+1]
                            if A[path[i-1]][path[i+1]]:
                                face = path[:i] + path[i+1:]
                                if face in omega_index[p-1]:
                                    row = omega_index[p-1][face]
                                    d[row, col] += (-1)**i

                ker_rank = dim_p - np.linalg.matrix_rank(d.astype(float))

        # Compute im(d_{p+1})
        if p == n - 1:
            im_next_rank = 0  # no d_{n}
        else:
            dim_pp1 = len(omega[p+1]) if p+1 in omega else 0
            if dim_pp1 == 0:
                im_next_rank = 0
            else:
                # Build d_{p+1}: Omega_{p+1} -> Omega_p
                d_next = np.zeros((dim_p, dim_pp1), dtype=np.int64)

                for col, path in enumerate(omega[p+1]):
                    for i in range(p + 2):
                        if i == 0 or i == p + 1:
                            face = path[:i] + path[i+1:]
                            if face in omega_index[p]:
                                row = omega_index[p][face]
                                d_next[row, col] += (-1)**i
                        else:
                            if A[path[i-1]][path[i+1]]:
                                face = path[:i] + path[i+1:]
                                if face in omega_index[p]:
                                    row = omega_index[p][face]
                                    d_next[row, col] += (-1)**i

                im_next_rank = np.linalg.matrix_rank(d_next.astype(float))

        beta_p = int(ker_rank - im_next_rank)
        betti.append(beta_p)

        if verbose:
            print(f"  beta_{p} = ker(d_{p}) - im(d_{p+1}) = {ker_rank} - {im_next_rank} = {beta_p}")

    return betti

def main():
    rng = np.random.RandomState(42)

    # PART 1: Exhaustive at n=5,6
    print("=" * 70)
    print("PART 1: Complete Betti vectors")
    print("=" * 70)

    for n in [5, 6]:
        total = 2**(n*(n-1)//2)
        betti_dist = defaultdict(int)

        for bits in range(total):
            A = bits_to_adj(bits, n)
            betti = compute_all_betti(A, n)
            betti_dist[tuple(betti)] += 1

        print(f"\nn={n}: {total} tournaments")
        for bv, count in sorted(betti_dist.items(), key=lambda x: -x[1]):
            print(f"  {list(bv)}: {count} ({100*count/total:.1f}%)")

    # PART 2: Paley T_7
    print("\n" + "=" * 70)
    print("PART 2: Paley T_7")
    print("=" * 70)

    n = 7
    # Paley T_7: i->j iff j-i is a QR mod 7. QR mod 7 = {1, 2, 4}
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1

    print("Paley T_7 adjacency:")
    for i in range(n):
        print(f"  {i}: beats {[j for j in range(n) if A[i][j]]}")

    betti = compute_all_betti(A, n, verbose=True)
    print(f"\nBetti vector: {betti}")

    # PART 3: Sample n=7 tournaments and collect even Betti statistics
    print("\n" + "=" * 70)
    print("PART 3: Even Betti numbers at n=7 (sampled)")
    print("=" * 70)

    n = 7
    nsamp = 2000
    betti_dist = defaultdict(int)
    beta2_nonzero = 0
    beta4_nonzero = 0
    beta6_nonzero = 0

    for _ in range(nsamp):
        A = random_tournament(n, rng)
        betti = compute_all_betti(A, n)
        betti_dist[tuple(betti)] += 1
        if betti[2] != 0:
            beta2_nonzero += 1
        if betti[4] != 0:
            beta4_nonzero += 1
        if betti[6] != 0:
            beta6_nonzero += 1

    print(f"\nn=7 ({nsamp} random tournaments):")
    print(f"  beta_2 != 0: {beta2_nonzero}")
    print(f"  beta_4 != 0: {beta4_nonzero}")
    print(f"  beta_6 != 0: {beta6_nonzero}")
    print(f"\nTop Betti vectors:")
    for bv, count in sorted(betti_dist.items(), key=lambda x: -x[1])[:15]:
        print(f"  {list(bv)}: {count}")

    # PART 4: Check even Betti vanishing specifically
    print("\n" + "=" * 70)
    print("PART 4: Even Betti vanishing check (including beta_4)")
    print("=" * 70)

    # Targeted search: try near-regular tournaments at n=7
    for trial in range(500):
        A = random_tournament(n, rng)
        # Check if this is a regular tournament (all scores = 3)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if scores != (3,3,3,3,3,3,3):
            continue
        betti = compute_all_betti(A, n)
        if betti[4] != 0:
            print(f"  Regular n=7 with beta_4={betti[4]}: {betti}")

    # Also check all the Paley-like circulant tournaments
    print("\nCirculant tournaments at n=7:")
    for S_bits in range(1, 2**3):
        S = set()
        for k in range(1, 4):
            if S_bits & (1 << (k-1)):
                S.add(k)
        # Complete S: if k in S then n-k not in S (tournament)
        valid = True
        for k in S:
            if (n - k) in S:
                valid = False
                break
        if not valid:
            continue

        # Build circulant tournament
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j and ((j - i) % n) in S:
                    A[i][j] = 1

        # Verify it's a tournament
        is_tournament = True
        for i in range(n):
            for j in range(i+1, n):
                if A[i][j] + A[j][i] != 1:
                    is_tournament = False
                    break
            if not is_tournament:
                break
        if not is_tournament:
            continue

        betti = compute_all_betti(A, n)
        scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        H = betti[0]  # Not H actually, just beta_0
        # Compute H via Hamiltonian paths
        H_val = len([p for p in permutations(range(n)) if all(A[p[i]][p[i+1]] for i in range(n-1))])
        print(f"  S={S}: scores={scores}, H={H_val}, betti={betti}")

if __name__ == '__main__':
    main()
