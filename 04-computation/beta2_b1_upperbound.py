#!/usr/bin/env python3
"""beta2_b1_upperbound.py - Prove b1(T) <= 1 for all tournaments

b1(T) = dim(ker d1|Om_1) - rk(d2|Om_2)
      = dim(Z_1) - rk(d_2)

For tournaments: dim(Z_1) = dim(Om_1) - rk(d_1) = |edges| - |TT_pairs|
where TT_pairs are (a,c) with a->c having at least one mediator b (a->b->c).

Wait, dim(Om_1) = |edges| for tournaments since all paths are allowed.
No: Om_1 is the space of 1-chains in ker(d_0) intersected with allowed.
d_0 is trivial at degree 1, so Om_1 = all allowed 1-paths.
For a tournament on n vertices, |edges| = C(n,2).

Actually, for path homology:
- Om_0 = R^n (vertex space)
- Om_1 = {f : sum_{e ending at v} f(e) = sum_{e starting at v} f(e) for each v}
  No, that's not right either. Let me re-derive.

Omega_p = ker(d_1) cap ... cap ker(d_{p-1}) on allowed p-paths.
At p=1: only inner face maps are d_i for 0 < i < 1, which is EMPTY.
So Om_1 = R^{A_1} = R^{edges} for tournaments. dim(Om_1) = C(n,2).

d_1: Om_1 -> Om_0 = R^n. d_1(a,b) = (b) - (a) (boundary map).
rk(d_1) = n-1 (if connected, which tournaments always are).
dim(Z_1) = C(n,2) - (n-1) = (n-1)(n-2)/2.

Now rk(d_2): d_2: Om_2 -> Om_1.
d_2(a,b,c) = (b,c) - (a,c) + (a,b).
rk(d_2) = dim(Om_2) - dim(ker d_2|Om_2).

So b1 = dim(Z_1) - rk(d_2) = (n-1)(n-2)/2 - rk(d_2).

We need b1 <= 1, i.e., rk(d_2) >= (n-1)(n-2)/2 - 1.

This means d_2 has AT MOST 1-dimensional cokernel in Z_1.

KEY: dim(Om_2) = |A_2| - #{constrained (a,c) pairs}...
but we showed this formula is wrong. Let me just compute directly.

Plan: For various n, compute dim(Om_2), dim(Z_2), rk(d_2), and verify
rk(d_2) >= (n-1)(n-2)/2 - 1.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def compute_dims(A, n):
    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)

    dim_Z1_expected = (n-1)*(n-2)//2

    if not paths1:
        return {'dim_O1': 0, 'dim_O2': 0, 'rk_d1': 0, 'rk_d2': 0,
                'dim_Z1': 0, 'b1': 0, 'dim_Z1_expected': dim_Z1_expected}

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    if dim_O1 > 0:
        D1_om = D1 @ omega1
        sv1 = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv1))
    else:
        rk_d1 = 0

    dim_Z1 = dim_O1 - rk_d1

    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    else:
        dim_O2 = 0
        omega2 = None

    if dim_O2 > 0:
        D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                        [tuple(p) for p in paths1])
        D2_om = D2 @ omega2
        sv2 = np.linalg.svd(D2_om, compute_uv=False)
        rk_d2 = int(sum(s > 1e-8 for s in sv2))
    else:
        rk_d2 = 0

    b1 = dim_O1 - rk_d1 - rk_d2

    return {
        'dim_O1': dim_O1, 'dim_O2': dim_O2,
        'rk_d1': rk_d1, 'rk_d2': rk_d2,
        'dim_Z1': dim_Z1, 'b1': b1,
        'dim_Z1_expected': dim_Z1_expected
    }


# ============================================================
# Part 1: Exhaustive at small n
# ============================================================
print("=" * 70)
print("b1 AND rk(d2) FOR ALL TOURNAMENTS")
print("=" * 70)

for n in [3, 4, 5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    b1_dist = Counter()
    rk_d2_dist = Counter()
    dim_O2_dist = Counter()
    t0 = time.time()

    if n <= 6:
        for bits in range(total):
            A = [[0] * n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> idx) & 1:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
            d = compute_dims(A, n)
            b1_dist[d['b1']] += 1
            rk_d2_dist[d['rk_d2']] += 1
            dim_O2_dist[d['dim_O2']] += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): total={total}, dim_Z1=(n-1)(n-2)/2={(n-1)*(n-2)//2}")
    print(f"  b1 distribution: {dict(sorted(b1_dist.items()))}")
    print(f"  rk(d2) distribution: {dict(sorted(rk_d2_dist.items()))}")
    print(f"  dim(Om_2) distribution: {dict(sorted(dim_O2_dist.items()))}")
    print(f"  Need rk(d2) >= {(n-1)*(n-2)//2 - 1} for b1 <= 1")


# ============================================================
# Part 2: Sampled at larger n
# ============================================================
print(f"\n{'=' * 70}")
print("SAMPLED b1 AT LARGER n")
print("=" * 70)

for n in [7, 8, 9, 10, 12, 15, 20]:
    random.seed(42)
    trials = max(10, min(200, 2000 // n))
    b1_dist = Counter()
    rk_d2_vals = []

    for trial in range(trials):
        A = random_tournament(n)
        d = compute_dims(A, n)
        b1_dist[d['b1']] += 1
        rk_d2_vals.append(d['rk_d2'])

    dim_Z1 = (n-1)*(n-2)//2
    min_rk = min(rk_d2_vals)
    max_rk = max(rk_d2_vals)
    print(f"\nn={n}: {trials} trials, dim_Z1={dim_Z1}")
    print(f"  b1 distribution: {dict(sorted(b1_dist.items()))}")
    print(f"  rk(d2): min={min_rk}, max={max_rk}, need >= {dim_Z1 - 1}")
    print(f"  rk(d2) deficit from dim_Z1: {dim_Z1 - min_rk}")


# ============================================================
# Part 3: Can we prove rk(d2) >= dim_Z1 - 1 algebraically?
# ============================================================
print(f"\n{'=' * 70}")
print("ALGEBRAIC ANALYSIS OF rk(d2)")
print("=" * 70)

# d_2: Om_2 -> Z_1 (actually d_2: Om_2 -> Om_1, but image is in Z_1)
# We need corank of d_2 in Z_1 to be at most 1.
#
# Key observation: at n=5, b1 in {0,1}. When b1=1, the cokernel has dim 1.
# The unique (up to scalar) element of coker is the "cycle around the tournament".
#
# For SC tournaments with b1=1: the cokernel generator z is a non-trivial 1-cycle
# that is NOT a boundary. Can we characterize z?
#
# At n=3 (3-cycle): Z_1 has dim 1 (one cycle: (0,1)-(0,2)+(1,2) or cyclic).
# Om_2 is empty (no allowed 2-paths in C_3). So rk(d_2) = 0, dim(Z_1) = 1, b1 = 1.
#
# At n=4: Z_1 has dim 3.
# For SC tournament 0->1->2->3->0, 0->2, 1->3:
# Om_2 should have at least some TT triples.

# Let me look at a specific n=4 tournament with b1=1.
print("\nn=4 example with b1=1:")
# Score (1,1,2,2): the near-regular
A = [[0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 0, 1, 0]]  # 0->1->2->3, 3->0, 3->2
# Actually let me make the 4-cycle: 0->1->2->3->0
A = [[0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 0, 0, 0]]
# This has only edges 0->1, 1->2, 2->3, 3->0 and need to fill in 0<->2 and 1<->3.
# Tournament: 0->2 and 1->3 (makes it score (2,2,2,2) - wait, n=4 can't be regular)
# Score: deg(0)=2, deg(1)=2, deg(2)=2, deg(3)=2. But sum=8, need C(4,2)=6. Error.
# n=4: 6 edges. Each vertex has out-degree summing to 6. Score sum = 6.
# Regular: each has out-deg 1.5. NOT integer. So no regular at n=4.

# Score (1,1,2,2) is the most balanced.
# 0->1, 1->2, 2->0, 2->3, 3->0, 0->3? Let me enumerate.
n = 4
b1_1_example = None
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    d = compute_dims(A, n)
    if d['b1'] == 1:
        b1_1_example = (bits, A, d)
        break

if b1_1_example:
    bits, A, d = b1_1_example
    print(f"  bits={bits}")
    print(f"  Adjacency: {A}")
    print(f"  dims: {d}")
    scores = [sum(A[i]) for i in range(n)]
    print(f"  Scores: {scores}")

    # Compute the cokernel generator
    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                    [tuple(p) for p in paths1])

    D1_om = D1 @ omega1
    D2_om = D2 @ omega2

    # Z_1 = ker(D1_om)
    U, S, Vt = np.linalg.svd(D1_om)
    nullity = sum(s < 1e-8 for s in S)
    # Null space of D1_om^T in Omega_1 coordinates
    # Actually: Z_1 in Om_1 coordinates = null space of D1_om
    # D1_om has shape (n, dim_O1). ker has dimension dim_O1 - rk = dim_Z1

    # Find ker(D1_om)
    if nullity > 0:
        null_vectors = Vt[-nullity:]  # rows of Vt
        print(f"\n  Z_1 basis (in Om_1 coords, dim={nullity}):")
        for i, v in enumerate(null_vectors):
            # Convert to path coordinates
            path_coords = omega1 @ v
            nonzero = [(paths1[j], path_coords[j]) for j in range(len(paths1)) if abs(path_coords[j]) > 1e-8]
            print(f"    z_{i}: {nonzero[:6]}...")

    # Image of d_2 in Z_1
    # Project D2_om into Om_1 coords, then into Z_1 coords
    # D2_om is in path1 coords. Project to Om_1 coords via omega1^+
    D2_in_om1 = np.linalg.lstsq(omega1, D2_om, rcond=None)[0]
    # Now project to Z_1 coords
    # Z_1 basis in Om_1 coords = null_vectors (each row is a basis vector)
    Z1_basis = null_vectors  # shape (nullity, dim_O1)
    D2_in_Z1 = Z1_basis @ D2_in_om1  # shape (nullity, dim_O2)
    print(f"\n  Image of d_2 in Z_1 (rank should be {d['rk_d2']}):")
    sv = np.linalg.svd(D2_in_Z1, compute_uv=False)
    rk = int(sum(s > 1e-8 for s in sv))
    print(f"    rank = {rk}, dim(Z_1) = {nullity}, corank = {nullity - rk}")

    # The cokernel generator
    if nullity - rk > 0:
        U2, S2, Vt2 = np.linalg.svd(D2_in_Z1.T)
        coker_in_Z1 = Vt2[-(nullity - rk):]  # shape (corank, nullity)
        coker_in_Om1 = coker_in_Z1 @ Z1_basis  # shape (corank, dim_O1)
        coker_in_paths = coker_in_Om1 @ omega1.T  # Wrong shape, let me fix
        # omega1 has shape (|paths1|, dim_O1). coker_in_Om1 has shape (corank, dim_O1).
        # To get path coordinates: multiply coker_in_Om1 @ omega1^T? No.
        # coker_in_Om1[i] is a vector in R^{dim_O1}. The path-space representation is omega1 @ coker_in_Om1[i].
        for i in range(coker_in_Z1.shape[0]):
            path_coords = omega1 @ coker_in_Om1[i]
            nonzero = [(paths1[j], round(path_coords[j], 4)) for j in range(len(paths1)) if abs(path_coords[j]) > 1e-8]
            print(f"\n  Cokernel generator {i} (1-cycle not in im(d_2)):")
            print(f"    {nonzero}")


# ============================================================
# Part 4: b1=1 always corresponds to a single cycle generator?
# ============================================================
print(f"\n{'=' * 70}")
print("b1=1 COKERNEL ANALYSIS (n=5)")
print("=" * 70)

n = 5
count_b1_1 = 0
coker_structures = Counter()

for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    d = compute_dims(A, n)
    if d['b1'] != 1:
        continue
    count_b1_1 += 1

    # Find the cokernel generator
    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1]

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    U, S, Vt = np.linalg.svd(D1_om)
    nullity = sum(s < 1e-8 for s in S)
    null_vectors = Vt[-nullity:]

    if not paths2:
        coker_structures['no_Om2'] += 1
        continue

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        coker_structures['dim_O2=0'] += 1
        continue

    D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                    [tuple(p) for p in paths1])
    D2_om = D2 @ omega2
    D2_in_om1 = np.linalg.lstsq(omega1, D2_om, rcond=None)[0]
    Z1_basis = null_vectors
    D2_in_Z1 = Z1_basis @ D2_in_om1
    sv = np.linalg.svd(D2_in_Z1, compute_uv=False)
    rk = int(sum(s > 1e-8 for s in sv))

    coker_structures[f'Z1={nullity},rk_d2={rk},corank={nullity-rk}'] += 1

    # Check: cokernel generator supported on how many edges?
    if nullity - rk == 1:
        U2, S2, Vt2 = np.linalg.svd(D2_in_Z1.T)
        coker_in_Z1 = Vt2[-1:]
        coker_in_Om1 = coker_in_Z1 @ Z1_basis
        path_coords = omega1 @ coker_in_Om1[0]
        num_edges = sum(abs(path_coords[j]) > 1e-8 for j in range(len(paths1)))
        coker_structures[f'edges_in_coker={num_edges}'] += 1

print(f"\nn=5: {count_b1_1} tournaments with b1=1")
print(f"  Cokernel structures:")
for k, v in sorted(coker_structures.items()):
    print(f"    {k}: {v}")


# ============================================================
# Part 5: What is dim(Om_2) as a function of tournament structure?
# ============================================================
print(f"\n{'=' * 70}")
print("dim(Om_2) vs TOURNAMENT INVARIANTS")
print("=" * 70)

n = 5
inv_data = Counter()
for bits in range(2 ** (n*(n-1)//2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    d = compute_dims(A, n)
    inv_data[(scores, d['dim_O2'], d['rk_d2'], d['b1'])] += 1

print(f"\nn=5: (score, dim_O2, rk_d2, b1) => count")
for (s, o, r, b), cnt in sorted(inv_data.items()):
    print(f"  score={s}, dim_O2={o}, rk_d2={r}, b1={b}: {cnt}")


print("\n\nDone.")
