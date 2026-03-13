"""
witness_rank3_decomp.py -- kind-pasteur-2026-03-13-S61

The diff matrix D = W_B - W_A under a Vitali atom ALWAYS has rank 3.
Hypothesis: the 3 rank-1 components correspond to the 3 outside vertices.

Each outside vertex k contributes a rank-1 change to the witness matrix:
  it changes which internal pairs it witnesses.

Also: connect this to the dc7 formula. Since dc7 is determined by
delta_sigma, and sigma changes, the rank-3 structure of D should
encode exactly the dc7 information.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter

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

def witness_matrix(A, n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    W = np.zeros((n, len(pairs)), dtype=int)
    for pi, (u, v) in enumerate(pairs):
        for k in range(n):
            if k == u or k == v:
                continue
            if (A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v]):
                W[k][pi] = 1
    return W, pairs

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
pair_idx = {p: i for i, p in enumerate(pairs)}

print("=" * 60)
print("RANK-3 DECOMPOSITION OF WITNESS DIFF")
print("=" * 60)

np.random.seed(42)
examples = []

for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        W_A, _ = witness_matrix(A, n)
        W_B, _ = witness_matrix(B, n)
        diff = W_B - W_A

        if np.count_nonzero(diff) == 0:
            continue

        c7_A = count_directed_k_cycles(A, n, 7)
        c7_B = count_directed_k_cycles(B, n, 7)
        dc7 = c7_B - c7_A

        S = set(subset)
        outside = [k for k in range(n) if k not in S]

        examples.append({
            'bits': bits,
            'subset': subset,
            'diff': diff.copy(),
            'dc7': dc7,
            'S': S,
            'outside': outside,
        })
        break

    if len(examples) >= 80:
        break

print(f"Collected {len(examples)} non-trivial Vitali examples")

# Decompose diff by outside vertex contribution
print(f"\n--- Per-outside-vertex decomposition ---")

for ex_idx in range(min(5, len(examples))):
    ex = examples[ex_idx]
    diff = ex['diff']
    S = ex['S']
    outside = ex['outside']

    print(f"\nExample {ex_idx}: S={ex['subset']}, outside={outside}, dc7={ex['dc7']}")

    # For each outside vertex k: extract the rows of diff corresponding to k
    # and also the changes in the k-row (what pairs does k's witnessing change for?)
    for k in outside:
        k_changes = []
        for pi in range(len(pairs)):
            if diff[k][pi] != 0:
                u, v = pairs[pi]
                k_changes.append((u, v, int(diff[k][pi])))
        if k_changes:
            print(f"  Outside k={k}: changes {k_changes}")

    # For each S vertex: extract its row changes
    for k in sorted(S):
        k_changes = []
        for pi in range(len(pairs)):
            if diff[k][pi] != 0:
                u, v = pairs[pi]
                k_changes.append((u, v, int(diff[k][pi])))
        if k_changes:
            print(f"  S-vertex k={k}: changes {k_changes}")

    # Now decompose diff into 3 components, one per outside vertex
    # Component for outside vertex k: restricted to columns involving at least one S vertex
    # Wait, the changes are: kX witnesses SS pairs, and kS witnesses SX pairs.
    # Let me think about the 3 "swap layers" differently.

    # The (1,1,2,2) sub-tournament on S has specific structure.
    # After reversal, all arcs within S flip.
    # This changes the 3-cycle structure:
    #   - For each outside vertex k, some (k, u, v) triples where u,v in S
    #     change their 3-cycle status (witness status).
    #   - For each pair (s, x) where s in S, x outside S,
    #     the internal vertices of S change their witnessing of (s, x).

    # The 3 outside vertices each contribute a "layer" of rank <= 1.
    # Let me extract these layers.

    layers = {}
    for k in outside:
        layer = np.zeros((n, len(pairs)), dtype=int)
        # Row k: direct changes to k's witnessing
        for pi in range(len(pairs)):
            if diff[k][pi] != 0:
                layer[k][pi] = diff[k][pi]

        # The corresponding column changes must balance.
        # For each column (pair) that k affects, the change in column sum must be 0.
        # So there must be compensating changes in S-vertex rows for the same pairs.
        # But we can't simply assign: the S-vertex changes correspond to MULTIPLE
        # outside vertices simultaneously.

        layers[k] = layer

    # Actually, the S-vertex changes for SX pairs involve one S vertex and one X vertex.
    # When the SX pair is (s, x), the witnesses are vertices OTHER than s and x.
    # So the witnesses can be: other S vertices, or other X vertices.
    # For diff[k'][(s,x)] where k' is an S-vertex, the change is because
    # k' witnessed (s,x) before but not after (or vice versa).

    # Let me try a different decomposition: by the PAIR type.
    # SS pairs: only witnessed by vertices not in the pair. At n=7, each SS pair
    #   has 5 potential witnesses: 2 in S, 3 outside S.
    #   Changes: only in k_outside rows (4 per outside vertex... let me check).

    print(f"\n  Column-wise analysis (which pairs change):")
    changed_cols = set()
    for pi in range(len(pairs)):
        col = diff[:, pi]
        if np.any(col != 0):
            u, v = pairs[pi]
            u_in = u in S
            v_in = v in S
            ptype = "SS" if (u_in and v_in) else ("SX" if (u_in or v_in) else "XX")
            plus = [k for k in range(n) if col[k] > 0]
            minus = [k for k in range(n) if col[k] < 0]
            changed_cols.add(pi)
            if ex_idx == 0:
                print(f"    pair ({u},{v}) [{ptype}]: +1 at k={plus}, -1 at k={minus}")

# Key question: does dc7 correlate with specific features of the diff?
print(f"\n{'='*60}")
print("dc7 CORRELATION WITH DIFF STRUCTURE")
print(f"{'='*60}")

dc7_vals = [ex['dc7'] for ex in examples]
dc7_dist = Counter(dc7_vals)
print(f"dc7 distribution: {dict(sorted(dc7_dist.items()))}")

# For each example, compute features of the diff matrix
for ex in examples:
    diff = ex['diff']
    S = ex['S']
    outside = ex['outside']

    # Feature: sum of diff entries restricted to outside rows and SS columns
    ss_pairs = [pi for pi in range(len(pairs)) if pairs[pi][0] in S and pairs[pi][1] in S]
    sx_pairs = [pi for pi in range(len(pairs))
                if (pairs[pi][0] in S) != (pairs[pi][1] in S)]

    # "Polarization" of outside vertices' changes
    polarizations = []
    for k in outside:
        pol = sum(diff[k][pi] for pi in ss_pairs)
        polarizations.append(pol)
    ex['polarizations'] = tuple(polarizations)

# Does polarization predict dc7?
print(f"\nPolarization pattern vs dc7:")
pol_dc7 = defaultdict(list)
for ex in examples:
    pol_dc7[ex['polarizations']].append(ex['dc7'])

for pol, dc7s in sorted(pol_dc7.items()):
    mean_dc7 = np.mean(dc7s)
    print(f"  pol={pol}: dc7 mean={mean_dc7:.2f}, values={Counter(dc7s)}")

# THE SIGMA CHANGE: compute delta_sigma for each example
print(f"\n{'='*60}")
print("SIGMA CHANGE VS dc7")
print(f"{'='*60}")

for ex in examples[:20]:
    A = bits_to_adj(ex['bits'], n)
    B = reverse_subtournament(A, n, list(ex['subset']))
    A2 = A @ A
    B2 = B @ B

    # sigma(u,v) = n-2 - A2[u][v] - A2[v][u]
    delta_sigma = []
    for u in range(n):
        for v in range(u+1, n):
            sig_A = n - 2 - int(A2[u][v]) - int(A2[v][u])
            sig_B = n - 2 - int(B2[u][v]) - int(B2[v][u])
            if sig_A != sig_B:
                delta_sigma.append((u, v, sig_B - sig_A))

    ex['delta_sigma'] = delta_sigma
    ex['sum_delta_sigma'] = sum(ds for _, _, ds in delta_sigma)
    ex['sum_abs_delta_sigma'] = sum(abs(ds) for _, _, ds in delta_sigma)

print(f"dc7 vs sigma change:")
print(f"{'dc7':>5} {'sum_dsig':>10} {'sum_|dsig|':>12} {'n_changed':>10}")
for ex in examples[:20]:
    print(f"{ex['dc7']:>5} {ex['sum_delta_sigma']:>10} {ex['sum_abs_delta_sigma']:>12} {len(ex['delta_sigma']):>10}")

# Is dc7 a LINEAR function of the sigma changes?
# Collect: delta_sigma for each pair, and dc7
# At n=7, there are C(7,2)=21 pairs.
# delta_sigma is sparse (few pairs change).

# Build design matrix
X_rows = []
y_rows = []
for ex in examples[:60]:
    A = bits_to_adj(ex['bits'], n)
    B = reverse_subtournament(A, n, list(ex['subset']))
    A2 = A @ A
    B2 = B @ B

    dsig = np.zeros(len(pairs))
    for u in range(n):
        for v in range(u+1, n):
            pi = pair_idx[(u, v)]
            sig_A = n - 2 - int(A2[u][v]) - int(A2[v][u])
            sig_B = n - 2 - int(B2[u][v]) - int(B2[v][u])
            dsig[pi] = sig_B - sig_A

    X_rows.append(dsig)
    y_rows.append(ex['dc7'])

X = np.array(X_rows, dtype=float)
y = np.array(y_rows, dtype=float)

from numpy.linalg import lstsq, matrix_rank

print(f"\nRegression: dc7 = linear function of delta_sigma")
print(f"  X shape: {X.shape}, rank: {matrix_rank(X)}")

coeffs, _, _, _ = lstsq(X, y, rcond=None)
err = np.max(np.abs(X @ coeffs - y))
print(f"  Max error: {err:.6f}")
print(f"  Perfect fit? {err < 0.001}")

if err < 0.001:
    print(f"  dc7 IS a linear function of delta_sigma!")
    # Find which pairs matter
    nonzero = [(i, float(coeffs[i])) for i in range(len(coeffs)) if abs(coeffs[i]) > 0.001]
    print(f"  Nonzero coefficients: {len(nonzero)}")
    for pi, c in nonzero[:10]:
        print(f"    pair {pairs[pi]}: coeff = {c:.4f}")
else:
    # Try sum of delta_sigma
    sum_dsig = np.array([sum(x) for x in X_rows]).reshape(-1, 1)
    c_sum, _, _, _ = lstsq(sum_dsig, y, rcond=None)
    err_sum = np.max(np.abs(sum_dsig * c_sum[0] - y))
    print(f"\n  dc7 = c * sum(delta_sigma)?")
    print(f"  c = {c_sum[0]:.6f}, max error = {err_sum:.6f}")

# What about dc7 = linear function of diff matrix entries?
X_diff = np.array([ex['diff'].flatten() for ex in examples[:60]], dtype=float)
coeffs_diff, _, _, _ = lstsq(X_diff, y, rcond=None)
err_diff = np.max(np.abs(X_diff @ coeffs_diff - y))
print(f"\n  dc7 = linear function of diff(W)?")
print(f"  X shape: {X_diff.shape}, rank: {matrix_rank(X_diff)}")
print(f"  Max error: {err_diff:.6f}")

print("\nDone.")
