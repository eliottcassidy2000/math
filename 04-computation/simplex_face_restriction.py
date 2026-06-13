"""
Simplex face restriction: How do Betti numbers and Omega dims
change under vertex deletion (= restricting to a face of the simplex)?

Key framing: Tournament T on n vertices = oriented (n-1)-simplex Delta_{n-1}.
Deleting vertex v = restricting to the face opposite v = (n-2)-simplex.

Questions:
1. Long exact sequence in path homology under face restriction?
2. How does beta_p(T) relate to beta_p(T-v)?
3. The "cone over T-v" construction: does it kill all homology?
4. Spectral sequence: Delta_{n-1} = colimit of faces?

For the GLMY path complex:
- Deleting v removes all paths through v
- But also changes which paths are "allowed" (boundary faces may change)
- The relationship is NON-TRIVIAL

Data collection: For each tournament T, compute betti(T) and betti(T-v)
for all v. Track: sum_v beta_p(T-v), max_v beta_p(T-v), etc.
"""
import numpy as np
from itertools import combinations
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

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    omega_dims = []
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        omega_dims.append(omega[p].shape[1] if omega[p].ndim == 2 else 0)
    betti = []
    for p in range(max_dim + 1):
        dim_om = omega_dims[p]
        if dim_om == 0:
            betti.append(0)
            continue
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else:
            rk = 0
        ker = dim_om - rk
        dim_om1 = omega_dims[p+1]
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else:
            im = 0
        betti.append(ker - im)
    return betti, omega_dims[:max_dim+1]

def main():
    print("=" * 70)
    print("SIMPLEX FACE RESTRICTION: Betti under vertex deletion")
    print("Tournament T on Delta_{n-1}; delete v = restrict to Delta_{n-2}")
    print("=" * 70)

    # Exhaustive at n=5,6
    for n in [5, 6]:
        d = n - 1
        total = 2 ** (n * (n-1) // 2)

        print(f"\n{'='*60}")
        print(f"d = {d} (n = {n})")
        print(f"{'='*60}")

        # Track: for each (betti(T), tuple of betti(T-v) for all v)
        transition_data = defaultdict(int)
        sum_deletion_betti = defaultdict(list)

        for bits in range(total):
            A = bits_to_adj(bits, n)
            betti_T, _ = path_betti_numbers(A, n)

            # Compute betti for all vertex deletions
            del_bettis = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                betti_Tv, _ = path_betti_numbers(B, n-1)
                del_bettis.append(tuple(int(b) for b in betti_Tv))

            # Sum of beta_p across deletions
            for p in range(len(betti_T)):
                s = sum(del_bettis[v][p] if p < len(del_bettis[v]) else 0 for v in range(n))
                sum_deletion_betti[p].append(s)

            betti_key = tuple(int(b) for b in betti_T)
            # Canonical deletion signature: sorted tuple of deletion Betti vectors
            del_sig = tuple(sorted(del_bettis))
            transition_data[(betti_key, del_sig)] += 1

        # Report
        print(f"\n  Dimension-change transitions (top 10):")
        for (bt, ds), cnt in sorted(transition_data.items(), key=lambda x: -x[1])[:10]:
            pct = 100 * cnt / total
            # Count distinct deletion types
            del_types = defaultdict(int)
            for db in ds:
                del_types[db] += 1
            del_summary = ", ".join(f"{db}x{c}" for db, c in sorted(del_types.items(), key=lambda x: -x[1]))
            print(f"    betti(T)={list(bt)} -> deletions: [{del_summary}] : {cnt} ({pct:.1f}%)")

        # Sum statistics
        print(f"\n  Sum_v beta_p(T-v) statistics:")
        for p in sorted(sum_deletion_betti.keys()):
            vals = sum_deletion_betti[p]
            if max(vals) > 0 or p <= 2:
                print(f"    p={p}: range [{min(vals)}, {max(vals)}], "
                      f"mean={np.mean(vals):.2f}, "
                      f"P(sum>0)={100*sum(1 for v in vals if v > 0)/len(vals):.1f}%")

    # Part 2: Dimensional reduction chains
    # For a specific tournament, trace Betti through successive vertex deletions
    print("\n" + "=" * 70)
    print("DIMENSIONAL REDUCTION CHAINS: Betti under successive deletions")
    print("=" * 70)

    # Paley T_7
    n = 7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1

    print(f"\nPaley T_7 reduction chain:")
    current = A.copy()
    current_n = n
    while current_n >= 3:
        betti, odims = path_betti_numbers(current, current_n)
        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        print(f"  d={current_n-1} (n={current_n}): betti={[int(b) for b in betti]}, "
              f"omega={odims}, chi={chi}")
        # Delete vertex 0
        current = delete_vertex(current, current_n, 0)
        current_n -= 1

    # Transitive tournament
    print(f"\nTransitive T_7 reduction chain:")
    A_trans = np.zeros((7, 7), dtype=int)
    for i in range(7):
        for j in range(i+1, 7):
            A_trans[i][j] = 1
    current = A_trans.copy()
    current_n = 7
    while current_n >= 3:
        betti, odims = path_betti_numbers(current, current_n)
        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        print(f"  d={current_n-1} (n={current_n}): betti={[int(b) for b in betti]}, "
              f"omega={odims}, chi={chi}")
        current = delete_vertex(current, current_n, 0)
        current_n -= 1

    # Random tournament
    rng = np.random.RandomState(42)
    A_rand = np.zeros((7, 7), dtype=int)
    for i in range(7):
        for j in range(i+1, 7):
            if rng.random() < 0.5:
                A_rand[i][j] = 1
            else:
                A_rand[j][i] = 1
    print(f"\nRandom T_7 (seed 42) reduction chain:")
    current = A_rand.copy()
    current_n = 7
    while current_n >= 3:
        betti, odims = path_betti_numbers(current, current_n)
        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        print(f"  d={current_n-1} (n={current_n}): betti={[int(b) for b in betti]}, "
              f"omega={odims}, chi={chi}")
        current = delete_vertex(current, current_n, 0)
        current_n -= 1

    # Part 3: Omega_p / C(n, p+1) filling ratio pattern
    print("\n" + "=" * 70)
    print("FILLING RATIO: dim(Omega_p) / C(n, p+1) across dimensions")
    print("For transitive tournaments: this should be 1.0")
    print("For Paley: should reflect vertex-transitive structure")
    print("=" * 70)

    from math import comb

    for label, make_tour in [
        ("Transitive", lambda n: np.array([[1 if j > i else 0 for j in range(n)] for i in range(n)])),
        ("Paley (nearest prime)", None),
    ]:
        if label == "Paley (nearest prime)":
            paley_primes = [3, 7, 11]
            for p in paley_primes:
                if p > 8: continue  # too expensive
                n = p
                A = np.zeros((n, n), dtype=int)
                qr_set = {k*k % n for k in range(1, n)}
                for i in range(n):
                    for j in range(n):
                        if i != j and ((j - i) % n) in qr_set:
                            A[i][j] = 1
                _, odims = path_betti_numbers(A, n)
                print(f"\n  Paley T_{p}:")
                for q in range(len(odims)):
                    c = comb(n, q+1)
                    ratio = odims[q] / c if c > 0 else 0
                    print(f"    p={q}: Omega={odims[q]}, C({n},{q+1})={c}, ratio={ratio:.4f}")
        else:
            for n in range(3, 8):
                A = make_tour(n)
                _, odims = path_betti_numbers(A, n)
                print(f"\n  Transitive T_{n}:")
                for q in range(len(odims)):
                    c = comb(n, q+1)
                    ratio = odims[q] / c if c > 0 else 0
                    print(f"    p={q}: Omega={odims[q]}, C({n},{q+1})={c}, ratio={ratio:.4f}")

if __name__ == '__main__':
    main()
