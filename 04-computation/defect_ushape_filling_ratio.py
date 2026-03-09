"""
defect_ushape_filling_ratio.py — Two remaining S44 recommendations

1. FILLING RATIO FORMULA: f_p = dim(Omega_p) / C(n, p+1).
   f_p > 1 for p >= 3 due to cyclic content. Can we find an exact formula?
   dim(Omega_p) = |A_p| - rank(constraint_matrix_p)
   |A_p| counts allowed p-paths. Need formula for both.

2. DEFECT RATE U-SHAPE: defect_p = dim(Omega_p) - rank(d_p) - rank(d_{p+1}).
   = dim(H_p) + (dim(ker d_p) - rank(d_{p+1})) = beta_p + "surplus"
   Actually defect_p = beta_p exactly.
   The U-shape is: P(beta_1>0) decreasing, P(beta_3>0) increasing.
   Why does beta_3 rate INCREASE with n while beta_1 rate DECREASES?

Author: kind-pasteur-S46 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

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

def compute_omega_dim(ap, p):
    if not ap.get(p, []):
        return 0
    if p == 0:
        return len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    sv = np.linalg.svd(P, compute_uv=False)
    rank = int(sum(s > 1e-10 for s in sv))
    return len(ap[p]) - rank


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] + A[j][i]) * (A[j][k] + A[k][j]) * (A[i][k] + A[k][i]) == 1:
                    # Check if it's a directed 3-cycle
                    if (A[i][j] == 1 and A[j][k] == 1 and A[k][i] == 1) or \
                       (A[i][k] == 1 and A[k][j] == 1 and A[j][i] == 1):
                        c3 += 1
    return c3


def main():
    print("=" * 70)
    print("DEFECT U-SHAPE AND FILLING RATIO ANALYSIS")
    print("=" * 70)

    # Part 1: Filling ratio f_p = dim(Omega_p) / C(n, p+1) at various n
    print("\n--- Part 1: Filling ratio by n ---")

    for n in [5, 6, 7]:
        if n <= 6:
            N_total = 2**(n*(n-1)//2)
            sample = False
        else:
            N_total = 500
            sample = True

        filling_sums = defaultdict(float)
        filling_counts = defaultdict(int)
        filling_min = defaultdict(lambda: float('inf'))
        filling_max = defaultdict(lambda: float('-inf'))

        rng = np.random.RandomState(42)

        count = 0
        for trial in range(N_total):
            if sample:
                A = random_tournament(n, rng)
            else:
                A = bits_to_adj(trial, n)

            ap = {}
            for p in range(min(n, 7)):
                ap[p] = enumerate_allowed_paths(A, n, p)

            for p in range(1, min(n-1, 6)):
                dim_o = compute_omega_dim(ap, p)
                simplex_dim = comb(n, p+1)
                if simplex_dim > 0:
                    f = dim_o / simplex_dim
                    filling_sums[p] += f
                    filling_counts[p] += 1
                    filling_min[p] = min(filling_min[p], f)
                    filling_max[p] = max(filling_max[p], f)

            count += 1
            if count % 5000 == 0 and not sample:
                print(f"  n={n}: {count}/{N_total} done", flush=True)
            if count % 100 == 0 and sample:
                print(f"  n={n}: {count}/{N_total} done", flush=True)

        print(f"\n  n={n} ({'exhaustive' if not sample else f'sampled {N_total}'}):")
        for p in sorted(filling_sums.keys()):
            avg = filling_sums[p] / filling_counts[p]
            print(f"    f_{p} = dim(Omega_{p})/C({n},{p+1}): avg={avg:.4f}, "
                  f"min={filling_min[p]:.4f}, max={filling_max[p]:.4f}")

    # Part 2: Filling ratio as function of c3
    print("\n--- Part 2: Filling ratio f_2 vs c3 at n=6 ---")
    n = 6
    c3_f2 = defaultdict(list)

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        c3 = count_3cycles(A, n)

        ap = {}
        for p in range(4):
            ap[p] = enumerate_allowed_paths(A, n, p)
        dim_o2 = compute_omega_dim(ap, 2)
        f2 = dim_o2 / comb(n, 3)
        c3_f2[c3].append(f2)

    print(f"  c3 -> avg(f_2):")
    for c3 in sorted(c3_f2.keys()):
        vals = c3_f2[c3]
        avg = sum(vals) / len(vals)
        mi = min(vals)
        ma = max(vals)
        print(f"    c3={c3}: avg={avg:.4f}, min={mi:.4f}, max={ma:.4f} ({len(vals)} tours)")

    # Part 3: dim(Omega_2) formula check: C(n,3) + 2*c3 - e_cyc - rank(constraints)
    # We know dim(A_2) = C(n,3) + 2*c3, and constraints come from non-allowed faces
    print("\n--- Part 3: Omega_2 = |A_2| - #constraints at n=5,6 ---")
    for n in [5, 6]:
        data = []
        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            c3 = count_3cycles(A, n)
            ap = {p: enumerate_allowed_paths(A, n, p) for p in range(3)}
            a2 = len(ap[2])
            dim_o2 = compute_omega_dim(ap, 2)
            constraints = a2 - dim_o2
            data.append((c3, a2, dim_o2, constraints))

        # Group by c3
        by_c3 = defaultdict(list)
        for c3, a2, o2, cons in data:
            by_c3[c3].append((a2, o2, cons))

        print(f"\n  n={n}:")
        for c3 in sorted(by_c3.keys()):
            vals = by_c3[c3]
            a2_vals = set(v[0] for v in vals)
            o2_vals = set(v[1] for v in vals)
            cons_vals = set(v[2] for v in vals)
            # Check |A_2| = C(n,3) + 2*c3
            expected_a2 = comb(n, 3) + 2 * c3
            a2_check = all(v[0] == expected_a2 for v in vals)
            print(f"    c3={c3}: |A_2|={a2_vals} (expect {expected_a2}, {'OK' if a2_check else 'FAIL'}), "
                  f"dim(Omega_2)={o2_vals}, #constraints={cons_vals}")

    # Part 4: Defect rate by n — the U-shape
    print("\n--- Part 4: Betti rate by n (the crossover / U-shape) ---")
    for n in [5, 6, 7, 8]:
        if n <= 6:
            rng = None
            N = 2**(n*(n-1)//2)
            sample = False
        else:
            rng = np.random.RandomState(42)
            N = 500 if n == 7 else 200
            sample = True

        beta1_pos = 0
        beta3_pos = 0
        total = 0

        for trial in range(N):
            if sample:
                A = random_tournament(n, rng)
            else:
                A = bits_to_adj(trial, n)

            ap = {}
            for p in range(min(6, n)):
                ap[p] = enumerate_allowed_paths(A, n, p)

            # beta_1: need rank(d_1) and rank(d_2)
            # rank(d_1) = n - 1 always
            # rank(d_2): compute from TT triples
            tts = []
            for i in range(n):
                for j in range(n):
                    if i == j or A[i][j] == 0: continue
                    for k in range(n):
                        if k == i or k == j: continue
                        if A[i][k] == 1 and A[j][k] == 1:
                            tts.append((i, j, k))

            edges = [(i, j) for i in range(n) for j in range(n) if A[i][j] == 1]
            edge_idx = {e: i for i, e in enumerate(edges)}
            bd = np.zeros((len(edges), len(tts)))
            for j, (a, b, c) in enumerate(tts):
                if (b, c) in edge_idx: bd[edge_idx[(b, c)], j] += 1
                if (a, c) in edge_idx: bd[edge_idx[(a, c)], j] -= 1
                if (a, b) in edge_idx: bd[edge_idx[(a, b)], j] += 1
            sv = np.linalg.svd(bd, compute_uv=False)
            rank_d2 = int(sum(s > 1e-8 for s in sv))
            dim_Z1 = comb(n, 2) - (n - 1)
            b1 = max(0, dim_Z1 - rank_d2)

            if b1 > 0:
                beta1_pos += 1

            # beta_3: simplified check
            from math import comb as C
            ob3_dim = compute_omega_dim(ap, 3)
            if ob3_dim > 0:
                ob3_basis, _ = compute_omega_basis_full(ap, 3)
                if ob3_basis is not None:
                    # Build d_3
                    bd3 = np.zeros((len(ap.get(2, [])), len(ap[3])))
                    idx2 = {path: i for i, path in enumerate(ap.get(2, []))}
                    for j, path in enumerate(ap[3]):
                        for sign, face in boundary_coeffs(path):
                            if face in idx2:
                                bd3[idx2[face], j] += sign
                    d3 = bd3 @ ob3_basis
                    sv3 = np.linalg.svd(d3, compute_uv=False)
                    rk3 = int(sum(s > 1e-8 for s in sv3))
                    ker3 = ob3_dim - rk3
                    if ker3 > 0:
                        beta3_pos += 1

            total += 1
            if total % 5000 == 0 and not sample:
                print(f"  n={n}: {total}/{N}", flush=True)

        print(f"  n={n}: P(beta_1>0) = {100*beta1_pos/total:.1f}%, P(beta_3>0) ~ {100*beta3_pos/total:.1f}%")

    print("\nDONE.")


def compute_omega_basis_full(ap, p):
    """Return full basis and dimension."""
    if not ap.get(p, []):
        return None, 0
    if p == 0:
        return np.eye(len(ap[p])), len(ap[p])

    apm1_set = set(ap.get(p-1, []))
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in apm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(len(ap[p])), len(ap[p])

    P = np.zeros((na_count, len(ap[p])))
    for j, path in enumerate(ap[p]):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    if ns.shape[1] > 0:
        return ns, ns.shape[1]
    else:
        return None, 0


if __name__ == '__main__':
    main()
