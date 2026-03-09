"""
seesaw_rank_constraints.py — Test if rank constraints force adjacent-odd seesaw.

At each level p, dim(Om_p) = rank(d_p) + ker(d_p) and ker(d_p) = im(d_{p+1}) + beta_p.
So: dim(Om_p) = rank(d_p) + im(d_{p+1}) + beta_p.

The key identities:
  rank(d_p) = im(d_p) (rank of the boundary map = dimension of image in Om_{p-1})
  im(d_p) = rank(d_p) <= ker(d_{p-1}) = dim(Om_{p-1}) - rank(d_{p-1})

Chain of inequalities:
  rank(d_p) + im(d_{p+1}) + beta_p = dim(Om_p)
  rank(d_p) <= dim(Om_{p-1}) - rank(d_{p-1})

For the seesaw at level 2k:
  rank(d_{2k}) + im(d_{2k+1}) + beta_{2k} = dim(Om_{2k})
  rank(d_{2k}) = im(d_{2k}) contributes to ker(d_{2k-1}) which affects beta_{2k-1}
  im(d_{2k+1}) contributes to ker(d_{2k+1}) which affects beta_{2k+1}

If beta_{2k-1} > 0 (deficiency at level 2k-1):
  rank(d_{2k}) = im(d_{2k}) < ker(d_{2k-1})
  => im(d_{2k}) = dim(Om_{2k-1}) - rank(d_{2k-1}) - beta_{2k-1}
  => im(d_{2k}) is SMALLER than ker(d_{2k-1}) by beta_{2k-1}

Does this force im(d_{2k+1}) to be large enough to prevent beta_{2k+1}>0?

Let's track all the rank data and see.

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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

def compute_chain_ranks(A, n, max_p=None):
    """Compute omega dimensions and boundary map ranks."""
    if max_p is None:
        max_p = n - 1

    ap = {}
    for p in range(max_p + 2):
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    omega_bases = {}
    omega_dims = {}
    for p in range(max_p + 2):
        if p > n - 1 or len(ap[p]) == 0:
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
            continue
        if p == 0:
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
            continue

        apm1_set = set(ap[p-1])
        non_allowed = {}
        na_count = 0
        for j, path in enumerate(ap[p]):
            for sign, face in boundary_coeffs(path):
                if len(set(face)) == len(face) and face not in apm1_set:
                    if face not in non_allowed:
                        non_allowed[face] = na_count
                        na_count += 1
        if na_count == 0:
            omega_bases[p] = np.eye(len(ap[p]))
            omega_dims[p] = len(ap[p])
        else:
            P = np.zeros((na_count, len(ap[p])))
            for j, path in enumerate(ap[p]):
                for sign, face in boundary_coeffs(path):
                    if face in non_allowed:
                        P[non_allowed[face], j] += sign
            U, S, Vt = np.linalg.svd(P, full_matrices=True)
            rank = int(sum(s > 1e-10 for s in S))
            ns = Vt[rank:].T
            if ns.shape[1] > 0:
                omega_bases[p] = ns
                omega_dims[p] = ns.shape[1]
            else:
                omega_bases[p] = np.zeros((len(ap[p]), 0))
                omega_dims[p] = 0

    ranks = {}
    for p in range(max_p + 2):
        if p == 0 or omega_dims.get(p, 0) == 0 or p > n - 1:
            ranks[p] = 0
            continue
        bd = np.zeros((len(ap[p-1]), len(ap[p])))
        idx_pm1 = {path: i for i, path in enumerate(ap[p-1])}
        for j, path in enumerate(ap[p]):
            for sign, face in boundary_coeffs(path):
                if face in idx_pm1:
                    bd[idx_pm1[face], j] += sign
        dp_om = bd @ omega_bases[p]
        sv = np.linalg.svd(dp_om, compute_uv=False)
        ranks[p] = int(sum(s > 1e-8 for s in sv))

    return omega_dims, ranks


def main():
    print("=" * 70)
    print("RANK CONSTRAINT ANALYSIS FOR ADJACENT-ODD SEESAW")
    print("=" * 70)

    # Part 1: At n=6, exhaustive: track (rank(d_2), rank(d_3), rank(d_4)) when beta_1>0 or beta_3>0
    print("\n--- Part 1: n=6 exhaustive rank data ---")
    n = 6
    N = 2**(n*(n-1)//2)

    rank_data = defaultdict(list)  # profile -> list of rank tuples

    for bits in range(N):
        A = bits_to_adj(bits, n)
        dims, rnks = compute_chain_ranks(A, n)

        betti = []
        for p in range(n):
            ker = dims.get(p, 0) - rnks.get(p, 0)
            im_next = rnks.get(p+1, 0)
            betti.append(ker - im_next)

        profile = tuple(betti)
        rank_tuple = tuple(rnks.get(p, 0) for p in range(n+1))
        dim_tuple = tuple(dims.get(p, 0) for p in range(n))
        rank_data[profile].append((rank_tuple, dim_tuple))

    for profile in sorted(rank_data.keys()):
        examples = rank_data[profile]
        # Extract unique (dim, rank) patterns
        rank_tuples = [r for r, d in examples]
        dim_tuples = [d for r, d in examples]
        unique_ranks = set(rank_tuples)
        unique_dims = set(dim_tuples)
        print(f"\n  Profile {list(profile)}: {len(examples)} tournaments")
        print(f"    # unique rank patterns: {len(unique_ranks)}")
        print(f"    # unique dim patterns: {len(unique_dims)}")

        # Show a few
        for i, (rt, dt) in enumerate(examples[:3]):
            print(f"    dims={list(dt)}, ranks={list(rt)}")

    # Part 2: At n=6, check: is rank(d_2) CONSTANT for all tournaments?
    print("\n--- Part 2: Is rank(d_p) constant across all n=6 tournaments? ---")
    all_ranks_by_p = defaultdict(set)
    for profile, examples in rank_data.items():
        for rt, dt in examples:
            for p in range(n+1):
                all_ranks_by_p[p].add(rt[p])

    for p in range(n+1):
        vals = sorted(all_ranks_by_p[p])
        print(f"  rank(d_{p}): {vals}")

    # Part 3: At n=7, sample: when beta_1>0, what are the ranks?
    print("\n--- Part 3: n=7 sampled rank data ---")
    n = 7
    rng = np.random.RandomState(12345)
    N = 500

    profile_ranks = defaultdict(list)
    for trial in range(N):
        A = random_tournament(n, rng)
        dims, rnks = compute_chain_ranks(A, n)

        betti = []
        for p in range(n):
            ker = dims.get(p, 0) - rnks.get(p, 0)
            im_next = rnks.get(p+1, 0)
            betti.append(ker - im_next)

        profile = tuple(betti)
        r2 = rnks.get(2, 0)
        r3 = rnks.get(3, 0)
        r4 = rnks.get(4, 0)
        d2 = dims.get(2, 0)
        d3 = dims.get(3, 0)
        d4 = dims.get(4, 0)
        profile_ranks[profile].append((r2, r3, r4, d2, d3, d4))

    for profile in sorted(profile_ranks.keys()):
        examples = profile_ranks[profile]
        r2_vals = sorted(set(e[0] for e in examples))
        r3_vals = sorted(set(e[1] for e in examples))
        r4_vals = sorted(set(e[2] for e in examples))
        print(f"\n  Profile {list(profile)} ({len(examples)} tournaments):")
        print(f"    rank(d_2) values: {r2_vals}")
        print(f"    rank(d_3) values: {r3_vals}")
        print(f"    rank(d_4) values: {r4_vals}")

        # Key: when beta_1=1, rank(d_2) = ker(d_1) - 1 = C(n,2)-n (one less)
        if profile[1] > 0:
            # rank(d_2) determines im(d_2) in Om_1
            print(f"    NOTE: beta_1>0 => rank(d_2) < ker(d_1) = {comb(n,2)-n+1}")
            print(f"    Expected: rank(d_2) = {comb(n,2)-n} (one less)")

        # Key: when beta_3=1, check if rank(d_4) is forced
        if profile[3] > 0:
            print(f"    NOTE: beta_3>0 => ker(d_3) > im(d_4)")
            for e in examples[:3]:
                r2, r3, r4, d2, d3, d4 = e
                ker3 = d3 - r3
                im4 = r4
                print(f"      ker(d_3)={ker3}, im(d_4)={im4}, gap={ker3-im4}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
