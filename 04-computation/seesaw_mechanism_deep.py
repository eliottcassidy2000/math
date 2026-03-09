"""
seesaw_mechanism_deep.py — Understand WHY adjacent-odd seesaw holds beyond beta_2=0.

THM-095 proves beta_1*beta_3=0 via:
  beta_2=0 => im(d_3) = ker(d_2) => im(d_2) + im(d_3) = dim(Omega_2)
  This coupling means if im(d_3) drops (beta_3>0), im(d_2) rises (beta_1 drops).

For beta_3*beta_5=0, we'd need beta_4=0 for the same argument.
But beta_4 CAN be nonzero! Yet beta_3*beta_5=0 still holds.

Theory: Maybe when beta_4>0, BOTH beta_3=0 AND beta_5=0.
Let's verify: does beta_4>0 => beta_3=0 AND beta_5=0?

Also: let's track the chain complex dimensions and ranks carefully
to understand the propagation of constraints.

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
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

def compute_full_chain_data(A, n, max_p=None):
    """Compute dimensions, ranks, and Betti numbers for the full chain complex."""
    if max_p is None:
        max_p = n - 1

    ap = {}
    for p in range(max_p + 2):
        if p <= n - 1:
            ap[p] = enumerate_allowed_paths(A, n, p)
        else:
            ap[p] = []

    # Compute Omega bases
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

    # Compute boundary map ranks
    ranks = {}  # ranks[p] = rank of d_p: Omega_p -> Omega_{p-1}
    for p in range(max_p + 2):
        if p == 0 or omega_dims[p] == 0 or p > n - 1:
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

    # Compute Betti numbers
    betti = []
    for p in range(max_p + 1):
        ker_dp = omega_dims[p] - ranks[p]
        im_dp1 = ranks[p+1] if p+1 in ranks else 0
        betti.append(ker_dp - im_dp1)

    return {
        'omega_dims': {p: omega_dims[p] for p in range(max_p + 2)},
        'ranks': {p: ranks[p] for p in range(max_p + 2)},
        'betti': betti,
        'Ap_sizes': {p: len(ap[p]) for p in range(max_p + 2)}
    }


def main():
    print("=" * 70)
    print("DEEP SEESAW MECHANISM ANALYSIS")
    print("=" * 70)

    # Part 1: Full chain data at n=7 for interesting profiles
    print("\n--- Part 1: Chain complex data for each Betti profile at n=7 ---")
    n = 7
    rng = np.random.RandomState(54321)
    N = 2000

    profile_examples = {}
    for trial in range(N):
        A = random_tournament(n, rng)
        data = compute_full_chain_data(A, n)
        bv = tuple(data['betti'])

        if bv not in profile_examples or len(profile_examples[bv]) < 3:
            if bv not in profile_examples:
                profile_examples[bv] = []
            profile_examples[bv].append(data)

    for profile in sorted(profile_examples.keys()):
        examples = profile_examples[profile]
        print(f"\n  Profile {list(profile)} ({len(examples)} examples):")
        for i, data in enumerate(examples[:2]):
            dims = [data['omega_dims'].get(p, 0) for p in range(n+1)]
            rnks = [data['ranks'].get(p, 0) for p in range(n+1)]
            aps = [data['Ap_sizes'].get(p, 0) for p in range(n+1)]
            print(f"    Example {i+1}:")
            print(f"      |A_p|:      {aps[:n]}")
            print(f"      dim(Om_p):   {dims[:n]}")
            print(f"      rank(d_p):  {rnks[:n]}")
            print(f"      ker(d_p):   {[dims[p]-rnks[p] for p in range(n)]}")
            print(f"      im(d_{'{p+1}'}):  {rnks[1:n+1]}")
            print(f"      beta_p:     {list(data['betti'])}")

            # Check conservation law: im(d_p) + im(d_{p+1}) vs dim(Omega_p)
            print(f"      Conservation: ", end="")
            for p in range(1, n-1):
                im_p = rnks[p]
                im_p1 = rnks[p+1]
                dim_p = dims[p]
                gap = dim_p - im_p - im_p1
                if gap != 0:
                    print(f"Om_{p}: gap={gap}(b_{p}={data['betti'][p]}), ", end="")
            print()

    # Part 2: The key question: when beta_4>0, what are the surrounding ranks?
    print("\n--- Part 2: Chain data when beta_4 > 0 ---")
    n = 7
    rng2 = np.random.RandomState(99999)
    b4_examples = []

    for trial in range(5000):
        A = random_tournament(n, rng2)
        data = compute_full_chain_data(A, n)
        if data['betti'][4] > 0:
            b4_examples.append(data)

    print(f"  Found {len(b4_examples)} tournaments with beta_4>0 out of 5000")
    for i, data in enumerate(b4_examples[:5]):
        dims = [data['omega_dims'].get(p, 0) for p in range(n)]
        rnks = [data['ranks'].get(p, 0) for p in range(n+1)]
        print(f"\n    Example {i+1}:")
        print(f"      dim(Om_p):  {dims}")
        print(f"      rank(d_p): {rnks[:n]}")
        print(f"      betti:     {data['betti']}")
        # Is beta_4>0 compatible with beta_3>0 or beta_5>0?
        print(f"      beta_3={data['betti'][3]}, beta_5={data['betti'][5]}")

    # Part 3: Generalized seesaw at n=8 - focusing on beta_3 and beta_5
    print("\n--- Part 3: Chain data at n=8 for beta_3>0 and beta_4>0 cases ---")
    n = 8
    rng3 = np.random.RandomState(77777)
    N = 300
    b3_examples = []
    b4_examples = []
    b5_examples = []

    for trial in range(N):
        A = random_tournament(n, rng3)
        data = compute_full_chain_data(A, n)
        bv = data['betti']
        if bv[3] > 0:
            b3_examples.append(data)
        if bv[4] > 0:
            b4_examples.append(data)
        if bv[5] > 0:
            b5_examples.append(data)

        if (trial + 1) % 100 == 0:
            print(f"  {trial+1}/{N}: b3+={len(b3_examples)}, b4+={len(b4_examples)}, b5+={len(b5_examples)}", flush=True)

    print(f"\n  n=8 summary: beta_3>0: {len(b3_examples)}, beta_4>0: {len(b4_examples)}, beta_5>0: {len(b5_examples)}")

    for label, examples in [("beta_3>0", b3_examples), ("beta_4>0", b4_examples), ("beta_5>0", b5_examples)]:
        if examples:
            print(f"\n  --- {label} ({len(examples)} examples) ---")
            for i, data in enumerate(examples[:3]):
                dims = [data['omega_dims'].get(p, 0) for p in range(n)]
                rnks = [data['ranks'].get(p, 0) for p in range(n+1)]
                print(f"    Example {i+1}: betti={data['betti']}")
                print(f"      dim(Om_p): {dims}")
                print(f"      rank(d_p):{rnks[:n]}")
                # Show the exact seesaw: ker-im at each level
                for p in range(1, n-1):
                    ker = dims[p] - rnks[p]
                    im_next = rnks[p+1]
                    bp = ker - im_next
                    if bp != 0:
                        print(f"      b_{p}={bp}: ker(d_{p})={ker}, im(d_{p+1})={im_next}, gap={bp}")

    # Part 4: Check if beta_3>0 forces rank(d_4) to be maximal
    # (which would prevent beta_3 and beta_5 from coexisting)
    print("\n--- Part 4: Rank saturation analysis ---")
    if b3_examples:
        print(f"  When beta_3>0:")
        for data in b3_examples[:5]:
            dims = [data['omega_dims'].get(p, 0) for p in range(n)]
            rnks = [data['ranks'].get(p, 0) for p in range(n+1)]
            # rank(d_4) = dim(Om_4) - ker(d_4) - how close is this to dim(Om_4)?
            dim4 = dims[4]
            rank4 = rnks[4]
            ker4 = dim4 - rank4
            im5 = rnks[5]
            b4 = ker4 - im5
            print(f"    dim(Om_4)={dim4}, rank(d_4)={rank4}, ker(d_4)={ker4}, im(d_5)={im5}, b_4={b4}")
            # Also check: is im(d_4) = ker(d_3)?
            dim3 = dims[3]
            rank3 = rnks[3]
            ker3 = dim3 - rank3
            im4 = rnks[4]
            b3 = ker3 - im4
            print(f"    dim(Om_3)={dim3}, rank(d_3)={rank3}, ker(d_3)={ker3}, im(d_4)={im4}, b_3={b3}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
