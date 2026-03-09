"""
beta3_seesaw_mechanism.py — Understand the seesaw beta_1*beta_3=0

Known: beta_1*beta_3=0 for ALL tournaments (THM-095).
Mechanism: im(d₂) acts as mediator. rank(d₂) takes exactly 2 values: n-1 and n.
  - rank(d₂) = n-1 <=> beta_1 = 1
  - rank(d₂) = n <=> beta_1 = 0

Question: How does rank(d₂) affect beta_3?
When rank(d₂) = n-1 (beta_1=1): does this FORCE beta_3=0?

From the chain complex:
  Omega_1 <- d_2 <- Omega_2 <- d_3 <- Omega_3 <- d_4 <- Omega_4

beta_2 = 0 means: im(d₃) = ker(d₂)

beta_3 = dim(ker(d₃)) - rank(d₄)

The seesaw works through im(d₂) dimensions:
  dim(ker(d₂)) = dim(Omega_2) - rank(d₂)

When beta_1=1: rank(d₂) = n-1 (one less), so ker(d₂) is one bigger.
Since beta_2=0: rank(d₃) = dim(ker(d₂)).
So rank(d₃) is one bigger when beta_1=1.

This means: rank(d₃) = dim(Omega_2) - (n-1) when beta_1=1
            rank(d₃) = dim(Omega_2) - n when beta_1=0

ker(d₃) = dim(Omega_3) - rank(d₃)

When beta_1=1: ker(d₃) = dim(Omega_3) - (dim(Omega_2) - (n-1))
         = dim(Omega_3) - dim(Omega_2) + n - 1
When beta_1=0: ker(d₃) = dim(Omega_3) - dim(Omega_2) + n

So ker(d₃) is ONE LESS when beta_1=1. This could be the reason beta_3=0 when beta_1=1:
ker(d₃) might drop to 0 (or = rank(d₄)), making beta_3=0.

Let's verify this mechanism quantitatively.

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

def compute_omega_basis(ap, p):
    if not ap.get(p, []):
        return np.zeros((0, 0)), 0
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
        return np.zeros((len(ap[p]), 0)), 0


def full_chain_complex_data(A, n):
    """Compute all chain complex dimensions and ranks."""
    ap = {}
    for p in range(min(6, n)):
        ap[p] = enumerate_allowed_paths(A, n, p)

    omega_bases = {}
    omega_dims = {}
    for p in range(min(6, n)):
        if not ap.get(p, []):
            omega_bases[p] = np.zeros((0, 0))
            omega_dims[p] = 0
        else:
            omega_bases[p], omega_dims[p] = compute_omega_basis(ap, p)

    # Build all boundary maps and compute ranks
    ranks = {}
    for p in range(1, min(6, n)):
        if omega_dims.get(p, 0) == 0 or omega_dims.get(p-1, 0) == 0:
            ranks[p] = 0
            continue

        bd = np.zeros((len(ap[p-1]), len(ap[p])))
        idx_prev = {path: i for i, path in enumerate(ap[p-1])}
        for j, path in enumerate(ap[p]):
            for sign, face in boundary_coeffs(path):
                if face in idx_prev:
                    bd[idx_prev[face], j] += sign

        Op = omega_bases[p]
        d_om = bd @ Op

        if omega_dims[p-1] < len(ap[p-1]):
            Oprev_pinv = np.linalg.pinv(omega_bases[p-1])
            d_omega = Oprev_pinv @ d_om
        else:
            d_omega = d_om

        sv = np.linalg.svd(d_omega, compute_uv=False)
        ranks[p] = int(sum(s > 1e-8 for s in sv))

    # Compute Betti numbers and ker dimensions
    kers = {}
    bettis = {}
    for p in range(min(6, n)):
        rk_p = ranks.get(p, 0)
        kers[p] = omega_dims.get(p, 0) - rk_p
        rk_next = ranks.get(p+1, 0)
        bettis[p] = max(0, kers[p] - rk_next)

    return {
        'omega_dims': omega_dims,
        'ranks': ranks,
        'kers': kers,
        'bettis': bettis,
    }


def main():
    print("=" * 70)
    print("BETA_3 SEESAW MECHANISM — Why beta_1*beta_3=0")
    print("=" * 70)

    # Part 1: Chain complex data grouped by beta_1 at n=6
    print("\n--- Part 1: Chain complex by beta_1 at n=6 (exhaustive) ---")
    n = 6
    by_b1 = defaultdict(list)

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        data = full_chain_complex_data(A, n)
        b1 = data['bettis'].get(1, 0)
        b3 = data['bettis'].get(3, 0)
        by_b1[b1].append(data)

    for b1_val in sorted(by_b1.keys()):
        entries = by_b1[b1_val]
        print(f"\n  beta_1 = {b1_val} ({len(entries)} tournaments):")

        # Omega dimensions
        for p in range(6):
            dims = [e['omega_dims'].get(p, 0) for e in entries]
            print(f"    Omega_{p}: min={min(dims)}, max={max(dims)}, avg={sum(dims)/len(dims):.1f}")

        # Ranks
        for p in range(1, 6):
            rks = [e['ranks'].get(p, 0) for e in entries]
            print(f"    rank(d_{p}): min={min(rks)}, max={max(rks)}, avg={sum(rks)/len(rks):.1f}")

        # ker(d_3) and beta_3
        kers3 = [e['kers'].get(3, 0) for e in entries]
        b3s = [e['bettis'].get(3, 0) for e in entries]
        rk4 = [e['ranks'].get(4, 0) for e in entries]
        print(f"    ker(d_3): min={min(kers3)}, max={max(kers3)}, avg={sum(kers3)/len(kers3):.1f}")
        print(f"    rank(d_4): min={min(rk4)}, max={max(rk4)}, avg={sum(rk4)/len(rk4):.1f}")
        print(f"    beta_3: {Counter(b3s)}")

    # Part 2: The seesaw equation
    print("\n--- Part 2: Seesaw equation verification ---")
    n = 6
    # From theory:
    # rank(d_2) = n-1+beta_1_offset where beta_1_offset = 0 when beta_1=1, 1 when beta_1=0
    # Wait, reversed: rank(d_2) = n-1 when beta_1=1, rank(d_2) = n when beta_1=0
    # ker(d_2) = dim(Omega_2) - rank(d_2)
    # rank(d_3) = ker(d_2) (since beta_2=0)
    # ker(d_3) = dim(Omega_3) - rank(d_3) = dim(Omega_3) - ker(d_2)
    #          = dim(Omega_3) - dim(Omega_2) + rank(d_2)
    # When beta_1=1: ker(d_3) = dim(Omega_3) - dim(Omega_2) + (n-1)
    # When beta_1=0: ker(d_3) = dim(Omega_3) - dim(Omega_2) + n

    seesaw_ok = 0
    seesaw_fail = 0

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        data = full_chain_complex_data(A, n)
        b1 = data['bettis'].get(1, 0)

        od2 = data['omega_dims'].get(2, 0)
        od3 = data['omega_dims'].get(3, 0)
        rk2 = data['ranks'].get(2, 0)

        # Verify rank(d_2) = n - beta_1
        expected_rk2 = n - b1
        if rk2 != expected_rk2:
            # Check n-1 formula instead
            if rk2 == n - 1 and b1 == 1:
                pass
            elif rk2 == n and b1 == 0:
                pass

        # Verify ker(d_3) formula
        ker_d2 = od2 - rk2
        rk3 = data['ranks'].get(3, 0)
        # beta_2=0 => rk3 = ker_d2
        if rk3 != ker_d2:
            pass  # Shouldn't happen

        ker_d3 = data['kers'].get(3, 0)
        expected_ker_d3 = od3 - ker_d2
        if ker_d3 == expected_ker_d3:
            seesaw_ok += 1
        else:
            seesaw_fail += 1

    print(f"  ker(d_3) = dim(Omega_3) - ker(d_2) formula: OK={seesaw_ok}, FAIL={seesaw_fail}")

    # Part 3: Critical computation — ker(d_3) difference between beta_1=0 and beta_1=1
    print("\n--- Part 3: ker(d_3) shift from beta_1 ---")
    # When beta_1=1: ker(d_3) = od3 - od2 + (n-1)
    # When beta_1=0: ker(d_3) = od3 - od2 + n
    # Difference: ker_d3(beta_1=0) - ker_d3(beta_1=1) = 1 (for same Omega dims)
    # But Omega dims also differ between beta_1=0 and beta_1=1 tournaments!

    # Let's compute the actual ker(d_3) - rank(d_4) = beta_3 value
    # and see if the beta_1=1 case systematically has smaller ker(d_3) - rank(d_4)

    b1_0_data = []
    b1_1_data = []

    for bits in range(2**(n*(n-1)//2)):
        A = bits_to_adj(bits, n)
        data = full_chain_complex_data(A, n)
        b1 = data['bettis'].get(1, 0)
        b3 = data['bettis'].get(3, 0)
        od3 = data['omega_dims'].get(3, 0)
        ker3 = data['kers'].get(3, 0)
        rk4 = data['ranks'].get(4, 0)
        od4 = data['omega_dims'].get(4, 0)

        if b1 == 0:
            b1_0_data.append((od3, ker3, rk4, od4, b3))
        else:
            b1_1_data.append((od3, ker3, rk4, od4, b3))

    print(f"\n  beta_1=0 ({len(b1_0_data)} tours):")
    print(f"    Omega_3: {min(d[0] for d in b1_0_data)}-{max(d[0] for d in b1_0_data)}")
    print(f"    ker(d_3): {min(d[1] for d in b1_0_data)}-{max(d[1] for d in b1_0_data)}")
    print(f"    rank(d_4): {min(d[2] for d in b1_0_data)}-{max(d[2] for d in b1_0_data)}")
    print(f"    Omega_4: {min(d[3] for d in b1_0_data)}-{max(d[3] for d in b1_0_data)}")
    print(f"    beta_3: {Counter(d[4] for d in b1_0_data)}")

    print(f"\n  beta_1=1 ({len(b1_1_data)} tours):")
    print(f"    Omega_3: {min(d[0] for d in b1_1_data)}-{max(d[0] for d in b1_1_data)}")
    print(f"    ker(d_3): {min(d[1] for d in b1_1_data)}-{max(d[1] for d in b1_1_data)}")
    print(f"    rank(d_4): {min(d[2] for d in b1_1_data)}-{max(d[2] for d in b1_1_data)}")
    print(f"    Omega_4: {min(d[3] for d in b1_1_data)}-{max(d[3] for d in b1_1_data)}")
    print(f"    beta_3: {Counter(d[4] for d in b1_1_data)}")

    # Part 4: Does beta_1=1 force ker(d_3)=0?
    print("\n--- Part 4: Does beta_1=1 force ker(d_3)=0? ---")
    ker3_when_b1_1 = Counter(d[1] for d in b1_1_data)
    print(f"  ker(d_3) when beta_1=1: {dict(sorted(ker3_when_b1_1.items()))}")

    # Part 5: Same analysis at n=7
    print("\n--- Part 5: Chain complex by beta_1 at n=7 (sampled) ---")
    n = 7
    rng = np.random.RandomState(42)
    N = 500
    b1_0_7 = []
    b1_1_7 = []

    for trial in range(N):
        A = random_tournament(n, rng)
        data = full_chain_complex_data(A, n)
        b1 = data['bettis'].get(1, 0)
        b3 = data['bettis'].get(3, 0)
        od3 = data['omega_dims'].get(3, 0)
        ker3 = data['kers'].get(3, 0)
        rk4 = data['ranks'].get(4, 0)
        od4 = data['omega_dims'].get(4, 0)

        if b1 == 0:
            b1_0_7.append((od3, ker3, rk4, od4, b3))
        else:
            b1_1_7.append((od3, ker3, rk4, od4, b3))

    print(f"  beta_1=0 ({len(b1_0_7)} tours):")
    if b1_0_7:
        print(f"    Omega_3: {min(d[0] for d in b1_0_7)}-{max(d[0] for d in b1_0_7)}")
        print(f"    ker(d_3): {min(d[1] for d in b1_0_7)}-{max(d[1] for d in b1_0_7)}")
        print(f"    beta_3: {Counter(d[4] for d in b1_0_7)}")

    print(f"  beta_1=1 ({len(b1_1_7)} tours):")
    if b1_1_7:
        print(f"    Omega_3: {min(d[0] for d in b1_1_7)}-{max(d[0] for d in b1_1_7)}")
        print(f"    ker(d_3): {min(d[1] for d in b1_1_7)}-{max(d[1] for d in b1_1_7)}")
        print(f"    beta_3: {Counter(d[4] for d in b1_1_7)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
