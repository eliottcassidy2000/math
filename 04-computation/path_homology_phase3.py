#!/usr/bin/env python3
"""
PHASE 3: Deep structure of β_3=1 tournaments and Fourier analysis

KEY DISCOVERIES SO FAR:
1. n<=5 tournaments: only β_0=1, β_1∈{0,1}
2. n=6 tournaments: STILL only β_0=1, β_1∈{0,1} (1000 samples)
3. n=7 tournaments: 8.4% have β_3=1! β_2 NEVER nonzero
4. β_2 is ALWAYS 0 for tournaments — β jumps from β_1 to β_3!
5. Circulant C_n^{[n-1]} (complete digraph): β_n-1 explodes (44, 265, ...)

QUESTIONS FOR THIS PHASE:
A. WHY does β_2=0 always for tournaments? Is there a structural reason?
B. What is the "shape" a β_3=1 tournament represents? (S^3?)
C. Does n=8 give β_4>0 or higher?
D. For β_3 tournaments: what is the explicit 3-cycle (homology generator)?
E. Connection to tournament invariants: does β_3 correlate with specific cycle counts?
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, circulant_digraph, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def count_kcycles(A, n, k):
    """Count directed k-cycles."""
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            is_cycle = True
            for i in range(k):
                if A[perm[i]][perm[(i+1)%k]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count // k

def transitivity_ratio(A, n):
    """Fraction of 3-vertex subsets that form transitive triples."""
    trans = 0
    total = 0
    for triple in combinations(range(n), 3):
        total += 1
        i, j, k = triple
        # Check if any orientation is a 3-cycle
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            pass  # 3-cycle
        else:
            trans += 1
    return trans / total if total > 0 else 0

# ===== A: WHY β_2=0 for tournaments? =====
print("=" * 70)
print("A: WHY β_2=0 FOR TOURNAMENTS — ALGEBRAIC ANALYSIS")
print("=" * 70)

print("""
For β_2>0, we need a 2-cycle in Ω_2 not in im(∂_3|Ω_3).

The boundary ∂_2: Ω_2 → Ω_1 has kernel = 2-cycles.
The boundary ∂_3: Ω_3 → Ω_2 has image = 2-boundaries.

For tournaments, Ω_2 = {2-paths (a,b,c) with a→b→c AND a→c (transitive triple)}.
Since ∂(a,b,c) = (b,c) - (a,c) + (a,b), all three edges exist iff a→b, b→c, a→c.

So Ω_2 = {transitive triples}. For β_2>0 we need 2-cycles among transitive triples.
""")

# Check Ω_2, ker∂_2, im∂_3 for some tournaments
for n in [5]:
    print(f"\nn={n}: detailed Ω_p structure")
    for trial in range(20):
        A = random_tournament(n)
        betti = path_betti_numbers(A, n, max_dim=n-1)
        if betti[1] > 0:
            # Compute dimensions of all chain spaces
            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
                om = compute_omega_basis(A, n, p, ap, apm1)
                dim_om = om.shape[1] if om.ndim == 2 else 0
                print(f"  p={p}: |A_p|={len(ap)}, dim(Ω_p)={dim_om}", end="")

                if p > 0 and dim_om > 0:
                    # Compute rank of ∂_p
                    bd = build_full_boundary_matrix(ap, apm1)
                    bd_omega = bd @ om
                    U, S, Vt = np.linalg.svd(bd_omega, full_matrices=False)
                    rank = sum(s > 1e-8 for s in S)
                    kernel_dim = dim_om - rank
                    print(f", rank(∂_p)={rank}, ker(∂_p)={kernel_dim}", end="")
                print()
            print(f"  β = {betti}")
            break

# ===== B: β_3=1 tournament detailed analysis =====
print("\n\n" + "=" * 70)
print("B: β_3=1 TOURNAMENT — FINDING THE 3-CYCLE GENERATOR")
print("=" * 70)

n = 7
found = False
for trial in range(2000):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    if betti[3] > 0:
        found = True
        break

if found:
    print(f"\nFound β_3=1 tournament at trial {trial}")
    t3 = count_3cycles(A, n)
    t4 = count_kcycles(A, n, 4)
    t5 = count_kcycles(A, n, 5)
    tr = transitivity_ratio(A, n)
    print(f"  t3={t3}, t4={t4}, t5={t5}")
    print(f"  transitivity ratio = {tr:.4f}")

    # Show chain complex dimensions
    print(f"\n  Chain complex Ω_p dimensions:")
    chain_dims = []
    omega_bases = []
    allowed_lists = []
    for p in range(n):
        ap = enumerate_allowed_paths(A, n, p)
        apm1 = enumerate_allowed_paths(A, n, p-1) if p > 0 else []
        om = compute_omega_basis(A, n, p, ap, apm1)
        dim_om = om.shape[1] if om.ndim == 2 else 0
        chain_dims.append(dim_om)
        omega_bases.append(om)
        allowed_lists.append(ap)
        print(f"    p={p}: |A_p|={len(ap)}, dim(Ω_p)={dim_om}")

    # Find the actual 3-cycle generator
    # Need ker(∂_3)/im(∂_4)
    if chain_dims[3] > 0:
        # ∂_3: Ω_3 → Ω_2
        ap3 = allowed_lists[3]
        ap2 = allowed_lists[2]
        om3 = omega_bases[3]
        om2 = omega_bases[2]

        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_omega = bd3 @ om3

        U, S, Vt = np.linalg.svd(bd3_omega, full_matrices=True)
        rank3 = sum(s > 1e-8 for s in S)
        ker3 = Vt[rank3:].T  # kernel vectors in Ω_3 coordinates

        print(f"\n    ∂_3: rank={rank3}, ker dim={chain_dims[3]-rank3}")

        # ∂_4: Ω_4 → Ω_3
        if chain_dims[4] > 0:
            ap4 = allowed_lists[4]
            om4 = omega_bases[4]
            bd4 = build_full_boundary_matrix(ap4, ap3)
            bd4_omega = bd4 @ om4
            im4 = bd4_omega  # image in Ω_3

            # Project im4 to A_3 coordinates, then find im(∂_4) ∩ ker(∂_3)
            print(f"    ∂_4: rank={np.linalg.matrix_rank(bd4_omega, tol=1e-8)}")
        else:
            print(f"    Ω_4 = 0, so im(∂_4) = 0")
            print(f"    → H_3 = ker(∂_3) directly")

        # Show the 3-cycle generator
        if ker3.shape[1] > 0:
            for k in range(min(ker3.shape[1], 3)):
                vec = om3 @ ker3[:, k]
                nonzero = [(ap3[i], round(vec[i], 4))
                           for i in range(len(ap3)) if abs(vec[i]) > 1e-8]
                print(f"\n    3-cycle generator #{k}: {len(nonzero)} nonzero paths")
                for path, coeff in nonzero[:20]:
                    print(f"      coeff={coeff:+.4f}: {path}")
                if len(nonzero) > 20:
                    print(f"      ... ({len(nonzero)-20} more)")

# ===== C: β_3 vs tournament invariants =====
print("\n\n" + "=" * 70)
print("C: β_3 CORRELATES — WHAT PREDICTS IT?")
print("=" * 70)

n = 7
data = []
for trial in range(500):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    t3 = count_3cycles(A, n)
    H = ham_path_count(A, n)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    tr = transitivity_ratio(A, n)
    data.append((betti[3], t3, H, scores, tr, betti[1]))

b3_1 = [d for d in data if d[0] == 1]
b3_0 = [d for d in data if d[0] == 0]

print(f"\n  β_3=1: {len(b3_1)} tournaments ({100*len(b3_1)/len(data):.1f}%)")
print(f"  β_3=0: {len(b3_0)} tournaments")

if b3_1:
    print(f"\n  β_3=1 tournaments:")
    print(f"    t3: mean={np.mean([d[1] for d in b3_1]):.1f}, range=[{min(d[1] for d in b3_1)},{max(d[1] for d in b3_1)}]")
    print(f"    H:  mean={np.mean([d[2] for d in b3_1]):.1f}, range=[{min(d[2] for d in b3_1)},{max(d[2] for d in b3_1)}]")
    print(f"    transitivity: mean={np.mean([d[4] for d in b3_1]):.4f}")
    print(f"    β_1: {Counter([d[5] for d in b3_1])}")

    print(f"\n  β_3=0 tournaments:")
    print(f"    t3: mean={np.mean([d[1] for d in b3_0]):.1f}")
    print(f"    H:  mean={np.mean([d[2] for d in b3_0]):.1f}")
    print(f"    transitivity: mean={np.mean([d[4] for d in b3_0]):.4f}")
    print(f"    β_1: {Counter([d[5] for d in b3_0])}")

    # Key: β_3=1 → β_1=0? Or can both be nonzero?
    print(f"\n  Critical: β_3=1 AND β_1=1 simultaneously? {sum(1 for d in data if d[0]==1 and d[5]==1)}")
    print(f"  β_3=1 AND β_1=0? {sum(1 for d in data if d[0]==1 and d[5]==0)}")

# ===== D: Complement duality for path homology =====
print("\n\n" + "=" * 70)
print("D: COMPLEMENT TOURNAMENT AND PATH HOMOLOGY")
print("=" * 70)

print("\nDoes β(T) = β(T^op)? (complement duality)")
n = 7
matches = 0
total_d = 0
for trial in range(200):
    A = random_tournament(n)
    # Build complement
    A_op = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_op[i][j] = 1 - A[i][j]

    b1 = path_betti_numbers(A, n, max_dim=6)
    b2 = path_betti_numbers(A_op, n, max_dim=6)
    total_d += 1
    if b1 == b2:
        matches += 1
    else:
        if total_d <= 5:
            print(f"  MISMATCH: β(T)={b1}, β(T^op)={b2}")

print(f"\n  β(T) = β(T^op): {matches}/{total_d} ({100*matches/total_d:.1f}%)")

# ===== E: Pattern in Euler characteristic for tournaments =====
print("\n\n" + "=" * 70)
print("E: EULER CHAR FOR n=7 TOURNAMENTS")
print("=" * 70)

n = 7
chi_dist = Counter()
chi_by_betti = defaultdict(list)
for trial in range(300):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    chi = sum((-1)**p * betti[p] for p in range(len(betti)))
    chi_dist[chi] += 1
    chi_by_betti[tuple(betti)].append(chi)

print(f"\n  Euler characteristic distribution:")
for chi in sorted(chi_dist.keys()):
    print(f"    χ={chi}: {chi_dist[chi]}")

print(f"\n  By Betti type:")
for bt in sorted(chi_by_betti.keys()):
    chi_vals = chi_by_betti[bt]
    print(f"    β={list(bt)}: χ={chi_vals[0]} ({len(chi_vals)} times)")

# ===== F: Poincaré polynomial =====
print("\n\n" + "=" * 70)
print("F: POINCARÉ POLYNOMIAL P(t) = Σ β_p t^p")
print("=" * 70)

print("\nCirculant digraphs — Poincaré polynomials:")
for n in [5, 7]:
    print(f"\n  n={n} (PRIME — all circulants connected):")
    for size in range(1, n):
        for S in combinations(range(1, n), size):
            A = circulant_digraph(n, list(S))
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 6))
            if any(b > 0 for b in betti[1:]):
                # Poincaré polynomial as string
                terms = []
                for p, b in enumerate(betti):
                    if b > 0:
                        if p == 0:
                            terms.append(str(b))
                        elif b == 1:
                            terms.append(f"t^{p}")
                        else:
                            terms.append(f"{b}t^{p}")
                poincare = " + ".join(terms)
                print(f"    C_{n}^{set(S)}: P(t) = {poincare}")

print("\n\nDone.")
