#!/usr/bin/env python3
"""
w_{n-5} DECOMPOSITION: Why does tr(c_2) at n=7 depend on more than t_3?

KEY INSIGHT from the pair-partition analysis:
  - overlap=0 (disjoint edges): IMPOSSIBLE at n=7 (need 8 vertices)
  - overlap=1 (one shared vertex): always contributes 0 (signed sum cancels)
  - overlap=2 (two shared vertices): depends on "complementary cycle structure"
  - overlap=3 (three shared vertices, fully consecutive): depends on 5-vertex transits

The overlap=1 cancellation happens because:
  When 4 edges split into a connected block + an isolated block,
  the isolated block's signed sum over its complement permutations = 0.

PROOF: For k vertices, sum over P(k,2) ordered pairs of (T(a,b)-1/2)
= sum T(a,b) - k(k-1)/2 = C(k,2) - k(k-1)/2 = 0.

So the non-trivial contributions come from overlap=2 (two chains) and
overlap=3 (one long chain). Let me compute what tournament invariant
each depends on.

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
from math import factorial, comb
import numpy as np
import random

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def count_3cycles_subset(A, S):
    """Count 3-cycles in induced sub-tournament on vertex set S."""
    count = 0
    S = list(S)
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for k in range(j+1, len(S)):
                a,b,c = S[i],S[j],S[k]
                count += (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a])
    return count

# =====================================================================
# VERIFY: overlap=1 always contributes 0
# =====================================================================
print("=" * 70)
print("VERIFY: overlap=1 contributions are always 0")
print("=" * 70)

n = 7
# overlap=1 positions: (0,1,3,5), (0,2,3,5), (0,2,4,5)
overlap1_positions = [(0,1,3,5), (0,2,3,5), (0,2,4,5)]

for trial in range(3):
    random.seed(trial * 999)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    for pos in overlap1_positions:
        total = 0
        for p in permutations(range(n)):
            s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
            total += s[pos[0]] * s[pos[1]] * s[pos[2]] * s[pos[3]]
        print(f"  Trial {trial}, pos={pos}: sum = {total:.6f} {'✓' if abs(total)<1e-10 else '✗'}")

# =====================================================================
# PROVE: why overlap=1 gives 0
# =====================================================================
print(f"\n  PROOF: For any 'disconnected' position quadruple (connected + isolated),")
print(f"  the isolated block sums to 0 because:")
print(f"  sum_{{ordered pairs (a,b) from k vertices}} (T(a,b)-1/2) = C(k,2)-k(k-1)/2 = 0")

for k in range(2, 7):
    random.seed(k*111)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    # Check sum over P(k,2) ordered pairs from first k vertices
    S = list(range(k))
    total = sum(A[a][b] - 0.5 for a in S for b in S if a != b)
    print(f"  k={k}: sum = {total:.6f} = {comb(k,2)} - {k*(k-1)//2} {'✓' if abs(total)<1e-10 else '✗'}")

# =====================================================================
# ANALYZE overlap=2 (two-chain pattern): positions like (0,1,3,4)
# =====================================================================
print("\n" + "=" * 70)
print("OVERLAP=2 ANALYSIS: Two-chain pattern")
print("=" * 70)

# Position (0,1,3,4): edges (v0,v1),(v1,v2),(v3,v4),(v4,v5)
# Two chains: {v0,v1,v2} and {v3,v4,v5}, using 6 of 7 vertices.
# Contribution = sum_{v6} G(V\{v6}) where G involves "complementary cycle structure".

# F(T) for a 3-element subset T:
# F(T) = sum over 6 orderings of (T(x,y)-1/2)(T(y,z)-1/2)
# = 3/2 if T is cyclic, -1/2 if T is transitive

print("  F(T) values:")
for trial in range(3):
    random.seed(trial * 555)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Check F for a few triples
    for T in [(0,1,2), (0,1,3), (2,3,4)]:
        F_val = 0
        for p in permutations(T):
            F_val += (A[p[0]][p[1]] - 0.5) * (A[p[1]][p[2]] - 0.5)
        is_cyclic = (A[T[0]][T[1]]*A[T[1]][T[2]]*A[T[2]][T[0]] +
                    A[T[0]][T[2]]*A[T[2]][T[1]]*A[T[1]][T[0]]) > 0
        expected = 1.5 if is_cyclic else -0.5
        if trial == 0:
            print(f"    T={T}: F={F_val:.4f}, cyclic={is_cyclic}, expected={expected:.1f} {'✓' if abs(F_val-expected)<0.01 else '✗'}")

# G(S) for 6-element subset S:
# G(S) = sum_{T⊂S, |T|=3} F(T)*F(S\T)
# = 2*c_3(S) - 4*c_2(S) + 5
# where c_3(S) = #cyclic triples, c_2(S) = #triples T with T cyclic and S\T transitive

print(f"\n  Computing G(S) = sum_T F(T)*F(S\\T) for 6-vertex subsets:")
random.seed(2026)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

for v_del in range(n):
    S = [v for v in range(n) if v != v_del]
    c3_S = count_3cycles_subset(A, S)

    # Compute G directly
    G_direct = 0
    for T in combinations(S, 3):
        comp = [v for v in S if v not in T]
        F_T = 0
        for p in permutations(T):
            F_T += (A[p[0]][p[1]] - 0.5) * (A[p[1]][p[2]] - 0.5)
        F_comp = 0
        for p in permutations(comp):
            F_comp += (A[p[0]][p[1]] - 0.5) * (A[p[1]][p[2]] - 0.5)
        G_direct += F_T * F_comp

    # Count c_2: triples where T cyclic and comp transitive
    c2_count = 0
    for T in combinations(S, 3):
        comp = tuple(v for v in S if v not in T)
        T_cyc = count_3cycles_subset(A, T)
        comp_cyc = count_3cycles_subset(A, comp)
        if T_cyc > 0 and comp_cyc == 0:
            c2_count += 1

    # Formula: G = 2*c3 - 4*c2 + 5
    G_formula = 2*c3_S - 4*c2_count + 5
    print(f"    v_del={v_del}: c3(S)={c3_S}, c2={c2_count}, G_direct={G_direct:.2f}, G_formula={G_formula:.2f} {'✓' if abs(G_direct-G_formula)<0.01 else '✗'}")

# =====================================================================
# ANALYZE overlap=3 (fully consecutive): positions like (0,1,2,3)
# =====================================================================
print("\n" + "=" * 70)
print("OVERLAP=3 ANALYSIS: Fully consecutive chain")
print("=" * 70)

# Position (0,1,2,3): edges (v0,v1),(v1,v2),(v2,v3),(v3,v4)
# All 4 edges consecutive, using 5 vertices.
# Contribution = sum over 7! perms of prod_{i=0}^{3} s_i
# where s_i = T(v_i,v_{i+1}) - 1/2.
# This equals (n-2)! * sum over P(n,5) ordered 5-tuples of
# (T(a,b)-1/2)(T(b,c)-1/2)(T(c,d)-1/2)(T(d,e)-1/2)

# Wait: positions 0,1,2,3 use v0,...,v4. Positions v5,v6 are free.
# Number of perms: P(7,5)*2! = 7!/(7-5)! * 2! = 2520*2 = 5040 = 7!. ✓
# More precisely: choose 5 vertices and order them as v0-v4: P(7,5) = 2520.
# Then v5,v6 are the remaining 2 in 2! orders. So each ordered 5-tuple
# appears with multiplier 2! = 2.

# The contribution is:
# 2! * sum over P(7,5) ordered 5-tuples of prod(T-1/2)

# What does sum over ordered 5-tuples of prod(T(a,b)-1/2)(T(b,c)-1/2)(T(c,d)-1/2)(T(d,e)-1/2) depend on?

print("  Computing 5-vertex chain sum for each 5-element subset:")
random.seed(2026)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

chain5_sums = {}
for S in combinations(range(n), 5):
    total = 0
    for p in permutations(S):
        prod = 1
        for i in range(4):
            prod *= (A[p[i]][p[i+1]] - 0.5)
        total += prod
    chain5_sums[S] = total
    c3 = count_3cycles_subset(A, S)
    # Also count 5-cycles
    c5 = 0
    for p in permutations(S):
        if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
            c5 += 1
    c5 //= 5
    H5 = sum(1 for p in permutations(S) if all(A[p[i]][p[i+1]] for i in range(4)))
    print(f"    S={S}: chain_sum={total:6.2f}, c3={c3}, c5={c5}, H5={H5}")

# Check: chain_sum = sum_P prod(f_i-1/2) where product is over all 4 edges
# This is related to H(T|S) by expansion:
# prod(f_i-1/2) = sum_{k=0}^4 (-1/2)^{4-k} * e_k(f)
# But f_i^2 = f_i, so the expansion simplifies.

# =====================================================================
# FULL w_{n-5} DECOMPOSITION
# =====================================================================
print("\n" + "=" * 70)
print("FULL w_{n-5} DECOMPOSITION AT n=7")
print("=" * 70)

# Total w_{n-5} = sum over 15 position quadruples of their contributions
# = 3 * overlap3_total + 9 * overlap2_total + 3 * overlap1_total (= 0)

# overlap3: 3 positions (0,1,2,3), (1,2,3,4), (2,3,4,5)
# overlap2: 9 positions in 3 groups
# overlap1: 3 positions, always = 0

random.seed(7777)
results = []
for trial in range(15):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3cycles(A)
    H = sum(1 for p in permutations(range(n)) if all(A[p[i]][p[i+1]] for i in range(n-1)))

    # Compute contributions by overlap type
    overlap_sums = {1: 0, 2: 0, 3: 0}
    w_total = 0

    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        for i in range(n-1):
            for j in range(i+1, n-1):
                for k in range(j+1, n-1):
                    for l in range(k+1, n-1):
                        val = s[i]*s[j]*s[k]*s[l]
                        w_total += val
                        # Count overlaps
                        edges = [(i,i+1),(j,j+1),(k,k+1),(l,l+1)]
                        ov = sum(1 for a in range(4) for b in range(a+1,4)
                                if set(edges[a]) & set(edges[b]))
                        overlap_sums[ov] += val

    # Compute invariant: sum_v c_3(T-v)
    sum_c3_del = sum(count_3cycles_subset(A, [w for w in range(n) if w != v]) for v in range(n))

    # Compute c2 invariant: sum_v #{T⊂V\{v}: T cyclic and (V\{v})\T transitive}
    sum_c2_del = 0
    for v in range(n):
        S = [w for w in range(n) if w != v]
        for T in combinations(S, 3):
            comp = [w for w in S if w not in T]
            if count_3cycles_subset(A, T) > 0 and count_3cycles_subset(A, comp) == 0:
                sum_c2_del += 1

    results.append({
        't3': t3, 'H': H, 'w': w_total,
        'ov1': overlap_sums[1], 'ov2': overlap_sums[2], 'ov3': overlap_sums[3],
        'sum_c3_del': sum_c3_del, 'sum_c2_del': sum_c2_del
    })

print(f"\n  {'t3':>3s} {'H':>4s} {'w':>8s} {'ov1':>6s} {'ov2':>8s} {'ov3':>8s} {'Σc3_del':>8s} {'Σc2_del':>8s}")
print("  " + "-" * 70)
for r in sorted(results, key=lambda x: x['t3']):
    print(f"  {r['t3']:3d} {r['H']:4d} {r['w']:8.0f} {r['ov1']:6.0f} {r['ov2']:8.0f} {r['ov3']:8.0f} {r['sum_c3_del']:8d} {r['sum_c2_del']:8d}")

# Try fitting: w = a*ov2 + b*ov3 (since ov1=0)
# Also: w = a*sum_c3_del + b*sum_c2_del + c
y = np.array([r['w'] for r in results])
X1 = np.array([[r['sum_c3_del'], r['sum_c2_del'], 1] for r in results])
c1, _, _, _ = np.linalg.lstsq(X1, y, rcond=None)
err1 = max(abs(y - X1 @ c1))
print(f"\n  w = {c1[0]:.4f}*Σc3_del + {c1[1]:.4f}*Σc2_del + {c1[2]:.4f}")
print(f"  Max error: {err1:.4f}")

# Also try: w = a*t3 + b*sum_c2_del + c
X2 = np.array([[r['t3'], r['sum_c2_del'], 1] for r in results])
c2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
err2 = max(abs(y - X2 @ c2))
print(f"  w = {c2[0]:.4f}*t3 + {c2[1]:.4f}*Σc2_del + {c2[2]:.4f}")
print(f"  Max error: {err2:.4f}")

# Just t3?
X3 = np.array([[r['t3'], 1] for r in results])
c3, _, _, _ = np.linalg.lstsq(X3, y, rcond=None)
err3 = max(abs(y - X3 @ c3))
print(f"  w = {c3[0]:.4f}*t3 + {c3[1]:.4f}")
print(f"  Max error: {err3:.4f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
