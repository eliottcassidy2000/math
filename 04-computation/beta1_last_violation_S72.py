#!/usr/bin/env python3
"""Find the last violation: (c3v, sc4) almost determines β₁ at n=6.
What's the 1 remaining case? What extra invariant closes the gap?"""
import numpy as np
from math import comb
from collections import defaultdict, Counter
from itertools import combinations
import sys
sys.stdout.reconfigure(line_buffering=True)

PRIME = 2**31 - 1

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def per_vertex_3cycles(A, n):
    c3v = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3v[i] += 1; c3v[j] += 1; c3v[k] += 1
    return tuple(sorted(c3v))

def count_sc4(A, n):
    sc4 = 0
    for verts in combinations(range(n), 4):
        scores = sorted(sum(A[vi, vj] for vj in verts if vj != vi) for vi in verts)
        if scores != [0, 1, 2, 3]:
            sc4 += 1
    return sc4

def count_sc5(A, n):
    """Count 5-vertex strongly-connected subtournaments."""
    if n < 5: return 0
    sc5 = 0
    for verts in combinations(range(n), 5):
        scores = sorted(sum(A[vi, vj] for vj in verts if vj != vi) for vi in verts)
        if scores != [0, 1, 2, 3, 4]:
            sc5 += 1
    return sc5

def count_4cycles(A, n):
    c4 = 0
    for i in range(n):
        for j in range(n):
            if j == i or A[i,j] == 0: continue
            for k in range(n):
                if k in (i,j) or A[j,k] == 0: continue
                for l in range(n):
                    if l in (i,j,k) or A[k,l] == 0 or A[l,i] == 0: continue
                    c4 += 1
    return c4 // 4

def compute_beta1(A, n):
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                edges.append((i, j))
                edge_idx[(i, j)] = len(edges) - 1
    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k in (i, j) or A[j][k] == 0: continue
                if A[i][k] == 1:
                    triples.append((i, j, k))
    if not triples:
        return comb(n, 2) - (n - 1)
    bd = np.zeros((len(edges), len(triples)), dtype=np.int64)
    for t, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            bd[edge_idx[(j, k)], t] = (bd[edge_idx[(j, k)], t] + 1) % PRIME
        if (i, k) in edge_idx:
            bd[edge_idx[(i, k)], t] = (bd[edge_idx[(i, k)], t] - 1) % PRIME
        if (i, j) in edge_idx:
            bd[edge_idx[(i, j)], t] = (bd[edge_idx[(i, j)], t] + 1) % PRIME
    M = bd % PRIME
    nrows, ncols = M.shape
    rank = 0
    for col in range(ncols):
        nonzero = np.where(M[rank:, col] != 0)[0]
        if len(nonzero) == 0: continue
        pivot = nonzero[0] + rank
        if pivot != rank: M[[rank, pivot]] = M[[pivot, rank]]
        inv = pow(int(M[rank, col]), PRIME - 2, PRIME)
        M[rank] = M[rank] * inv % PRIME
        factors = M[:, col].copy()
        factors[rank] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank])) % PRIME
        rank += 1
    return max(0, comb(n, 2) - (n - 1) - rank)

print("=" * 72)
print("LAST VIOLATION ANALYSIS")
print("=" * 72)

n = 6
total = 2**(n*(n-1)//2)

# Collect all data
groups = defaultdict(lambda: defaultdict(list))
import time
t0 = time.time()

for bits in range(total):
    A = bits_to_adj(bits, n)
    b1 = compute_beta1(A, n)
    c3v = per_vertex_3cycles(A, n)
    sc4 = count_sc4(A, n)
    groups[(c3v, sc4)][b1].append(bits)

print(f"Collected in {time.time()-t0:.1f}s")

# Find violations
print("\nViolations where (c3v, sc4) gives multiple β₁ values:")
for key, b1_dict in sorted(groups.items()):
    if len(b1_dict) > 1:
        c3v, sc4 = key
        print(f"\n  c3v={c3v}, sc4={sc4}:")
        for b1_val, bits_list in sorted(b1_dict.items()):
            print(f"    β₁={b1_val}: {len(bits_list)} tournaments")

            # Examine these tournaments for distinguishing features
            for bits in bits_list[:3]:
                A = bits_to_adj(bits, n)
                c4 = count_4cycles(A, n)
                scores = tuple(sorted(int(np.sum(A[i])) for i in range(n)))
                print(f"      bits={bits}: scores={scores}, c4={c4}")

        # Test additional discriminants on this violation
        b0_bits = b1_dict.get(0, [])
        b1_bits = b1_dict.get(1, [])

        if b0_bits and b1_bits:
            # c4 discriminant?
            c4_b0 = set()
            c4_b1 = set()
            for bits in b0_bits:
                A = bits_to_adj(bits, n)
                c4_b0.add(count_4cycles(A, n))
            for bits in b1_bits:
                A = bits_to_adj(bits, n)
                c4_b1.add(count_4cycles(A, n))

            overlap = c4_b0 & c4_b1
            if not overlap:
                print(f"    ★ c4 RESOLVES: β₁=0 has c4∈{c4_b0}, β₁=1 has c4∈{c4_b1}")
            else:
                print(f"    ✗ c4 overlap: {overlap}")

            # Score seq discriminant?
            ss_b0 = set()
            ss_b1 = set()
            for bits in b0_bits:
                A = bits_to_adj(bits, n)
                ss_b0.add(tuple(sorted(int(np.sum(A[i])) for i in range(n))))
            for bits in b1_bits:
                A = bits_to_adj(bits, n)
                ss_b1.add(tuple(sorted(int(np.sum(A[i])) for i in range(n))))

            ss_overlap = ss_b0 & ss_b1
            if not ss_overlap:
                print(f"    ★ score_seq RESOLVES: β₁=0 has {ss_b0}, β₁=1 has {ss_b1}")
            else:
                print(f"    score_seq overlap: {ss_overlap}")

            # sc5 discriminant? (if n=6)
            if n >= 5:
                sc5_b0 = set()
                sc5_b1 = set()
                for bits in b0_bits[:20]:  # sample for speed
                    A = bits_to_adj(bits, n)
                    sc5_b0.add(count_sc5(A, n))
                for bits in b1_bits[:20]:
                    A = bits_to_adj(bits, n)
                    sc5_b1.add(count_sc5(A, n))

                sc5_overlap = sc5_b0 & sc5_b1
                if not sc5_overlap:
                    print(f"    ★ sc5 RESOLVES: β₁=0 has sc5∈{sc5_b0}, β₁=1 has sc5∈{sc5_b1}")
                else:
                    print(f"    sc5 overlap: {sc5_overlap}")

# Final: test (c3v, sc4, c4) at n=6
print("\n" + "=" * 72)
print("FINAL: Does (c3v, sc4, c4) determine β₁ at n=6?")
print("=" * 72)

groups_final = defaultdict(set)
for bits in range(total):
    A = bits_to_adj(bits, n)
    b1 = compute_beta1(A, n)
    c3v = per_vertex_3cycles(A, n)
    sc4 = count_sc4(A, n)
    c4 = count_4cycles(A, n)
    groups_final[(c3v, sc4, c4)].add(b1)

violations = sum(1 for v in groups_final.values() if len(v) > 1)
print(f"  (c3v, sc4, c4) determines β₁: {'YES!' if violations==0 else f'NO ({violations} violations)'}")

if violations == 0:
    print("\n  ★★★ FORMULA FOUND: β₁ is determined by (c3v, sc4, c4) ★★★")
    print("  This gives β₁ in O(n⁴) instead of O(n³·rank) — MASSIVE speedup!")
    print("  NEXT: prove this algebraically → unlock instant β₁ at any n")

print("\nDONE")
