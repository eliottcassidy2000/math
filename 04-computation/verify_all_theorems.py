#!/usr/bin/env python3
"""
Comprehensive verification script for ALL theorems in the repository.

Runs targeted computational checks for each theorem that has a computable claim.
Designed to be re-run from the repo root: python3 04-computation/verify_all_theorems.py

Author: opus-2026-03-06-S17
"""

import sys
import os
import itertools
import random
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph, hamiltonian_path_count)

random.seed(42)


def indep_poly(adj):
    """Compute independence polynomial coefficients of a graph given by adjacency matrix."""
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))


def transpose_tournament(T):
    n = len(T)
    return [[T[j][i] for j in range(n)] for i in range(n)]


def delete_vertex(T, v):
    n = len(T)
    return [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]


def score_seq(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))


PASS = 0
FAIL = 0
SKIP = 0


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  *** FAIL: {name} *** {detail}")


def skip(name, reason=""):
    global SKIP
    SKIP += 1
    print(f"  SKIP: {name} ({reason})")


# ============================================================
# THM-001: Redei's theorem — H(T) is always odd
# ============================================================
print("\n" + "=" * 70)
print("THM-001: Redei's theorem — H(T) always odd")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n - 1) // 2
    fails = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        if H % 2 == 0:
            fails += 1
    check(f"n={n}: all {1 << m} tournaments have odd H", fails == 0, f"{fails} failures")

# n=7 sampled
fails = 0
for _ in range(1000):
    T = random_tournament(7)
    if hamiltonian_path_count(T) % 2 == 0:
        fails += 1
check("n=7 (1000 random): H always odd", fails == 0)

# ============================================================
# THM-002/CONJ-001: OCF — H(T) = I(Omega(T), 2)
# ============================================================
print("\n" + "=" * 70)
print("THM-002/CONJ-001: OCF — H(T) = I(Omega(T), 2)")
print("=" * 70)

for n in [3, 4, 5]:
    m = n * (n - 1) // 2
    fails = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles) if cycles else []
        ip = indep_poly(cg)
        I2 = eval_poly(ip, 2)
        if H != I2:
            fails += 1
    check(f"n={n} exhaustive: H = I(Omega,2)", fails == 0, f"{fails} failures")

for n in [6, 7]:
    fails = 0
    num = 500 if n == 6 else 200
    for _ in range(num):
        T = random_tournament(n)
        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        cg = conflict_graph(cycles) if cycles else []
        ip = indep_poly(cg)
        I2 = eval_poly(ip, 2)
        if H != I2:
            fails += 1
    check(f"n={n} ({num} random): H = I(Omega,2)", fails == 0, f"{fails} failures")

# ============================================================
# THM-004: inshat algebraic identity (inshat-1)/2 = #TypeII
# ============================================================
print("\n" + "=" * 70)
print("THM-004: (inshat-1)/2 = #TypeII positions")
print("=" * 70)

def compute_inshat_typeII(T, v):
    """Compute inshat and TypeII count for vertex v in tournament T."""
    n = len(T)
    others = [i for i in range(n) if i != v]
    # Find all Ham paths of T-v
    Tv = delete_vertex(T, v)
    nv = n - 1

    count_inshat = 0
    count_typeII = 0
    count_paths = 0

    # Enumerate Ham paths of T-v by brute force
    for perm in itertools.permutations(range(nv)):
        # Check if this is a valid Ham path
        valid = True
        for k in range(nv - 1):
            if Tv[perm[k]][perm[k + 1]] != 1:
                valid = False
                break
        if not valid:
            continue
        count_paths += 1

        # Map back to original vertices
        path = [others[perm[k]] for k in range(nv)]

        # Compute signature
        sig = [1 if T[v][path[j]] else 0 for j in range(nv)]

        # Boundary term
        b = sig[0] + (1 - sig[-1])

        # Type I and Type II
        typeI = sum(1 for j in range(nv - 1) if sig[j] == 0 and sig[j + 1] == 1)
        typeII = sum(1 for j in range(nv - 1) if sig[j] == 1 and sig[j + 1] == 0)

        inshat = b + typeI + typeII
        count_inshat += inshat
        count_typeII += typeII

        # Check (inshat-1)/2 = typeII for THIS path
        if (inshat - 1) != 2 * typeII:
            return False, count_paths

    return True, count_paths

for n in [4, 5]:
    m = n * (n - 1) // 2
    fails = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        for v in range(n):
            ok, _ = compute_inshat_typeII(T, v)
            if not ok:
                fails += 1
    check(f"n={n} exhaustive: (inshat-1)/2 = #TypeII", fails == 0, f"{fails} failures")

# ============================================================
# THM-022 Thm 1: T^op preserves I(Omega(T), x)
# ============================================================
print("\n" + "=" * 70)
print("THM-022 Thm 1: I(Omega(T),x) = I(Omega(T^op),x)")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    fails = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        Top = transpose_tournament(T)
        c1 = find_odd_cycles(T)
        c2 = find_odd_cycles(Top)
        ip1 = indep_poly(conflict_graph(c1) if c1 else [])
        ip2 = indep_poly(conflict_graph(c2) if c2 else [])
        if ip1 != ip2:
            fails += 1
    check(f"n={n} exhaustive: I(Omega(T)) = I(Omega(T^op))", fails == 0, f"{fails} failures")

# ============================================================
# THM-019: Omega perfectness / C5 structure
# ============================================================
print("\n" + "=" * 70)
print("THM-019: Omega(T) perfectness — holds n<=7, fails n=8")
print("=" * 70)

def has_odd_hole(adj, k=5):
    """Check if graph has an induced C_k (odd hole of length k)."""
    m = len(adj)
    if m < k:
        return False
    for combo in itertools.combinations(range(m), k):
        # Check if induced subgraph on combo is a cycle
        sub = [[adj[combo[i]][combo[j]] for j in range(k)] for i in range(k)]
        # Each vertex must have degree exactly 2
        degs = [sum(sub[i]) for i in range(k)]
        if not all(d == 2 for d in degs):
            continue
        # Check it's a single cycle of length k
        visited = [False] * k
        visited[0] = True
        cur = 0
        for step in range(k - 1):
            found = False
            for nxt in range(k):
                if sub[cur][nxt] and not visited[nxt]:
                    visited[nxt] = True
                    cur = nxt
                    found = True
                    break
            if not found:
                break
        if all(visited) and sub[cur][0]:
            return True
    return False

# n=5 exhaustive check (small enough)
c5_found = 0
for bits in range(1 << 10):
    T = tournament_from_bits(5, bits)
    cycles = find_odd_cycles(T)
    if len(cycles) < 5:
        continue
    cg = conflict_graph(cycles)
    if has_odd_hole(cg, 5):
        c5_found += 1
check("n=5: no C5 in Omega(T) (exhaustive)", c5_found == 0, f"{c5_found} found")

# n=8 sampled — expect C5
c5_found = 0
for _ in range(100):
    T = random_tournament(8)
    cycles = [c for c in find_odd_cycles(T) if len(c) == 3]  # 3-cycle subgraph
    if len(cycles) < 5:
        continue
    cg = conflict_graph(cycles)
    if has_odd_hole(cg, 5):
        c5_found += 1
check("n=8: C5 found in some Omega_3(T) (expected)", c5_found > 0, f"found in {c5_found}/100")

# ============================================================
# THM-021: Discriminant real-rootedness for n<=8
# ============================================================
print("\n" + "=" * 70)
print("THM-021: Real roots of I(Omega(T),x) via discriminant")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    fails = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cycles = find_odd_cycles(T)
        if not cycles:
            continue
        cg = conflict_graph(cycles)
        ip = indep_poly(cg)
        if len(ip) == 3:  # degree 2
            a1, a2 = ip[1], ip[2]
            if a1 * a1 < 4 * a2:
                fails += 1
    check(f"n={n} exhaustive: disc >= 0 for I(Omega)", fails == 0, f"{fails} failures")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print(f"VERIFICATION SUMMARY: {PASS} passed, {FAIL} failed, {SKIP} skipped")
print("=" * 70)
if FAIL > 0:
    print("*** THERE WERE FAILURES — INVESTIGATE ***")
else:
    print("All checks passed.")
