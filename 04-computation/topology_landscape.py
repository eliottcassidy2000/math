#!/usr/bin/env python3
"""
TOPOLOGICAL LANDSCAPE of circulant digraphs and tournaments.

Maps Betti numbers to geometric shapes and studies:
1. What topological spaces appear at each n?
2. The "odd-only" phenomenon for tournaments (β_2 = 0 always)
3. How connection set structure determines topology
4. Circulant topology at larger n via Fourier v3
"""
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_fourier_v3 import fourier_betti_v3
from path_homology_v2 import path_betti_numbers, circulant_digraph

# ===== GEOMETRIC DICTIONARY =====
def classify_topology(betti):
    """Map Betti numbers to known topological space."""
    b = tuple(betti)
    # Trim trailing zeros
    while len(b) > 1 and b[-1] == 0:
        b = b[:-1]

    if b == (1,):
        return "point (contractible)"
    elif b == (1, 1):
        return "S^1 (circle)"
    elif b == (1, 2, 1):
        return "T^2 (torus)"
    elif b == (1, 0, 1):
        return "S^2 (2-sphere)"
    elif b == (1, 0, 0, 1):
        return "S^3 (3-sphere)"
    elif b == (1, 0, 0, 0, 1):
        return "S^4 (4-sphere)"
    elif b[0] == 1 and b[1] == 0 and all(x == 0 for x in b[2:-1]) and b[-1] > 0:
        dim = len(b) - 1
        mult = b[-1]
        if mult == 1:
            return f"S^{dim} ({dim}-sphere)"
        else:
            return f"{mult}×S^{dim}"
    elif b == (1, 1, 0, 1):
        return "S^1 ∨ S^3"  # wedge
    else:
        # Try to decompose
        parts = []
        for i in range(1, len(b)):
            if b[i] > 0:
                parts.append(f"β_{i}={b[i]}")
        return "mixed: " + ", ".join(parts)


# ===== 1. Complete topology census for small circulants =====
print("=" * 70)
print("TOPOLOGICAL CENSUS: ALL CIRCULANTS C_n^S")
print("=" * 70)

for n in [5, 7, 11]:
    print(f"\nn={n} (prime):", flush=True)
    topology_count = Counter()
    examples = defaultdict(list)

    for size in range(1, min(n, 5)):
        for S in combinations(range(1, n), size):
            S_set = set(S)
            if n <= 7:
                A = circulant_digraph(n, list(S))
                betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
            else:
                betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 5))

            topo = classify_topology(betti)
            topology_count[topo] += 1
            if len(examples[topo]) < 3:
                examples[topo].append((S_set, betti))

    for topo, count in topology_count.most_common():
        exs = examples[topo]
        ex_str = "; ".join(f"S={e[0]}" for e in exs[:2])
        print(f"  {topo}: {count} ({ex_str})")

# ===== 2. How |S| determines topology =====
print("\n\n" + "=" * 70)
print("|S| vs TOPOLOGY for prime n")
print("=" * 70)

for n in [7, 11, 13]:
    print(f"\nn={n}:", flush=True)
    for size in range(1, min(n, 6)):
        topo_by_size = Counter()
        for S in combinations(range(1, n), size):
            S_set = set(S)
            if n <= 7:
                A = circulant_digraph(n, list(S))
                betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
            else:
                betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 5))
            topo = classify_topology(betti)
            topo_by_size[topo] += 1
        print(f"  |S|={size}: {dict(topo_by_size.most_common(5))}", flush=True)

# ===== 3. The "tournament angle" — near-complete circulants =====
print("\n\n" + "=" * 70)
print("NEAR-COMPLETE CIRCULANTS (tournament-like)")
print("=" * 70)

print("\nA tournament on n vertices has ~n(n-1)/2 directed edges.")
print("C_n^S has n*|S| edges. For tournament-like density, |S| ≈ (n-1)/2.\n")

for n in [7, 11, 13]:
    target_size = (n - 1) // 2
    print(f"\nn={n}, target |S|={target_size}:", flush=True)
    topo_count = Counter()
    for S in combinations(range(1, n), target_size):
        S_set = set(S)
        if n <= 7:
            A = circulant_digraph(n, list(S))
            betti = path_betti_numbers(A, n, max_dim=min(n-1, 5))
        else:
            betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 5))
        topo = classify_topology(betti)
        topo_count[topo] += 1
    for topo, count in topo_count.most_common():
        print(f"  {topo}: {count}")

# ===== 4. Complement pairs: S vs [n-1]\S =====
print("\n\n" + "=" * 70)
print("COMPLEMENT PAIRS: C_n^S vs C_n^{complement(S)}")
print("=" * 70)

for n in [7, 11]:
    print(f"\nn={n}:", flush=True)
    full = set(range(1, n))
    seen = set()
    for size in range(1, n // 2 + 1):
        for S in combinations(range(1, n), size):
            S_set = frozenset(S)
            comp = frozenset(full - S_set)
            pair = (min(S_set, comp), max(S_set, comp))
            if pair in seen:
                continue
            seen.add(pair)

            if n <= 7:
                A1 = circulant_digraph(n, list(S_set))
                b1 = path_betti_numbers(A1, n, max_dim=min(n-1, 5))
                A2 = circulant_digraph(n, list(comp))
                b2 = path_betti_numbers(A2, n, max_dim=min(n-1, 5))
            else:
                b1 = fourier_betti_v3(set(S_set), n, max_dim=min(n-1, 5))
                b2 = fourier_betti_v3(set(comp), n, max_dim=min(n-1, 5))

            t1 = classify_topology(b1)
            t2 = classify_topology(b2)

            if t1 != t2 or any(b > 0 for b in b1[2:]):
                print(f"  S={set(S_set)}: {t1} (β={b1})")
                print(f"  S'={set(comp)}: {t2} (β={b2})")
                print()

# ===== 5. The critical question: β_2 for circulants vs tournaments =====
print("\n\n" + "=" * 70)
print("β_2 LANDSCAPE: When is β_2 > 0?")
print("=" * 70)

print("\nFor tournaments: β_2 = 0 ALWAYS (HYP-301).")
print("For circulants: β_2 can be nonzero.\n")
print("Which connection sets give β_2 > 0?")

for n in [5, 6, 7, 8, 9, 11, 13]:
    beta2_cases = []
    for size in range(1, min(n, 5)):
        for S in combinations(range(1, n), size):
            S_set = set(S)
            if n <= 8:
                A = circulant_digraph(n, list(S))
                betti = path_betti_numbers(A, n, max_dim=min(n-1, 4))
            else:
                betti = fourier_betti_v3(S_set, n, max_dim=min(n-1, 4))
            if len(betti) > 2 and betti[2] > 0:
                beta2_cases.append((S_set, betti))

    if beta2_cases:
        print(f"\nn={n}: {len(beta2_cases)} sets with β_2 > 0")
        for S_set, betti in beta2_cases[:5]:
            # Check if S+S contains all of S (partial closure)
            sums = set()
            for a in S_set:
                for b in S_set:
                    if a != b:
                        sums.add((a + b) % n)
            closure = sums & S_set
            print(f"  S={S_set}: β={betti}, |S∩(S+S)|={len(closure)}/{len(S_set)}")
    else:
        print(f"n={n}: NO β_2 > 0 (up to |S|=3)")

print("\nDone.")
