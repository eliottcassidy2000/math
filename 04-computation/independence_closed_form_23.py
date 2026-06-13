#!/usr/bin/env python3
"""
Independence polynomial closed forms and tournament H-spectrum.
opus-2026-03-14-S85

CROWN JEWELS FROM forbidden_H_23.py:
- I(C_k, 2) = 2^k + (-1)^k (Lucas-like!)
- I(P_k, 2) = (2^{k+2} - (-1)^k)/3 (Fibonacci-like!)
- Missing H values {7, 21, 63} have specific structural origins

This script:
1. Verify and extend the closed forms
2. Characterize ALL possible I(G, 2) values for Ω components
3. Determine which connected graphs can/cannot appear in Ω
4. Connect to generating functions and zeta-type objects
"""

import math
from collections import defaultdict, Counter
from itertools import combinations

# ============================================================
# Part 1: Verify I(C_k, 2) = 2^k + (-1)^k
# ============================================================
print("=" * 70)
print("PART 1: I(C_k, 2) = 2^k + (-1)^k (Lucas at x=2)")
print("=" * 70)

def I_cycle(k, x=2):
    """Independence polynomial of C_k evaluated at x, by enumeration."""
    total = 0
    for mask in range(1 << k):
        indep = True
        for i in range(k):
            if (mask & (1 << i)) and (mask & (1 << ((i+1) % k))):
                indep = False
                break
        if indep:
            total += x ** bin(mask).count('1')
    return total

def I_path(k, x=2):
    """Independence polynomial of P_k evaluated at x, by enumeration."""
    if k == 0:
        return 1
    total = 0
    for mask in range(1 << k):
        indep = True
        for i in range(k-1):
            if (mask & (1 << i)) and (mask & (1 << (i+1))):
                indep = False
                break
        if indep:
            total += x ** bin(mask).count('1')
    return total

print("\nVerification — I(C_k, 2) vs 2^k + (-1)^k:")
for k in range(3, 16):
    computed = I_cycle(k) if k <= 12 else None
    formula = 2**k + (-1)**k
    match = "✓" if computed == formula else "✗" if computed is not None else "?"
    if computed is not None:
        print(f"  k={k:2d}: I(C_k,2)={computed:6d}, 2^k+(-1)^k={formula:6d} {match}")
    else:
        print(f"  k={k:2d}: formula={formula:6d}")

print("\nVerification — I(P_k, 2) vs (2^{k+2} - (-1)^k)/3:")
for k in range(1, 16):
    computed = I_path(k) if k <= 12 else None
    formula = (2**(k+2) - (-1)**k) // 3
    match = "✓" if computed == formula else "✗" if computed is not None else "?"
    if computed is not None:
        print(f"  k={k:2d}: I(P_k,2)={computed:6d}, formula={formula:6d} {match}")
    else:
        print(f"  k={k:2d}: formula={formula:6d}")

# ============================================================
# Part 2: General I(G, 2) via transfer matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 2: LUCAS / FIBONACCI CONNECTION")
print("=" * 70)

print("""
Key relations (x = 2):
  I(P_k, x) satisfies: f(k) = f(k-1) + x·f(k-2), f(0)=1, f(1)=1+x
  At x=2: f(k) = f(k-1) + 2·f(k-2), f(0)=1, f(1)=3
  Char eq: t² = t + 2 → t = 2, -1
  Solution: f(k) = (2^{k+2} - (-1)^k) / 3 (Jacobsthal-related!)

  I(C_k, x) = I(P_k, x) - x·I(P_{k-2}, x) [standard identity]
  Wait, let's verify:
""")

# Verify: I(C_k, x) = I(P_k, x) - x·I(P_{k-2}, x)? No...
# Standard: I(C_k) = I(P_{k-1}) + x·I(P_{k-3}) for k≥4? Let me check directly.
# Actually I(C_n, x) = f_n(x) where f satisfies same recurrence as paths
# but with different initial: I(C_n) = L_n(-x) where L_n are Lucas polynomials
# Direct check:
for k in range(3, 10):
    c = I_cycle(k)
    p = I_path(k)
    p1 = I_path(k-1)
    p2 = I_path(k-2) if k >= 2 else 0
    # Common identity: I(C_n, x) = I(P_{n-1}, x) + x * I(P_{n-3}, x) ?
    # Or: I(C_n) = I(P_n) - 2*I(P_{n-2}) ?
    diff = p - c
    print(f"  k={k}: I(C_k)={c}, I(P_k)={p}, I(P_k)-I(C_k)={diff}, 2*I(P_{k-2})={2*p2}")

print()
# Let's try: I(C_n) = I(P_{n-1}) + x·I(P_{n-3})
for k in range(3, 10):
    c = I_cycle(k)
    rhs = I_path(k-1) + 2 * (I_path(k-3) if k >= 3 else 0)
    print(f"  k={k}: I(C_k)={c}, I(P_{k-1})+2·I(P_{k-3})={rhs}, match={c==rhs}")

# ============================================================
# Part 3: Jacobsthal Connection
# ============================================================
print("\n" + "=" * 70)
print("PART 3: JACOBSTHAL NUMBERS AND H SPECTRUM")
print("=" * 70)

# Jacobsthal numbers: J(n) = J(n-1) + 2·J(n-2), J(0)=0, J(1)=1
# = (2^n - (-1)^n) / 3
# So I(P_k, 2) = 4 * J(k) + 1? Let's check:
# J(k) = (2^k - (-1)^k)/3
# 4*J(k) = (2^{k+2} - 4·(-1)^k)/3
# I(P_k) = (2^{k+2} - (-1)^k)/3
# Not the same (factor of 4 vs 1 on (-1)^k). But related!

print("Jacobsthal: J(n) = (2^n - (-1)^n)/3")
print("Path:  I(P_k, 2) = (2^{k+2} - (-1)^k)/3 = 4·J(k) + (-1)^k·(4-1)/3")
print()

for k in range(0, 12):
    J = (2**k - (-1)**k) // 3
    Ip = (2**(k+2) - (-1)**k) // 3
    Ic = 2**k + (-1)**k if k >= 3 else None
    print(f"  k={k:2d}: J(k)={J:5d}, I(P_k,2)={Ip:5d}, I(C_k,2)={'N/A' if Ic is None else Ic:>5}")

# Lucas numbers: L(n) = 2^n + (-1)^n (at specific argument)
# Wait, I(C_k, 2) = 2^k + (-1)^k. This IS the Lucas sequence at α=2, β=-1!
# Standard Lucas: L(n) = α^n + β^n with α+β=1, αβ=-2
# So α=2, β=-1 gives α+β=1, αβ=-2. YES!
# These are Lucas numbers of the second kind for the recurrence t²-t-2=0.

print("\n--- IDENTIFICATION ---")
print("I(C_k, 2) = L_k where L_k = 2^k + (-1)^k")
print("  = Lucas numbers for char eq t² - t - 2 = 0 (roots 2, -1)")
print("I(P_k, 2) = U_k where U_k = (2^{k+1} - (-1)^{k+1})/3")
print("  = Chebyshev-like second kind for same characteristic")
print()
print("Both satisfy SAME recurrence: a(k) = a(k-1) + 2·a(k-2)")
print("Cycle initial: a(1)=1(?), a(2)=... hmm, let's just use closed form")

# ============================================================
# Part 4: Multiplicative Structure of Achievable H
# ============================================================
print("\n" + "=" * 70)
print("PART 4: MULTIPLICATIVE SEMIGROUP OF ACHIEVABLE H VALUES")
print("=" * 70)

# H = product of I(C_i, 2) over components of Ω
# Where each C_i is a connected graph that CAN appear as Ω component.
# Key: which I values are "atoms" (come from connected realizable components)?

# From n≤5 exhaustive: component sizes are 1,2,3,4,5 (cycles in Ω)
# Component I values we observed:
# n=4: empty→1, (1,)→3, (2,)→5
# n=5: empty→1, (1,)→3, (2,)→5, (3,)→9, (4,)→{11,13,15}, (5,)→15

# The (1,) means 1 cycle = single vertex in Ω → I=3
# (2,) means 2 cycles sharing a vertex = edge in Ω → I=5
# etc.

# Let's compute the actual I(component, 2) for each Ω component at n=5

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset({i, j, k}))
                elif adj[j][i] and adj[i][k] and adj[k][j]:
                    cycles.append(frozenset({i, j, k}))
    return cycles

# At n=5, each tournament has a specific Ω and we can verify H = I(Ω, 2)
n = 5
m = n * (n - 1) // 2
N = 1 << m

print(f"\nn=5: Verifying H = I(Ω, 2) for all {N} tournaments:")
component_I_dist = Counter()  # distribution of component-I tuples

verified = 0
failed = 0
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)

    # Find 3-cycles (the odd cycles; 5-cycles also exist but for n=5,
    # all odd cycles are 3-cycles or the full 5-cycle)
    cycles_3 = find_3cycles(adj, n)

    # Build conflict graph
    c = len(cycles_3)
    adj_conflict = defaultdict(set)
    for a in range(c):
        for b in range(a+1, c):
            if cycles_3[a] & cycles_3[b]:
                adj_conflict[a].add(b)
                adj_conflict[b].add(a)

    # Compute I(Ω, 2) directly
    I_omega = 0
    for mask in range(1 << c):
        indep = True
        bits_list = [i for i in range(c) if mask & (1 << i)]
        for a in bits_list:
            for b in bits_list:
                if a < b and b in adj_conflict[a]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            I_omega += 2 ** len(bits_list)

    if I_omega == H:
        verified += 1
    else:
        failed += 1
        if failed <= 5:
            print(f"  MISMATCH at bits={bits}: H={H}, I(Ω_3,2)={I_omega}, #3-cycles={c}")

print(f"\n  Verified: {verified}/{N}")
print(f"  Failed: {failed}/{N}")
if failed > 0:
    print("  NOTE: Failures mean 5-cycles contribute to Ω beyond 3-cycles!")

# ============================================================
# Part 5: Including 5-cycles at n=5
# ============================================================
print("\n" + "=" * 70)
print("PART 5: FULL Ω WITH 5-CYCLES AT n=5")
print("=" * 70)

from itertools import permutations

def find_all_odd_cycles(adj, n, max_len=None):
    """Find all directed odd cycles up to length max_len."""
    if max_len is None:
        max_len = n
    cycles = set()
    for length in range(3, max_len + 1, 2):  # odd lengths only
        for combo in combinations(range(n), length):
            # Check all cyclic orderings
            # Fix smallest element first
            min_v = min(combo)
            rest = [v for v in combo if v != min_v]
            for perm in permutations(rest):
                path = (min_v,) + perm
                is_cycle = True
                for i in range(length):
                    if not adj[path[i]][path[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(path)
                    break  # only need one canonical form per vertex set + direction
    return list(cycles)

verified_full = 0
failed_full = 0
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)

    cycles = find_all_odd_cycles(adj, n)
    c = len(cycles)

    # Build conflict graph
    adj_conflict = defaultdict(set)
    for a in range(c):
        for b in range(a+1, c):
            if set(cycles[a]) & set(cycles[b]):
                adj_conflict[a].add(b)
                adj_conflict[b].add(a)

    I_omega = 0
    for mask in range(1 << c):
        indep = True
        bits_list = [i for i in range(c) if mask & (1 << i)]
        for a in bits_list:
            for b in bits_list:
                if a < b and b in adj_conflict[a]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            I_omega += 2 ** len(bits_list)

    if I_omega == H:
        verified_full += 1
    else:
        failed_full += 1
        if failed_full <= 5:
            print(f"  MISMATCH: bits={bits}, H={H}, I(Ω_full,2)={I_omega}")
            print(f"    #cycles={c}, cycles={cycles}")

print(f"\nFull Ω (3+5 cycles) at n=5:")
print(f"  Verified: {verified_full}/{N}")
print(f"  Failed: {failed_full}/{N}")

# ============================================================
# Part 6: Generating Function for Achievable H
# ============================================================
print("\n" + "=" * 70)
print("PART 6: GENERATING FUNCTION FOR H SPECTRUM")
print("=" * 70)

# The "H-zeta function" ζ_H(s) = Σ_{T} H(T)^{-s}
# Or the "H-generating function" F_n(z) = Σ_T z^{H(T)}

for n in [4, 5]:
    m = n * (n-1) // 2
    N = 1 << m

    H_dist = Counter()
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        H_dist[H] += 1

    print(f"\nn={n}: H-generating function F_n(z) = Σ N_H · z^H:")
    for H in sorted(H_dist.keys()):
        print(f"  z^{H}: coefficient = {H_dist[H]}")

    # Check: F_n(1) = 2^m (total tournaments)
    total = sum(H_dist.values())
    print(f"  F_n(1) = {total} = 2^{m} = {2**m}")

    # F_n'(1) / F_n(1) = <H> = mean H
    mean_H = sum(H * count for H, count in H_dist.items()) / total
    print(f"  <H> = {mean_H:.4f}, n!/2^{n-1} = {math.factorial(n)/2**(n-1):.4f}")

    # Product structure: if H = I(Ω, 2) and components are independent,
    # the generating function should factor!
    # F_n(z) = ??? Is there a nice product/Euler product?

# ============================================================
# Part 7: The OCF bridge — I(Ω, 2) = H is the OCF evaluation
# ============================================================
print("\n" + "=" * 70)
print("PART 7: OCF BRIDGE — I(Ω, 2) = H")
print("=" * 70)

print("""
The Odd-Cycle Collection Formula (OCF):
  H(T) = Σ_{S independent in Ω(T)} 2^|S|

This IS the independence polynomial I(Ω(T), x) evaluated at x = 2!

Key consequences:
1. H is ALWAYS odd (Rédei) because I(G, 2) = 2^k + (-1)^k for cycles,
   and products of odd numbers are odd.

2. H = 1 iff Ω = empty (no odd cycles) iff T is transitive.

3. H = 3 iff Ω has exactly one vertex (one odd cycle, isolated).

4. H = 5 iff Ω has exactly one edge (two cycles sharing a vertex).

5. H = 7 iff Ω = K₃ (triangle in conflict graph).
   But K₃ as Ω component requires three pairwise-sharing cycles
   that don't share a common vertex — IMPOSSIBLE by THM-201!

6. H = 9 = 3² iff Ω has two isolated vertices (two independent cycles).

7. The closed form I(C_k, 2) = 2^k + (-1)^k = Lucas(k, 2, -1)
   means cycle-components of Ω of size k contribute factor 2^k + (-1)^k.

DEEP CONNECTION: The evaluation point x=2 corresponds to:
- Binary tournament choice (each arc has 2 orientations)
- OCF weight 2^|S| per independent set of size |S|
- Jacobsthal recurrence a(k) = a(k-1) + 2·a(k-2)
""")

# ============================================================
# Part 8: What IS Ω for max-H tournaments?
# ============================================================
print("=" * 70)
print("PART 8: CONFLICT GRAPH Ω OF MAX-H TOURNAMENTS")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n-1) // 2
    N = 1 << m

    if n > 6:
        break

    max_H = 0
    max_bits_list = []
    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)
        if H > max_H:
            max_H = H
            max_bits_list = [bits]
        elif H == max_H:
            max_bits_list.append(bits)

    print(f"\nn={n}: max H = {max_H}, achieved by {len(max_bits_list)} tournaments")

    # Analyze Ω for first max-H tournament
    bits = max_bits_list[0]
    adj = get_tournament(n, bits)

    if n <= 6:
        cycles = find_all_odd_cycles(adj, n)
        c = len(cycles)
        print(f"  #odd cycles: {c}")

        # Conflict graph edges
        edges = []
        for a in range(c):
            for b in range(a+1, c):
                if set(cycles[a]) & set(cycles[b]):
                    edges.append((a, b))
        print(f"  #conflict edges: {len(edges)}")

        # Independence number
        max_indep = 0
        for mask in range(1 << c):
            indep = True
            bits_list = [i for i in range(c) if mask & (1 << i)]
            for a in bits_list:
                for b in bits_list:
                    if a < b and (a, b) in [(e[0], e[1]) for e in edges]:
                        indep = False
                        break
                if not indep:
                    break
            if indep:
                size = len(bits_list)
                if size > max_indep:
                    max_indep = size

        print(f"  Independence number α(Ω): {max_indep}")
        print(f"  I(Ω, 2) = {max_H}")
        print(f"  log₂(H) ≈ {math.log2(max_H):.4f}")

        # Chromatic number of Ω (clique cover = indep sets of complement)
        # Actually, for max H we want large independent sets.
        # H = Σ 2^|S| is maximized when Ω has many large independent sets.

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — INDEPENDENCE POLYNOMIAL CLOSED FORMS")
print("=" * 70)
print("""
CROWN JEWELS:

1. I(C_k, 2) = 2^k + (-1)^k [Lucas sequence for t²-t-2=0]
   - I(C_3)=7, I(C_5)=31, I(C_7)=127 (Mersenne primes!)
   - I(C_4)=17, I(C_6)=65

2. I(P_k, 2) = (2^{k+2} - (-1)^k)/3 [Jacobsthal-related]
   - I(P_1)=3, I(P_2)=5, I(P_3)=11, I(P_4)=21, I(P_5)=43

3. H(T) = I(Ω(T), 2) where Ω = odd-cycle conflict graph
   - Product formula: H = ∏ I(Cᵢ, 2) over components
   - x=2 evaluation ↔ binary arc choice ↔ OCF weight

4. Forbidden H: 7 = I(K₃, 2) forbidden because K₃ cannot be Ω-component
   - 21 forbidden despite P₄ being connected with I=21
   - Deeper constraints on which graphs CAN appear as Ω components

5. JACOBSTHAL BRIDGE: The recurrence a(k) = a(k-1) + 2·a(k-2)
   governing both path and cycle I-values is EXACTLY the
   Jacobsthal recurrence. The "2" in "+2·a(k-2)" is the
   OCF evaluation point = number of arc orientations!
""")
