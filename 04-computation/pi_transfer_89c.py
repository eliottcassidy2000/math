#!/usr/bin/env python3
"""
pi_transfer_89c.py — Transfer matrix analysis for Paley tournaments
opus-2026-03-14-S89c

For a tournament T on n vertices, define the transfer matrix approach:
The Hamiltonian path count can be written using inclusion-exclusion
on the "next vertex" at each step. But for circulant tournaments,
the transition only depends on the difference (j-i) mod p.

Key insight: For P_p, the arc set is {(i,j) : j-i ∈ QR mod p}.
So the adjacency matrix A is a CIRCULANT matrix.

H(T) = sum over all orderings π of vertices such that
π(i) → π(i+1) for all i.

For a circulant tournament, by symmetry we can fix π(0) = 0.
Then H = p × (number of Hamiltonian paths starting at 0).

Let h_0 = number of Ham paths starting at vertex 0.
Then H = p × h_0 (by circulant symmetry).

We can compute h_0 using the transfer matrix on subsets.

But more interestingly: what about the CHARACTERISTIC POLYNOMIAL
of the transfer matrix?
"""

import numpy as np
from itertools import permutations
from math import factorial
from fractions import Fraction

def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = [[] for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].append(j)
    return adj, qr

print("=" * 70)
print("PART 1: Circulant symmetry of H — paths from vertex 0")
print("=" * 70)

for p in [3, 7, 11, 19]:
    adj, qr = paley_tournament(p)

    # Count H paths starting from vertex 0
    dp = [dict() for _ in range(1 << p)]
    dp[1][0] = 1  # mask={0}, endpoint=0

    for mask in range(1, 1 << p):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]

    full = (1 << p) - 1
    h0 = sum(dp[full].values())

    # Count total H
    H = 0
    for start in range(p):
        # Actually just use the full DP
        pass

    # Recompute full H using the optimized method
    dp_full = [dict() for _ in range(1 << p)]
    for v in range(p):
        dp_full[1 << v][v] = 1
    for mask in range(1, 1 << p):
        for v in dp_full[mask]:
            if dp_full[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp_full[new_mask][u] = dp_full[new_mask].get(u, 0) + dp_full[mask][v]
    H = sum(dp_full[full].values())

    print(f"\n  P_{p}: H = {H}, h₀ = {h0}, H/h₀ = {H/h0}, p = {p}")
    print(f"    H = p × h₀? {'✓' if H == p * h0 else '✗'}")

print()
print("=" * 70)
print("PART 2: Distribution of h_v (paths from each vertex)")
print("=" * 70)

for p in [3, 7, 11]:
    adj, qr = paley_tournament(p)

    # For each starting vertex, count paths
    h_v = []
    for start in range(p):
        dp = [dict() for _ in range(1 << p)]
        dp[1 << start][start] = 1
        for mask in range(1, 1 << p):
            for v in dp[mask]:
                if dp[mask][v] == 0:
                    continue
                for u in adj[v]:
                    if mask & (1 << u) == 0:
                        new_mask = mask | (1 << u)
                        dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]
        full = (1 << p) - 1
        h_v.append(sum(dp[full].values()))

    print(f"\n  P_{p}: h_v = {h_v}")
    print(f"    All equal? {'✓' if len(set(h_v)) == 1 else '✗'}")
    print(f"    Sum = {sum(h_v)} = H = {sum(h_v)}")

print()
print("=" * 70)
print("PART 3: Distribution of paths by ending vertex")
print("=" * 70)

# For each start vertex 0, how many paths end at each vertex?
for p in [7, 11]:
    adj, qr = paley_tournament(p)

    dp = [dict() for _ in range(1 << p)]
    dp[1][0] = 1
    for mask in range(1, 1 << p):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]

    full = (1 << p) - 1
    end_counts = {}
    for v in dp[full]:
        end_counts[v] = dp[full][v]

    print(f"\n  P_{p}: paths from 0, by ending vertex:")
    for v in range(p):
        c = end_counts.get(v, 0)
        is_qr = (v % p) in qr if v > 0 else False
        print(f"    end at {v}: {c} {'(QR)' if is_qr else '(NQR)' if v > 0 else '(start)'}")

    # Is there a QR/NQR pattern?
    qr_sum = sum(end_counts.get(v, 0) for v in range(1, p) if v in qr)
    nqr_sum = sum(end_counts.get(v, 0) for v in range(1, p) if v not in qr)
    print(f"    Total ending at QR vertices: {qr_sum}")
    print(f"    Total ending at NQR vertices: {nqr_sum}")
    print(f"    QR/NQR ratio: {Fraction(qr_sum, nqr_sum)}")

print()
print("=" * 70)
print("PART 4: The 'reduced' transfer matrix (QR classes)")
print("=" * 70)

# For P_p, the vertices can be partitioned into two classes:
# QR and NQR (plus 0). The transition matrix between classes
# has a specific structure.
#
# From a QR vertex, we go to: (QR vertex - current) ∈ QR?
# This depends on whether the DIFFERENCE is a QR.
#
# Key: if current is QR and next is QR, then next-current is QR iff
# both are QR and their difference is QR. This is controlled by
# the Jacobi symbol or more precisely by the number of pairs (a,b)
# with a, b, a-b all QR.

for p in [7, 11, 19]:
    adj, qr = paley_tournament(p)

    # Count transitions between classes
    # Classes: QR = {q ∈ QR}, NQR = {n ∈ NQR}
    qr_set = qr
    nqr_set = set(range(1, p)) - qr_set

    # Transition counts from QR to QR, QR to NQR, etc.
    # An arc from i to j exists iff j-i ∈ QR.
    # If i ∈ QR and j ∈ QR: arc iff j-i ∈ QR
    # j-i = j × (1 - i/j) mod p. Since i,j ∈ QR, i/j ∈ QR.
    # So j-i ∈ QR iff j × (1 - i/j) ∈ QR iff (1 - i/j) ∈ QR (since j ∈ QR).
    # Let r = i/j mod p (r ∈ QR \ {0}). Then arc iff (1-r) ∈ QR.

    # Count: how many r ∈ QR \ {0, 1} have (1-r) ∈ QR?
    count_qq = 0  # QR r with QR (1-r)
    count_qn = 0  # QR r with NQR (1-r)
    for r in qr_set:
        if r == 0:
            continue
        one_minus_r = (1 - r) % p
        if one_minus_r == 0:  # r = 1
            continue  # same vertex, skip
        if one_minus_r in qr_set:
            count_qq += 1
        else:
            count_qn += 1

    # r=1 gives 1-r=0 (same vertex)
    # Total QR \ {0}: (p-1)/2 elements. Minus r=0 (not possible) and r=1.
    # So (p-1)/2 - 1 = (p-3)/2 values of r.
    # count_qq + count_qn should be (p-3)/2.

    print(f"\n  P_{p}:")
    print(f"    QR r with QR (1-r): {count_qq}")
    print(f"    QR r with NQR (1-r): {count_qn}")
    print(f"    Total: {count_qq + count_qn} = (p-3)/2 = {(p-3)//2}")

    # For the transition from a QR vertex to another QR vertex:
    # Of the (p-1)/2 - 1 other QR vertices, how many have QR difference?
    # This is count_qq (up to multiplicative factor from j).
    # Actually: from QR vertex i to QR vertex j (j ≠ i):
    # arc from i→j iff j-i ∈ QR.
    # Number of such pairs: ?

    # Direct count:
    arcs_QQ = 0
    arcs_QN = 0
    arcs_NQ = 0
    arcs_NN = 0
    for i in range(p):
        for j in range(p):
            if i == j:
                continue
            if (j - i) % p not in qr_set:
                continue
            # arc i→j exists
            i_class = 'Q' if i in qr_set else ('N' if i in nqr_set else '0')
            j_class = 'Q' if j in qr_set else ('N' if j in nqr_set else '0')
            if i_class == 'Q' and j_class == 'Q':
                arcs_QQ += 1
            elif i_class == 'Q' and j_class == 'N':
                arcs_QN += 1
            elif i_class == 'N' and j_class == 'Q':
                arcs_NQ += 1
            elif i_class == 'N' and j_class == 'N':
                arcs_NN += 1

    nQ = len(qr_set)
    nN = len(nqr_set)
    print(f"    |QR| = {nQ}, |NQR| = {nN}")
    print(f"    Arcs Q→Q: {arcs_QQ} (per QR vertex: {arcs_QQ/nQ:.1f})")
    print(f"    Arcs Q→N: {arcs_QN} (per QR vertex: {arcs_QN/nQ:.1f})")
    print(f"    Arcs N→Q: {arcs_NQ} (per NQR vertex: {arcs_NQ/nN:.1f})")
    print(f"    Arcs N→N: {arcs_NN} (per NQR vertex: {arcs_NN/nN:.1f})")

    # Including vertex 0:
    arcs_from_0 = sum(1 for j in range(1, p) if (j % p) in qr_set)
    arcs_to_0 = sum(1 for i in range(1, p) if ((-i) % p) in qr_set)
    print(f"    Arcs 0→Q: {sum(1 for j in qr_set if j in adj[0])}")
    print(f"    Arcs 0→N: {sum(1 for j in nqr_set if j in adj[0])}")

print()
print("=" * 70)
print("PART 5: H(P_p) via eigenvalue structure?")
print("=" * 70)

# The adjacency matrix A of P_p has eigenvalues:
# λ_0 = (p-1)/2, λ_k = (-1 + χ(k)i√p)/2 for k=1,...,p-1
# where χ(k) = (k/p) = Legendre symbol.
#
# Can we express H in terms of these eigenvalues?
# H = Σ_{all permutations σ ∈ S_p} Π_{i=0}^{p-2} A(σ(i), σ(i+1))
#
# This is NOT the permanent, and NOT a simple function of eigenvalues.
# It's the "path permanent" or "serial correlation."
#
# For CIRCULANT matrices, there might be a simplification.
# Let ω = e^{2πi/p}. Then A = Σ_k λ_k P_k where P_k is the
# projection onto the k-th Fourier mode.
#
# P_k = (1/p) × (ω^{jk})_{j=0,...,p-1} column × (ω^{-jk})_{j=0,...,p-1} row
#
# But H involves a SUM over all orderings, which doesn't simplify
# through circulant structure in an obvious way.

# However, there's a beautiful result for PERMANENTS of circulant matrices!
# The permanent of a circulant matrix C with first row (c_0,...,c_{p-1}) is:
# perm(C) = Σ_{S ⊆ {0,...,p-1}, |S|=p} Π_{k∈S} (Σ_j c_j ω^{jk})
# = Σ over multisets of eigenvalues...
#
# This is Minc's formula / the permanent of a circulant.
# But H is not the permanent.

# Let's just check: is there a numerical pattern?
print("\n  Checking H vs eigenvalue expressions:")
for p in [3, 7, 11, 19, 23]:
    d = (p-1)//2  # = λ_0
    lam_sq = (1 + p) / 4  # |λ_k|² for k≠0

    # Product of all eigenvalues: det(A) = d × (lam_sq)^{d}
    det_A = d * ((p+1)/4)**d

    # Sum of eigenvalues^k
    # Σ λ_k = trace(A) = 0
    # Σ λ_k² = trace(A²) = d(d+1)/2 × p... no
    # trace(A²) = Σ_{i,j} A_{ij}A_{ji} = number of pairs (i,j) with i→j
    # Since A_{ij} + A_{ji} = 1 for i≠j (tournament), A_{ij}A_{ji} = 0
    # So trace(A²) = 0! (No 2-cycles in a tournament.)

    # trace(A³) = number of directed 3-cycles × 3 (?)
    # Actually trace(A³) = Σ_{i,j,k} A_{ij}A_{jk}A_{ki} = 3 × t_3
    # where t_3 is the number of directed 3-cycles.

    if p <= 19:
        adj, qr = paley_tournament(p)
        # Count 3-cycles
        t3 = 0
        for i in range(p):
            for j in adj[i]:
                for k in adj[j]:
                    if i in adj[k]:
                        t3 += 1
        t3 //= 3  # each cycle counted 3 times
    else:
        t3 = None

    H = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}[p]

    print(f"\n  p={p}:")
    print(f"    det(A) = {det_A:.0f}")
    print(f"    H = {H}")
    print(f"    H / det(A) = {H / det_A:.6f}")
    if t3:
        print(f"    t_3 = {t3}")
        print(f"    t_3 / C(p,3) = {t3 / (p*(p-1)*(p-2)//6):.6f}")
        # For Paley: t_3 = p(p-1)(p-5)/24 + p/4 = p(p²-6p+5+6)/24 = p(p²-6p+11)/24
        # Actually t_3 = (1/4)(C(p,3) - p(p-1)/4) for regular tournament
        # = p(p-1)(p-2)/24 - p(p-1)/16
        # Hmm, for doubly regular:
        # t_3 = C(p,3) - (p choose 2) × d/... I should just verify.
        expected_t3 = p * (p*p - 4*p + 7) // 24 if p >= 3 else 0
        # For Paley: every triple has exactly 1 or 3 directed 3-cycles
        # Actually: a triple forms a 3-cycle with prob 1/4 in random tournament
        # For Paley: t_3 = p(p-1)(p-5)/24 when p ≡ 3 mod 4
        # No wait: for a regular tournament, t_3 = p(p-1)(p-5)/24 + p(p-1)/24
        # = p(p-1)(p-4)/24... let me just compute.
        print(f"    t_3 formula check: p(p-1)(2p-1-3(p-1))/... too messy")

print()
print("=" * 70)
print("PART 6: The Hamiltonian count ratio H(P_p) / H_max")
print("=" * 70)

# For Paley primes, is P_p the tournament with max H?
# We know: YES for p=3, 7, 11. NO for p=19 (cyclic does better).
# What about p=23?

# We can't compute max over ALL tournaments on 23 vertices.
# But we can compare with specific competitors.

# The cyclic tournament C_p: vertex set Z/pZ, arc i→j iff
# j-i ∈ {1, 2, ..., (p-1)/2}. This is the "standard" circulant.

for p in [7, 11, 19, 23]:
    # Build cyclic tournament
    adj_cyc = [[] for _ in range(p)]
    half = (p-1) // 2
    for i in range(p):
        for d in range(1, half+1):
            adj_cyc[i].append((i + d) % p)

    # Count H for cyclic tournament
    dp = [dict() for _ in range(1 << p)]
    for v in range(p):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << p):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj_cyc[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]
    full = (1 << p) - 1
    H_cyc = sum(dp[full].values())

    H_paley = {7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}[p]

    print(f"\n  p={p}:")
    print(f"    H(P_p) = {H_paley}")
    print(f"    H(C_p) = {H_cyc}")
    print(f"    Ratio H(P_p)/H(C_p) = {H_paley/H_cyc:.6f}")
    winner = "Paley" if H_paley >= H_cyc else "Cyclic"
    print(f"    Winner: {winner}")

print("\n\nDone!")
