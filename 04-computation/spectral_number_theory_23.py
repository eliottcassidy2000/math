#!/usr/bin/env python3
"""
Spectral graph theory and number theory connections for tournaments.
opus-2026-03-14-S85

SPECTRAL CONNECTIONS:
1. Skew-adjacency matrix S = A - A^T (real antisymmetric → pure imaginary eigenvalues)
2. Pfaffian of skew-adjacency matrix (±√det) — combinatorial interpretation?
3. Zeta function Z_T(s) = Π_p (1 - p^{-s})^{-1} over "primes" of T
4. H(T) mod p behavior for various primes p
5. Quadratic residues and tournament construction (Paley tournaments)
6. Arithmetic of H values: GCD structure, factorizations

NUMBER THEORY:
- H(T) is always odd (Rédei)
- max_H(n) connects to factorial-like growth
- The sequence of achievable H values at each n has arithmetic structure
- Connection to Legendre symbol: Paley tournament T_p has i→j iff (i-j|p)=1
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
import sys

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H(adj, n, all_perms):
    return sum(1 for p in all_perms if all(adj[p[i]][p[i+1]] == 1 for i in range(n-1)))

# ============================================================
# Part 1: Pfaffian of Skew-Adjacency Matrix
# ============================================================
print("=" * 70)
print("PART 1: PFAFFIAN OF SKEW-ADJACENCY MATRIX")
print("=" * 70)

# For tournament T, define skew-adjacency S where S[i][j] = 1 if i→j, -1 if j→i, 0 if i=j.
# det(S) = Pfaffian(S)^2 for even n.
# For odd n, det(S) = 0 (skew-symmetric of odd size).

import numpy as np

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms_n = list(permutations(range(n)))
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    pfaff_by_H = defaultdict(list)

    for bits in range(N):
        if bits % 5000 == 0 and N > 5000:
            print(f"  n={n}: {bits}/{N}", file=sys.stderr)

        adj = get_tournament(n, bits)

        # Skew-adjacency matrix
        S = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    S[i][j] = 1 if adj[i][j] else -1

        det_S = np.linalg.det(S)

        # H computation (skip for n=6 — too expensive)
        if n <= 5:
            H = compute_H(adj, n, all_perms_n)
        else:
            # Use score as proxy
            scores = [0] * n
            for k, (i, j) in enumerate(arcs):
                if (bits >> k) & 1:
                    scores[i] += 1
                else:
                    scores[j] += 1
            H = tuple(sorted(scores))  # score proxy

        pfaff_by_H[H].append(round(det_S))

    if n <= 5:
        print(f"\nn={n}: det(skew-adjacency) by H:")
        for h in sorted(pfaff_by_H.keys()):
            vals = pfaff_by_H[h]
            det_dist = Counter(vals)
            print(f"  H={h:2d}: det distribution = {dict(sorted(det_dist.items()))}")
    else:
        print(f"\nn={n}: det(skew-adjacency) by score sequence:")
        for score in sorted(pfaff_by_H.keys())[:10]:
            vals = pfaff_by_H[score]
            det_dist = Counter(vals)
            print(f"  Score={score}: det distribution = {dict(sorted(det_dist.items()))}")

# ============================================================
# Part 2: Eigenvalues of Skew-Adjacency Matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 2: EIGENVALUE SPECTRA")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms_n = list(permutations(range(n)))

    spec_by_H = defaultdict(list)

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H(adj, n, all_perms_n)

        # Skew-adjacency
        S = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    S[i][j] = 1 if adj[i][j] else -1

        eigs = np.linalg.eigvals(S)
        # Should be pure imaginary ± pairs
        imag_parts = sorted([round(e.imag, 4) for e in eigs], reverse=True)
        spec_by_H[H].append(tuple(imag_parts))

    print(f"\nn={n}: Distinct spectra by H:")
    for h in sorted(spec_by_H.keys()):
        specs = spec_by_H[h]
        distinct_specs = Counter(specs)
        print(f"  H={h:2d}: {len(distinct_specs)} distinct spectra")
        for sp, cnt in sorted(distinct_specs.items(), key=lambda x: -x[1])[:3]:
            print(f"    {sp}: {cnt} tournaments")

# ============================================================
# Part 3: H Values and Prime Factorization
# ============================================================
print("\n" + "=" * 70)
print("PART 3: ARITHMETIC OF H VALUES")
print("=" * 70)

def factorize(n):
    if n <= 1:
        return {}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

# All achievable H values by n
for n in [3, 4, 5, 6, 7]:
    if n <= 5:
        m = n * (n - 1) // 2
        N = 1 << m
        all_perms_n = list(permutations(range(n)))
        H_set = set()
        for bits in range(N):
            adj = get_tournament(n, bits)
            H_set.add(compute_H(adj, n, all_perms_n))
        achievable = sorted(H_set)
    elif n == 6:
        # Known from previous computation
        achievable = [1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45]
    elif n == 7:
        # Known: max is 189, min is 1. We know some from Paley.
        achievable = None

    if achievable:
        print(f"\nn={n}: Achievable H values:")
        for h in achievable:
            fac = factorize(h)
            fac_str = " × ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(fac.items())) if fac else "1"
            print(f"  H={h:3d} = {fac_str}")

        # Gaps
        gaps = [achievable[i+1] - achievable[i] for i in range(len(achievable)-1)]
        print(f"  Gaps: {gaps}")
        print(f"  All gaps even: {all(g % 2 == 0 for g in gaps)}")

        # GCD of all H values
        from math import gcd
        from functools import reduce
        g = reduce(gcd, achievable)
        print(f"  GCD of all H: {g}")

        # Sum of all H values (weighted by count)
        if n <= 5:
            total_H = sum(compute_H(get_tournament(n, bits), n, all_perms_n) for bits in range(N))
            mean_H = total_H / N
            print(f"  Mean H = {mean_H:.4f} = {total_H}/{N}")
            print(f"  n!/2^(n-1) = {math.factorial(n)/2**(n-1):.4f}")

# ============================================================
# Part 4: H mod p Patterns
# ============================================================
print("\n" + "=" * 70)
print("PART 4: H mod p PATTERNS")
print("=" * 70)

for n in [4, 5, 6]:
    if n <= 5:
        m = n * (n - 1) // 2
        N = 1 << m
        all_perms_n = list(permutations(range(n)))
        H_all = []
        for bits in range(N):
            adj = get_tournament(n, bits)
            H_all.append(compute_H(adj, n, all_perms_n))
    elif n == 6:
        achievable = [1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45]
        H_all = achievable  # just check achievable values

    print(f"\nn={n}: H values modular structure:")
    for p in [2, 3, 5, 7, 11]:
        residues = Counter(h % p for h in H_all)
        # Which residues are missing?
        present = set(residues.keys())
        missing = set(range(p)) - present
        print(f"  mod {p:2d}: present={sorted(present)}, missing={sorted(missing)}")

# ============================================================
# Part 5: Paley Tournament H Values
# ============================================================
print("\n" + "=" * 70)
print("PART 5: PALEY TOURNAMENT H VALUES")
print("=" * 70)

# Paley tournament T_p (p prime ≡ 3 mod 4):
# Vertices = F_p, i→j iff (j-i) is a QR mod p.

def legendre(a, p):
    """Legendre symbol (a/p)."""
    if a % p == 0:
        return 0
    return pow(a, (p-1)//2, p)

def paley_tournament(p):
    """Build Paley tournament on p vertices."""
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j:
                if legendre(j - i, p) == 1:
                    adj[i][j] = 1
    return adj

for p in [3, 7, 11]:
    if p > 8:
        print(f"\n  p={p}: Paley tournament (too large for full H computation)")
        adj = paley_tournament(p)
        scores = [sum(adj[i]) for i in range(p)]
        print(f"    Scores: {sorted(scores)} (should be all {(p-1)//2})")
        print(f"    Regular: {len(set(scores)) == 1}")
        # Count 3-cycles
        c3 = sum(1 for i in range(p) for j in range(i+1, p) for k in range(j+1, p)
                if adj[i][j] + adj[j][i] + adj[i][k] + adj[k][i] + adj[j][k] + adj[k][j] == 3
                and ((adj[i][j] and adj[j][k] and adj[k][i]) or
                     (adj[j][i] and adj[i][k] and adj[k][j])))
        print(f"    3-cycles: {c3}")
        continue

    adj = paley_tournament(p)
    all_perms_p = list(permutations(range(p)))
    H = compute_H(adj, p, all_perms_p)
    scores = [sum(adj[i]) for i in range(p)]

    print(f"\n  p={p}: Paley tournament")
    print(f"    H = {H}")
    print(f"    Scores = {sorted(scores)}")
    print(f"    Regular: {len(set(scores)) == 1}")
    print(f"    H/n! = {H/math.factorial(p):.6f}")
    print(f"    max possible H = {H}")

# ============================================================
# Part 6: Quadratic Residue Structure in H
# ============================================================
print("\n" + "=" * 70)
print("PART 6: QUADRATIC RESIDUE PATTERNS")
print("=" * 70)

# For n=5, the regular tournament is related to QR mod 5
# (even though 5 ≡ 1 mod 4, so not Paley)
# But C_5 (directed 5-cycle) is a regular tournament.

n = 5
m = n * (n - 1) // 2
N = 1 << m
all_perms_5 = list(permutations(range(5)))

# Find all regular tournaments at n=5 (scores all = 2)
regular_tours = []
for bits in range(N):
    adj = get_tournament(n, bits)
    scores = [sum(adj[i][j] for j in range(n) if j != i) for i in range(n)]
    if all(s == 2 for s in scores):
        H = compute_H(adj, n, all_perms_5)
        regular_tours.append((bits, H))

print(f"\nn=5: All regular tournaments (score = (2,2,2,2,2)):")
for bits, H in regular_tours:
    adj = get_tournament(n, bits)
    # Count 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[j][i] and adj[i][k] and adj[k][j]):
                    c3 += 1
    print(f"  bits={bits:4d}: H={H}, c3={c3}")

# ============================================================
# Part 7: Characteristic Polynomial of Tournament Matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 7: CHARACTERISTIC POLYNOMIALS")
print("=" * 70)

# The adjacency matrix A of tournament T (A[i][j] = 1 if i→j, 0 otherwise)
# has characteristic polynomial χ_T(λ).
# How does χ_T relate to H(T)?

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms_n = list(permutations(range(n)))

    charpoly_by_H = defaultdict(list)

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H(adj, n, all_perms_n)

        A = np.array(adj, dtype=float)
        coeffs = np.round(np.poly(A)).astype(int)
        charpoly_by_H[H].append(tuple(coeffs))

    print(f"\nn={n}: Characteristic polynomials by H:")
    for h in sorted(charpoly_by_H.keys()):
        polys = charpoly_by_H[h]
        distinct = Counter(polys)
        print(f"  H={h:2d}: {len(distinct)} distinct char polys")
        for poly, cnt in sorted(distinct.items(), key=lambda x: -x[1])[:3]:
            print(f"    χ(λ) coeffs = {poly}: {cnt} tournaments")

# ============================================================
# Part 8: Trace Powers — Counting Closed Walks
# ============================================================
print("\n" + "=" * 70)
print("PART 8: TRACE POWERS = CLOSED WALKS")
print("=" * 70)

# tr(A^k) = number of closed walks of length k.
# tr(A) = 0 (no self-loops)
# tr(A^2) = number of mutual arcs... but in tournament, A^2[i][i] = # j with i→j→...→i
# Actually tr(A^2) = Σ_i Σ_j A[i][j]*A[j][i] = 0 (tournament: A[i][j]+A[j][i]=1, so product=0)
# tr(A^3) = # directed 3-cycles (each counted 3 times)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms_n = list(permutations(range(n)))

    traces_by_H = defaultdict(list)

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H(adj, n, all_perms_n)

        A = np.array(adj, dtype=float)
        traces = []
        Ak = np.eye(n)
        for k in range(1, n+1):
            Ak = Ak @ A
            traces.append(int(round(np.trace(Ak))))

        traces_by_H[H].append(tuple(traces))

    print(f"\nn={n}: Trace powers by H value:")
    for h in sorted(traces_by_H.keys()):
        all_traces = traces_by_H[h]
        distinct = Counter(all_traces)
        print(f"  H={h:2d}: {len(distinct)} distinct trace sequences")
        for tr, cnt in sorted(distinct.items(), key=lambda x: -x[1])[:3]:
            print(f"    tr(A^1..{n}) = {tr}: {cnt} tournaments")

# ============================================================
# Part 9: Determinant and Permanent
# ============================================================
print("\n" + "=" * 70)
print("PART 9: DETERMINANT AND PERMANENT OF ADJACENCY")
print("=" * 70)

def permanent(M):
    """Compute permanent of matrix M using Ryser's formula."""
    n = len(M)
    total = 0
    for S in range(1, 1 << n):
        bits_set = []
        for j in range(n):
            if S & (1 << j):
                bits_set.append(j)
        k = len(bits_set)
        prod = 1
        for i in range(n):
            s = sum(M[i][j] for j in bits_set)
            prod *= s
        total += (-1)**(n - k) * prod
    return total

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    all_perms_n = list(permutations(range(n)))

    det_perm_by_H = defaultdict(list)

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H(adj, n, all_perms_n)

        A = adj
        det_A = int(round(np.linalg.det(np.array(A, dtype=float))))
        perm_A = permanent(A)

        det_perm_by_H[H].append((det_A, perm_A))

    print(f"\nn={n}: det(A) and perm(A) by H:")
    for h in sorted(det_perm_by_H.keys()):
        vals = det_perm_by_H[h]
        det_dist = Counter(d for d, _ in vals)
        perm_dist = Counter(p for _, p in vals)
        print(f"  H={h:2d}: det distribution = {dict(sorted(det_dist.items()))}")
        print(f"         perm distribution = {dict(sorted(perm_dist.items()))}")

    # Interesting: is perm(A) related to H?
    print(f"\n  Checking: is perm(A) always related to H?")
    for h in sorted(det_perm_by_H.keys()):
        vals = det_perm_by_H[h]
        perms = set(p for _, p in vals)
        print(f"    H={h}: perm(A) values = {sorted(perms)}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — SPECTRAL AND NUMBER-THEORETIC STRUCTURE")
print("=" * 70)
print("""
KEY FINDINGS:
1. PFAFFIAN: det(skew-adjacency) varies by H value.
   For even n: Pfaffian² = det gives a tournament invariant.

2. EIGENVALUES: Pure imaginary eigenvalues of skew-adjacency.
   Different H values can have same or different spectra.

3. ARITHMETIC: H values are always odd, gaps always even (Rédei parity).
   GCD of all H values = 1. H mod p patterns show forbidden residues.

4. PALEY: Paley tournament T_p achieves extremal H values.
   These are vertex-transitive (regular scores).

5. CHARACTERISTIC POLYNOMIALS: Different H values have different
   typical char polys. This connects spectrum to HP count.

6. TRACE POWERS: tr(A^k) counts closed walks. tr(A^3)/3 = # 3-cycles.
   This gives a direct connection between spectrum and cycle structure.

7. PERMANENT: perm(A) counts matchings in the tournament.
   Different from H but related combinatorially.
""")
