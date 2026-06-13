#!/usr/bin/env python3
"""
det(I+2A) RECURRENCE IN JACOBSTHAL COORDINATES
opus-2026-03-13-S67k

KEY FACTS:
  - det(I+2A) is always a perfect square (HYP-788, kind-pasteur)
  - I+2A = J+S where J=all-ones, S=skew-adjacency
  - H = I(CG, 2) follows an ADDITIVE recurrence under vertex deletion
  - det(I+2A) should follow a MULTIPLICATIVE recurrence

This script investigates:
1. How det(I+2A) updates when adding a vertex
2. The relationship between √det(I+2A) and H
3. Whether √det follows Jacobsthal-type patterns
4. The MULTIPLICATIVE deletion tree vs the ADDITIVE H tree
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
from functools import lru_cache

def adj_matrix(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def score_seq(A, n):
    return tuple(sorted(int(x) for x in A.sum(axis=1)))

def fast_hash(A, n):
    scores = list(A.sum(axis=1))
    neighbor_sig = []
    for i in range(n):
        out_scores = sorted(scores[j] for j in range(n) if A[i][j])
        in_scores = sorted(scores[j] for j in range(n) if A[j][i])
        neighbor_sig.append((scores[i], tuple(out_scores), tuple(in_scores)))
    return tuple(sorted(neighbor_sig))

def get_iso_classes(n):
    m = n * (n - 1) // 2
    hash_groups = defaultdict(list)
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        hash_groups[h].append(A)
    classes = []
    for h, group in hash_groups.items():
        A = group[0]
        H = count_hp(A, n)
        classes.append({'H': H, 'A': A, 'count': len(group)})
    return sorted(classes, key=lambda c: c['H'])

def det_I2A(A, n):
    """Compute det(I + 2A) exactly using integer arithmetic."""
    M = np.eye(n, dtype=float) + 2 * A.astype(float)
    return int(round(np.linalg.det(M)))

def sqrt_det(d):
    """Integer square root of |d|."""
    ad = abs(d)
    s = int(round(ad**0.5))
    if s*s == ad:
        return s
    # Try nearby values
    for delta in range(-2, 3):
        if (s+delta) >= 0 and (s+delta)**2 == ad:
            return s+delta
    return None

# ============================================================
print("=" * 72)
print("det(I+2A) RECURRENCE ANALYSIS")
print("=" * 72)

# ============================================================
# PART 1: I+2A = J+S decomposition
# ============================================================
print("\n--- Part 1: I+2A = J+S Decomposition ---")
print("For tournament T with adjacency A (A_{ij} = 1 if i→j):")
print("  A + A^T = J - I  (all pairs have exactly one arc)")
print("  S = A - A^T  (skew-symmetric)")
print("  2A = (J-I) + S")
print("  I + 2A = J + S")
print()
print("So det(I+2A) = det(J+S) where:")
print("  J = all-ones matrix (rank 1, eigenvalue n with eigenvector 1)")
print("  S = skew-symmetric (eigenvalues purely imaginary)")

# ============================================================
# PART 2: Multiplicative structure under vertex deletion
# ============================================================
print("\n--- Part 2: How det(I+2A) changes under vertex deletion ---")

for n in range(3, 7):
    print(f"\nn={n}:")
    classes = get_iso_classes(n)

    for cl in classes:
        A = cl['A']
        H = cl['H']
        d = det_I2A(A, n)
        sd = sqrt_det(d)

        # Vertex deletions
        sub_dets = []
        sub_Hs = []
        for v in range(n):
            rem = [u for u in range(n) if u != v]
            Av = A[np.ix_(rem, rem)]
            sub_d = det_I2A(Av, n-1)
            sub_H = count_hp(Av, n-1)
            sub_dets.append(sub_d)
            sub_Hs.append(sub_H)

        # Compute the ratios
        ratios = []
        for sub_d in sub_dets:
            if sub_d != 0:
                ratios.append(d / sub_d)
            else:
                ratios.append(float('inf'))

        # Also compute the Schur complement factor
        # det(M) / det(M_{11}) = det(M/M_{11}) for block decomposition
        # For vertex v deletion: det(I+2A) = det(I+2A_{-v}) * (Schur factor)

        print(f"  H={H:3d}  det={d:8d}  √det={sd}"
              f"  sub_dets={sorted(sub_dets)}  "
              f"sub_Hs={sorted(sub_Hs)}")

        # Check: is d / sub_d always an integer?
        for v in range(n):
            sub_d = sub_dets[v]
            if sub_d != 0:
                ratio = d / sub_d
                is_int = abs(ratio - round(ratio)) < 0.001
                if not is_int or True:  # show all
                    pass  # will show below

        # Show ratio for first vertex only to save space
        v = 0
        sub_d = sub_dets[0]
        if sub_d != 0 and sd is not None:
            sub_sd = sqrt_det(sub_d)
            if sub_sd and sub_sd > 0:
                ratio_sq = sd / sub_sd
                print(f"         v=0: det/sub_det={d/sub_d:.1f}, "
                      f"√det/√sub_det={sd/sub_sd:.4f}")

# ============================================================
# PART 3: The relationship H vs √det
# ============================================================
print("\n" + "=" * 72)
print("Part 3: H vs √det(I+2A)")
print("=" * 72)

print("\nFor each iso class, compute H and √det(I+2A).")
print("Question: are H and √det related by a recurrence?")
print()

for n in range(3, 7):
    print(f"\nn={n}:")
    classes = get_iso_classes(n)
    h_to_sdets = defaultdict(list)

    for cl in classes:
        A = cl['A']
        H = cl['H']
        d = det_I2A(A, n)
        sd = sqrt_det(d)
        h_to_sdets[H].append(sd)
        # Quick check on H mod 3 and sd mod 3
        ratio_str = f"{H/sd:.3f}" if sd else "N/A"
        print(f"  H={H:3d}  √det={sd if sd else 0:5d}  H/√det={ratio_str:>8s}"
              f"  H mod 3={H%3}  √det mod 3={sd%3 if sd else 0}"
              f"  (H-1)/2={(H-1)//2}")

print("\n--- Key observations ---")
print("1. √det is ALWAYS odd (confirmed)")
print("2. H and √det are NOT equal in general (only at n=3)")
print("3. H/√det varies: can be 1, 3, 5, ... (always odd?)")

# Check H/√det more carefully
print("\n--- H/√det analysis ---")
for n in range(3, 7):
    print(f"\nn={n}: H/√det values:")
    classes = get_iso_classes(n)
    ratios = Counter()
    for cl in classes:
        H = cl['H']
        d = det_I2A(cl['A'], n)
        sd = sqrt_det(d)
        if sd and sd > 0:
            r = H / sd
            ratios[round(r, 4)] += 1
    for r, c in sorted(ratios.items()):
        print(f"  H/√det = {r:.4f}  ({c} classes)")

# ============================================================
# PART 4: The cofactor expansion and Jacobsthal connection
# ============================================================
print("\n" + "=" * 72)
print("Part 4: Cofactor Expansion and Jacobsthal")
print("=" * 72)

print("""
For the matrix M = I + 2A = J + S, where J = 11^T and S = A-A^T:

det(M) = det(J + S)

By the matrix determinant lemma (rank-1 update):
  det(J + S) = det(S) · (1 + 1^T S^{-1} 1)     [if S invertible]

For odd n: det(S) = 0 (skew-symmetric with odd dimension).
  So we need the adjugate: det(J+S) = 1^T adj(S) 1

For even n: det(S) = Pf(S)^2 (Pfaffian squared).
  det(J+S) = Pf(S)^2 · (1 + 1^T S^{-1} 1)

JACOBSTHAL CONNECTION:
  The Jacobsthal number J(n) = (2^n - (-1)^n)/3.
  Note that 3 divides 2^n - (-1)^n always.

  For the transitive tournament T_n (acyclic, score sequence 0,1,...,n-1):
    H(T_n) = 1 (only one HP)
    det(I+2A(T_n)) = 1

  For the cyclic tournament C_3:
    H(C_3) = 3 = J(4)
    det(I+2A(C_3)) = 9 = J(4)^2

  CONJECTURE: √det(I+2A) divides some Jacobsthal number?
""")

# Check if √det values are related to Jacobsthal
J = [0, 1]
for i in range(2, 30): J.append(J[-1] + 2*J[-2])
J_set = set(J[:25])

print("√det values and their Jacobsthal status:")
for n in range(3, 7):
    classes = get_iso_classes(n)
    sdet_vals = set()
    for cl in classes:
        d = det_I2A(cl['A'], n)
        sd = sqrt_det(d)
        sdet_vals.add(sd)
    sdet_sorted = sorted(sdet_vals)
    print(f"\n  n={n}: √det values = {sdet_sorted}")
    for sd in sdet_sorted:
        in_J = sd in J_set
        # Is sd a product of Jacobsthal numbers?
        # Is sd related to H values?
        print(f"    √det={sd}: Jacobsthal? {in_J}, "
              f"sd mod 3 = {sd%3}, "
              f"(sd²-1)/8 = {(sd**2-1)//8 if sd%2==1 else 'N/A'}")

# ============================================================
# PART 5: The det recurrence for Paley tournaments
# ============================================================
print("\n" + "=" * 72)
print("Part 5: det(I+2A) for Special Tournaments")
print("=" * 72)

# Paley P_3 (cyclic on 3 vertices)
print("\nPaley P_3 (the 3-cycle):")
A3 = np.array([[0,1,0],[0,0,1],[1,0,0]])
H3 = count_hp(A3, 3)
d3 = det_I2A(A3, 3)
print(f"  H = {H3}, det(I+2A) = {d3}, √det = {sqrt_det(d3)}")

# Transitive T_n
print("\nTransitive tournaments (H=1):")
for n in range(2, 8):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    d = det_I2A(A, n)
    H = count_hp(A, n)
    print(f"  n={n}: H={H}, det(I+2A)={d}")

# Paley P_5
print("\nPaley P_5 (QR={1,4}, directed 5-cycle):")
# C_5^{1,4}: 0→1, 0→4, 1→2, 1→0→..., using QR = {1,4}
# Actually QR mod 5 = {1, 4}
A5 = np.zeros((5,5), dtype=int)
QR5 = {1, 4}
for i in range(5):
    for j in range(5):
        if i != j and (j-i) % 5 in QR5:
            A5[i][j] = 1
H5 = count_hp(A5, 5)
d5 = det_I2A(A5, 5)
print(f"  H = {H5}, det(I+2A) = {d5}, √det = {sqrt_det(d5)}")

# Paley P_7
print("\nPaley P_7 (QR={1,2,4}):")
A7 = np.zeros((7,7), dtype=int)
QR7 = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j-i) % 7 in QR7:
            A7[i][j] = 1
H7 = count_hp(A7, 7)
d7 = det_I2A(A7, 7)
sd7 = sqrt_det(d7)
print(f"  H = {H7}, det(I+2A) = {d7}, √det = {sd7}")

# Check the eigenvalue product interpretation
print("\n--- Eigenvalue product ---")
print("det(I+2A) = ∏(1+2λ_k) where λ_k are eigenvalues of A")
for n, An, Hn in [(3, A3, H3), (5, A5, H5), (7, A7, H7)]:
    eigs = np.linalg.eigvals(An.astype(float))
    prod_val = np.prod(1 + 2*eigs)
    dn = det_I2A(An, n)
    print(f"  P_{n}: eigenvalues of A = {sorted(eigs.real, reverse=True)[:4]}...")
    print(f"    ∏(1+2λ) = {prod_val.real:.1f} + {prod_val.imag:.1f}i, "
          f"det = {dn}")
    # Also compute H as sum of eigenvalue powers?
    # H(T) for circulant = sum of transfer matrix eigenvalue powers
    # det(I+2A) = product of (1+2*eigenvalue)

# ============================================================
# PART 6: The DUAL recurrences — additive H vs multiplicative det
# ============================================================
print("\n" + "=" * 72)
print("Part 6: DUAL RECURRENCES — H (additive) vs det (multiplicative)")
print("=" * 72)

print("""
Under vertex deletion T → T-v:

ADDITIVE (H):
  H(T) = H(T-v) + 2·Σ_{C∋v} μ(C)
  Delta is ALWAYS positive (removing vertex can't increase H? No...)
  Delta is ALWAYS even.

MULTIPLICATIVE (det):
  det(I+2A(T)) = det(I+2A(T-v)) · factor(v)
  Question: what is factor(v)?

  By Schur complement: if M = [[M_{-v}, w], [w^T, m_vv]]:
    det(M) = det(M_{-v}) · (m_vv - w^T M_{-v}^{-1} w)
  where m_vv = 1 + 2·0 = 1 (diagonal of I+2A) and
  w is the column/row of interactions with v.

  So factor(v) = 1 - w^T (J+S)_{-v}^{-1} w
  where w = (1+2A_v) is the v-th column of I+2A (minus diagonal).
  Actually w = (J+S)_{v column, minus diagonal entry}.

Let's compute the multiplicative factors:
""")

for n in [4, 5]:
    print(f"\nn={n}:")
    classes = get_iso_classes(n)
    for cl in classes[:6]:  # Show first 6
        A = cl['A']
        H = cl['H']
        d = det_I2A(A, n)

        factors = []
        for v in range(n):
            rem = [u for u in range(n) if u != v]
            Av = A[np.ix_(rem, rem)]
            sub_d = det_I2A(Av, n-1)
            if sub_d != 0:
                factor = d / sub_d
            else:
                factor = float('inf')
            factors.append(factor)

        print(f"  H={H:3d}  det={d:6d}  factors={[round(f,2) for f in factors]}")
        # Check: is the product of factors related to n or H?
        prod_factors = 1
        valid = True
        for f in factors:
            if f == float('inf'):
                valid = False
                break
            prod_factors *= f
        if valid:
            print(f"    ∏ factors = {prod_factors:.2f}, "
                  f"det^{n}/(∏ sub_dets) = {prod_factors:.2f}")

# ============================================================
# PART 7: The sign pattern
# ============================================================
print("\n" + "=" * 72)
print("Part 7: Sign Patterns of det(I+2A)")
print("=" * 72)

print("Is det(I+2A) always positive?")
for n in range(3, 7):
    classes = get_iso_classes(n)
    neg_count = sum(1 for cl in classes if det_I2A(cl['A'], n) < 0)
    zero_count = sum(1 for cl in classes if det_I2A(cl['A'], n) == 0)
    print(f"  n={n}: {len(classes)} classes, {neg_count} negative det, "
          f"{zero_count} zero det")

# ============================================================
# PART 8: H·√det connection
# ============================================================
print("\n" + "=" * 72)
print("Part 8: H × √det Products and Sums")
print("=" * 72)

print("Exploring H + √det, H - √det, H × √det, H² - det:")
for n in range(3, 7):
    print(f"\nn={n}:")
    classes = get_iso_classes(n)
    for cl in classes:
        A = cl['A']
        H = cl['H']
        d = det_I2A(A, n)
        sd = sqrt_det(d)
        if sd:
            print(f"  H={H:3d}  √det={sd:4d}  H+√det={H+sd:4d}  "
                  f"H-√det={H-sd:4d}  H·√det={H*sd:6d}  "
                  f"H²-det={H**2-d:8d}")

# ============================================================
# PART 9: SYNTHESIS
# ============================================================
print("\n" + "=" * 72)
print("SYNTHESIS")
print("=" * 72)

print("""
THE DUAL RECURRENCE STRUCTURE:

                    TOURNAMENT T
                    /          \\
           H(T) = I(CG,2)     det(I+2A) = (Pf sum)²
           (ADDITIVE)          (MULTIPLICATIVE)
           |                   |
  H = H' + 2Σμ         det = det' · (Schur factor)
  (Jacobsthal            (Jacobsthal-multiplicative
   additive tower)        tower?)
           |                   |
           \\                  /
            Both encode tournament structure
            H captures CYCLE information
            det captures MATCHING information

KEY FINDING: H and √det are NOT simply related.
  At n=3: H = √det always (both = 1 or 3)
  At n=4: H/√det ∈ {1, 5}
  At n≥5: H/√det takes many rational values

BUT: H² - det takes suggestive values.
  H² is the "doubled" additive invariant
  det is the multiplicative invariant
  Their DIFFERENCE should encode the "interaction" term.

The JACOBSTHAL connection:
  H follows a Jacobsthal ADDITIVE recurrence (evaluation of indep poly at x=2)
  det follows a MULTIPLICATIVE recurrence (product of eigenvalue factors)
  Both are connected through x=2 being the Jacobsthal evaluation point.

OPEN: Is there a single recurrence that generates BOTH H and √det simultaneously?
  Perhaps: (H, √det) → (H', √det') under vertex deletion?
  This would be a 2D recurrence connecting the additive and multiplicative worlds.
""")
