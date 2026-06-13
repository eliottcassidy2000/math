#!/usr/bin/env python3
"""
cartan_commutator_bch_s95.py — opus-2026-03-21-S95

CARTAN COMMUTATOR AND BCH ANALYSIS

For the flip-chain transition matrix T, decomposed as T = A + S_tl + σI:
  [A, S_tl] = A·S_tl - S_tl·A

This commutator measures how "entangled" the tournament and cooperation
modes are. If [A,S] = 0, the sectors are independent.

Key questions:
1. Does ||[A,S]||_F have a closed form?
2. How does it depend on n?
3. What does it tell us about the BCH expansion of exp(-t·T)?
4. For individual tournament adjacency matrices (not flip chain):
   how does the commutator [B, B^2_sym] relate to invariants?

Also: exact BCH analysis at n=3,4.

Author: opus-2026-03-21-S95
"""

import sys
import numpy as np
from itertools import combinations, permutations
from fractions import Fraction
from collections import defaultdict
from scipy.linalg import expm, logm
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

def score_sequence(A, n):
    return tuple(sorted(A.sum(axis=1)))

def cartan_decompose(M):
    n = M.shape[0]
    M = M.astype(float)
    anti = (M - M.T) / 2.0
    sym = (M + M.T) / 2.0
    sigma = np.trace(M) / n
    scalar = sigma * np.eye(n)
    sym_tl = sym - scalar
    return anti, sym_tl, sigma

print("=" * 72)
print("  CARTAN COMMUTATOR AND BCH ANALYSIS")
print("  opus-2026-03-21-S95")
print("=" * 72)

# ========================================================================
# PART 1: Commutator for individual tournament adjacency matrices
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: COMMUTATOR [A_anti, A_sym] FOR TOURNAMENT ADJACENCY")
print("=" * 72)
print("""
  For tournament adjacency A (entries 0,1):
    A_anti = B/2 = (A - A^T)/2
    A_sym_tl = (A + A^T)/2 - (n-1)/(2n)*I = (J - I)/2 - (n-1)/(2n)*I
             = J/2 - (2n-1)/(2n)*I   (CONSTANT for all tournaments!)

  So [A_anti, A_sym_tl] = [B/2, (J-I)/2 - (n-1)/(2n)*I]
                        = [B/2, (J-I)/2]   (scalar commutes with everything)
                        = (1/4)[B, J-I]
                        = (1/4)[B, J]       (since [B, I] = 0)
                        = (1/4)(BJ - JB)

  BJ[i,j] = sum_k B[i,k]*1 = sum_k B[i,k] = sum_{k!=i} B[i,k]
  But B[i,k] = +1 if i->k, -1 if k->i. So sum_k B[i,k] = s_i - (n-1-s_i)
  = 2*s_i - (n-1), where s_i = out-degree of vertex i.

  So BJ[i,j] = 2*s_i - (n-1) for ALL j. (Constant row.)
  Similarly JB[i,j] = sum_k B[k,j] = 2*pred(j) - (n-1) = -(2*s_j - (n-1)) = (n-1) - 2*s_j.

  Wait: sum_k B[k,j] = sum_{k!=j} (A[k,j] - A[j,k]) = in_j - out_j = (n-1-s_j) - s_j
  = (n-1) - 2*s_j.

  So JB[i,j] = (n-1) - 2*s_j.

  Therefore:
    [B, J][i,j] = BJ[i,j] - JB[i,j] = (2*s_i - (n-1)) - ((n-1) - 2*s_j)
                = 2*s_i + 2*s_j - 2*(n-1)
                = 2*(s_i + s_j - (n-1))

  And [A_anti, A_sym_tl] = (1/4) * [B, J] has entries:
    [A_anti, A_sym_tl][i,j] = (s_i + s_j - (n-1)) / 2

  THIS IS A RANK-2 MATRIX! It's (s + c)/2 where s = score vector and c is a constant.
  Specifically: [A_anti, A_sym_tl] = (1/2)(s*1^T + 1*s^T - (n-1)*J)
  where s = (s_1, ..., s_n) is the score vector.

  Wait, let me be more careful:
  [A_anti, A_sym_tl][i,j] = (s_i + s_j - (n-1)) / 2

  Define d_i = s_i - (n-1)/2 (the deviation from the mean score).
  Then [A_anti, A_sym_tl][i,j] = (d_i + d_j) / ... hmm

  s_i + s_j - (n-1) = (s_i - (n-1)/2) + (s_j - (n-1)/2) = d_i + d_j

  So [A_anti, A_sym_tl][i,j] = (d_i + d_j) / 2

  This is a RANK-2 MATRIX (it's d*1^T + 1*d^T, divided by 2).
  Its Frobenius norm squared:
  ||[A_anti, A_sym_tl]||_F^2 = (1/4) * sum_{i,j} (d_i + d_j)^2
  = (1/4) * (n * sum d_i^2 + 2*(sum d_i)^2 + n * sum d_j^2)
  = (1/4) * (2n * Var + 2*(sum d_i)^2 ...

  Actually sum d_i = 0 (since sum s_i = C(n,2) and sum (n-1)/2 = n(n-1)/2 = C(n,2)).

  So sum_{i,j} (d_i + d_j)^2 = sum_{i,j} d_i^2 + 2*sum_{i,j} d_i*d_j + sum_{i,j} d_j^2
  = n * sum d_i^2 + 2*(sum d_i)*(sum d_j) + n * sum d_j^2
  = 2n * sum d_i^2 + 0 = 2n * sigma_s^2 * n = 2n * S2

  where S2 = sum d_i^2 = Landau's second-order statistic.

  So ||[A_anti, A_sym_tl]||_F^2 = (1/4) * 2n * S2 = n*S2/2

  where S2 = sum_i (s_i - (n-1)/2)^2 = sum s_i^2 - n*(n-1)^2/4.

  Let me verify this.
""")

for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        anti, sym_tl, sigma = cartan_decompose(A)
        comm = anti @ sym_tl - sym_tl @ anti
        comm_norm2 = np.sum(comm**2)

        # Predicted: n * S2 / 2 where S2 = sum (s_i - (n-1)/2)^2
        scores = A.sum(axis=1)
        S2 = np.sum((scores - (n-1)/2)**2)
        predicted = n * S2 / 2

        results.append({
            'H': H, 'c3': c3, 'sc': sc,
            'comm_norm2': comm_norm2,
            'predicted': predicted,
            'S2': S2,
            'match': abs(comm_norm2 - predicted) < 1e-10
        })

    all_match = all(r['match'] for r in results)
    print(f"  ||[A_anti, A_sym_tl]||² = n·S₂/2 verified: "
          f"{'ALL MATCH ✓' if all_match else 'SOME FAIL ✗'} "
          f"({len(results)} tournaments)")

    if n <= 5:
        groups = defaultdict(list)
        for r in results:
            groups[(r['H'], r['c3'], r['sc'])].append(r)

        print(f"  {'H':>4s} {'c3':>4s} {'score':>20s} {'||[A,S]||²':>12s} "
              f"{'n·S₂/2':>12s} {'S₂':>8s}")
        for key in sorted(groups.keys()):
            H, c3, sc = key
            r = groups[key][0]
            print(f"  {H:4d} {c3:4d} {str(sc):>20s} {r['comm_norm2']:12.4f} "
                  f"{r['predicted']:12.4f} {r['S2']:8.2f}")

# ========================================================================
# PART 2: THE COMMUTATOR THEOREM
# ========================================================================
print("\n\n" + "=" * 72)
print("  THEOREM: COMMUTATOR OF CARTAN SECTORS FOR TOURNAMENT ADJACENCY")
print("=" * 72)
print("""
  THEOREM (Cartan Commutator Formula):

  Let T be a tournament on n vertices with adjacency matrix A and
  score sequence (s_1, ..., s_n). Let d_i = s_i - (n-1)/2.

  Define S₂ = Σᵢ dᵢ² (the Landau irregularity statistic).

  Then in the Cartan decomposition A = A_anti + A_sym_tl + σI:

    [A_anti, A_sym_tl][i,j] = (dᵢ + dⱼ) / 2

    ||[A_anti, A_sym_tl]||²_F = n · S₂ / 2

  Proof:
    A_anti = B/2 where B = A - A^T (skew-adjacency).
    A_sym_tl = (J-I)/2 - σI (constant for all tournaments).
    [B/2, (J-I)/2 - σI] = (1/4)[B, J] since [B, I] = [B, σI] = 0.
    BJ has (i,j)-entry = Σ_k B[i,k] = 2sᵢ - (n-1).
    JB has (i,j)-entry = Σ_k B[k,j] = (n-1) - 2sⱼ.
    [B,J][i,j] = 2(sᵢ + sⱼ - (n-1)) = 2(dᵢ + dⱼ).
    [A_anti, A_sym_tl][i,j] = (dᵢ + dⱼ)/2. ✓

    ||·||² = (1/4)Σᵢⱼ(dᵢ+dⱼ)² = (1/4)(2n·Σdᵢ² + 2(Σdᵢ)²)
    Since Σdᵢ = 0: = (1/4)(2n·S₂) = n·S₂/2. ∎

  COROLLARIES:

  1. [A_anti, A_sym_tl] = 0 iff S₂ = 0 iff the tournament is REGULAR
     (all scores equal = (n-1)/2, only possible at odd n).
     THE CARTAN SECTORS COMMUTE IFF THE TOURNAMENT IS REGULAR.

  2. The commutator is MAXIMIZED for the transitive tournament where
     S₂ = n(n²-1)/12 (the maximum of Landau's statistic).

  3. The commutator has rank at most 2 (it's d·1^T + 1·d^T, scaled).
     Rank = 0 iff regular; rank = 1 iff d is proportional to 1 (impossible
     since Σdᵢ=0, so impossible if d≠0); rank = 2 generically.

  4. For EVEN n, no tournament is regular, so the commutator is ALWAYS
     nonzero: the Cartan sectors are always entangled.

  5. For ODD n, Paley tournaments T_p are regular (all scores = (p-1)/2),
     so their Cartan sectors COMMUTE. This is consistent with Paley being
     "maximally symmetric."
""")

# Verify Corollary 5: Paley tournament at n=7
print("  Verification — Paley T_7:")
n = 7
QR = set()
for a in range(1, n):
    QR.add((a*a) % n)
print(f"  Quadratic residues mod 7: {sorted(QR)}")

A_paley = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j and ((j - i) % n) in QR:
            A_paley[i][j] = 1

scores = A_paley.sum(axis=1)
S2 = np.sum((scores - (n-1)/2)**2)
print(f"  Scores: {scores}")
print(f"  S₂ = {S2}")
print(f"  Regular: {S2 == 0}")

anti, sym_tl, sigma = cartan_decompose(A_paley)
comm = anti @ sym_tl - sym_tl @ anti
print(f"  ||[A_anti, A_sym_tl]||² = {np.sum(comm**2):.10f}")
print(f"  Sectors commute: {np.sum(comm**2) < 1e-20}")

# ========================================================================
# PART 3: BCH ANALYSIS — exact expansion at n=3
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: BCH EXPANSION AT n=3")
print("=" * 72)
print("""
  For the cyclic tournament C₃ (1->2->3->1), B is:
    B = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]

  This is the generator of rotations around (1,1,1)/√3 in R³.
  eigenvalues: 0, ±i√3.

  exp(tB) = I + sin(√3·t)/√3 · B + (1-cos(√3·t))/3 · B²

  (Rodrigues' rotation formula for so(3).)

  For the transitive tournament T₃ (1->2, 1->3, 2->3):
    B = [[0, 1, 1], [-1, 0, 1], [-1, -1, 0]]
  eigenvalues: 0, ±i√3 (SAME eigenvalues!)

  In fact, ALL tournaments on 3 vertices have the SAME eigenvalue spectrum
  for B: {0, ±i√3}. This is because ||B||² = n(n-1) = 6, and for n=3
  with eigenvalues 0, ±iθ: 2θ² = 6, so θ = √3.
""")

# Cyclic tournament C_3
A_cyc = np.array([[0,1,0],[0,0,1],[1,0,0]])
B_cyc = (A_cyc - A_cyc.T).astype(float)
print("  Cyclic C₃: B =")
for row in B_cyc.astype(int):
    print(f"    {row}")

# Transitive tournament T_3
A_trans = np.array([[0,1,1],[0,0,1],[0,0,0]])
B_trans = (A_trans - A_trans.T).astype(float)
print("\n  Transitive T₃: B =")
for row in B_trans.astype(int):
    print(f"    {row}")

# Check eigenvalues
eigs_cyc = np.linalg.eigvals(B_cyc)
eigs_trans = np.linalg.eigvals(B_trans)
print(f"\n  C₃ eigenvalues: {sorted(eigs_cyc.imag, reverse=True)}")
print(f"  T₃ eigenvalues: {sorted(eigs_trans.imag, reverse=True)}")

# Verify Rodrigues for C₃
for t in [0.1, 0.5, 1.0, np.log(2), np.log(3)]:
    E_exact = expm(t * B_cyc)
    sqrt3 = np.sqrt(3)
    B2_cyc = B_cyc @ B_cyc
    E_rodrigues = (np.eye(3) + np.sin(sqrt3*t)/sqrt3 * B_cyc
                   + (1 - np.cos(sqrt3*t))/3 * B2_cyc)
    err = np.max(np.abs(E_exact - E_rodrigues))
    print(f"  t={t:.4f}: Rodrigues error = {err:.2e}")

# ========================================================================
# PART 4: BCH for the FULL adjacency matrix A = B/2 + (J-I)/2
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: BCH FOR FULL ADJACENCY exp(A) DECOMPOSITION")
print("=" * 72)
print("""
  A = B/2 + (J-I)/2  where B is skew and (J-I)/2 is symmetric.

  BCH: exp(X+Y) = exp(X) exp(Y) exp(-[X,Y]/2) exp(higher order...)

  Let X = t·B/2 (tournament) and Y = t·(J-I)/2 (cooperation).

  [X, Y] = (t²/4)[B, J-I] = (t²/4)[B, J]  (since [B,I]=0)

  From PART 2: [B, J][i,j] = 2(sᵢ + sⱼ - (n-1)).

  So [X, Y][i,j] = (t²/2)(sᵢ + sⱼ - (n-1)) = t² · (dᵢ + dⱼ).

  For REGULAR tournaments (dᵢ = 0 ∀i): [X, Y] = 0.
  Then exp(A) = exp(X) · exp(Y) EXACTLY.
  The tournament and cooperation modes are COMPLETELY DECOUPLED.

  For non-regular tournaments: [X, Y] ≠ 0 and higher BCH terms matter.

  Let's verify by comparing exp(tA) with exp(tX)·exp(tY) at n=3,5.
""")

for n, name in [(3, "cyclic C₃"), (3, "transitive T₃"), (5, "regular (Paley-like)")]:
    if n == 3 and "cyclic" in name:
        A = np.array([[0,1,0],[0,0,1],[1,0,0]])
    elif n == 3 and "transitive" in name:
        A = np.array([[0,1,1],[0,0,1],[0,0,0]])
    else:
        # n=5 regular: Paley T_5
        n = 5
        QR5 = set()
        for a in range(1, 5):
            QR5.add((a*a) % 5)
        A = np.zeros((5, 5), dtype=int)
        for i in range(5):
            for j in range(5):
                if i != j and ((j - i) % 5) in QR5:
                    A[i][j] = 1
        name = f"Paley T₅"

    Af = A.astype(float)
    B = (Af - Af.T) / 2  # tournament part (= A_anti)
    S = (Af + Af.T) / 2  # symmetric part

    print(f"\n  {name} (n={n}):")
    scores = A.sum(axis=1)
    d = scores - (n-1)/2
    S2 = np.sum(d**2)
    print(f"    Scores: {scores}, S₂={S2:.1f}, Regular: {S2 < 1e-10}")

    for t in [0.1, 0.5, 1.0]:
        E_full = expm(t * Af)
        E_bch0 = expm(t * B) @ expm(t * S)  # zero-order BCH
        E_bch1 = E_bch0 @ expm(-t*t/2 * (B @ S - S @ B))  # first correction

        err0 = np.max(np.abs(E_full - E_bch0))
        err1 = np.max(np.abs(E_full - E_bch1))

        print(f"    t={t:.1f}: BCH₀ err={err0:.2e}, BCH₁ err={err1:.2e}"
              + (" ← EXACT" if err0 < 1e-14 else ""))

# ========================================================================
# PART 5: The Cartan entanglement ratio κ
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: CARTAN ENTANGLEMENT RATIO κ = ||[A,S]||/(||A||·||S||)")
print("=" * 72)
print("""
  For tournament adjacency:
    ||[A_anti, A_sym_tl]||² = n·S₂/2
    ||A_anti||² = n(n-1)/4
    ||A_sym_tl||² = ... let's compute this.

  A_sym_tl = (J-I)/2 - σ·I where σ = (n-1)/(2n).
  A_sym_tl[i,j] = 1/2 for i≠j, and (1/2 - σ) = 1/(2n) for i=j.

  Wait: A_sym = (A+A^T)/2 = (J-I)/2. So A_sym[i,j] = 1/2 for i≠j, 0 diagonal.
  σ = (n-1)/(2n).
  A_sym_tl = A_sym - σI.
  A_sym_tl[i,j] = 1/2 for i≠j, -(n-1)/(2n) for i=j.
  Wait: A_sym_tl[i,i] = 0 - σ = -(n-1)/(2n).
  A_sym_tl[i,j] = 1/2 for i≠j.

  ||A_sym_tl||² = n*(n-1)²/(4n²) + n(n-1)/4 = (n-1)²/(4n) + n(n-1)/4
                = (n-1)/(4) * ((n-1)/n + n) = (n-1)/(4) * (n²+n-1)/n
                Hmm, let me just compute directly.
""")

for n in [3, 4, 5, 6, 7, 8, 9, 10]:
    sigma = (n-1) / (2*n)

    # ||A_sym_tl||²
    # diagonal: n * sigma² = n * (n-1)² / (4n²) = (n-1)² / (4n)
    diag_norm2 = n * sigma**2
    # off-diagonal: n*(n-1) * (1/2)² = n*(n-1)/4
    off_norm2 = n * (n - 1) / 4

    # But wait: A_sym_tl[i,i] = -sigma, A_sym_tl[i,j] = 1/2 for i!=j
    norm2_stl = diag_norm2 + off_norm2

    # ||A_anti||² = n(n-1)/4
    norm2_anti = n * (n - 1) / 4

    # S₂ for transitive: s_i = i for i=0,...,n-1. d_i = i - (n-1)/2.
    # S₂ = sum (i-(n-1)/2)² = n(n²-1)/12.
    S2_trans = n * (n*n - 1) / 12
    # S₂ for regular: 0
    S2_reg = 0

    # ||[A,S]||² = n*S₂/2
    comm_trans = n * S2_trans / 2
    comm_reg = 0

    # κ = ||[A,S]||² / (||A||² · ||S||²)
    kappa_trans = comm_trans / (norm2_anti * norm2_stl) if norm2_anti * norm2_stl > 0 else 0

    # max κ is for transitive
    kappa_max = kappa_trans

    print(f"  n={n:2d}: ||A||²={(norm2_anti):.2f}, ||S_tl||²={norm2_stl:.2f}, "
          f"κ_max={kappa_max:.6f}, "
          f"κ_max·simple = {kappa_max * 12:.6f}")

print("""
  EXACT FORMULA for κ:

  κ = ||[A_anti, A_sym_tl]||² / (||A_anti||² · ||A_sym_tl||²)

  Numerator: n·S₂/2
  ||A_anti||² = n(n-1)/4
  ||A_sym_tl||² = n*(n-1)²/(4n²) + n(n-1)/4 = (n-1)/(4n) * (n-1+n²) = (n-1)(n²+n-1)/(4n)

  Wait, let me recompute:
  ||A_sym_tl||² = n * ((n-1)/(2n))² + n(n-1) * (1/2)²
                = (n-1)²/(4n) + n(n-1)/4
                = (n-1)/(4) * ((n-1)/n + n)
                = (n-1)/(4n) * (n-1 + n²)
                = (n-1)(n²+n-1)/(4n)

  For the transitive tournament: S₂ = n(n²-1)/12
  κ_trans = (n · n(n²-1)/12 / 2) / (n(n-1)/4 · (n-1)(n²+n-1)/(4n))
          = n²(n²-1) / 24 / (n(n-1)²(n²+n-1) / (16n))
          = n²(n²-1) / 24 · 16n / (n(n-1)²(n²+n-1))
          = 16n²(n+1)(n-1) / (24 · n · (n-1)² · (n²+n-1))
          = 16n(n+1) / (24(n-1)(n²+n-1))
          = 2n(n+1) / (3(n-1)(n²+n-1))
""")

for n in range(3, 15):
    predicted = 2*n*(n+1) / (3*(n-1)*(n*n+n-1))
    S2_trans = n * (n*n - 1) / 12
    norm2_anti = n * (n - 1) / 4
    sigma = (n-1) / (2*n)
    norm2_stl = n * sigma**2 + n * (n - 1) / 4
    actual = (n * S2_trans / 2) / (norm2_anti * norm2_stl) if norm2_anti * norm2_stl > 0 else 0

    print(f"  n={n:2d}: κ_predicted = {predicted:.8f}, κ_actual = {actual:.8f}, "
          f"match: {'✓' if abs(predicted-actual) < 1e-10 else '✗'}")

# ========================================================================
# PART 6: What determines the Cartan entanglement of A²?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: CARTAN STRUCTURE OF A² — DEEPER ANALYSIS")
print("=" * 72)
print("""
  A² has zero trace and zero scalar. Its Cartan split:
    (A²)_anti = (A² - (A²)^T)/2   [varies by tournament]
    (A²)_sym  = (A² + (A²)^T)/2   [varies by tournament]

  Q1: Is (A²)_sym determined by the score sequence?
  Q2: Is (A²)_anti related to c₃?
  Q3: What's the commutator [(A²)_anti, (A²)_sym]?

  Key relation: A² + (A^T)² = A² + (J-I-A)² ... let me expand.
  A^T = J - I - A (since A + A^T = J - I).
  (A^T)² = (J-I-A)² = (J-I)² - (J-I)A - A(J-I) + A²

  (J-I)² = J² - 2J + I = nJ - 2J + I = (n-2)J + I
  (J-I)A = (nJ-2J+I)... hmm, (J-I)A[i,j] = sum_k (J-I)[i,k]*A[k,j]
         = sum_{k!=i} A[k,j] = (n-1-s_j) - (1-A[i,j]) = (n-2-s_j+A[i,j])
  Wait: sum_{k!=i} A[k,j] = (in-degree of j) - (1 if i->j, 0 otherwise... no)

  sum_{k=0}^{n-1} A[k,j] = in_j = n-1-s_j (where s_j = out-degree).
  So sum_{k!=i} A[k,j] = (n-1-s_j) - A[i,j].
  (J-I)A[i,j] = (n-1-s_j) - A[i,j].

  Similarly: A(J-I)[i,j] = sum_k A[i,k]*(J-I)[k,j] = sum_{k!=j} A[i,k]
  = s_i - A[i,j].

  So (A^T)² = (n-2)J + I - ((n-1-s_j)-A[i,j]) - (s_i-A[i,j]) + A²
            = (n-2)J + I - n + 1 + s_j + A[i,j] - s_i + A[i,j] + A²
  Hmm, this is getting messy with mixed indices. Let me just work matrix-wise.

  (A²)_sym = (A² + (A^T)²)/2

  Key identity: A + A^T = J - I
  So A^T = J - I - A
  (A^T)² = (J-I)² - (J-I)A - A(J-I) + A²
          = (n-2)J + I - (J-I)A - A(J-I) + A²

  (A²)_sym = (A² + (A^T)²)/2 = A² + ((n-2)J + I)/2 - ((J-I)A + A(J-I))/2
  Hmm, let me just verify numerically what (A²)_sym depends on.
""")

for n in [5]:
    print(f"\n--- n = {n}: checking what determines (A²)_sym ---")

    by_score = defaultdict(list)
    for A in all_tournaments(n):
        sc = score_sequence(A, n)
        A2 = A @ A
        A2_sym = (A2 + A2.T) / 2.0
        by_score[sc].append(A2_sym)

    for sc in sorted(by_score.keys()):
        mats = by_score[sc]
        # Check if all are the same
        all_same = all(np.allclose(m, mats[0]) for m in mats[1:])
        if not all_same:
            # How many distinct?
            distinct = []
            for m in mats:
                found = False
                for d in distinct:
                    if np.allclose(m, d):
                        found = True
                        break
                if not found:
                    distinct.append(m)
            print(f"  Score {sc}: {len(distinct)} distinct (A²)_sym matrices "
                  f"(out of {len(mats)} tournaments)")
        else:
            print(f"  Score {sc}: (A²)_sym is UNIQUE ({len(mats)} tournaments)")

# ========================================================================
# PART 7: Exact (A²)_sym formula
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: EXACT FORMULA FOR (A²)_sym AND (A²)_anti")
print("=" * 72)
print("""
  (A²)[i,j] = #{k : i->k->j}   for i≠j
  ((A^T)²)[i,j] = #{k : i<-k<-j} = #{k : j->k->i}   for i≠j

  (A²)_sym[i,j] = (#{k : i->k->j} + #{k : j->k->i}) / 2  for i≠j
  (A²)_anti[i,j] = (#{k : i->k->j} - #{k : j->k->i}) / 2  for i≠j

  For i≠j, define:
    w(i,j) = #{k≠i,j : i->k AND k->j}

  Then A²[i,j] = w(i,j) and (A^T)²[i,j] = w(j,i).

  (A²)_sym[i,j] = (w(i,j) + w(j,i)) / 2
  (A²)_anti[i,j] = (w(i,j) - w(j,i)) / 2

  Now, the total w(i,j) + w(j,i) counts k≠i,j such that:
    (i->k AND k->j) OR (j->k AND k->i)
  i.e., k is a "2-witness" that i and j are connected through k.

  The difference w(i,j) - w(j,i) = #{k: i->k->j} - #{k: j->k->i}
  = the net directional dominance of i over j through 2-step paths.
""")

# Check: is w(i,j) + w(j,i) determined by score sequence?
for n in [5]:
    print(f"\n--- n = {n}: checking w(i,j)+w(j,i) ---")

    by_score_pair = defaultdict(set)
    for A in all_tournaments(n):
        sc = score_sequence(A, n)
        scores = A.sum(axis=1)
        # For each pair (i,j), compute w(i,j)+w(j,i)
        A2 = A @ A
        for i in range(n):
            for j in range(i+1, n):
                w_sum = A2[i,j] + A2[j,i]
                si, sj = int(scores[i]), int(scores[j])
                by_score_pair[(sc, min(si,sj), max(si,sj))].add(int(w_sum))

    constant_count = 0
    variable_count = 0
    for key, vals in sorted(by_score_pair.items()):
        if len(vals) > 1:
            variable_count += 1
        else:
            constant_count += 1

    print(f"  {constant_count} constant, {variable_count} variable "
          f"(w(i,j)+w(j,i) by score pair)")

    # Actually, let's check the exact formula:
    # w(i,j) = #{k≠i,j : i->k AND k->j}
    # For each k, the probability that i->k AND k->j depends on the
    # specific arcs, not just scores.
    # But: w(i,j) = s_i - A[i,j] - #{k≠j : i->k AND j->k}
    # Hmm, this uses the number of common successors...

    # The TOTAL for all k: w(i,j) + w(j,i) = total number of 2-paths
    # from i to j in either direction.

    # Key insight: for k≠i,j, the arc between i and k is independent
    # of the arc between k and j (in a tournament). But they're NOT
    # independent because they share the tournament structure.

    # Let me check a simpler formula:
    # w(i,j) = s_i·(1-A[j,i]/(n-2))... no, this doesn't work.

    # Let's verify: is w(i,j)+w(j,i) = n-2 when i->j?
    print(f"\n  Checking w(i,j)+w(j,i) vs arcs:")
    for A in [np.array([[0,1,1,1,0],[0,0,0,1,1],[0,1,0,0,1],[0,0,1,0,1],[1,0,0,0,0]])]:
        A2 = A @ A
        scores = A.sum(axis=1)
        print(f"  Scores: {scores}")
        for i in range(5):
            for j in range(i+1, 5):
                w_sum = A2[i,j] + A2[j,i]
                arc = "i->j" if A[i,j] else "j->i"
                print(f"    ({i},{j}) {arc}: w+w'={w_sum}, "
                      f"s_i+s_j={int(scores[i]+scores[j])}")

print("\n\nDone. opus-2026-03-21-S95")
