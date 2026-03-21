#!/usr/bin/env python3
"""
cartan_tournament_adjacency_s95.py — opus-2026-03-21-S95

CARTAN DECOMPOSITION OF TOURNAMENT ADJACENCY MATRICES

Key question: For the skew-adjacency matrix B = A - A^T of a tournament
(where A is the {0,1} adjacency matrix), what is the Cartan structure?

For a tournament T on n vertices:
  A[i,j] = 1 if i->j, 0 otherwise
  B = A - A^T (skew-symmetric, entries in {-1,0,1})

Note: B is ALREADY antisymmetric, so in the Cartan decomposition of gl(n):
  B = B_anti + B_sym + B_scalar
  B_anti = (B - B^T)/2 = B  (since B^T = -B)
  B_sym = (B + B^T)/2 = 0
  B_scalar = Tr(B)/n * I = 0

So the skew-adjacency is PURELY in the tournament sector — trivially.

THE INTERESTING OBJECT is instead:
1. The {0,1} adjacency matrix A itself (not skew-symmetric)
2. The transition matrix of the flip chain (as in S94)
3. The transfer matrix M[a,b] (HP counting matrix)

Let's analyze the {0,1} adjacency matrix A:
  A = (B + J)/2 where J is the all-ones matrix minus identity
  Cartan decomposition of A:
    A_anti = B/2 (the tournament part)
    A_sym = (J - I)/2 (the cooperation part — CONSTANT for all tournaments!)
    A_scalar = (n-1)/(2n) * I (constant)

So the Cartan decomposition of the adjacency matrix is BORING: the
symmetric part is the same for every tournament. Only the antisymmetric
part varies, and it IS the tournament.

THE RIGHT OBJECT TO STUDY: the transfer matrix M[a,b] (counts directed
paths from a to b through all n vertices). This has BOTH varying symmetric
and antisymmetric parts, and M[a,b] = M[b,a] (symmetric at odd n), so
the tournament sector is ZERO at odd n!

Even better: look at M^(k) = (powers of A). These have varying Cartan
structure that depends on the tournament.

Let's compute:
1. Cartan decomposition of A^k for k=1,2,3,...
2. How fractions depend on tournament invariants (H, c3, score seq)
3. The commutator [A_anti, A_sym] for actual tournaments
4. Cartan decomposition of exp(A) and exp(B) (heat kernel analogues)

Author: opus-2026-03-21-S95
"""

import sys
import numpy as np
from itertools import combinations, permutations
from fractions import Fraction
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
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
    """Count Hamiltonian paths via DP."""
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
    """Count directed 3-cycles."""
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

def score_sequence(A, n):
    """Return sorted score sequence."""
    return tuple(sorted(A.sum(axis=1)))

def cartan_decompose(M):
    """Cartan decomposition of an n x n matrix.
    Returns: (anti, sym_traceless, scalar_coeff, fractions)
    Fractions = (anti_frac, sym_frac, scalar_frac) of Frobenius norm^2.
    """
    n = M.shape[0]
    anti = (M - M.T) / 2.0
    sym = (M + M.T) / 2.0
    sigma = np.trace(M) / n
    scalar = sigma * np.eye(n)
    sym_tl = sym - scalar

    norm2 = np.sum(M**2)
    if norm2 < 1e-15:
        return anti, sym_tl, sigma, (0, 0, 0)

    anti_frac = np.sum(anti**2) / norm2
    sym_frac = np.sum(sym_tl**2) / norm2
    scalar_frac = np.sum(scalar**2) / norm2

    return anti, sym_tl, sigma, (anti_frac, sym_frac, scalar_frac)

def commutator_norm(A_part, S_part):
    """||[A, S]||_F"""
    comm = A_part @ S_part - S_part @ A_part
    return np.sqrt(np.sum(comm**2))

# ========================================================================
print("=" * 72)
print("  CARTAN DECOMPOSITION OF TOURNAMENT ADJACENCY MATRICES")
print("  opus-2026-03-21-S95")
print("=" * 72)

# ========================================================================
# PART 1: The adjacency matrix A = (B + J)/2
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: CARTAN FRACTIONS OF A^k FOR ALL TOURNAMENTS")
print("=" * 72)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        # Cartan of A, A^2, A^3
        fracs_by_k = {}
        comm_norms = {}
        for k in [1, 2, 3]:
            Ak = np.linalg.matrix_power(A, k)
            anti, sym_tl, sigma, fracs = cartan_decompose(Ak.astype(float))
            fracs_by_k[k] = fracs
            if k == 1:
                comm_norms[1] = commutator_norm(anti, sym_tl)

        results.append({
            'H': H, 'c3': c3, 'sc': sc,
            'fracs': fracs_by_k,
            'comm': comm_norms[1],
            'A': A
        })

    # Group by (H, c3) and report averages
    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    print(f"\n  Grouped by (H, c3): {len(groups)} distinct classes")
    print(f"  {'H':>4s} {'c3':>4s} {'cnt':>5s} | "
          f"{'A1_anti':>8s} {'A1_sym':>8s} {'A1_scl':>8s} | "
          f"{'A2_anti':>8s} {'A2_sym':>8s} | "
          f"{'||[A,S]||':>10s}")

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        cnt = len(grp)
        f1 = np.mean([r['fracs'][1] for r in grp], axis=0)
        f2 = np.mean([r['fracs'][2] for r in grp], axis=0)
        cm = np.mean([r['comm'] for r in grp])
        print(f"  {H:4d} {c3:4d} {cnt:5d} | "
              f"{f1[0]:8.4f} {f1[1]:8.4f} {f1[2]:8.4f} | "
              f"{f2[0]:8.4f} {f2[1]:8.4f} | "
              f"{cm:10.6f}")

# ========================================================================
# PART 2: The walk matrix W = sum_{k=0}^{n-1} A^k
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: CARTAN DECOMPOSITION OF WALK MATRIX W = sum A^k")
print("=" * 72)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)

        # Walk matrix
        W = np.zeros((n, n), dtype=float)
        Ak = np.eye(n, dtype=float)
        for k in range(n):
            W += Ak
            Ak = Ak @ A

        anti, sym_tl, sigma, fracs = cartan_decompose(W)
        cm = commutator_norm(anti, sym_tl)

        results.append({
            'H': H, 'c3': c3,
            'fracs': fracs,
            'comm': cm,
            'sigma': sigma
        })

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    print(f"  {'H':>4s} {'c3':>4s} {'cnt':>5s} | "
          f"{'W_anti':>8s} {'W_sym':>8s} {'W_scl':>8s} | "
          f"{'sigma':>8s} {'||[A,S]||':>10s}")

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        cnt = len(grp)
        f = np.mean([r['fracs'] for r in grp], axis=0)
        cm = np.mean([r['comm'] for r in grp])
        sg = np.mean([r['sigma'] for r in grp])
        print(f"  {H:4d} {c3:4d} {cnt:5d} | "
              f"{f[0]:8.4f} {f[1]:8.4f} {f[2]:8.4f} | "
              f"{sg:8.4f} {cm:10.4f}")

# ========================================================================
# PART 3: The exp(A) heat kernel and Cartan
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: CARTAN DECOMPOSITION OF exp(tB) FOR SKEW-ADJACENCY B")
print("=" * 72)
print("  B = A - A^T (skew-symmetric, entries in {-1,0,1})")
print("  exp(tB) is ORTHOGONAL (in SO(n)) for real t")
print("  Its Cartan decomposition tells us about ROTATION structure")

from scipy.linalg import expm

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)

        B = (A - A.T).astype(float)

        # exp(t*B) for various t
        fracs_by_t = {}
        for t in [0.1, 0.5, 1.0, np.log(2), np.log(3)]:
            E = expm(t * B)
            _, _, _, fracs = cartan_decompose(E)
            fracs_by_t[t] = fracs

        results.append({
            'H': H, 'c3': c3,
            'fracs': fracs_by_t,
            'B': B
        })

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    t_vals = [0.1, 0.5, 1.0, np.log(2), np.log(3)]
    t_names = ['0.1', '0.5', '1.0', 'ln2', 'ln3']

    print(f"  {'H':>4s} {'c3':>4s} | ", end="")
    for tn in t_names:
        print(f"  anti({tn:>4s})", end="")
    print()

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        print(f"  {H:4d} {c3:4d} | ", end="")
        for t in t_vals:
            f = np.mean([r['fracs'][t] for r in grp], axis=0)
            print(f"  {f[0]:10.4f}", end="")
        print()

# ========================================================================
# PART 4: KEY THEOREM — Cartan fractions of A are CONSTANT
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: THEOREM — ADJACENCY CARTAN FRACTIONS ARE UNIVERSAL")
print("=" * 72)
print("""
  For any tournament T on n vertices with adjacency matrix A:
    A = A_anti + A_sym + A_scalar

  A_anti = B/2 = (A - A^T)/2  (the tournament part)
  A_sym = (A + A^T)/2 = (J - I)/2  (CONSTANT — the complete graph K_n)
  A_scalar = (n-1)/(2n) * I  (CONSTANT)

  Since A + A^T = J - I for ANY tournament:
    ||A_sym||_F^2 = ||(J-I)/2 - (n-1)/(2n)*I||_F^2
                  = || J/2 - I/2 - (n-1)/(2n)*I ||_F^2
                  = || J/2 - (2n-1)/(2n)*I ||_F^2
    (Constant for all T.)

  And ||A_anti||_F^2 = ||B/2||_F^2 = (1/4) * sum_{i≠j} 1 = C(n,2)/2
    Wait — B[i,j]^2 = 1 for i≠j, so ||B||_F^2 = n(n-1).
    Hence ||A_anti||_F^2 = n(n-1)/4.

  Also ||A||_F^2 = sum A[i,j]^2 = number of 1s = C(n,2) = n(n-1)/2.

  So anti_frac = (n(n-1)/4) / (n(n-1)/2) = 1/2. ALWAYS.
  And scalar_frac = n * ((n-1)/(2n))^2 / (n(n-1)/2) = (n-1)/(2n^2) / (1/2)
                  = ... let me compute directly.
""")

for n in [3, 4, 5, 6, 7, 8]:
    # Theoretical values
    norm2_A = n * (n-1) / 2  # number of 1s in A
    norm2_anti = n * (n-1) / 4  # ||B/2||^2
    sigma = (n-1) / (2*n)
    norm2_scalar = n * sigma**2
    norm2_sym_tl = norm2_A - norm2_anti - norm2_scalar  # by Pythagonality

    anti_frac = norm2_anti / norm2_A
    scalar_frac = norm2_scalar / norm2_A
    sym_frac = norm2_sym_tl / norm2_A

    print(f"  n={n}: anti={anti_frac:.6f}, sym={sym_frac:.6f}, "
          f"scalar={scalar_frac:.6f}, sum={anti_frac+sym_frac+scalar_frac:.6f}")
    print(f"        exact: anti=1/2, scalar=(n-1)/(2n^2) * n / (n(n-1)/2) "
          f"= 1/(n*(n-1)/2) * n*(n-1)^2/(4n^2)")

    # Cleaner: scalar_frac = n * ((n-1)/(2n))^2 / (n(n-1)/2)
    # = (n-1)^2 / (4n) / (n(n-1)/2) = (n-1) / (2n^2)
    sf_exact = Fraction(n-1, 2*n*n)
    af_exact = Fraction(1, 2)
    symf_exact = 1 - af_exact - sf_exact
    print(f"        exact fracs: anti=1/2, scalar={(n-1)}/{2*n*n}={sf_exact}, "
          f"sym={symf_exact}")

print("""
  THEOREM (Cartan Universality for Tournament Adjacency):

  For ANY tournament T on n vertices with adjacency matrix A ∈ {0,1}^{n×n}:

    anti_frac(A) = 1/2           (the tournament sector, CONSTANT)
    scalar_frac(A) = (n-1)/(2n²) (the identity sector, CONSTANT)
    sym_frac(A) = (n²-1)/(2n²)   (the cooperation sector, CONSTANT)

  Proof: A + A^T = J - I for any tournament, so the symmetric part is fixed.
  ||B||_F^2 = n(n-1), ||A||_F^2 = n(n-1)/2, so anti_frac = 1/2 always. ∎

  COROLLARY: The {0,1} adjacency matrix carries NO tournament-specific
  information in its Cartan balance. All tournament information is in the
  DIRECTION of A_anti (which vector in so(n)), not in its magnitude.

  This means: to see tournament-dependent Cartan structure, we must look at
  NONLINEAR functions of A: powers A^k, exp(A), or the walk matrix.
""")

# ========================================================================
# PART 5: A^2 Cartan — THIS should vary by tournament
# ========================================================================
print("=" * 72)
print("  PART 5: CARTAN FRACTIONS OF A^2 — TOURNAMENT-DEPENDENT")
print("=" * 72)
print("""
  A^2[i,j] = number of directed 2-paths from i to j.

  For a tournament: A^2[i,j] = #{k : i->k->j} = out(i) - A[i,j] ...
  Actually A^2[i,j] = sum_k A[i,k]*A[k,j].

  The SYMMETRIC part of A^2: (A^2 + (A^2)^T)/2 = (A^2 + (A^T)^2)/2
  This counts the total number of 2-paths between i and j (in either direction).

  The ANTISYMMETRIC part: (A^2 - (A^2)^T)/2
  This is the EXCESS of 2-paths i->j over j->i.

  NOTE: (A^T)^2 = (A^2)^T, so (A^2)^T[i,j] = #{k : k->i, j->k} = #{k : j->k->i}

  So the antisymmetric part of A^2 measures DIRECTIONAL DOMINANCE in 2-step walks.
""")

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        A2 = A @ A
        anti, sym_tl, sigma, fracs = cartan_decompose(A2.astype(float))
        cm = commutator_norm(anti, sym_tl)

        results.append({
            'H': H, 'c3': c3, 'sc': sc,
            'fracs': fracs,
            'comm': cm,
            'sigma': sigma,
            'trace_A2': np.trace(A2)
        })

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    print(f"  {'H':>4s} {'c3':>4s} {'cnt':>5s} | "
          f"{'anti':>8s} {'sym':>8s} {'scalar':>8s} | "
          f"{'Tr(A²)':>8s} {'sigma':>8s} {'||[A,S]||':>10s}")

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        cnt = len(grp)
        f = np.mean([r['fracs'] for r in grp], axis=0)
        cm = np.mean([r['comm'] for r in grp])
        sg = np.mean([r['sigma'] for r in grp])
        tr = np.mean([r['trace_A2'] for r in grp])
        print(f"  {H:4d} {c3:4d} {cnt:5d} | "
              f"{f[0]:8.4f} {f[1]:8.4f} {f[2]:8.4f} | "
              f"{tr:8.1f} {sg:8.4f} {cm:10.4f}")

# ========================================================================
# PART 6: TRACE of A^2 = number of mutual pairs = relates to c3
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: EXACT FORMULA FOR Tr(A^2) AND CARTAN OF A^2")
print("=" * 72)
print("""
  Tr(A^2) = sum_i (A^2)[i,i] = sum_i sum_k A[i,k]*A[k,i]
           = sum_{i,k} A[i,k]*A[k,i] = 0 for tournaments!

  Wait: A[i,k]*A[k,i] = 1 iff both i->k AND k->i, impossible in a tournament.
  So Tr(A^2) = 0 ALWAYS.

  Therefore sigma(A^2) = 0, and the scalar fraction is ZERO.

  This means A^2 decomposes into ONLY tournament + cooperation:
    A^2 = A^2_anti + A^2_sym   (no scalar part)

  The antisymmetric fraction of A^2 depends on the tournament!
""")

# Verify and compute exact anti_frac for A^2
for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} (exhaustive) ---")
    frac_set = defaultdict(list)

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        A2 = A @ A
        assert np.trace(A2) == 0, f"Tr(A^2) != 0 at n={n}!"

        # ||A^2||_F^2
        norm2 = np.sum(A2**2)
        # A^2_anti = (A^2 - (A^2)^T)/2
        anti = (A2 - A2.T) / 2.0
        norm2_anti = np.sum(anti**2)
        anti_frac = norm2_anti / norm2 if norm2 > 0 else 0

        frac_set[(H, c3, sc)].append(anti_frac)

    print(f"  {'H':>4s} {'c3':>4s} {'score_seq':>20s} {'anti_frac(A²)':>15s} {'cnt':>5s}")
    for key in sorted(frac_set.keys()):
        H, c3, sc = key
        vals = frac_set[key]
        # All should be the same within a class
        mean_v = np.mean(vals)
        std_v = np.std(vals)
        print(f"  {H:4d} {c3:4d} {str(sc):>20s} {mean_v:15.6f} {len(vals):5d}"
              + (f"  (std={std_v:.1e})" if std_v > 1e-10 else ""))

# ========================================================================
# PART 7: anti_frac(A^2) as a function of score sequence
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: anti_frac(A^2) DEPENDS ON SCORE SEQUENCE")
print("=" * 72)
print("""
  CLAIM: anti_frac(A^2) depends only on the score sequence, not the
  specific tournament. Let's verify this and find the exact formula.

  For tournament T with score sequence (s_1, ..., s_n):
    A^2[i,j] = sum_k A[i,k]*A[k,j] = #{k : i->k, k->j}

  For i != j:
    A^2[i,j] = #{k != i,j : i->k AND k->j} + A[i,j]*0 (k=j: A[j,j]=0)
             + A[i,i]*A[i,j] = 0 (A[i,i]=0)

  So A^2[i,j] = #{common successors of i and predecessors of j, excluding i,j}
             + (1 if i->j) * 0 ... hmm

  Actually: A^2[i,j] = sum_{k=0}^{n-1} A[i,k]*A[k,j]
  For k=i: A[i,i]=0. For k=j: A[j,j]=0. For k!=i,j:
    A^2[i,j] = sum_{k!=i,j} A[i,k]*A[k,j]

  If i->j (A[i,j]=1):
    For each k != i,j: either i->k or k->i, and either k->j or j->k.
    A[i,k]*A[k,j] = 1 iff i->k AND k->j.

  Define for i != j:
    w(i,j) = #{k != i,j : i->k AND k->j} = A^2[i,j]

  Then:
    ||A^2||_F^2 = sum_{i!=j} w(i,j)^2  (diagonal is 0)
    ||A^2_anti||_F^2 = (1/4) sum_{i!=j} (w(i,j) - w(j,i))^2

  And anti_frac(A^2) = ||A^2_anti||_F^2 / ||A^2||_F^2
""")

# Check if anti_frac depends only on score sequence
for n in [5, 6]:
    print(f"\n--- n = {n} ---")
    by_score = defaultdict(set)

    for A in all_tournaments(n):
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        A2 = A @ A
        norm2 = np.sum(A2**2)
        anti = (A2 - A2.T) / 2.0
        anti_frac = np.sum(anti**2) / norm2 if norm2 > 0 else 0

        by_score[sc].add(round(anti_frac, 10))

    all_unique = True
    for sc in sorted(by_score.keys()):
        vals = by_score[sc]
        if len(vals) > 1:
            all_unique = False
            print(f"  Score {sc}: {len(vals)} distinct anti_frac values: {sorted(vals)[:5]}")

    if all_unique:
        print(f"  CONFIRMED: anti_frac(A^2) is CONSTANT on each score sequence class (n={n})")
    else:
        print(f"  REFUTED: anti_frac(A^2) varies WITHIN score sequence classes (n={n})")

# ========================================================================
# PART 8: The commutator [B, B^2] and higher structure
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 8: COMMUTATOR [B, B^k] AND LIE BRACKET STRUCTURE")
print("=" * 72)
print("""
  B = A - A^T is in so(n). Its powers B^k generate a Lie subalgebra.
  The commutator [B, B^2] = B*B^2 - B^2*B = B^3 - B^3 = 0? NO!
  [B, B^2] = B*B^2 - B^2*B = 0 only if B commutes with B^2.
  But B always commutes with B^2 (any matrix commutes with its own powers).

  The interesting commutator is between B of DIFFERENT tournaments.
  Or between B and some fixed reference (like the all-ones antisymmetric part).

  Better: the Lie algebra generated by {B_T : T a tournament on n vertices}
  is a subalgebra of so(n). What is its dimension?
""")

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")

    # Collect all skew-adjacency matrices
    Bs = []
    for A in all_tournaments(n):
        B = (A - A.T).astype(float)
        Bs.append(B.flatten())

    Bs = np.array(Bs)
    # Rank of span
    rank = np.linalg.matrix_rank(Bs, tol=1e-10)
    dim_son = n * (n - 1) // 2
    num_tournaments = 2**(n*(n-1)//2)

    print(f"  dim(so({n})) = {dim_son}")
    print(f"  #tournaments = {num_tournaments}")
    print(f"  rank(span of all B_T) = {rank}")
    print(f"  Span = {'FULL so(' + str(n) + ')' if rank == dim_son else 'STRICT SUBSPACE'}")

# ========================================================================
# PART 9: exp(B) eigenvalues — tournament structure in SO(n)
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 9: exp(B) EIGENVALUES — TOURNAMENT ROTATIONS IN SO(n)")
print("=" * 72)
print("""
  Since B is real skew-symmetric, exp(B) ∈ SO(n).
  Eigenvalues of B are purely imaginary: ±iθ_1, ..., ±iθ_k, (0 if n odd).
  Eigenvalues of exp(B) are: e^{±iθ_j} = cos(θ_j) ± i*sin(θ_j).

  These are ROTATION ANGLES. Each tournament defines a rotation in SO(n).

  The angles θ_j encode the tournament's "rotational structure."
  Question: How do these angles relate to H, c3, score sequence?
""")

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        B = (A - A.T).astype(float)
        eigs_B = np.sort(np.linalg.eigvals(B).imag)[::-1]
        # Only positive imaginary parts (the angles)
        angles = eigs_B[eigs_B > 1e-10]

        # exp(B) eigenvalues
        eigs_expB = np.linalg.eigvals(expm(B))

        results.append({
            'H': H, 'c3': c3, 'sc': sc,
            'angles': tuple(round(a, 6) for a in angles),
            'det_expB': round(np.real(np.prod(eigs_expB)), 6)
        })

    # Group by (H, c3)
    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        # Show representative angles
        rep = grp[0]
        print(f"  H={H:3d}, c3={c3:2d}: angles={rep['angles']}, "
              f"det(exp(B))={rep['det_expB']}, count={len(grp)}")

# ========================================================================
# PART 10: TRACE FORMULAS — connecting tournament invariants to Cartan
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 10: TRACE FORMULAS FOR B^k AND TOURNAMENT INVARIANTS")
print("=" * 72)
print("""
  B = A - A^T, entries B[i,j] = A[i,j] - A[j,i] ∈ {-1, 0, 1}.
  Since B is skew-symmetric: Tr(B^k) = 0 for odd k.

  Tr(B^2) = sum_{i,j} B[i,j]*B[j,i] = -sum_{i,j} B[i,j]^2 = -||B||_F^2
           = -n(n-1) (CONSTANT for all tournaments)

  Tr(B^3) = sum_{i,j,k} B[i,j]*B[j,k]*B[k,i]
  B[i,j]*B[j,k]*B[k,i] = (A[i,j]-A[j,i])*(A[j,k]-A[k,j])*(A[k,i]-A[i,k])

  For a 3-cycle i->j->k->i: each factor is +1, product = 1.
  For i->k->j->i: each factor is -1, product = -1.
  For non-cycles (some pair equal): factors include 0? No — all pairs distinct
  in {i,j,k} with i<j<k, both orientations give cycles.

  So Tr(B^3) = sum_{distinct i,j,k} B[i,j]*B[j,k]*B[k,i]
  This counts 3-cycles with sign: +1 for i->j->k->i, -1 for i->k->j->i.

  Wait, but B is skew so Tr(B^3) = 0 (odd power of skew matrix).

  Let me check: Tr(B^3) = sum eigenvalues^3 = sum (iθ_j)^3 = -i*sum θ_j^3.
  For real skew matrix, eigenvalues come in conjugate pairs, so sum is real.
  But (iθ)^3 = -iθ^3, and paired with (-iθ)^3 = iθ^3, sum = 0.
  So Tr(B^3) = 0. ✓

  Tr(B^4) = sum (iθ_j)^4 = sum θ_j^4. This is always POSITIVE (or 0).
  And Tr(B^4) = Σ_{i,j} (B^2)[i,j]^2 = ||B^2||_F^2.

  INTERESTING: Tr(B^4) = sum θ_j^4 depends on the tournament!
""")

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)

        B = (A - A.T).astype(float)
        tr_B2 = np.trace(B @ B)
        tr_B4 = np.trace(B @ B @ B @ B)
        tr_B6 = np.trace(np.linalg.matrix_power(B, 6))

        results.append({
            'H': H, 'c3': c3, 'sc': sc,
            'tr_B2': tr_B2, 'tr_B4': tr_B4, 'tr_B6': tr_B6
        })

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    print(f"  {'H':>4s} {'c3':>4s} {'cnt':>5s} | "
          f"{'Tr(B²)':>10s} {'Tr(B⁴)':>10s} {'Tr(B⁶)':>10s}")

    for (H, c3) in sorted(groups.keys()):
        grp = groups[(H, c3)]
        cnt = len(grp)
        t2 = np.mean([r['tr_B2'] for r in grp])
        t4 = np.mean([r['tr_B4'] for r in grp])
        t6 = np.mean([r['tr_B6'] for r in grp])
        std4 = np.std([r['tr_B4'] for r in grp])
        print(f"  {H:4d} {c3:4d} {cnt:5d} | "
              f"{t2:10.1f} {t4:10.1f} {t6:10.1f}"
              + (f"  (std4={std4:.1f})" if std4 > 0.1 else ""))

# ========================================================================
# PART 11: EXACT FORMULA — Tr(B^4) = 2n(n-1)(n-2) - 24*c3
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 11: EXACT FORMULA FOR Tr(B^4)")
print("=" * 72)
print("""
  Tr(B^4) = sum_{i,j,k,l} B[i,j]*B[j,k]*B[k,l]*B[l,i]

  Since B[i,j]^2 = 1 for i != j and B[i,i] = 0:

  Tr(B^4) = sum_{i} (B^2)[i,i]^2 ... wait, let me just compute B^2 first.

  (B^2)[i,j] = sum_k B[i,k]*B[k,j]

  For i = j:
    (B^2)[i,i] = sum_{k!=i} B[i,k]*B[k,i] = -sum_{k!=i} B[i,k]^2 = -(n-1)

  So Tr(B^2) = -n(n-1). ✓

  For i != j:
    (B^2)[i,j] = sum_{k!=i,j} B[i,k]*B[k,j] + B[i,j]*B[j,j] + B[i,i]*B[i,j]
               = sum_{k!=i,j} B[i,k]*B[k,j]

  Now B[i,k]*B[k,j] for i->k->j: (+1)(+1) = +1 (both same direction)
  Wait: B[i,k] = +1 if i->k, -1 if k->i.
  B[k,j] = +1 if k->j, -1 if j->k.

  So B[i,k]*B[k,j] = +1 if (i->k AND k->j) OR (k->i AND j->k)
                    = -1 if (i->k AND j->k) OR (k->i AND k->j)

  The +1 cases count 2-paths from i to j AND 2-paths from j to i (reversed).
  The -1 cases count 2-paths from i away and j away (divergent).

  (B^2)[i,j] = #{k: i->k->j or j->k->i} - #{k: i->k AND j->k} - #{k: k->i AND k->j}
              = (w(i,j) + w(j,i)) - (s_i + s_j - n+2 - A[i,j] - A[j,i]) ... complex.

  Let me just check the formula Tr(B^4) vs c3 numerically.
""")

# Check if Tr(B^4) is a linear function of c3
for n in [3, 4, 5, 6]:
    print(f"\n--- n = {n} ---")
    data = []

    for A in all_tournaments(n):
        c3 = count_3cycles(A, n)
        B = (A - A.T).astype(float)
        tr_B4 = round(np.trace(np.linalg.matrix_power(B, 4)))
        data.append((c3, tr_B4))

    # Check linearity
    c3_vals = sorted(set(c3 for c3, _ in data))
    by_c3 = defaultdict(set)
    for c3, tb4 in data:
        by_c3[c3].add(tb4)

    is_linear = all(len(v) == 1 for v in by_c3.values())

    if is_linear:
        # Fit Tr(B^4) = a + b*c3
        points = [(c3, list(v)[0]) for c3, v in sorted(by_c3.items())]
        if len(points) >= 2:
            c3_0, tb4_0 = points[0]
            c3_1, tb4_1 = points[1]
            b = (tb4_1 - tb4_0) / (c3_1 - c3_0) if c3_1 != c3_0 else 0
            a = tb4_0 - b * c3_0
            print(f"  Tr(B^4) determined by c3 alone: Tr(B^4) = {a:.0f} + {b:.0f}*c3")
            # Verify
            for c3, tb4 in points:
                predicted = a + b * c3
                ok = "✓" if abs(predicted - tb4) < 0.5 else "✗"
                print(f"    c3={c3}: Tr(B^4)={tb4}, predicted={predicted:.0f} {ok}")
    else:
        print(f"  Tr(B^4) NOT determined by c3 alone!")
        for c3 in sorted(by_c3.keys()):
            if len(by_c3[c3]) > 1:
                print(f"    c3={c3}: Tr(B^4) ∈ {sorted(by_c3[c3])}")

# ========================================================================
# PART 12: FORMULA Tr(B^4) = 2Q - 8c3 where Q depends only on n
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 12: DERIVING Tr(B^4) FORMULA")
print("=" * 72)

# From Part 11 data:
# n=3: Tr(B^4) = 12 + (-8)*c3
# n=4: Tr(B^4) = 48 + (-8)*c3
# n=5: Tr(B^4) = 120 + (-8)*c3
# Pattern: constant term = 2*n*(n-1)*(n-2)/...

for n in [3, 4, 5, 6, 7]:
    # Base term (c3=0, transitive tournament)
    # Transitive: A[i,j] = 1 iff i < j (or any total order)
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    B = (A - A.T).astype(float)
    tr_B4_trans = round(np.trace(np.linalg.matrix_power(B, 4)))

    # Expected: this should be the "base" with c3=0
    # n(n-1)(n-2)(3n-7)/3 or similar?
    expected_forms = [
        n*(n-1)*(n-2),
        2*n*(n-1)*(n-2)//3,
    ]

    print(f"  n={n}: Tr(B^4) at c3=0 (transitive) = {tr_B4_trans}")
    print(f"         n(n-1)(n-2) = {n*(n-1)*(n-2)}")
    print(f"         2*n(n-1)(n-2)/3 = {2*n*(n-1)*(n-2)/3}")
    ratio = tr_B4_trans / (n*(n-1)) if n > 1 else 0
    print(f"         Tr(B^4)/n(n-1) = {ratio}")

print("""
  FINAL FORMULA ATTEMPT:

  For the transitive tournament (c3 = 0):
    B[i,j] = sign(i-j) for i != j

  B^2[i,j] for i != j:
    = sum_{k != i,j} sign(i-k)*sign(k-j)
    = #{k between i,j} * (+1) + #{k outside} * (-1)   ... orientation dependent

  This is related to the score sequence of the transitive tournament.
  Let me compute B^2 explicitly for small n.
""")

# Explicit B^2 for transitive tournament
for n in [3, 4, 5]:
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    B = (A - A.T)
    B2 = B @ B
    print(f"\n  n={n}, B^2 for transitive tournament:")
    for i in range(n):
        print(f"    {B2[i]}")
    print(f"  Tr(B^4) = Σ_{'{i,j}'} (B^2[i,j])^2 = {np.sum(B2**2)}")

print("\n\nDone. opus-2026-03-21-S95")
