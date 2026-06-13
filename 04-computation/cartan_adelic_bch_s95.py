#!/usr/bin/env python3
"""
cartan_adelic_bch_s95.py — opus-2026-03-21-S95

ADELIC PRODUCT FORMULA AND HEAT KERNEL BCH ANALYSIS

Part A: Adelic Product Formula
  R_ad(T) = r_∞ × Π_{p|D} r_p where r is the tournament fraction
  Does this satisfy a normalization constraint?

Part B: Exact BCH expansion for exp(tA) = exp(t·tournament) × exp(t·cooperation)
  For n=3,4: exact symbolic computation.
  For regular tournaments: BCH is exact (sectors commute).
  For non-regular: compute correction terms.

Author: opus-2026-03-21-S95
"""

import sys
import numpy as np
from fractions import Fraction
from collections import defaultdict
from itertools import combinations
from scipy.linalg import expm
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

print("=" * 72)
print("  ADELIC PRODUCT FORMULA AND HEAT KERNEL BCH")
print("  opus-2026-03-21-S95")
print("=" * 72)

# ========================================================================
# PART A: ADELIC PRODUCT FORMULA
# ========================================================================
print("\n" + "=" * 72)
print("  PART A: ADELIC CARTAN PRODUCT FORMULA")
print("=" * 72)
print("""
  For the tournament adjacency A, the Cartan fractions are:
    anti_frac = 1/2 (CONSTANT)
    sym_frac = (n²-1)/(2n²)
    scalar_frac = (n-1)/(2n²)

  These are UNIVERSAL — so the adelic product is trivially constant.
  There's no variation between tournaments at the adjacency level.

  For the A^2 matrix (2-step walks), the fractions vary.
  Let's define the adelic product for A^2:

  At the archimedean place:
    r_∞ = anti_frac(A^2) = ||antisym(A²)||² / ||A²||²

  At each prime p dividing the conductor D = odd_part(sum of A^2 entries):
    Reduce A^2 mod p, compute Frobenius norms of parts.

  Question: does r_∞ × Π r_p = constant?

  Actually, the "adelic product formula" |x|_∞ · Π |x|_p = 1
  applies to ELEMENTS of Q, not to matrices. The Cartan fraction
  is not an element of Q — it's a ratio of norms.

  Let's instead look at whether the DETERMINANT of A^2 satisfies
  a product formula, or whether Tr(A^{2k}) mod p gives useful info.
""")

# Check det(A^2) and its factorization
for n in [3, 4, 5]:
    print(f"\n--- n = {n}: determinants of A² ---")
    det_set = defaultdict(list)

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        A2 = A @ A
        det_A2 = round(np.linalg.det(A2))
        det_set[(H, c3)].append(det_A2)

    groups = defaultdict(set)
    for (H, c3), dets in det_set.items():
        for d in dets:
            groups[(H, c3)].add(d)

    print(f"  {'H':>4s} {'c3':>4s} {'det(A²)':>30s}")
    for key in sorted(groups.keys()):
        H, c3 = key
        vals = sorted(groups[key])
        print(f"  {H:4d} {c3:4d} {str(vals):>30s}")

# ========================================================================
# PART A2: Alternative adelic invariant — characteristic polynomial mod p
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART A2: CHARACTERISTIC POLYNOMIAL OF B² MOD PRIMES")
print("=" * 72)
print("""
  B² is symmetric with integer entries. Its characteristic polynomial
  χ(B², λ) has integer coefficients.

  Reducing mod p gives χ(B², λ) mod p ∈ F_p[λ].
  The factorization pattern of χ mod p is an invariant.

  For the transitive tournament T_n, B² has a specific pattern.
  We can check if different tournaments give different χ mod p.
""")

for n in [5]:
    print(f"\n--- n = {n}: char poly of B² mod primes ---")
    results = defaultdict(list)

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        B = (A - A.T).astype(float)
        B2 = B @ B

        # Characteristic polynomial coefficients
        eigs = np.linalg.eigvalsh(B2)
        charpoly = np.round(np.poly(eigs)).astype(int)

        results[(H, c3)].append(tuple(charpoly))

    groups = defaultdict(set)
    for (H, c3), polys in results.items():
        for p in polys:
            groups[(H, c3)].add(p)

    print(f"  {'H':>4s} {'c3':>4s} {'#distinct char polys':>20s}")
    for key in sorted(groups.keys()):
        H, c3 = key
        n_distinct = len(groups[key])
        rep = list(groups[key])[0]
        print(f"  {H:4d} {c3:4d} {n_distinct:20d}  (rep: {rep})")

# ========================================================================
# PART B: EXACT BCH EXPANSION
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART B: EXACT BCH EXPANSION FOR TOURNAMENT ADJACENCY")
print("=" * 72)
print("""
  A = B/2 + (J-I)/2 where B is skew-symmetric and (J-I)/2 is symmetric.

  For REGULAR tournaments (Paley), [B, J] = 0, so
    exp(tA) = exp(t·B/2) · exp(t·(J-I)/2)
  exactly.

  For non-regular tournaments, we need BCH corrections:
    exp(t(X+Y)) = exp(tX) · exp(tY) · exp(Z₂) · exp(Z₃) · ...
  where Z₂ = -t²/2 · [X,Y], Z₃ involves [X,[X,Y]] and [Y,[X,Y]], etc.

  We proved that [B/2, (J-I)/2] = (1/4)[B,J] has entries
  [B,J][i,j] = 2(d_i + d_j) where d_i = s_i - (n-1)/2.

  So Z₂ = -t²/2 · (1/4)[B,J] = -t²/8 · [B,J].

  The question is: how quickly do the BCH corrections converge?
  And do they have tournament-theoretic meaning?
""")

# Compute BCH quality for all tournaments at n=3,4,5
for n in [3, 4, 5]:
    print(f"\n--- n = {n}: BCH convergence ---")
    results = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        c3 = count_3cycles(A, n)
        scores = A.sum(axis=1)
        S2 = np.sum((scores - (n-1)/2)**2)

        Af = A.astype(float)
        B = (Af - Af.T) / 2   # A_anti = B/2
        S = (Af + Af.T) / 2   # A_sym (includes scalar)

        # BCH at various orders, t=1
        E_exact = expm(Af)
        E_bch0 = expm(B) @ expm(S)
        err0 = np.max(np.abs(E_exact - E_bch0))

        # Z2 correction
        comm = B @ S - S @ B
        E_bch1 = E_bch0 @ expm(-0.5 * comm)
        err1 = np.max(np.abs(E_exact - E_bch1))

        # Z3 correction (approximate)
        comm_XY = comm
        comm_X_XY = B @ comm_XY - comm_XY @ B
        comm_Y_XY = S @ comm_XY - comm_XY @ S
        Z3 = (1.0/12) * (comm_X_XY + comm_Y_XY)
        E_bch2 = E_bch1 @ expm(Z3)
        err2 = np.max(np.abs(E_exact - E_bch2))

        results.append({
            'H': H, 'c3': c3, 'S2': S2,
            'err0': err0, 'err1': err1, 'err2': err2
        })

    groups = defaultdict(list)
    for r in results:
        groups[(r['H'], r['c3'])].append(r)

    print(f"  {'H':>4s} {'c3':>4s} {'S₂':>6s} | "
          f"{'BCH₀':>10s} {'BCH₁':>10s} {'BCH₂':>10s}")

    for key in sorted(groups.keys()):
        H, c3 = key
        grp = groups[key]
        S2 = grp[0]['S2']
        e0 = np.mean([r['err0'] for r in grp])
        e1 = np.mean([r['err1'] for r in grp])
        e2 = np.mean([r['err2'] for r in grp])
        print(f"  {H:4d} {c3:4d} {S2:6.1f} | "
              f"{e0:10.2e} {e1:10.2e} {e2:10.2e}"
              + ("  ← regular" if S2 < 0.01 else ""))

# ========================================================================
# PART B2: BCH CORRECTIONS AND TOURNAMENT INVARIANTS
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART B2: BCH CORRECTION NORMS AS TOURNAMENT INVARIANTS")
print("=" * 72)
print("""
  THEOREM: The BCH₀ error ||exp(A) - exp(B/2)·exp((J-I)/2)|| depends
  on the score sequence through S₂.

  Specifically:
  ||Z₂||_F = (1/8)||[B,J]||_F = (1/8)·√(8n·S₂) = √(n·S₂/8)

  This follows from ||[B,J]||_F² = 8n·S₂ (proved earlier: factor 4
  from the 1/4 scaling cancels in the double application).

  Wait, let me recompute: [B/2, (J-I)/2] has ||·||² = n·S₂/2.
  Z₂ = -t²/2 · [B/2, (J-I)/2], so ||Z₂||_F = t²/2 · √(n·S₂/2).

  At t=1: ||Z₂||_F = √(n·S₂/2) / 2 = √(n·S₂) / (2√2).

  So the first BCH correction is PROPORTIONAL to √(n·S₂).
  For the transitive tournament: S₂ = n(n²-1)/12.
  ||Z₂||_trans ∝ √(n²(n²-1)/12) = n√((n²-1)/12).

  For the regular tournament: S₂ = 0, so Z₂ = 0 and BCH is exact. ✓
""")

print("  Verification: ||Z₂||_F vs √(n·S₂/2)/2")
for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    for A in all_tournaments(n):
        scores = A.sum(axis=1)
        S2 = np.sum((scores - (n-1)/2)**2)
        if S2 < 0.01:
            continue

        Af = A.astype(float)
        B = (Af - Af.T) / 2
        S = (Af + Af.T) / 2
        comm = B @ S - S @ B
        Z2_norm = np.sqrt(np.sum(comm**2)) / 2  # ||Z₂||_F

        predicted = np.sqrt(n * S2 / 2) / 2
        match = abs(Z2_norm - predicted) < 1e-10

        H = count_hp(A, n)
        if H == 1 or (n == 5 and H == 15):
            print(f"    H={H}: ||Z₂||={Z2_norm:.6f}, "
                  f"predicted={predicted:.6f}, match={match}")
        break

# ========================================================================
# PART B3: The heat kernel K(t) = Tr(exp(tA)) through Cartan
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART B3: HEAT KERNEL Tr(exp(tA)) AND CARTAN DECOMPOSITION")
print("=" * 72)
print("""
  K(t) = Tr(exp(tA)) = Σ_i exp(t·λ_i) where λ_i are eigenvalues of A.

  For the tournament adjacency A (integer matrix):
    Eigenvalues are algebraic (but generally not rational).
    A = (B + J - I)/2, so λ_i(A) = (λ_i(B) + λ_i(J-I))/2?
    NO — A = B/2 + (J-I)/2 but B and (J-I) don't commute in general.

  For REGULAR tournaments (Paley), they DO commute, so:
    Eigenvalues of A = eigenvalues of B/2 + eigenvalues of (J-I)/2
    = (±iθ_k/2) + ((n-1)/2 or -1/2)

  Wait — for regular tournaments, the eigenvectors of B align with
  those of J-I? Let me check.

  For B = A - A^T (regular tournament): B·1 = 0 (row sums zero)
  And (J-I)·1 = (n-1)·1. So 1 is an eigenvector of both.

  For the other eigenvectors of J-I (which are orthogonal to 1,
  with eigenvalue -1): do they align with B's eigenvectors?

  For CIRCULANT regular tournaments (like Paley T_p), YES — they
  share Fourier eigenvectors.
""")

# Verify: eigenvalue structure for Paley T_5 and T_7
for n, name in [(5, "Paley T_5"), (7, "Paley T_7")]:
    QR = set()
    for a in range(1, n):
        QR.add((a * a) % n)

    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in QR:
                A[i][j] = 1.0

    B = (A - A.T) / 2.0
    sym = (A + A.T) / 2.0
    sigma = np.trace(A) / n
    sym_tl = sym - sigma * np.eye(n)

    eigs_A = np.sort(np.linalg.eigvals(A))[::-1]
    eigs_B = np.sort(np.linalg.eigvals(B))[::-1]
    eigs_sym = np.sort(np.linalg.eigvals(sym))[::-1]

    print(f"\n  {name}:")
    print(f"    A eigenvalues: {[f'{e.real:.4f}{e.imag:+.4f}j' for e in eigs_A]}")
    print(f"    B/2 eigenvalues: {[f'{e.real:.4f}{e.imag:+.4f}j' for e in eigs_B/2]}")
    print(f"    sym eigenvalues: {[f'{e:.4f}' for e in np.real(eigs_sym)]}")

    # Check if A eigenvalues = B/2 + sym eigenvalues
    # For circulant matrices, sort by DFT index
    from numpy.fft import fft

    row = A[0]  # first row of circulant
    fft_A = fft(row)
    print(f"    DFT of first row of A: {[f'{e.real:.4f}{e.imag:+.4f}j' for e in fft_A]}")

    row_B = (A - A.T)[0] / 2.0
    fft_B = fft(row_B)
    print(f"    DFT of first row of B/2: {[f'{e.real:.4f}{e.imag:+.4f}j' for e in fft_B]}")

    row_sym = (A + A.T)[0] / 2.0
    fft_sym = fft(row_sym)
    print(f"    DFT of first row of sym: {[f'{e.real:.4f}' for e in np.real(fft_sym)]}")

    # Check: DFT(A) = DFT(B/2) + DFT(sym)
    sum_check = fft_B + fft_sym
    print(f"    B/2 + sym DFT: {[f'{e.real:.4f}{e.imag:+.4f}j' for e in sum_check]}")
    match = np.allclose(fft_A, sum_check)
    print(f"    DFT(A) = DFT(B/2) + DFT(sym): {match}")

    # Heat kernel at logarithmic rationals
    print(f"\n    Heat kernel K(t) = Tr(exp(tA)):")
    for t in [0.1, 0.5, 1.0, np.log(2), np.log(3)]:
        K = np.real(np.trace(expm(t * A)))
        K_B = np.real(np.trace(expm(t * B / 2)))
        K_sym = np.real(np.trace(expm(t * sym)))
        # For commuting matrices: K = K_B * K_sym / n (not quite)
        print(f"      t={t:.4f}: K={K:.4f}, K_B={K_B:.4f}, K_sym={K_sym:.4f}, "
              f"K_B*K_sym/n={K_B*K_sym/n:.4f}")

# ========================================================================
# PART B4: Exact heat kernel for circulant Paley tournament
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART B4: EXACT HEAT KERNEL FOR PALEY TOURNAMENTS")
print("=" * 72)
print("""
  For a circulant tournament T_p with adjacency A, the eigenvalues are:
    λ_k = Σ_j A[0,j] · ω^{jk}  where ω = e^{2πi/p}

  For Paley T_p: A[0,j] = 1 iff j is a QR mod p.
  So λ_k = Σ_{j ∈ QR} ω^{jk} = Gauss sum type expression.

  For k ≠ 0: λ_k = (χ(k)·√p - 1) / 2 where χ = Legendre symbol,
  if we're careful about conventions.

  Actually: Σ_{j ∈ QR} ω^{jk} is a partial Gauss sum.
  For Paley T_p: QR = {1², 2², ..., ((p-1)/2)²} mod p.

  By the standard result:
    Σ_{j=1}^{p-1} (j/p) · ω^{jk} = (k/p) · g(p)
  where g(p) = Σ_{j=1}^{p-1} (j/p) · ω^j is the Gauss sum.
  g(p) = √p for p ≡ 1 (mod 4), i√p for p ≡ 3 (mod 4).

  And Σ_{j ∈ QR} ω^{jk} = (1/2)(Σ_{j=1}^{p-1}(1+(j/p))·ω^{jk} - 1)
  Hmm, this needs care. Let me just verify numerically.
""")

for p in [5, 7, 11, 13]:
    QR = set()
    for a in range(1, p):
        QR.add((a * a) % p)

    omega = np.exp(2j * np.pi / p)
    eigs = []
    for k in range(p):
        lk = sum(omega**(j*k) for j in range(p) if ((j - 0) % p) in QR and j != 0)
        eigs.append(lk)

    # Sort by real part
    eigs_sorted = sorted(eigs, key=lambda x: -x.real)
    print(f"\n  T_{p}: eigenvalues of A:")
    for k, e in enumerate(eigs):
        print(f"    k={k}: λ = {e.real:.6f} {e.imag:+.6f}j")

    # Heat kernel K(t) = Σ exp(t·λ_k)
    print(f"  K(ln2) = Σ exp(ln2·λ_k) = Σ 2^λ_k:")
    K_ln2 = sum(np.exp(np.log(2) * e) for e in eigs)
    print(f"    K(ln2) = {K_ln2.real:.6f}")

    # Check if K(ln2) is algebraic
    # For rational eigenvalues, 2^λ would be algebraic.
    # But eigenvalues involve Gauss sums (√p), so they're algebraic.
    # 2^(a+b√p) = 2^a · 2^{b√p} — the second factor is transcendental!
    # So K(t) is generally transcendental for t ≠ 0.

    # But K(ln2) for T_5 should be checkable
    if p == 5:
        print(f"    T_5 eigenvalues: λ_0 = (p-1)/2 = 2,")
        print(f"    λ_k for k=1,...,4 involve Gauss sum g(5) = √5")
        print(f"    K(ln2) = 2^2 + 4·2^{eigs[1].real:.4f} · cos(...)")

print("\n\nDone. opus-2026-03-21-S95")
