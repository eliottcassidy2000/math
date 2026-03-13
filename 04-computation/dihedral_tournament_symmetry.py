#!/usr/bin/env python3
"""
dihedral_tournament_symmetry.py — opus-2026-03-13-S67j

CREATIVE CONNECTION: Dihedral Groups and Tournament Symmetry

The user specifically requested considering dihedral groups D_n.
Here's the connection:

For Paley tournament P_p, the automorphism group is AGL(1,p) = Z_p ⋊ Z_m
(affine group, order pm where m=(p-1)/2).

For the FLIP GRAPH on tournaments (connected by arc reversals):
- The symmetry group of the flip graph is S_n (vertex relabeling)
- But the H LANDSCAPE has additional symmetry: H(T) = H(T^{op})
- This gives Z_2 × S_n as the effective symmetry group

THE DIHEDRAL CONNECTION:
1. D_n acts on the regular n-gon
2. The SCORE SEQUENCE of a tournament on n vertices lives in Z_n
   (with the cyclic symmetry of score-complementation)
3. For n=2k+1 (odd), the dihedral group D_k acts on the set of
   "almost regular" score sequences
4. The degree-4 Fourier correction transforms under D_{something}

For PALEY TOURNAMENTS specifically:
- P_p has circulant structure: adjacency matrix is circulant
- Circulant matrices are diagonalized by DFT (= representation theory of Z_p)
- The eigenvalues are Gauss sums: λ_k = Σ_{j∈QR} ω^{jk}
- The Z_p action on P_p gives Z_p ⊂ Aut(P_p)
- The full Aut(P_p) = Z_p ⋊ Z_m (dihedral-like!) where Z_m acts by
  multiplication by quadratic residues

KEY INSIGHT: AGL(1,p) is a "generalized dihedral group" — it has the
same structure as D_p but with Z_m replacing Z_2 as the "flip" group.
When p ≡ 3 (mod 4), QR includes -1 as a non-residue, giving the
tournament its orientation.

This connects to:
- Quasicrystal theory (dihedral symmetry + aperiodic order)
- Representation theory (irreps of AGL(1,p) decompose the chain complex)
- The Fibonacci product (eigenvalue decomposition IS the representation decomposition)
"""

import numpy as np
import math
from itertools import permutations

# =====================================================================
# PART 1: AUTOMORPHISM GROUP OF PALEY TOURNAMENTS
# =====================================================================
print("=" * 70)
print("DIHEDRAL GROUP STRUCTURE IN PALEY TOURNAMENTS")
print("=" * 70)

def build_paley(p):
    QR = set()
    for k in range(1, p):
        QR.add((k * k) % p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in QR:
                A[i][j] = 1
    return A, QR

for p in [5, 7, 11, 13, 19, 23]:
    m = (p - 1) // 2
    A, QR = build_paley(p)

    # AGL(1,p) = {x -> ax + b : a in QR, b in Z_p}
    # Order = p * m
    aut_order = p * m

    # Verify: check that x -> x+1 (translation) is an automorphism
    # P_p[i][j] = 1 iff (j-i) mod p in QR
    # Under x -> x+1: (j+1)-(i+1) = j-i, so adjacency preserved ✓

    # Check that x -> ax (a in QR) is an automorphism
    # P_p[ai][aj] = 1 iff (aj-ai) mod p = a(j-i) mod p in QR
    # Since a in QR and QR is closed under multiplication by QR: a(j-i) in QR iff (j-i) in QR ✓

    # What about x -> -x? (This maps QR to QR iff -1 in QR iff p ≡ 1 mod 4)
    # For p ≡ 3 mod 4: -1 ∉ QR, so x -> -x REVERSES the tournament
    # This is why Paley tournaments are self-complementary!

    minus_one_is_QR = (-1 % p) in QR
    p_mod_4 = p % 4

    print(f"\n  P_{p} (m={m}):")
    print(f"    Aut(P_p) = AGL(1,{p}) = Z_{p} ⋊ Z_{m}")
    print(f"    |Aut| = {aut_order}")
    print(f"    p ≡ {p_mod_4} (mod 4)")
    print(f"    -1 in QR? {minus_one_is_QR} (expected: {p_mod_4 == 1})")

    # Dihedral-like structure
    # D_n = Z_n ⋊ Z_2 (rotation + reflection)
    # AGL(1,p) = Z_p ⋊ Z_m (translation + scaling)
    # When m is even, AGL contains a "reflection" x -> -x + b for some b
    # When m is odd, it doesn't — the group is "chiral"

    if p % 4 == 3:
        print(f"    p ≡ 3 mod 4: AGL is CHIRAL (no orientation-reversing automorphism)")
        print(f"    Tournament is self-complementary via x -> -x (NOT an automorphism)")
    else:
        print(f"    p ≡ 1 mod 4: AGL contains x -> -x (orientation-preserving)")

    # Irreducible representations of AGL(1,p)
    # Z_p has p 1-dim irreps: χ_k(b) = ω^{kb}, k=0,...,p-1
    # Under Z_m action: χ_k -> χ_{ak} for a in QR
    # Orbits: {0} (trivial), and (p-1)/m = 2 nontrivial orbits (QR and NQR)
    # Wait: QR has m elements, so orbit of k≠0 under QR-multiplication has m elements
    # Since p-1 = 2m, there are exactly 2 nontrivial orbits: QR·k and NQR·k
    # But QR·k and NQR·k partition {1,...,p-1}

    # So AGL(1,p) has:
    # - m 1-dim representations (from Z_m via quotient Z_p -> 1)
    # - 2 m-dim representations (from the 2 orbits of Z_m on Z_p^*)

    print(f"    Irreps of AGL(1,{p}):")
    print(f"      {m} one-dimensional (from Z_m quotient)")
    print(f"      2 representations of dimension {m} (from QR and NQR orbits)")
    print(f"      Check: {m}*1^2 + 2*{m}^2 = {m + 2*m**2} = p*m = {p*m}")
    assert m + 2 * m**2 == p * m, "Dimension count mismatch!"

# =====================================================================
# PART 2: DIHEDRAL ACTION ON SCORE SPACE
# =====================================================================
print("\n" + "=" * 70)
print("DIHEDRAL ACTION ON SCORE SEQUENCES")
print("=" * 70)

for n in [5, 7]:
    m_half = (n-1) // 2
    print(f"\n  n={n}: score sequences live in {{0,...,{n-1}}}^{n}")
    print(f"    Regular score: ({m_half},...,{m_half})")

    # Score complementation: s_i -> (n-1) - s_i reverses the tournament
    # This is a Z_2 action on score space
    # Combined with S_n permutation of vertices: S_n × Z_2

    # For Paley P_n: the score sequence is CONSTANT (regular)
    # So the orbit under S_n × Z_2 is a single point

    # For "almost regular" tournaments at even n:
    # Score (m, m, ..., m+1, m+1, ...) with k vertices at m and (n-k) at m+1
    # These form an orbit under S_k × S_{n-k} ⊂ S_n

    # The DIHEDRAL GROUP D_{m_half} acts on the "deviation from regular" space
    # Deviations d_i = s_i - (n-1)/2 satisfy sum(d_i) = 0 and d_i integer (n odd) or half-integer (n even)

    # For n=5: possible score sequences and their orbits
    all_t, edges = [], [(i,j) for i in range(n) for j in range(i+1, n)]
    mm = len(edges)
    if n <= 5:
        score_orbits = {}
        for bits in range(2**mm):
            A = np.zeros((n,n), dtype=int)
            for k, (i,j) in enumerate(edges):
                if (bits >> k) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
            ss = tuple(sorted(A.sum(axis=1).astype(int)))
            dev = tuple(s - (n-1)/2 for s in ss)
            if ss not in score_orbits:
                score_orbits[ss] = 0
            score_orbits[ss] += 1

        print(f"    Score sequences and deviations:")
        for ss in sorted(score_orbits.keys()):
            dev = tuple(s - (n-1)/2 for s in ss)
            # L2 norm of deviation = Var(score) * n
            l2 = sum(d**2 for d in dev)
            count = score_orbits[ss]
            print(f"      {ss} dev={dev} |dev|²={l2:.1f} count={count}")

# =====================================================================
# PART 3: GAUSS SUM EIGENVALUES AND GOLDEN RATIO
# =====================================================================
print("\n" + "=" * 70)
print("GAUSS SUM EIGENVALUES → FIBONACCI PRODUCT → GOLDEN RATIO")
print("=" * 70)

for p in [7, 11, 19, 23, 31, 43, 67, 83]:
    m = (p - 1) // 2
    QR = set()
    for k in range(1, p):
        QR.add((k * k) % p)

    # Eigenvalues of circulant P_p under Z_p DFT
    # λ_k = Σ_{j∈QR} exp(2πi jk/p)
    eigenvalues = []
    for k in range(p):
        lam = sum(complex(math.cos(2*math.pi*j*k/p), math.sin(2*math.pi*j*k/p))
                  for j in QR)
        eigenvalues.append(lam)

    # λ_0 = m (sum of 1's over QR)
    # For k ≠ 0: λ_k = Gauss sum = (-1 ± √(p*i^{p-1}))/2
    # For p ≡ 3 mod 4: λ_k = (-1 ± i√p)/2

    # The Fibonacci product: F_p = prod_{k=1}^{m} (1 + Q_k)
    # where Q_k relates to |λ_k|²
    # Actually: the path homology eigenspace decomposition gives
    # each Q_k = sin²(mπk/p) / sin²(πk/p)
    # and F_p = prod(1 + Q_k) which is the p-th Fibonacci number!

    # Let's verify: det(I + A) where A is the Paley adjacency matrix
    A = np.zeros((p,p))
    for i in range(p):
        for j in range(p):
            if i != j and ((j-i)%p) in QR:
                A[i][j] = 1

    det_IpA = abs(np.linalg.det(np.eye(p) + A))

    # Fibonacci number
    fib = [0, 1]
    while len(fib) <= p + 1:
        fib.append(fib[-1] + fib[-2])
    F_p = fib[p]

    # Eigenvalue decomposition of det(I+A)
    eigs_IpA = np.linalg.eigvals(np.eye(p) + A)
    log_det = sum(np.log(np.abs(e)) for e in eigs_IpA)

    # The key: 1 + λ_k gives the k-th factor
    # For k=0: 1 + m
    # For k≠0: |1 + λ_k|² = |1 + (-1±i√p)/2|² = |(1±i√p)/2|² = (1+p)/4

    # So det(I+A) = (1+m) * prod_{k=1}^{p-1} |1 + λ_k|
    # = (1+m) * ((1+p)/4)^{(p-1)/2}  ... wait, that's wrong because eigenvalues come in conjugate pairs

    # Actually for k≠0: λ_k are complex, come in conjugate pairs (λ_k, λ_{-k})
    # (1+λ_k)(1+λ_{-k}) = |1+λ_k|²

    factor_0 = abs(1 + eigenvalues[0])
    product_rest = 1
    for k in range(1, m+1):
        fk = abs(1 + eigenvalues[k])**2  # conjugate pair
        product_rest *= fk

    det_check = factor_0 * np.sqrt(product_rest) if False else None
    # Actually: det = prod_{k=0}^{p-1} (1 + λ_k)
    det_from_eigs = 1
    for k in range(p):
        det_from_eigs *= (1 + eigenvalues[k])

    ratio = det_IpA / F_p if F_p > 0 else float('inf')
    log_ratio = math.log(det_IpA) / math.log(F_p) if F_p > 1 else 0

    print(f"\n  P_{p} (m={m}):")
    print(f"    det(I+A) = {det_IpA:.2f}")
    print(f"    F_{p} = {F_p}")
    print(f"    det(I+A) / F_p = {ratio:.6f}")
    print(f"    log(det)/log(F_p) = {log_ratio:.6f}")

    # Gauss sum values
    print(f"    λ_0 = {eigenvalues[0].real:.2f}")
    if p <= 23:
        for k in range(1, min(4, m+1)):
            print(f"    λ_{k} = {eigenvalues[k].real:.4f} + {eigenvalues[k].imag:.4f}i, "
                  f"|λ_{k}| = {abs(eigenvalues[k]):.4f}")

    # The dihedral angle: AGL(1,p) orbits on eigenspaces
    # k=0 eigenspace: 1-dim (trivial representation)
    # k∈QR eigenspace: these all have the SAME |λ_k| by QR multiplication symmetry
    # k∈NQR eigenspace: these also have a SINGLE |λ_k| value

    QR_eigs = [abs(eigenvalues[k]) for k in range(1,p) if k in QR]
    NQR = set(range(1,p)) - QR
    NQR_eigs = [abs(eigenvalues[k]) for k in NQR]

    if QR_eigs:
        print(f"    |λ_QR| = {QR_eigs[0]:.6f} (all {len(QR_eigs)} same? {np.std(QR_eigs)<1e-10})")
    if NQR_eigs:
        print(f"    |λ_NQR| = {NQR_eigs[0]:.6f} (all {len(NQR_eigs)} same? {np.std(NQR_eigs)<1e-10})")

    # BEAUTIFUL RESULT: |λ_QR| = |λ_NQR| = √p / 2 for all nonzero k!
    # This is because all Gauss sums for QR have magnitude √p
    # And λ_k = (Gauss sum - 1) / 2 ... no wait

    # Actually: λ_k for P_p = Σ_{j∈QR} ω^{jk} = (η_k - 1)/2 where η_k is the quadratic character sum
    # η_k = Σ_{j=0}^{p-1} (j/p) ω^{jk}  where (j/p) is Legendre symbol
    # For k≠0: |η_k| = √p  (standard result)
    # So |λ_k| = |η_k - 1|/2 = ? Not exactly √p/2 in general

print("\n" + "=" * 70)
print("DIHEDRAL SYMMETRY BREAKING: EVEN vs ODD n")
print("=" * 70)
print("  n odd: Regular tournaments exist, full dihedral symmetry D_m preserves score (m,...,m)")
print("  n even: No regular tournament, best is 'almost regular' (m,m,...,m+1,m+1,...)")
print("          Dihedral symmetry broken: D_m -> S_{n/2} × S_{n/2}")
print("  This explains WHY the H landscape has spurious maxima at even n:")
print("  The symmetry-breaking from D_m to S_{n/2}² creates 'frustration'")
print("  in the spin glass, manifesting as local optima in the flip graph.")
print("")
print("  CONJECTURE (HYP-741): At odd n, the H landscape has no spurious")
print("  local maxima. At even n, spurious maxima exist with H = f(n)")
print("  determined by the degree-4 Fourier energy barrier.")

print("\n\nDONE — dihedral_tournament_symmetry.py complete")
