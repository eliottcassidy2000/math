#!/usr/bin/env python3
"""
polygon_helix_knacci_s90.py — Deep connections: polygons, helices, k-nacci
opus-2026-03-15-S90

GRAND QUESTION: The transfer matrix char poly λ³-λ²-xλ-x at x=1 gives
the tribonacci equation. But what if we PARAMETRIZE x and trace the
eigenvalues? And what is the k-nacci generalization?

Key insight: The NUD transfer matrix for k-block clusters would have
char poly of degree k+1, giving k-nacci at x=1. The 3-state (k=2)
case from S112 gives tribonacci. The 2-state (k=1) Fibonacci case
would be the "no adjacent descents" count.

CONNECTION TO TOURNAMENTS ON POLYGONS:
- n-gon vertices at roots of unity
- k-nacci = constraint on k consecutive 1s
- Independent sets on path graph P_n with gap ≥ k
- This IS the Lucas number / Zeckendorf constraint
- The circulant independent set count on C_n = Lucas analog
"""

import numpy as np
from math import factorial, comb, log, pi, sqrt, gcd
from itertools import permutations, combinations
from fractions import Fraction

# ======================================================================
# PART 1: THE k-NACCI TRANSFER MATRICES
# ======================================================================
print("=" * 70)
print("PART 1: k-NACCI TRANSFER MATRICES AND THEIR MEANING")
print("=" * 70)

print("""
For each k ≥ 1, there is a CONSTRAINT: "no k consecutive 1s."

This constraint is enforced by a (k+1)-state transfer matrix M_k.
The states track the current "run length" of consecutive 1s:
  State 0: just saw a 0
  State 1: saw exactly 1 consecutive 1
  ...
  State k-1: saw exactly k-1 consecutive 1s (next 1 is FORBIDDEN)

The transfer matrix M_k has size k × k:
  M_k[i][0] = 1 for all i (can always place 0, resetting run to 0)
  M_k[i][i+1] = 1 for i < k-1 (placing 1, extending run)
  M_k[k-1][k] is FORBIDDEN (would give k consecutive 1s)
""")

def knacci_transfer_matrix(k):
    """Construct the k-nacci transfer matrix (size k × k).
    Entry M[i][j] = 1 if state j can follow state i.
    State i = "run of i consecutive 1s ending here".
    """
    M = np.zeros((k, k))
    for i in range(k):
        M[i][0] = 1  # Can always place 0 (go to state 0)
        if i < k-1:
            M[i][i+1] = 1  # Can place 1 (extend run)
    return M

for k in range(2, 8):
    M = knacci_transfer_matrix(k)
    eigenvalues = np.linalg.eigvals(M)
    eigenvalues = sorted(eigenvalues, key=lambda z: -abs(z))
    dom = eigenvalues[0].real

    names = {2: "Fibonacci", 3: "tribonacci", 4: "tetranacci",
             5: "pentanacci", 6: "hexanacci", 7: "heptanacci"}
    name = names.get(k, f"{k}-nacci")

    # Count binary strings of length n with no k consecutive 1s
    # This is the (n+k)-th k-nacci number (approximately)
    v = np.ones(k)  # initial state: one way to be in each state
    counts = []
    for n in range(1, 15):
        v = M.T @ v  # note: we need M transpose for state evolution
        counts.append(int(round(sum(v))))

    # The characteristic polynomial of M
    charpoly = np.round(np.poly(M), 10)

    print(f"\nk={k} ({name}):")
    print(f"  Transfer matrix M ({k}×{k}):")
    for row in M:
        print(f"    {[int(x) for x in row]}")
    print(f"  Char poly coefficients: {charpoly}")
    print(f"  Dominant eigenvalue: {dom:.10f}")
    print(f"  # binary strings of length n with no {k} consecutive 1s:")
    print(f"    n=1..14: {counts}")

# ======================================================================
# PART 2: COUNTING INDEPENDENT SETS ON PATH AND CYCLE GRAPHS
# ======================================================================
print()
print("=" * 70)
print("PART 2: INDEPENDENT SETS ON PATHS AND CYCLES — LUCAS NUMBERS")
print("=" * 70)

print("""
The number of independent sets of size k on the path graph P_n
(vertices 1,...,n, edges {i,i+1}) is:
  i_k(P_n) = C(n-k+1, k)

The total number of independent sets (including empty set) is:
  I(P_n) = Σ_k C(n-k+1, k) = F(n+2) (Fibonacci!)

For the CYCLE graph C_n:
  I(C_n) = L(n) (Lucas numbers!) for n ≥ 3.

These connect to the Zeckendorf representation:
  The independent sets of P_n biject with binary strings of length n
  with no two adjacent 1s. The total count is F(n+2).

For k-GAP independent sets (no two selected vertices within distance k):
  On P_n: count = k-nacci(n + something)
  On C_n: count = k-Lucas analog
""")

# Compute I(P_n, x) and I(C_n, x) for small n
def independence_poly_path(n):
    """Independence polynomial of path graph P_n."""
    if n == 0:
        return [1]
    if n == 1:
        return [1, 1]  # 1 + x
    # Recurrence: I(P_n, x) = I(P_{n-1}, x) + x * I(P_{n-2}, x)
    prev2 = [1]         # P_0
    prev1 = [1, 1]      # P_1
    for i in range(2, n+1):
        # Add: prev1 + x * prev2
        new = [0] * (max(len(prev1), len(prev2)+1))
        for j in range(len(prev1)):
            new[j] += prev1[j]
        for j in range(len(prev2)):
            new[j+1] += prev2[j]
        prev2, prev1 = prev1, new
    return prev1

def independence_poly_cycle(n):
    """Independence polynomial of cycle graph C_n."""
    # I(C_n, x) = I(P_{n-1}, x) + x * I(P_{n-3}, x)  for n ≥ 4
    # Or: I(C_n, x) = I(P_n, x) - x^2 * I(P_{n-4}, x) ... various formulas
    # Better: I(C_n, x) = I(P_{n-1}, x) + x * I(P_{n-3}, x) for n ≥ 3
    if n == 3:
        return [1, 3]  # 1 + 3x (three independent singletons)
    p_nm1 = independence_poly_path(n-1)
    p_nm3 = independence_poly_path(n-3)
    result = [0] * max(len(p_nm1), len(p_nm3)+1)
    for j in range(len(p_nm1)):
        result[j] += p_nm1[j]
    for j in range(len(p_nm3)):
        result[j+1] += p_nm3[j]
    return result

print("Path graph independence polynomials I(P_n, x):")
for n in range(1, 11):
    poly = independence_poly_path(n)
    at_1 = sum(poly)  # Should be Fibonacci
    at_2 = sum(c * 2**k for k, c in enumerate(poly))
    print(f"  P_{n}: {poly}  I(P_{n},1)={at_1}  I(P_{n},2)={at_2}")

print("\nCycle graph independence polynomials I(C_n, x):")
for n in range(3, 11):
    poly = independence_poly_cycle(n)
    at_1 = sum(poly)  # Should be Lucas
    at_2 = sum(c * 2**k for k, c in enumerate(poly))
    print(f"  C_{n}: {poly}  I(C_{n},1)={at_1}  I(C_{n},2)={at_2}")

# ======================================================================
# PART 3: THE FIBONACCI-TOURNAMENT CONNECTION
# ======================================================================
print()
print("=" * 70)
print("PART 3: FIBONACCI NUMBERS IN TOURNAMENTS")
print("=" * 70)

print("""
The tournament-Fibonacci connection runs deep:

1. INDEPENDENCE POLYNOMIAL: I(Ω(T), x) gives H(T) at x=2.
   For the conflict graph Ω(T), the independence polynomial at x=2
   counts Hamiltonian paths (OCF/Grinberg-Stanley).

2. FIBONACCI IN TRANSFERS: The 2-state transfer matrix for "no unit descent"
   gives Fibonacci-like counts. The 3-state matrix gives tribonacci.

3. ZECKENDORF AND TOURNAMENTS:
   H(T) = I(Ω(T), 2) = Σ_k α_k · 2^k
   where α_k = #{independent sets of size k in Ω(T)}.
   This is a "base-2 expansion" of H with Fibonacci-like constraints
   on the α_k (from the structure of Ω).

4. CYCLIC TOURNAMENTS AND LUCAS:
   For circulant tournaments on Z_n, the conflict graph has circular
   structure. The independence polynomial of circular graphs involves
   Lucas numbers.

5. PALEY = QUADRATIC RESIDUE CODE:
   The Paley tournament T_p uses QR(p) as connection set.
   The eigenvalues are Gauss sums at p-th roots of unity.
   The "base" is the Gauss sum modulus √p.
""")

# Compute: for circulant tournaments, is Ω(T) related to a cycle graph?
print("Independence polynomials of Ω(T) for circulant tournaments:")

def circulant_tournament(n, S):
    """Build circulant tournament on Z_n with connection set S."""
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            T[i][j] = 1
    return T

def conflict_graph(T, n):
    """Build conflict graph Ω(T): vertices = odd cycles, edges = share vertex."""
    # For small n, find all directed 3-cycles
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] + T[j][k] + T[k][i] == 3:
                    cycles.append(frozenset([i,j,k]))
                elif T[i][k] + T[k][j] + T[j][i] == 3:
                    cycles.append(frozenset([i,j,k]))
    # Remove duplicates
    cycles = list(set(cycles))

    # Build adjacency
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for a in range(nc):
        for b in range(a+1, nc):
            if cycles[a] & cycles[b]:  # share a vertex
                adj[a][b] = adj[b][a] = 1

    return cycles, adj

def independence_poly_graph(adj, nc):
    """Compute independence polynomial of a graph by brute force."""
    poly = [0] * (nc + 1)
    for mask in range(2**nc):
        # Check if selected vertices form an independent set
        selected = [i for i in range(nc) if (mask >> i) & 1]
        is_indep = True
        for a in range(len(selected)):
            for b in range(a+1, len(selected)):
                if adj[selected[a]][selected[b]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            poly[len(selected)] += 1
    return poly

for p in [5, 7]:
    print(f"\n--- Circulant tournaments on Z_{p} ---")

    if p % 4 == 3:
        # Paley
        qr = set()
        for a in range(1, p):
            qr.add((a*a) % p)
        S_paley = sorted(qr)
        T = circulant_tournament(p, S_paley)
        cycles, adj = conflict_graph(T, p)
        ipoly = independence_poly_graph(adj, len(cycles))
        H = sum(c * 2**k for k, c in enumerate(ipoly))
        print(f"  Paley T_{p}: {len(cycles)} 3-cycles")
        print(f"  I(Ω, x) = {ipoly}")
        print(f"  H = I(Ω, 2) = {H}")

    # Interval
    S_int = list(range(1, (p+1)//2))
    T = circulant_tournament(p, S_int)
    cycles, adj = conflict_graph(T, p)
    ipoly = independence_poly_graph(adj, len(cycles))
    H = sum(c * 2**k for k, c in enumerate(ipoly))
    print(f"  Interval C_{p}: {len(cycles)} 3-cycles")
    print(f"  I(Ω, x) = {ipoly}")
    print(f"  H = I(Ω, 2) = {H}")

# ======================================================================
# PART 4: THE HELIX AS EIGENVALUE EVOLUTION
# ======================================================================
print()
print("=" * 70)
print("PART 4: THE GRAND HELIX — EIGENVALUE EVOLUTION ACROSS k")
print("=" * 70)

print("""
As k increases (Fibonacci → tribonacci → ... → k-nacci → ...):
  - The dominant eigenvalue τ_k increases toward 2
  - More complex conjugate pairs appear
  - Each pair traces a helix as we raise to the n-th power
  - The helices NEST: inner helices (earlier k) inside outer helices (later k)

The LIMIT k → ∞:
  - τ_∞ = 2 (every binary string allowed = 2^n total)
  - All complex roots collapse to the unit circle
  - The "helix" flattens into a circle

This is the SPECTRAL GEOMETRY of the k-nacci family:
  k=2 (Fibonacci): 1 helix with frequency arg(ψ) = π (alternating)
  k=3 (tribonacci): 1 helix with frequency ≈ 124.7°
  k=4 (tetranacci): 2 helices with different frequencies
  ...
  k→∞: all helices crowd the unit circle
""")

print("Eigenvalue table for k-nacci transfer matrices:")
print(f"{'k':>3} | {'τ_k':>12} | {'|λ_c|':>8} | {'arg°':>8} | {'period':>8} | {'decays?':>8}")
print("-" * 65)

for k in range(2, 15):
    M = knacci_transfer_matrix(k)
    eigs = np.linalg.eigvals(M)
    eigs = sorted(eigs, key=lambda z: -abs(z))
    tau_k = eigs[0].real

    # Find the largest complex eigenvalue
    complex_eigs = [e for e in eigs if abs(e.imag) > 1e-10]
    if complex_eigs:
        lam_c = max(complex_eigs, key=abs)
        mod = abs(lam_c)
        arg = np.angle(lam_c) * 180 / pi
        period = 360 / abs(np.angle(lam_c) * 180 / pi) if abs(np.angle(lam_c)) > 1e-10 else float('inf')
        decays = "yes" if mod < 1 else "NO"
    else:
        mod = 0
        arg = 0
        period = float('inf')
        decays = "N/A"

    print(f"{k:3d} | {tau_k:12.8f} | {mod:8.5f} | {arg:8.2f} | {period:8.2f} | {decays:>8}")

# ======================================================================
# PART 5: THE TOURNAMENT POLYGON — WHERE EIGENVECTORS LIVE
# ======================================================================
print()
print("=" * 70)
print("PART 5: EIGENVECTORS ON THE POLYGON")
print("=" * 70)

print("""
For a circulant tournament on Z_n with connection set S:
  - Vertices sit at ω^j = e^{2πij/n} on the unit circle
  - The adjacency matrix is diagonalized by the DFT
  - Eigenvector k: χ_k(j) = ω^{jk} (the k-th character)
  - Eigenvalue k: λ_k = Σ_{s∈S} ω^{ks}

The eigenvectors ARE the vertices of the polygon (rotated by k)!

For the Paley tournament T_p:
  λ_k = (χ(k)·√p - 1)/2 where χ = Legendre symbol
  - For k = QR: λ_k = (√p - 1)/2 > 0
  - For k = NR: λ_k = (-√p - 1)/2 < 0
  - The eigenvalues split into TWO values (not n-1 distinct ones!)
  - This 2-level structure is why Paley is "almost Fibonacci" (2 eigenvalues)

For the interval tournament C_p^{1,...,(p-1)/2}:
  λ_k = (sin(kπ(p+1)/(2p)) / sin(kπ/p) - 1) / 2
  - These give (p-1)/2 distinct eigenvalue magnitudes
  - The dominant λ_1 ~ p/(2π) as p → ∞
  - The subdominant eigenvalues decay like 1/sin(kπ/p)
""")

# Show eigenvalue distribution for Paley vs Interval
for p in [7, 11, 13, 19, 23]:
    print(f"\n--- p = {p} ---")
    omega = np.exp(2j * pi / p)

    # Paley eigenvalues (for p ≡ 3 mod 4)
    if p % 4 == 3:
        qr = set()
        for a in range(1, p):
            qr.add((a * a) % p)
        S = sorted(qr)
        eigs_paley = []
        for k in range(p):
            lam = sum(omega**(k*s) for s in S)
            eigs_paley.append(lam)

        # Should be: λ_0 = |S| = (p-1)/2, all others = (±√p - 1)/2
        distinct = set(round(e.real, 4) for e in eigs_paley)
        print(f"  Paley T_{p}: eigenvalues = {sorted(distinct, reverse=True)}")
        print(f"    (p-1)/2 = {(p-1)/2}, (√p-1)/2 = {(sqrt(p)-1)/2:.4f}, (-√p-1)/2 = {(-sqrt(p)-1)/2:.4f}")

    # Interval eigenvalues
    S_int = list(range(1, (p+1)//2))
    eigs_int = []
    for k in range(p):
        lam = sum(omega**(k*s) for s in S_int)
        eigs_int.append(lam)
    mods = sorted(set(round(abs(e), 4) for e in eigs_int), reverse=True)
    print(f"  Interval: |eigenvalues| = {mods[:6]}{'...' if len(mods) > 6 else ''}")

# ======================================================================
# PART 6: THE HELICAL PERMANENT — H(T) AS SUM OF HELICES
# ======================================================================
print()
print("=" * 70)
print("PART 6: THE HELICAL PERMANENT — H(T) DECOMPOSED")
print("=" * 70)

print("""
H(T) = permanent(A) for the adjacency matrix A.

For circulant T with eigenvalues λ_0, ..., λ_{n-1}:
  The permanent is NOT simply Π λ_k (that's the determinant).
  But the permanent CAN be expressed via the eigenvalues using
  MacMahon's master theorem or the Ryser formula.

The key insight: H(T) = I(Ω(T), 2), and I(Ω, x) is a polynomial
in x whose coefficients α_k count independent sets of size k.

For Paley T_p with only 2 distinct eigenvalue magnitudes:
  The permanent concentrates on a few terms.
  H(T_p) ≈ C · ((p-1)/2)^{p-1} · correction
  The correction involves the Gauss sum eigenvalues.

For the interval tournament:
  H = Σ over all perfect matchings of arc weights
  The dominant eigenvalue λ_1 ~ p/(2π) gives the main term
  The subdominant eigenvalues give helical corrections

ON THE RIEMANN SURFACE:
  log H(T) = Σ_k log(contribution from eigenspace k)
  Each eigenspace k contributes a term that, on the Riemann surface,
  is a straight line with slope (log|λ_k|, arg(λ_k)).
  The total log H is a SUPERPOSITION of these lines.
  In the flat complex plane, this superposition appears as a
  sum of spirals = a complex helix pattern.
""")

# Compute H for interval tournaments via direct enumeration at small p
for p in [5, 7]:
    S = list(range(1, (p+1)//2))
    T = circulant_tournament(p, S)

    # Count Hamiltonian paths
    H = 0
    for perm in permutations(range(p)):
        ok = True
        for i in range(p-1):
            if T[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            H += 1

    # Eigenvalues
    omega = np.exp(2j * pi / p)
    eigs = [sum(omega**(k*s) for s in S) for k in range(p)]

    print(f"\nInterval T on Z_{p}:")
    print(f"  H = {H}")
    print(f"  Eigenvalues: {[f'{e.real:.3f}+{e.imag:.3f}i' for e in eigs[:4]]}")
    print(f"  |λ_0| = {abs(eigs[0]):.4f} = (p-1)/2 = {(p-1)/2}")
    print(f"  |λ_1| = {abs(eigs[1]):.4f}")
    if p == 7:
        # H(interval on Z_7) should be same as H(Paley T_7) since they're isomorphic
        qr = set()
        for a in range(1, p):
            qr.add((a*a) % p)
        S_p = sorted(qr)
        T_p = circulant_tournament(p, S_p)
        H_p = sum(1 for perm in permutations(range(p))
                  if all(T_p[perm[i]][perm[i+1]] for i in range(p-1)))
        print(f"  H(Paley T_7) = {H_p}")
        print(f"  Same? {H == H_p}")

# ======================================================================
# PART 7: BRINGING IT ALL TOGETHER
# ======================================================================
print()
print("=" * 70)
print("PART 7: THE FULL PICTURE")
print("=" * 70)

print("""
╔══════════════════════════════════════════════════════════════════════╗
║                    THE HELICAL SIMPLICIAL RÉDEI                    ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  POLYGON (roots of unity)                                          ║
║     │                                                              ║
║     ├── Circulant tournaments ←→ Connection sets S                 ║
║     ├── DFT diagonalizes adjacency matrix                          ║
║     ├── Eigenvalues = character sums over S                        ║
║     └── Eigenvectors = the polygon vertices themselves!            ║
║                                                                    ║
║  DIHEDRAL SYMMETRY                                                 ║
║     │                                                              ║
║     ├── Rotations ↔ Circulant automorphisms                        ║
║     ├── Reflections ↔ Anti-automorphisms (self-complementary)      ║
║     └── Paley = full D_p symmetry                                  ║
║                                                                    ║
║  HELIX (Riemann surface of log)                                    ║
║     │                                                              ║
║     ├── Complex eigenvalues trace spirals under M^n                ║
║     ├── On covering space, spirals = straight lines                ║
║     ├── Dominant eigenvalue = flat growth (Pisot root)             ║
║     └── Subdominant = helical correction (complex conjugates)      ║
║                                                                    ║
║  k-NACCI FAMILY                                                    ║
║     │                                                              ║
║     ├── k=2: Fibonacci, φ, no adjacent 1s, 1 helix (degenerate)   ║
║     ├── k=3: Tribonacci, τ₃, no 3 consecutive 1s, 1 helix        ║
║     ├── k=4: Tetranacci, τ₄, no 4 consecutive 1s, 2 helices      ║
║     ├── ...                                                        ║
║     └── k→∞: τ_∞ = 2, all binary strings, helices → unit circle   ║
║                                                                    ║
║  ZECKENDORF REPRESENTATION                                         ║
║     │                                                              ║
║     ├── Base φ: digits {0,1}, no 11, uniqueness = Zeckendorf      ║
║     ├── Base τ₃: digits {0,1}, no 111, tribonacci Zeckendorf      ║
║     ├── 131 ≈ τ₃^8 (OPEN-Q-029), 504 = T(13)                     ║
║     └── Transfer matrix enforces the digit constraint              ║
║                                                                    ║
║  SIMPLICIAL RÉDEI (THM-220)                                        ║
║     │                                                              ║
║     ├── sim_H ∈ {0,1} for all tournaments (n ≥ 4)                 ║
║     ├── 2·n! have sim_H = 1 (n! transitive + n! near-transitive)  ║
║     ├── Near-transitive: single reversed arc, H = 2^{n-2}+1       ║
║     ├── Monodromy group (Z/2)^{n-2} → covering space structure    ║
║     └── Transitive core = universal cover (remove 3-cycle loops)   ║
║                                                                    ║
║  THE CIRCLE CLOSES:                                                ║
║     Polygon → DFT → Eigenvalues → Helix → k-nacci →               ║
║     Zeckendorf → Transfer matrix → Simplicial → Tournament →      ║
║     Polygon                                                        ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝
""")

print("opus-2026-03-15-S90: THE HELICAL SIMPLICIAL RÉDEI EXPLORATION")
