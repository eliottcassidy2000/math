#!/usr/bin/env python3
"""
dihedral_tournament_geometry.py — Dihedral group structure of tournament H-maximization

KEY GEOMETRIC INSIGHT (from user):
  - Tournament on n vertices = oriented n-gon (vertices on unit circle)
  - D_{2n} = symmetries of the regular n-gon (|D_{2n}| = 2n, always even)
  - Tournaments live at ODD n; dihedral groups at EVEN order 2n
  - Groups with ODD order (C_3, C_5, C_7, ...) interlace between tournament sizes

REPRESENTATION-THEORETIC FRAMEWORK:
  For a circulant tournament T(S) on Z_p (p odd prime):

  1. CHORD DECOMPOSITION: The p-gon has (p-1)/2 chord types.
     Chord type k connects vertices distance k apart (k=1,...,(p-1)/2).
     D_{2p} acts transitively on each chord type.

  2. ORIENTATION CHOICE: A circulant tournament assigns each chord type
     an "orientation" — clockwise (k in S) or counterclockwise (k not in S).
     This is an element of {+1,-1}^{(p-1)/2}.

  3. IRREP DECOMPOSITION: Under D_{2p}, the adjacency matrix decomposes:
     - Trivial irrep: eigenvalue lambda_0 = m = (p-1)/2
     - 2D irreps rho_k (k=1,...,(p-1)/2): eigenvalue pair (lambda_k, conj(lambda_k))

     The AMPLITUDE |lambda_k| and PHASE arg(lambda_k) encode the chord mode.

  4. PARSEVAL CONSTRAINT: sum_{k=1}^{m} |lambda_k|^2 = m*p/4

  5. H-MAXIMIZATION = finding the optimal energy distribution across chord modes

For Paley: all |lambda_k| = sqrt(p)/2 (uniform energy distribution)
For Interval: |lambda_1| ~ p/pi >> |lambda_2| >> ... (concentrated energy)

The CROSSOVER at p=19 occurs because concentrated energy beats uniform energy
for Hamiltonian path counting once p is large enough.

Author: opus-2026-03-12-S62
"""

import numpy as np
import cmath
import math
from itertools import combinations

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

def eigenvalues_circulant(S, p):
    """Compute eigenvalues of circulant tournament with connection set S."""
    omega = cmath.exp(2j * cmath.pi / p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

def chord_orientation_vector(S, p):
    """Return the orientation vector sigma in {+1,-1}^m.
    sigma_k = +1 if k in S (clockwise), -1 if k not in S (counterclockwise).
    """
    m = (p - 1) // 2
    return [1 if k in S else -1 for k in range(1, m+1)]

def chord_mode_spectrum(S, p):
    """Decompose into (amplitude, phase) pairs for each chord mode.

    For chord mode k: the 2D irrep rho_k has
      amplitude = |lambda_k|
      phase = arg(lambda_k)
    """
    eigs = eigenvalues_circulant(S, p)
    m = (p - 1) // 2
    result = []
    for k in range(1, m+1):
        lam = eigs[k]
        result.append({
            'k': k,
            'amplitude': abs(lam),
            'phase': cmath.phase(lam),
            'lambda': lam,
            'lambda_conj': eigs[p-k]  # = conj(lambda_k) for real A
        })
    return result

def geometric_chord_length(k, n):
    """Physical chord length on unit circle: 2*sin(pi*k/n)."""
    return 2 * math.sin(math.pi * k / n)


# ============================================================
# SECTION 1: CHORD DECOMPOSITION AND ORIENTATION VECTORS
# ============================================================

print("=" * 70)
print("DIHEDRAL TOURNAMENT GEOMETRY")
print("=" * 70)

print("\n--- SECTION 1: CHORD TYPES AND ORIENTATIONS ---")
print()

for p in [7, 11, 19]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_interval = frozenset(range(1, m+1))

    sigma_paley = chord_orientation_vector(S_paley, p)
    sigma_interval = chord_orientation_vector(S_interval, p)

    print(f"p={p} ({m} chord types):")
    print(f"  Chord types:  k = {list(range(1, m+1))}")
    print(f"  Chord lengths: {['%.3f' % geometric_chord_length(k, p) for k in range(1, m+1)]}")
    print(f"  Paley QR:     S = {sorted(S_paley)}")
    print(f"  Paley sigma:  {sigma_paley}")
    print(f"  Interval:     S = {sorted(S_interval)}")
    print(f"  Interval sigma: {sigma_interval}")

    # Inner product of orientation vectors
    dot = sum(a*b for a, b in zip(sigma_paley, sigma_interval))
    print(f"  <sigma_P, sigma_I> = {dot} (of {m})")
    print(f"  Agreement fraction: {(dot + m) / (2*m):.4f}")

    # Hamming distance
    hamming = sum(1 for a, b in zip(sigma_paley, sigma_interval) if a != b)
    print(f"  Hamming distance: {hamming}")
    print()


# ============================================================
# SECTION 2: IRREP AMPLITUDE SPECTRA
# ============================================================

print("--- SECTION 2: CHORD MODE SPECTRA (D_{2p} irreps) ---")
print()

for p in [7, 11, 19]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_interval = frozenset(range(1, m+1))

    spec_P = chord_mode_spectrum(S_paley, p)
    spec_I = chord_mode_spectrum(S_interval, p)

    print(f"p={p}:")
    print(f"  {'k':>3} {'chord_len':>10} {'|lam|_P':>10} {'|lam|_I':>10} {'ratio':>8} "
          f"{'phase_P/pi':>11} {'phase_I/pi':>11}")

    total_E_P = sum(s['amplitude']**2 for s in spec_P)
    total_E_I = sum(s['amplitude']**2 for s in spec_I)

    for sP, sI in zip(spec_P, spec_I):
        k = sP['k']
        chord = geometric_chord_length(k, p)
        ratio = sI['amplitude'] / sP['amplitude'] if sP['amplitude'] > 0 else float('inf')
        print(f"  {k:>3} {chord:>10.4f} {sP['amplitude']:>10.4f} {sI['amplitude']:>10.4f} "
              f"{ratio:>8.4f} {sP['phase']/math.pi:>11.4f} {sI['phase']/math.pi:>11.4f}")

    print(f"  Total energy: Paley={total_E_P:.4f}, Interval={total_E_I:.4f}, "
          f"expected={m*p/4:.4f}")

    # Gini coefficient of amplitude distribution
    amps_P = sorted(s['amplitude'] for s in spec_P)
    amps_I = sorted(s['amplitude'] for s in spec_I)

    def gini(values):
        n = len(values)
        total = sum(values)
        if total == 0: return 0
        cumsum = 0
        area = 0
        for i, v in enumerate(sorted(values)):
            cumsum += v
            area += cumsum
        return 1 - 2 * area / (n * total) + 1/n

    print(f"  Gini(|lam|): Paley={gini(amps_P):.4f}, Interval={gini(amps_I):.4f}")

    # Fraction of energy in mode 1
    frac1_P = spec_P[0]['amplitude']**2 / total_E_P
    frac1_I = spec_I[0]['amplitude']**2 / total_E_I
    print(f"  Energy in mode 1: Paley={frac1_P:.4f}, Interval={frac1_I:.4f}")
    print()


# ============================================================
# SECTION 3: GEOMETRIC MEANING — POLYGON WINDING
# ============================================================

print("--- SECTION 3: POLYGON WINDING AND DIRECTED ARCS ---")
print()

print("""
GEOMETRIC PICTURE:
  Place vertices v_0, v_1, ..., v_{p-1} at p-th roots of unity on the unit circle.

  A CIRCULANT tournament with connection set S orients each edge:
    v_j -> v_k  iff  (k-j) mod p in S

  CHORD TYPE k: all edges {v_j, v_{j+k}} for j = 0,...,p-1.
  The ORIENTATION of chord type k:
    - If k in S: v_j -> v_{j+k} (clockwise jump of k positions)
    - If k not in S: v_{j+k} -> v_j (counterclockwise jump)

  INTERVAL TOURNAMENT S={1,...,m}:
    ALL chords oriented CLOCKWISE. Every vertex beats its m
    nearest clockwise neighbors. This is a "one-way flow" around the circle.

    Geometrically: the tournament is a FLOW on the polygon,
    where material moves consistently in one direction.

  PALEY TOURNAMENT S=QR_p:
    Chord orientations follow the Legendre symbol chi(k).
    Short chords (k=1) and scattered longer chords go clockwise,
    while other chords go counterclockwise.

    Geometrically: a TURBULENT flow with both clockwise and
    counterclockwise components at different length scales.
""")

# Visualize the orientation pattern
for p in [7, 11, 19]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))

    print(f"p={p}: Orientation pattern by chord length:")
    for k in range(1, m+1):
        chord_len = geometric_chord_length(k, p)
        direction_P = "CW" if k in S_paley else "CCW"
        direction_I = "CW"  # always clockwise for interval
        print(f"  k={k:2d} (len={chord_len:.3f}): Paley={direction_P:>3}, Interval={direction_I:>3}")
    print()


# ============================================================
# SECTION 4: DIHEDRAL GROUP ACTION ON TOURNAMENTS
# ============================================================

print("--- SECTION 4: D_{2p} ACTION ON TOURNAMENTS ---")
print()

print("""
The DIHEDRAL GROUP D_{2p} has 2p elements:
  - Rotations: r^0, r^1, ..., r^{p-1}  (cyclic subgroup C_p)
  - Reflections: s, sr, sr^2, ..., sr^{p-1}

Action on tournaments T(S) on Z_p:
  1. ROTATION r: v_j -> v_{j+1 mod p}.
     This maps T(S) to T(S) (circulants are rotation-invariant).

  2. REFLECTION s: v_j -> v_{-j mod p}.
     This maps T(S) to T(-S) where -S = {p-k : k in S}.
     For tournaments, -S = S^c (complement), so s maps T to T^op.

  KEY: No circulant tournament can be D_{2p}-invariant,
  because reflection always gives T^op (the reverse tournament).

  But D_{2p} acts on the SPACE of all circulant tournaments!
  The rotation subgroup C_p acts trivially (preserves all circulants).
  The reflection s acts as T -> T^op, i.e., sigma -> -sigma.

  So the relevant symmetry for H-maximization is:
  H(T) = H(T^op) (since reversing all arcs preserves HP count)
  This means H(sigma) = H(-sigma) for orientation vectors.

  The D_{2p} action on orientation space {+1,-1}^m:
    Rotation: trivial (doesn't change sigma)
    Reflection: sigma -> -sigma

  Orbits: {sigma, -sigma} pairs. The orbit of sigma_interval = (+1,...,+1)
  is {(+1,...,+1), (-1,...,-1)} = {interval, reverse interval}.
  Both have the same H (since H(T) = H(T^op) for odd n).
""")

# Verify H(T) = H(T^op) for circulants at small p
for p in [7, 11]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_paley_comp = frozenset(range(1, p)) - S_paley
    S_interval = frozenset(range(1, m+1))
    S_interval_comp = frozenset(range(m+1, p))

    eP = eigenvalues_circulant(S_paley, p)
    ePc = eigenvalues_circulant(S_paley_comp, p)
    eI = eigenvalues_circulant(S_interval, p)
    eIc = eigenvalues_circulant(S_interval_comp, p)

    print(f"p={p}: sigma -> -sigma preserves |lambda|:")
    for k in range(1, m+1):
        print(f"  k={k}: |lam_P|={abs(eP[k]):.4f} vs |lam_Pc|={abs(ePc[k]):.4f}  "
              f"|lam_I|={abs(eI[k]):.4f} vs |lam_Ic|={abs(eIc[k]):.4f}")
    print()


# ============================================================
# SECTION 5: INTERLACING WITH ODD-ORDER GROUPS
# ============================================================

print("--- SECTION 5: ODD-ORDER GROUPS AND INTERLACING ---")
print()

print("""
The user's key observation: dihedral groups D_{2n} have EVEN order 2n.
Groups with ODD order (C_1, C_3, C_5, C_7, ...) are INTERLACED
between tournament sizes.

Tournament sizes:     n = 3   4   5   6   7   8   9   ...
D_{2n} order:            6   8  10  12  14  16  18   ...
Interlaced C_k:       3   5   7   9  11  13  15  17   ...

The ODD-ORDER cyclic groups C_k act on the tournament vertex set
differently from D_{2n}:
  - C_n acts as the rotation subgroup of D_{2n}
  - C_{n+1} cannot act on an n-gon at all!

But there IS a deep connection: the CIRCULANT structure.

For a circulant tournament on Z_p, the rotation group C_p acts on vertices.
The EIGENVALUES of the adjacency matrix are the DFT coefficients,
which live in the GROUP ALGEBRA C[C_p].

The representation ring of C_p: all irreps are 1-dimensional,
  rho_k: j -> omega^{jk} for k = 0, ..., p-1.

Under D_{2p}, these pair up into 2D irreps:
  rho_k + rho_{p-k} for k = 1, ..., (p-1)/2.

The PALEY tournament's special property: its eigenvalues are
determined by the GAUSS SUM, which is a character sum over
the multiplicative group (Z/pZ)^*.

(Z/pZ)^* is cyclic of order p-1. When p = 3 mod 4, p-1 = 2 mod 4.
The QR subgroup has index 2, order (p-1)/2.

So the QR structure involves:
  C_{(p-1)/2} < C_{p-1} = (Z/pZ)^*

This is an ODD-order group C_{(p-1)/2} acting on the multiplicative side!
The interlacing the user notices is between:
  - Additive structure: C_p acts on vertices (order p, odd)
  - Multiplicative structure: C_{(p-1)/2} acts on QR (order (p-1)/2)

The interplay between these two cyclic actions (additive and multiplicative)
is what creates the Gauss sum and determines the Paley spectrum.
""")

# Display the additive/multiplicative group orders
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    m = (p - 1) // 2
    pm4 = p % 4
    print(f"  p={p:2d}: |C_p|={p}, |(Z/pZ)*|={p-1}, |QR|={m}, "
          f"p mod 4 = {pm4}, D_{{2p}} order = {2*p}")


# ============================================================
# SECTION 6: THE POLYGON FLOW PICTURE
# ============================================================

print(f"\n--- SECTION 6: POLYGON FLOW AND H-MAXIMIZATION ---")
print()

print("""
FLOW ON THE n-GON:
  Think of a circulant tournament as a FLOW on the regular n-gon.
  Each chord type k contributes a "flow mode" at angular frequency k.

  The total flow field at position theta on the circle:
    F(theta) = sum_{k=1}^{m} sigma_k * f_k(theta)
  where f_k is the k-th harmonic flow mode and sigma_k = +/-1.

  INTERVAL tournament: F = sum_{k=1}^m f_k (all modes constructive)
    -> MAXIMUM DIRECTIONAL BIAS (strongest "wind")
    -> Creates a dominant "highway" for Hamiltonian paths

  PALEY tournament: F = sum chi(k) * f_k (mixed signs)
    -> Each mode has EQUAL energy but phases partially cancel
    -> More "turbulent" but more "isotropic"

H-MAXIMIZATION TRADEOFF:
  - Directional bias (interval): good for LONG paths (the HP count)
    because paths can ride the "wind" all the way around
  - Isotropy (Paley): good for SHORT cycles (k=5 structure)
    because 5-cycles can form in any orientation

  At small p: there aren't many chord modes, so isotropy helps more.
  At large p: the "wind" effect compounds over many steps, dominating.
""")

# Compute the "directional bias" for each tournament type
for p in [7, 11, 19, 23]:
    if p % 4 != 3:
        continue
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_interval = frozenset(range(1, m+1))

    # The "bias" = |sum of eigenvalues| = |sum sigma_k * lambda_k(interval)|
    eigs_I = eigenvalues_circulant(S_interval, p)
    eigs_P = eigenvalues_circulant(S_paley, p)

    # Directional bias: magnitude of lambda_1 (dominant mode)
    bias_I = abs(eigs_I[1])
    bias_P = abs(eigs_P[1])

    # Total "flow strength": sum of |lambda_k|
    flow_I = sum(abs(eigs_I[k]) for k in range(1, p))
    flow_P = sum(abs(eigs_P[k]) for k in range(1, p))

    # "Peakedness": ratio of max |lambda| to rms |lambda|
    rms_I = math.sqrt(sum(abs(eigs_I[k])**2 for k in range(1, p)) / (p-1))
    rms_P = math.sqrt(sum(abs(eigs_P[k])**2 for k in range(1, p)) / (p-1))
    peak_I = max(abs(eigs_I[k]) for k in range(1, p)) / rms_I
    peak_P = max(abs(eigs_P[k]) for k in range(1, p)) / rms_P

    print(f"p={p}:")
    print(f"  Dominant mode |lam_1|: Interval={bias_I:.4f}, Paley={bias_P:.4f}")
    print(f"  Total flow:           Interval={flow_I:.4f}, Paley={flow_P:.4f}")
    print(f"  Peakedness (max/rms): Interval={peak_I:.4f}, Paley={peak_P:.4f} (Paley always 1.0)")
    print(f"  Dirichlet peak ratio: |lam_1|/sqrt(p)/2 = {bias_I/(math.sqrt(p)/2):.4f}")
    print()


# ============================================================
# SECTION 7: CONNECTION TO SPECIFIC DIHEDRAL REPRESENTATIONS
# ============================================================

print("--- SECTION 7: D_{2p} CHARACTER TABLE AND TOURNAMENT MODES ---")
print()

print("""
D_{2n} irreducible representations (n odd):
  - 1D trivial: chi_0(r^k) = 1, chi_0(s*r^k) = 1
  - 1D sign:    chi_1(r^k) = 1, chi_1(s*r^k) = -1
  - 2D rho_j (j=1,...,(n-1)/2):
      chi_j(r^k) = 2*cos(2*pi*j*k/n)
      chi_j(s*r^k) = 0

The adjacency matrix A of a circulant tournament decomposes:
  A = (m/p)*J + sum_{j=1}^{m} [lambda_j P_j + conj(lambda_j) P_{p-j}]

where P_j is the projection onto the j-th Fourier mode.

Under D_{2p}, the pair (P_j, P_{p-j}) transforms as the 2D irrep rho_j.
The eigenvalue PAIR (lambda_j, conj(lambda_j)) is the "D_{2p} eigenvalue"
of the tournament in representation rho_j.

The REFLECTION s maps lambda_j -> conj(lambda_j) (complex conjugation).
For the interval: lambda_j has phase ~ pi*j/p (Dirichlet kernel).
For Paley: lambda_j has phase ~ +/- theta_0 (two possible phases).

D_{2p} EQUIVARIANCE OF H:
  H is invariant under the FULL D_{2p} action (since H(T) = H(T^op)).
  So H is a function on D_{2p}-orbits of circulant tournaments.
  These orbits are parametrized by:
    (|lambda_1|, |lambda_2|, ..., |lambda_m|)
  since the reflection fixes |lambda_j| but conjugates the phase.

  But H also depends on the PHASES! So the orbit structure alone
  doesn't determine H. The phases encode the "chirality" of the flow.
""")

# Compute the character decomposition for specific tournaments
for p in [7, 11]:
    m = (p - 1) // 2
    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_interval = frozenset(range(1, m+1))

    eigs_P = eigenvalues_circulant(S_paley, p)
    eigs_I = eigenvalues_circulant(S_interval, p)

    print(f"p={p}: D_{2*p} character decomposition")
    print(f"  {'mode':>5} {'|lam|_P':>10} {'|lam|_I':>10} {'phase_P/pi':>12} {'phase_I/pi':>12}")
    for j in range(1, m+1):
        print(f"  rho_{j:>2} {abs(eigs_P[j]):>10.4f} {abs(eigs_I[j]):>10.4f} "
              f"{cmath.phase(eigs_P[j])/math.pi:>12.4f} {cmath.phase(eigs_I[j])/math.pi:>12.4f}")
    print()


# ============================================================
# SECTION 8: SCANNING ALL ORIENTATION VECTORS AT p=7
# ============================================================

print("--- SECTION 8: EXHAUSTIVE SCAN AT p=7 (2^3 = 8 orientations) ---")
print()

p = 7
m = 3

# Count HP for each orientation vector
# An orientation vector sigma in {+1,-1}^3 determines S:
#   S = {k : sigma_k = +1}
# Actually we also need the complement: if k not in S, then p-k must be in S.
# For a valid tournament: S and -S partition {1,...,p-1}.
# sigma_k = +1 means k in S and p-k not in S.

def sigma_to_S(sigma, p):
    m = (p - 1) // 2
    S = set()
    for k in range(1, m+1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return frozenset(S)

# For each of 2^m = 8 orientations, compute eigenvalues and estimate H
# For p=7, we can compute H exactly via Held-Karp
def count_hp_exact(A):
    """Held-Karp DP for Hamiltonian path count."""
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c
    return int(np.sum(dp[full]))

def adjacency_matrix(S, p):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S:
                A[i][j] = 1
    return A

print(f"p=7: All 2^{m} = {2**m} orientation vectors")
print(f"{'sigma':>12} {'S':>20} {'H':>8} {'note':>20}")

results_7 = []
S_paley_7 = frozenset(j for j in range(1, p) if is_qr(j, p))

for bits in range(2**m):
    sigma = [1 if (bits >> i) & 1 else -1 for i in range(m)]
    S = sigma_to_S(sigma, p)
    A = adjacency_matrix(S, p)
    H = count_hp_exact(A)

    note = ""
    if S == S_paley_7:
        note = "PALEY"
    elif S == frozenset(range(1, p)) - S_paley_7:
        note = "PALEY complement"
    elif S == frozenset(range(1, m+1)):
        note = "INTERVAL"
    elif S == frozenset(range(m+1, p)):
        note = "INTERVAL complement"

    results_7.append((sigma, S, H, note))
    print(f"  {str(sigma):>12} {str(sorted(S)):>20} {H:>8} {note:>20}")

# Sort by H
results_7.sort(key=lambda x: -x[2])
print(f"\nRanked by H:")
for rank, (sigma, S, H, note) in enumerate(results_7, 1):
    print(f"  #{rank}: H={H}, S={sorted(S)} {note}")


# ============================================================
# SECTION 9: PHASE COHERENCE AND POLYGON WINDING NUMBER
# ============================================================

print(f"\n--- SECTION 9: WINDING NUMBER INTERPRETATION ---")
print()

print("""
A Hamiltonian PATH on the p-gon visits all vertices in some order.
The PATH can be characterized by its WINDING NUMBER around the circle:
  w(P) = (1/2pi) * sum of angular displacements along the path.

For the interval tournament (all clockwise):
  Most Hamiltonian paths have HIGH winding number (they go mostly
  in one direction around the circle).

For the Paley tournament (mixed directions):
  Paths have LOWER average winding number (forced to change direction).

CONJECTURE: H-maximization at large p is equivalent to maximizing
the number of paths with high winding number, which is achieved by
the maximally coherent orientation (the interval tournament).
""")

# For p=7, compute the "net angular displacement" of each Hamiltonian path
# and see if high-winding paths are more common in interval vs Paley

def enumerate_hp(A, p):
    """Enumerate all Hamiltonian paths and their step sequences."""
    n = p
    paths = []

    def backtrack(path, visited):
        if len(path) == n:
            paths.append(list(path))
            return
        v = path[-1]
        for u in range(n):
            if u not in visited and A[v][u]:
                visited.add(u)
                path.append(u)
                backtrack(path, visited)
                path.pop()
                visited.remove(u)

    for start in range(n):
        backtrack([start], {start})
    return paths

def path_winding(path, p):
    """Compute the net angular displacement of a path on the p-gon.
    Each step v_i -> v_{i+1} has angular displacement 2*pi*(v_{i+1}-v_i)/p.
    Normalize by p to get a "step direction" in {-(p-1)/2,...,(p-1)/2}.
    """
    m = (p - 1) // 2
    total = 0
    for i in range(len(path) - 1):
        diff = (path[i+1] - path[i]) % p
        if diff > m:
            diff -= p  # normalize to [-m, m]
        total += diff
    return total

p = 7
m = 3
S_paley_7 = frozenset(j for j in range(1, p) if is_qr(j, p))
S_interval_7 = frozenset(range(1, m+1))

A_paley = adjacency_matrix(S_paley_7, p)
A_interval = adjacency_matrix(S_interval_7, p)

paths_P = enumerate_hp(A_paley, p)
paths_I = enumerate_hp(A_interval, p)

windings_P = [path_winding(path, p) for path in paths_P]
windings_I = [path_winding(path, p) for path in paths_I]

print(f"p={p}: Winding number distribution")
print(f"  Paley:    H={len(paths_P)}, mean_winding={np.mean(windings_P):.2f}, "
      f"std={np.std(windings_P):.2f}")
print(f"  Interval: H={len(paths_I)}, mean_winding={np.mean(windings_I):.2f}, "
      f"std={np.std(windings_I):.2f}")

# Histogram of winding numbers
from collections import Counter
wind_hist_P = Counter(windings_P)
wind_hist_I = Counter(windings_I)
all_winds = sorted(set(list(wind_hist_P.keys()) + list(wind_hist_I.keys())))

print(f"\n  {'winding':>8} {'Paley':>8} {'Interval':>8}")
for w in all_winds:
    print(f"  {w:>8} {wind_hist_P.get(w,0):>8} {wind_hist_I.get(w,0):>8}")

# Also check p=5
print()
p5 = 5
m5 = 2
S_p5 = frozenset(j for j in range(1, p5) if is_qr(j, p5))
S_i5 = frozenset(range(1, m5+1))
A_p5 = adjacency_matrix(S_p5, p5)
A_i5 = adjacency_matrix(S_i5, p5)
paths_P5 = enumerate_hp(A_p5, p5)
paths_I5 = enumerate_hp(A_i5, p5)
wind_P5 = [path_winding(p_path, p5) for p_path in paths_P5]
wind_I5 = [path_winding(p_path, p5) for p_path in paths_I5]
print(f"p={p5}: Paley H={len(paths_P5)}, winding mean={np.mean(wind_P5):.2f}")
print(f"p={p5}: Interval H={len(paths_I5)}, winding mean={np.mean(wind_I5):.2f}")


# ============================================================
# SECTION 10: DEEP CONNECTION — QR AS DIHEDRAL EIGENVECTOR
# ============================================================

print(f"\n--- SECTION 10: QR STRUCTURE IN DIHEDRAL FRAMEWORK ---")
print()

print("""
DEEP CONNECTION: The Paley tournament's connection set S = QR_p
is the UNIQUE set (up to complement) that produces a FLAT eigenvalue
spectrum: |lambda_k| = sqrt(p)/2 for all k.

In the D_{2p} representation theory:
  Paley = ISOTROPIC tournament (equal energy in all chord modes)
  Interval = ANISOTROPIC tournament (energy concentrated in mode 1)

This is like a crystal vs. an amorphous material:
  - Crystal (Paley): uniform diffraction pattern (flat spectrum)
  - Amorphous (Interval): peaked diffraction pattern (Dirichlet kernel)

H-MAXIMIZATION:
  At small p (few modes): isotropy wins, like how a perfect crystal
  has special electronic properties.

  At large p (many modes): anisotropy wins, like how directing all
  energy into one beam is better than spreading it uniformly.

THE CROSSOVER is when the "dimensionless coupling"
  g = 2*sqrt(p)/pi
exceeds a critical threshold:
  g < g_c: Paley (isotropic) wins
  g > g_c: Interval (anisotropic) wins

From the data: g_c is between g(11) = 2.11 and g(19) = 2.78.
""")

# Compute the dimensionless coupling for each prime
print(f"Dimensionless coupling g = 2*sqrt(p)/pi:")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    g = 2 * math.sqrt(p) / math.pi
    pm4 = p % 4
    paley = "Paley exists" if pm4 == 3 else "no Paley"
    print(f"  p={p:2d}: g={g:.4f}, p mod 4 = {pm4} ({paley})")


# ============================================================
# SECTION 11: THE CRITICAL QUESTION — DOES p=13 HAVE A PALEY?
# ============================================================

print(f"\n--- SECTION 11: CROSSOVER ANALYSIS INCLUDING p=1 mod 4 ---")
print()

print("""
For p = 1 mod 4 (p = 5, 13, 17, 29, ...):
  There is NO Paley tournament (since -1 is a QR, so QR_p is symmetric
  and doesn't define a valid tournament).

  But the INTERVAL tournament still exists for any odd p!
  So the question becomes: among ALL circulant tournaments on Z_p,
  which maximizes H?

For p = 3 mod 4: Paley vs Interval is the main competition.
For p = 1 mod 4: Interval may be the universal winner.

The dihedral framework suggests:
  p = 1 mod 4: -1 is QR, so reflection (j -> -j) maps QR to QR.
    The QR set is "palindromic" (symmetric under reflection).
    No valid tournament can use QR as connection set.

  p = 3 mod 4: -1 is NQR, so reflection maps QR to NQR.
    QR and NQR are "chiral" — they differ by a reflection.
    The Paley tournament T_p has the MAXIMAL symmetry group
    among circulant tournaments (size p*(p-1)/2 instead of just p).
""")

# Summarize what we know
print("CROSSOVER SUMMARY:")
print(f"{'p':>4} {'p mod 4':>8} {'g=2sqrt(p)/pi':>14} {'winner':>10}")
known_winners = {
    3: "TIE", 5: "TIE", 7: "Paley", 11: "Paley",
    13: "?", 17: "?", 19: "Interval", 23: "Interval?"
}
for p in [3, 5, 7, 11, 13, 17, 19, 23]:
    g = 2 * math.sqrt(p) / math.pi
    winner = known_winners.get(p, "?")
    print(f"  {p:>3} {p%4:>8} {g:>14.4f} {winner:>10}")


# ============================================================
# SECTION 12: FOURIER ANALYSIS OF H AS FUNCTION ON {+1,-1}^m
# ============================================================

print(f"\n--- SECTION 12: H AS BOOLEAN FUNCTION ON ORIENTATION CUBE ---")
print()

print("""
H is a function on the hypercube {+1,-1}^m where m = (p-1)/2.
Each coordinate sigma_k encodes the orientation of chord type k.

As a function on {+1,-1}^m, H has a FOURIER-WALSH expansion:
  H(sigma) = sum_{S subset [m]} hat{H}(S) * prod_{k in S} sigma_k

This is different from the Walsh expansion of H over tournament space!
Here sigma lives in the ORIENTATION CUBE (dimension m),
not the TOURNAMENT CUBE (dimension C(p,2)).

The D_{2p} symmetry H(sigma) = H(-sigma) means:
  hat{H}(S) = 0 for |S| odd
  Only EVEN-degree Walsh coefficients survive.

Let's compute this for p=7 (m=3, 8 orientations).
""")

# Compute the Walsh-Fourier expansion of H on {+1,-1}^3 at p=7
p = 7
m = 3

# Build H values for all 8 orientations
H_vals = {}
for bits in range(2**m):
    sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
    S = sigma_to_S(list(sigma), p)
    A = adjacency_matrix(S, p)
    H = count_hp_exact(A)
    H_vals[sigma] = H

# Walsh-Fourier transform
print(f"p=7: Walsh-Fourier expansion of H on orientation cube")
from itertools import chain

def subsets(n):
    """All subsets of {0,...,n-1}."""
    for r in range(n+1):
        for combo in combinations(range(n), r):
            yield combo

walsh_coeffs = {}
for S in subsets(m):
    coeff = 0
    for bits in range(2**m):
        sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
        prod = 1
        for k in S:
            prod *= sigma[k]
        coeff += H_vals[sigma] * prod
    coeff /= 2**m
    walsh_coeffs[S] = coeff
    if abs(coeff) > 0.001:
        degree = len(S)
        print(f"  hat{{H}}({S}) = {coeff:.4f}  (degree {degree})")

print(f"\n  Verification: H(sigma) = H(-sigma)?")
for bits in range(2**(m-1)):
    sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
    neg_sigma = tuple(-s for s in sigma)
    print(f"    H{sigma} = {H_vals[sigma]}, H{neg_sigma} = {H_vals[neg_sigma]}, "
          f"equal? {H_vals[sigma] == H_vals[neg_sigma]}")


# Final theoretical summary
print(f"\n{'=' * 70}")
print(f"THEORETICAL SYNTHESIS")
print(f"{'=' * 70}")
print("""
THE DIHEDRAL TOURNAMENT FRAMEWORK:

1. CHORD DECOMPOSITION: A circulant tournament on Z_p decomposes into
   m = (p-1)/2 chord modes, each either CW or CCW oriented.
   This is an element sigma of the orientation cube {+1,-1}^m.

2. DIHEDRAL ACTION: D_{2p} acts on the orientation cube:
   - Rotations: trivial (sigma -> sigma)
   - Reflection: sigma -> -sigma (reverse all orientations = T^op)
   H is invariant under this action.

3. IRREP SPECTRUM: Each chord mode k carries energy |lambda_k|^2,
   with Parseval constraint sum |lambda_k|^2 = mp/4.
   Paley = uniform distribution; Interval = concentrated in mode 1.

4. CROSSOVER CRITERION: The dimensionless coupling g = 2*sqrt(p)/pi
   determines whether isotropy (Paley) or concentration (interval) wins.
   Crossover at g ~ 2.4 (between p=11 and p=19).

5. TRACE ALTERNATION: Paley wins at k = 1 mod 4 (short-cycle advantage),
   Interval wins at k = 3 mod 4 (long-cycle advantage).
   The balance tips toward interval at large p.

6. BOOLEAN FUNCTION: H on {+1,-1}^m has only even-degree Walsh coefficients
   (consequence of D_{2p} symmetry H(sigma) = H(-sigma)).
   The degree-0 coefficient is the MEAN H over all orientations.
   The degree-2 coefficients encode pairwise chord interactions.

OPEN QUESTIONS:
  - Is the crossover monotone? (Once interval beats Paley, does it always win?)
  - What is the exact critical g_c?
  - Can this framework explain the Satake NDRT results?
  - Does the boolean function structure reveal why c_3 is always equal?
""")
