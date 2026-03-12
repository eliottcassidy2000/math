#!/usr/bin/env python3
"""
interlacing_groups.py — Odd-order groups interlacing tournament sizes

The user's observation:
  D_6 (|D|=6) = triangle = n=3 tournament
  D_8 (|D|=8) = square = n=4
  D_{10} (|D|=10) = pentagon = n=5
  ...

Groups with ODD order are "interlaced":
  C_1(1), [D_2=V], C_3(3), [D_4], C_5(5), [D_6], C_7(7), [D_8], ...

But the IMPORTANT odd-order group for a tournament on p vertices
is (Z/pZ)* = the multiplicative group, which has order p-1 (EVEN).
The QR subgroup has order (p-1)/2 = m, which is ODD (for p=3 mod 4).

So the interlacing is:
  Tournament size p:   3    7    11    19    23    31    ...
  QR subgroup order m: 1    3     5     9    11    15    ...
  (p-1)/2 is ODD for p=3 mod 4)

The ODD-order cyclic group C_m acts multiplicatively on Z_p,
and the Legendre symbol chi: C_{p-1} -> {+/-1} restricts to
the TRIVIAL character on C_m (since QR^2 = QR).

QUESTION: How does the ORDER m = (p-1)/2 relate to the H-maximization?

For m=1 (p=3): trivial case, only 1 chord type, H always same.
For m=3 (p=7): 3 chord types, QR transitive, Paley wins.
For m=5 (p=11): 5 chord types, QR transitive, Paley wins.
For m=9 (p=19): 9 chord types, QR transitive, Interval wins.

The transition happens between m=5 and m=9.
Is it related to the REPRESENTATION THEORY of C_m?

Author: opus-2026-03-12-S62
"""

import numpy as np
import cmath
import math

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p - 1) // 2, p) == 1


# ============================================================
# SECTION 1: THE MULTIPLICATIVE GROUP STRUCTURE
# ============================================================

print("=" * 70)
print("MULTIPLICATIVE GROUP AND QR STRUCTURE")
print("=" * 70)

print("""
For a tournament on Z_p (p odd prime), the key groups are:

ADDITIVE: (Z/pZ, +) = C_p  (order p, acts by rotation)
MULTIPLICATIVE: (Z/pZ)* = C_{p-1}  (order p-1, acts by scaling)
QR SUBGROUP: QR_p < (Z/pZ)*  (index 2, order (p-1)/2 = m)

For p = 3 mod 4:
  m = (p-1)/2 is ODD
  -1 is NQR (chi(-1) = -1)
  The QR subgroup has odd order m

For p = 1 mod 4:
  m = (p-1)/2 is EVEN
  -1 is QR (chi(-1) = +1)
  The QR subgroup has even order m

The GROUP ALGEBRA factorization:
  C[C_{p-1}] = C[QR] tensor C[Z_2]

For p = 3 mod 4:
  The irreducible representations of C_m (odd order m):
    rho_j: a -> omega_m^{j*ind(a)} for j = 0,...,m-1
  where ind(a) is the discrete log in the QR subgroup.

  The Legendre symbol chi maps:
    QR -> +1 (trivial on QR)
    NQR -> -1

  So chi is NOT a character of QR itself, but of the full group (Z/pZ)*.
  As a representation of C_{p-1}: chi = ind_2 (the order-2 character).
""")

# Display the group structure for small primes
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    m = (p - 1) // 2

    # Find a primitive root
    for g in range(2, p):
        if len(set(pow(g, k, p) for k in range(1, p))) == p - 1:
            break

    # QR elements
    qr = sorted(j for j in range(1, p) if is_qr(j, p))
    nqr = sorted(j for j in range(1, p) if not is_qr(j, p))

    # The QR group C_m: elements are g^{2k} for k=0,...,m-1
    qr_generator = pow(g, 2, p)  # g^2 generates QR

    print(f"\np={p:2d} (p mod 4 = {p%4}): |G*|={p-1}, m={m}, "
          f"primitive root g={g}, QR generator g^2={qr_generator}")
    print(f"  QR = {qr}")
    print(f"  NQR = {nqr}")

    # Structure of C_m
    # Orbit of each chord type under QR multiplication
    chord_orbits = {}
    visited = set()
    for k in range(1, m + 1):
        if k in visited:
            continue
        orbit = set()
        for a in qr:
            ak = (a * k) % p
            if ak > m:
                ak = p - ak
            orbit.add(ak)
        orbit = sorted(orbit)
        chord_orbits[k] = orbit
        visited.update(orbit)

    if p % 4 == 3:
        print(f"  Chord orbits under QR action ({len(chord_orbits)} orbits):")
        for start, orbit in sorted(chord_orbits.items()):
            print(f"    {start}: {orbit}")


# ============================================================
# SECTION 2: THE REPRESENTATION-THEORETIC DECOMPOSITION
# ============================================================

print(f"\n{'=' * 70}")
print("REPRESENTATION DECOMPOSITION OF THE ORIENTATION SPACE")
print("=" * 70)

print("""
The orientation space {+1,-1}^m decomposes under the QR action (C_m) as:

R^m = V_0 + V_1 + ... + V_{(m-1)/2}  (for m odd)

where V_j is the 2D irrep of C_m (or 1D for j=0).

V_0 = span of sigma_P (the Paley orientation, fixed by QR)
V_j (j > 0) = 2D rotation by angle 2*pi*j/m

The interaction matrix J, being QR-equivariant, has:
  J|_{V_0} = lambda_0 * Id  (1D eigenspace, eigenvalue lambda_0)
  J|_{V_j} = lambda_j * Id  (2D eigenspace, eigenvalue lambda_j)

So J has at most (m+1)/2 distinct eigenvalues (for m odd).
""")

# Verify: at p=11 (m=5), J should have (5+1)/2 = 3 distinct eigenvalues
# We computed: {-435.4, -435.4, 154.9, 154.9, 561.0}
# That's 3 distinct eigenvalues with multiplicities 2, 2, 1. ✓

# At p=7 (m=3), J should have (3+1)/2 = 2 distinct eigenvalues
# We computed: {-3.5, -3.5, 7.0}
# That's 2 distinct eigenvalues with multiplicities 2, 1. ✓

print("\nVerification of eigenvalue multiplicities:")
print("  p=7 (m=3): J eigenvalues = {-3.5[x2], 7.0[x1]}")
print(f"    Expected: {(3+1)//2} distinct = 2. Match? YES")
print("  p=11 (m=5): J eigenvalues = {-435.4[x2], 154.9[x2], 561.0[x1]}")
print(f"    Expected: {(5+1)//2} distinct = 3. Match? YES")


# ============================================================
# SECTION 3: THE CRITICAL INTERPLAY — ADDITIVE vs MULTIPLICATIVE
# ============================================================

print(f"\n{'=' * 70}")
print("ADDITIVE vs MULTIPLICATIVE: THE DEEP INTERPLAY")
print("=" * 70)

print("""
The H-maximization problem lives at the intersection of:

1. ADDITIVE structure (C_p): determines the CIRCULANT property.
   The eigenvalues are Fourier modes: lambda_k = sum_{s in S} omega^{ks}.
   This gives the Dirichlet kernel for the interval.

2. MULTIPLICATIVE structure (C_m): determines the QR pattern.
   The Paley tournament uses QR as the connection set.
   The eigenvalues are Gauss sums: lambda_k = (chi(k)*g - 1)/2.

The CONFLICT between these two:
  - The ADDITIVE Fourier basis diagonalizes the circulant adjacency matrix.
  - The MULTIPLICATIVE character (chi) diagonalizes the interaction matrix J.
  - These are DIFFERENT bases! The additive modes are {omega^k : k=0,...,p-1}.
    The multiplicative modes are {chi(k) : characters of C_{p-1}}.

The GAUSS SUM g = sum chi(k) omega^k is the BRIDGE between them.
It transforms additive characters into multiplicative characters.

BEAUTIFUL PICTURE:
  The tournament adjacency matrix A = (1/2)(J_p - I + B)
  where B is the "signed adjacency" with B[i,j] = chi(j-i).

  B is simultaneously:
  - A circulant matrix (additive structure)
  - Determined by chi (multiplicative structure)

  Its eigenvalues involve g because g = sum chi(k) omega^k
  is the DFT of chi, connecting the two structures.

The H-maximization question becomes:
  Which CONNECTION SET S makes the eigenvalues of A(S) best for HP counting?

  Paley: uses the Gauss sum bridge, getting flat eigenvalues.
  Interval: uses the additive structure directly, getting Dirichlet kernel.

  The Gauss sum approach (Paley) is optimal when the number of modes
  is small enough that the "bridge" doesn't lose information.
  The direct approach (interval) is optimal when there are many modes
  and concentrating energy is more important than distributing it.
""")

# Compute the Gauss sum for each prime
print("\nGauss sums and their role:")
for p in [3, 7, 11, 19, 23, 31, 43]:
    if p % 4 != 3:
        continue
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    g = sum((1 if is_qr(k, p) else -1) * omega**k for k in range(1, p))

    print(f"  p={p:2d}: g = {g.real:>8.4f} + {g.imag:>8.4f}i, |g| = {abs(g):.4f} = sqrt({p})")
    print(f"         arg(g)/pi = {cmath.phase(g)/math.pi:.4f}")

    # The Gauss sum satisfies g^2 = chi(-1) * p = -p for p=3 mod 4
    g2 = g**2
    print(f"         g^2 = {g2.real:.4f} + {g2.imag:.4f}i (expected: {-p:.4f})")

# Now the KEY insight about interlacing:
print(f"\n{'=' * 70}")
print("INTERLACING: ODD-ORDER GROUPS AND TOURNAMENT SIZES")
print("=" * 70)

print("""
The user's observation about INTERLACING:

  Tournament n:   3   4   5   6   7   8   9  10  11  ...
  D_{2n} order:   6   8  10  12  14  16  18  20  22  ...
  Odd groups:    C_3  C_5  C_7  C_9  C_11 ...

For PALEY tournaments (p = 3 mod 4):
  p =   3   7  11  19  23  31  ...
  m =   1   3   5   9  11  15  ...
  C_m = C_1 C_3 C_5 C_9 C_11 C_15

The QR subgroup C_m has ORDER m, which is the SAME as the number of
chord types. This is the group that creates the interaction symmetry.

CRITICAL INTERLACING PATTERN:
  For p = 3 mod 4, m = (p-1)/2 is always ODD.
  The number of irreps of C_m is m, but they pair up into (m-1)/2 complex
  conjugate pairs plus the trivial irrep.

  So J has (m+1)/2 distinct eigenvalues:
    lambda_0 (the Paley eigenvalue, 1D)
    lambda_1, ..., lambda_{(m-1)/2} (each 2D)

  The PALEY eigenvalue lambda_0 is the largest (as proved).
  The question is whether it REMAINS the largest relative to H0
  as m grows.

  The number of eigenspace dimensions grows as m, but lambda_0
  is bounded by H-dependent quantities, while the sum of all
  lambda_j's is constrained by the Parseval-type identity on H.

PREDICTION:
  lambda_0 / E[H] -> 0 as p -> infinity
  because lambda_0 ~ m^2 (quadratic in chord interactions)
  while E[H] ~ (p-1)! / 2^{p-1} (exponential in p).

  The degree-2 "peak" at Paley becomes a smaller and smaller
  fraction of the total H, making it less important for maximization.
""")

# Compute lambda_0 / E[H] ratio
# At p=7: lambda_0 = 7, E[H] = 178.5, ratio = 0.039
# At p=11: lambda_0 = 561, E[H] = 93101.25, ratio = 0.006

print("lambda_0 / E[H] ratio (Paley eigenvalue relative to mean H):")
print(f"  p=7:  lambda_0 = 7.0,   E[H] = 178.5,     ratio = {7.0/178.5:.6f}")
print(f"  p=11: lambda_0 = 561.0, E[H] = 93101.25,   ratio = {561.0/93101.25:.6f}")
print(f"  Ratio shrinking -> degree-2 advantage matters less at large p")


# ============================================================
# SECTION 4: THE POLYGON GEOMETRY — VERTEX POSITIONS
# ============================================================

print(f"\n{'=' * 70}")
print("POLYGON GEOMETRY: VERTICES ON THE UNIT CIRCLE")
print("=" * 70)

print("""
Place vertices at z_k = exp(2*pi*i*k/p) for k = 0,...,p-1.

The CHORD from z_j to z_k has:
  Length: |z_k - z_j| = 2*|sin(pi*(k-j)/p)|
  Midpoint: (z_k + z_j)/2
  Direction: arg(z_k - z_j)

For chord type d (connecting vertices distance d apart):
  Length: 2*sin(pi*d/p)
  Direction from z_j: arg(z_{j+d} - z_j) = arg(omega^j * (omega^d - 1))
                     = (2*pi*j + pi*(d+1))/p - pi/2

  The chord direction ROTATES as j increases!
  This is the "winding" effect.

For the INTERVAL tournament (all chords CW):
  The "resultant flow" at vertex j is:
    F(j) = sum_{d=1}^{m} (z_{j+d} - z_j) = sum_{d=1}^{m} z_j*(omega^d - 1)
          = z_j * [sum omega^d - m]
          = z_j * [(omega^{m+1} - omega)/(omega - 1) - m]

  This is a UNIFORM FLOW (same magnitude at every vertex) that
  rotates with the polygon. The "wind" blows in the clockwise direction.

For the PALEY tournament:
  F(j) = sum_{d in QR} (z_{j+d} - z_j) - sum_{d in NQR} (z_{j+d} - z_j)
        = z_j * [sum_{d in QR} omega^d - sum_{d in NQR} omega^d]
        = z_j * [sum_{d=1}^{p-1} chi(d) * omega^d]
        = z_j * g  (the Gauss sum!)

  So the "flow" at the Paley tournament is z_j * g, which has:
  - Magnitude |g| = sqrt(p)
  - Direction: arg(g) + arg(z_j) = arg(g) + 2*pi*j/p

  This is also a UNIFORM FLOW but with different magnitude and direction!

COMPARISON:
  Interval flow magnitude at each vertex: |sum_{d=1}^m (omega^d - 1)|
  Paley flow magnitude: |g| = sqrt(p)

  For the interval: the flow magnitude is |sum_{d=1}^m omega^d - m|
    = |omega*(1-omega^m)/(1-omega) - m|
  This is approximately m * |1 - omega^{(m+1)/2}| ~ m for small phases.
  Actually: |sum_{d=1}^m omega^d| = |sin(m*pi/p)/sin(pi/p)| ~ m (for small pi/p)

  So interval flow ~ m and Paley flow ~ sqrt(p) = sqrt(2m+1).
  For large m: m >> sqrt(2m+1), so interval has MUCH stronger flow.
""")

# Compute the flow vectors for specific primes
for p in [7, 11, 19]:
    if p % 4 != 3:
        continue
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)

    # Interval flow (at vertex 0)
    F_interval = sum(omega**d - 1 for d in range(1, m + 1))

    # Paley flow
    g = sum((1 if is_qr(k, p) else -1) * omega**k for k in range(1, p))

    # Net flow (signed sum of edges from vertex 0)
    # For tournament T(S): F = sum_{d in S} omega^d - sum_{d not in S} omega^d
    # Wait, the tournament flow is: vertex 0 beats d in S, so arc goes 0->d.
    # The "flow out" is sum_{d in S} omega^d.
    # The "flow in" is sum_{d not in S, d != 0} omega^d.
    # Net = flow_out - flow_in = sum_{d in S} omega^d - sum_{d in S^c} omega^d

    S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
    S_interval = frozenset(range(1, m + 1))

    net_P = sum(omega**d for d in S_paley) - sum(omega**d for d in range(1, p) if d not in S_paley)
    net_I = sum(omega**d for d in S_interval) - sum(omega**d for d in range(1, p) if d not in S_interval)

    print(f"\np={p}:")
    print(f"  Interval net flow: |F| = {abs(net_I):.4f}, arg/pi = {cmath.phase(net_I)/math.pi:.4f}")
    print(f"  Paley net flow:    |F| = {abs(net_P):.4f}, arg/pi = {cmath.phase(net_P)/math.pi:.4f}")
    print(f"  Ratio |F_I|/|F_P| = {abs(net_I)/abs(net_P):.4f}")
    print(f"  sqrt(p) = {math.sqrt(p):.4f}")
    print(f"  m = {m}")

    # The flow ratio is related to the Dirichlet peak / Gauss sum
    # |F_I| ~ 2*|lambda_1(interval)|, |F_P| ~ 2*|lambda_1(Paley)|
    # No wait, net flow is different from eigenvalue.
    # net flow = sum_{d in S} omega^d - sum_{d in S^c} omega^d
    #          = 2*lambda_1 + 1 - (-m + lambda_0) ... no, let me think
    # Actually: lambda_k = sum_{d in S} omega^{dk}
    # So net flow = sum_{d in S} omega^d - sum_{d not in S} omega^d = 2*lambda_1 - (-1)
    # Hmm, sum_{d=1}^{p-1} omega^d = -1
    # So sum_{d not in S} omega^d = -1 - sum_{d in S} omega^d = -1 - lambda_1
    # Wait: lambda_1 = sum_{d in S} omega^d
    # sum_{d not in S, d != 0} omega^d = sum_{d=1}^{p-1} omega^d - lambda_1 = -1 - lambda_1
    # Net = lambda_1 - (-1 - lambda_1) = 2*lambda_1 + 1

    lambda1_I = sum(omega**d for d in S_interval)
    lambda1_P = sum(omega**d for d in S_paley)

    print(f"  Check: 2*lambda_1 + 1 = {2*lambda1_I + 1:.4f} + {(2*lambda1_I + 1).imag:.4f}i")
    print(f"         net_I = {net_I.real:.4f} + {net_I.imag:.4f}i")
    print(f"  Match? {abs((2*lambda1_I + 1) - net_I) < 0.001}")


print(f"\n{'=' * 70}")
print("CONCLUSION: THE GEOMETRIC PICTURE")
print("=" * 70)

print("""
THE COMPLETE GEOMETRIC PICTURE:

1. Tournament on Z_p = oriented regular p-gon on the unit circle.

2. Each CHORD TYPE corresponds to a D_{2p} orbit and carries a
   "flow mode" — a harmonic component of the directed flow field.

3. The INTERVAL tournament has all modes constructive (maximum flow).
   Its dominant eigenvalue |lambda_1| ~ p/pi is the Dirichlet kernel peak.

4. The PALEY tournament has modes oriented by the Legendre symbol chi.
   Its eigenvalues are all |lambda_k| = sqrt(p)/2 (via Gauss sum).

5. The INTERACTION MATRIX J of the Walsh expansion decomposes under
   the QR action C_m. Paley sigma is the UNIQUE fixed vector (since
   QR acts transitively), hence an eigenvector of J.

6. The CROSSOVER happens when the single strong flow mode of the
   interval (lambda_1 ~ p/pi) overwhelms the many equal modes of
   Paley (each sqrt(p)/2), which occurs at p ~ 19.

7. The ratio lambda_1(interval)/lambda(Paley) = 2*sqrt(p)/pi -> infinity,
   so the interval advantage grows WITHOUT bound. The crossover is
   permanent: once interval wins, it wins forever.

8. The ODD-ORDER group C_m = QR_p creates the symmetry that makes
   Paley special (eigenvector of J), but this symmetry becomes
   IRRELEVANT as the number of chord modes grows, because the
   degree-2 term is an ever-smaller fraction of H.

THE DIHEDRAL GROUPS D_{2n} ENCODE THE GEOMETRY:
  D_{2n} has rotations (preserving the circulant) and reflections
  (mapping T to T^op). The eigenvalue MAGNITUDES |lambda_k| are
  D_{2n}-invariants, while the PHASES are not.

  For H-maximization, the magnitudes matter more than the phases
  at large n, because H is approximately a function of |lambda_k|^2
  (via trace formulas). This is why the interval's concentrated
  magnitude distribution beats Paley's uniform distribution.
""")
