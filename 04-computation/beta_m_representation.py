"""
Representation-theoretic analysis of β_m^{(0)} = m(m-3)/2.

KEY IDEA: The k=0 eigenspace with Z_m = QR* acting freely has:
  β_m^{(0)} = m · β_m^{orb}   where β_m^{orb} = (m-3)/2

So the orbit complex has orbit Betti number (m-3)/2 at degree m.

The orbit complex is a quotient of the diff-seq complex by Z_m.
It has:
  - 1 vertex
  - 1 edge-orbit (Omega_1/m = 1)
  - (m-1) face-orbits (Omega_2/m = m-1)
  - etc.

The orbit Euler characteristic is:
  chi_orb = 1 + sum_{d≥1} (-1)^d Omega_d/m
          = 1 + (chi - 1)/m = 1 + 0/m = 1

And the orbit Betti is:
  β_m^{orb} = β_{m+1}^{orb} = (m-3)/2

For m=5: 1 orbit cycle at degree 5, 1 at degree 6.
For m=9: 3 orbit cycles at each.

QUESTION: Is (m-3)/2 related to the representation theory of Z_m on Omega_m?

Since Z_m acts freely, Omega_m decomposes into Omega_m/m copies of the regular rep.
The boundary map commutes with Z_m, so it respects the regular-rep decomposition.

In the regular representation, every irreducible appears once. So:
  Omega_m = (Omega_m/m) copies of Z_m-reg
  The boundary maps decompose accordingly.

The orbit complex is the "trivial isotypic component" of the chain complex.
But since the action is free, the trivial component has dim = Omega_d/m.

Actually, for a free Z_m action on a chain complex:
  The orbit complex has the SAME homotopy type as each eigenspace complex.
  So β^{orb} = β^{(eigenspace)} for each eigenspace.

  For k=0: β_m^{(0)} = m · β_m^{orb} (since k=0 is the whole Z_m-invariant part,
  which has m times the orbit complex dimensions... wait, that's not right.

  For a FREE action: the k=0 eigenspace IS the Z_m-invariant chains.
  Since the action is free, dim(Omega_d^{Z_m-inv}) = Omega_d/m.
  And the homology of the invariant chains is H^{Z_m}_d, which equals β_d^{orb}.

  BUT: β_m^{(0)} = m(m-3)/2, and β_m^{orb} = (m-3)/2.
  So β_m^{(0)} = m · β_m^{orb}. This means the k=0 eigenspace homology
  is m times the orbit complex homology.

  WHY? Because the k=0 eigenspace has Omega_d^{(0)} = Omega_d (all diff-seq chains),
  NOT just the Z_m-invariant ones. The Z_m action decomposes Omega_d into
  m eigenspaces (for characters χ_0, ..., χ_{m-1} of Z_m), each of dim Omega_d/m.

  The k=0 eigenspace of Z_p is ALL of Omega_d. The Z_m action is ADDITIONAL.

  The k=0 eigenspace homology β_m^{(0)} is the homology of the full Omega_• complex
  with the standard (k=0) boundary. The Z_m action on this complex is free,
  so β_m^{(0)} = m · β_m^{orb}.

So the question reduces to: why is β_m^{orb} = (m-3)/2?

The orbit complex has dimensions:
  Omega_d^{orb} = Omega_d / m for d ≥ 1
  Omega_0^{orb} = 1

For P_7 (m=3): Omega^{orb} = [1, 1, 2, 3, 3, 2, 1]
For P_11 (m=5): Omega^{orb} = [1, 1, 4, 14, 41, 92, 140, 138, 90, 36, 6]

These are SMALL complexes. Let me compute the orbit complex Betti numbers directly.

Actually, the orbit complex boundary map is obtained from the Z_m quotient.
Each orbit is represented by one diff-seq, and the boundary of an orbit
is the orbit of the boundary. So the orbit boundary map is the same as
restricting to Z_m-invariant chains.

For a Z_m-equivariant chain complex C_• with free action:
  The orbit complex C_•/Z_m has boundary ∂^{orb} where
  ∂^{orb}([σ]) = [∂(σ)]

  This is well-defined since ∂ commutes with Z_m.
  H_*(C_•/Z_m) ≅ H_*(C_•)^{Z_m} for free actions (by Shapiro's lemma or direct).
  Wait, that's not quite right.

For a free action of a FINITE group G on a chain complex over a field of char 0:
  H_*(C_•/G) ≅ H_*(C_•)^G ≅ H_*(C_•)_G
  (invariants ≅ coinvariants since G is finite and char 0)

  dim H_*(C_•)^G = (1/|G|) sum_{g∈G} tr(g | H_*(C_•))

  By Hopf trace: tr(g | H_*(C_•)) via Lefschetz formula.
  For g ≠ id: tr(g | Omega_d) = 0 (free action). So L(g) = 1 (just from Omega_0).
  L(g) = sum (-1)^d tr(g | H_d)
  = tr(g | H_0) + sum_{d≥1} (-1)^d tr(g | H_d)
  = 1 + sum_{d≥1} (-1)^d tr(g | H_d) = 1

  So sum_{d≥1} (-1)^d tr(g | H_d) = 0 for all g ≠ id.

  For the fixed-point formula:
  dim H_d^G = (1/m) sum_{g∈Z_m} tr(g | H_d)
            = (1/m) [tr(id | H_d) + sum_{g≠id} tr(g | H_d)]
            = (1/m) [β_d + sum_{g≠id} tr(g | H_d)]

  Since H_d is a free Z_m-module (m | β_d), tr(g | H_d) = 0 for g ≠ id.
  So dim H_d^G = β_d / m. ✓

  This confirms β_d^{orb} = β_d^{(0)} / m.

So β_m^{orb} = m(m-3)/(2m) = (m-3)/2.

The ORBIT complex has Betti numbers:
  β_0^{orb} = 1
  β_m^{orb} = β_{m+1}^{orb} = (m-3)/2

For m=3: β^{orb} = [1,0,0,0,0,0,0] (contractible)
For m=5: β^{orb} = [1,0,0,0,0,1,1,0,0,0,0]
For m=9: β^{orb} = [1,0,...,0,3,3,0,...,0]

The m=5 orbit complex is interesting: it has the same Betti numbers as S^0 ∨ S^5 ∨ S^6.
A contractible complex with a 5-sphere and 6-sphere attached.

Can I compute the orbit complex for P_11 EXPLICITLY?
Orbit complex: 1 vertex, 1 edge, 4 faces, 14 tetrahedra, ...
This is small enough to analyze!
"""
from fractions import Fraction

# Known Omega dimensions
omega_7 = [1, 3, 6, 9, 9, 6, 3]
omega_11 = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

print("=== ORBIT COMPLEX DIMENSIONS ===\n")

for name, O, m_val in [("P_7", omega_7, 3), ("P_11", omega_11, 5)]:
    p = 2 * m_val + 1
    orb = [1] + [O[d] // m_val for d in range(1, len(O))]
    print(f"{name} (m={m_val}):")
    print(f"  Omega = {O}")
    print(f"  Orbit Omega = {orb}")
    chi_orb = sum((-1)**d * o for d, o in enumerate(orb))
    print(f"  chi_orb = {chi_orb}")

    # Boundary rank recursion for orbit complex
    R = [0, 0]  # R_0^orb = 0, R_1^orb = 0 (since β_0^orb = 1)
    for d in range(1, m_val):
        R.append(orb[d] - R[-1])
    print(f"  R_orb (bottom) = {R}")

    R_top = [0]  # R_{2m+1}^orb = 0
    for d in range(2*m_val, m_val + 1, -1):
        R_top.insert(0, orb[d] - R_top[0])
    print(f"  R_orb (top) = {R_top} (starting from d={m_val+2})")

    R_m_orb = R[-1]
    R_m2_orb = R_top[0]
    budget = orb[m_val] - R_m_orb
    beta_orb = (m_val - 3) // 2
    print(f"  R_m^orb = {R_m_orb}")
    print(f"  Budget = Omega_m^orb - R_m^orb = {budget}")
    print(f"  β_m^orb = (m-3)/2 = {beta_orb}")
    print(f"  R_{{m+1}}^orb = {budget - beta_orb}")
    print()

# For P_11: orbit complex has 1+1+4+14+41+92+140+138+90+36+6 = 563 cells
# This is small enough to compute directly!

# The alternating sum:
# For P_11 orbit: 1-1+4-14+41-92+140-138+90-36+6 = 1
# β_5^orb = 1, β_6^orb = 1

# For P_7 orbit: [1,1,2,3,3,2,1]
# 1-1+2-3+3-2+1 = 1
# β^orb = 0 everywhere except β_0 = 1

# KEY OBSERVATION: the orbit complex for P_7 is CONTRACTIBLE.
# 1 vertex, 1 edge (star). But it also has 2 faces, etc.
# The 1 edge connects vertex 0 to itself (loop), since the orbit
# of the edge (s_1) under Z_m gives m edges, all starting from 0 and
# going to partial sum s_1. In the orbit, there's just 1 edge.

# Actually, the orbit complex is NOT a simplicial complex -- it's a CW complex.
# The "vertex" is the orbit {0, 1, ..., p-1} under Z_p.
# The "edge" orbit {(s_1)} for s_1 ∈ QR/Z_m has one representative per orbit.

# For the boundary of the orbit edge: ∂(s_1) = face_0 - face_1 = () - () = 0
# (both faces are the same vertex orbit). So all edges are cycles!
# But there's 1 vertex and 1 edge, with ∂ = 0. So H_0 = Z, H_1 = Z.
# Wait, that gives β_1 = 1, but we computed β_1^orb = 0!

# The issue: the boundary in the k=0 eigenspace IS the standard boundary,
# which maps (s_1) to () - () = 0. So R_1 = 0, ker(∂_1) = Omega_1 = 1.
# Then β_1 = ker - im(∂_2) = 1 - R_2.
# R_2: the boundary ∂_2(s_1, s_2) has face_0 = (s_2), face_1 = (s_1+s_2), face_2 = (s_1).
# In the orbit complex, if there's only 1 edge orbit, then face_0 = face_2 = edge,
# and face_1 is... well, s_1 + s_2 might be in a different QR orbit.

# Hmm, but Z_m acts by multiplication, not addition. The orbit complex under Z_m
# doesn't collapse by addition.

# Let me reconsider. The orbit complex under Z_m:
# Vertices: 1 (the Z_m orbit of (), i.e., the single 0-cell)
# Edge orbits: QR / Z_m = {1} (since Z_m = QR acts transitively on QR by mult)
#   Actually QR = Z_m as a set, so QR/Z_m has 1 orbit.
# So there's 1 edge orbit.
#
# But Omega_1/m = m/m = 1 ✓
# Omega_2/m = m(m-1)/m = m-1 ✓
#
# 2-cell orbits: A_2 diff-seqs (s_1, s_2) with Z_m acting by (s_1,s_2) → (qs_1, qs_2).
# The orbits have size m (free action). Number of orbits = Omega_2/m = m-1.
# For m=5: 4 orbits of 2-cells.
#
# The boundary of a 2-cell (s_1, s_2):
#   ∂(s_1, s_2) = (s_2) - (s_1 + s_2) + (s_1)  [for k=0]
#                 wait, this is the standard boundary:
#   face_0 = (s_2), face_1 = (s_1 + s_2), face_2 = (s_1)
#   In Omega_1: (s_2) and (s_1) are in Omega_1 (they're in QR).
#   But (s_1 + s_2) might NOT be in Omega_1 if s_1 + s_2 ∉ QR.
#   That's exactly the constraint: (s_1, s_2) ∈ Omega_2 means
#   the junk faces cancel, so ∂ only has allowed faces.
#
#   In the orbit complex: (s_2) and (s_1) are both in the single edge orbit.
#   And (s_1 + s_2) is also in the edge orbit (since it's in QR for valid diff-seqs,
#   and Z_m acts transitively).
#
#   So ∂_orb(s_1, s_2) = () - () + () = ()? No:
#   ∂_orb([s_1, s_2]) = [s_2] - [s_1 + s_2] + [s_1]
#   In the orbit, [s_2] = [s_1] = [s_1+s_2] = [e] (the single edge class).
#   So ∂_orb = e - e + e = e. Non-zero!
#   R_2^orb = 1, so β_1^orb = 1 - 1 = 0. ✓ (since dim(ker ∂_1^orb) = 1, im ∂_2^orb = 1)

# Wait, that assumes all edges land in the same orbit class, which is true
# since Z_m acts transitively on QR.

# For the boundary map on orbits, we need to track the ORBIT of each face:
# ∂_orb([σ]) = sum of [face_i(σ)] with signs
# Since face_i(σ) is some element of A_{d-1}, we map it to its Z_m orbit class.
# The orbit class of (s_1) is {(qs_1) : q ∈ QR} = QR = single orbit.
# So [s_1] = [s_2] = [any QR element].

# More generally: (s_1, ..., s_{d-1}) and (qs_1, ..., qs_{d-1}) are in the same orbit.
# The number of orbits = Omega_{d-1}/m.

# The orbit boundary matrix has size (Omega_{d-1}/m) × (Omega_d/m).

# For P_11, m=5:
# d=2: ∂_orb: 4 orbits → 1 orbit. Each 2-cell orbit maps to 3× the edge orbit.
# Actually, ∂_orb([s_1,s_2]) = [s_2] - [s_1+s_2] + [s_1] = 3 × [e]? No.
# Not all face orbits are the same for higher d!

# For d=2: face_0(s_1,s_2) = (s_2), face_1(s_1,s_2) = (s_1+s_2), face_2(s_1,s_2) = (s_1).
# The orbits of (s_2), (s_1+s_2), (s_1) are all the single QR orbit.
# So ∂_orb([s_1,s_2]) = [e] - [e] + [e] = [e].
# All 4 orbit reps give [e], so rank = 1.
# β_1 = dim(ker ∂_1^orb) - rank(∂_2^orb) = 1 - 1 = 0 ✓

# For d=3: face_i gives 2-tuples. The orbit of (s_1, s_2) under Z_m has m elements.
# Number of 2-cell orbits = (m-1) = 4.
# ∂_orb([s_1,s_2,s_3]) involves 4 faces, each mapped to one of the 4 orbits of 2-cells.
# This is more interesting...

# I should compute this explicitly for P_11.

print("=== EXPLICIT ORBIT COMPLEX COMPUTATION ===")
print("(Not implemented yet - needs orbit identification)")
print()
print("KEY THEORETICAL RESULT:")
print("  β_m^{(0)} = m · β_m^{orb}")
print("  β_m^{orb} = (m-3)/2")
print("  For m=3: β_orb = 0 (contractible)")
print("  For m=5: β_orb = 1 (one cycle pair at degrees 5,6)")
print("  For m=9: β_orb = 3 (predicted)")
print()
print("  The orbit complex sequence (m-3)/2 = 0, 1, 2, 3, 4, ...")
print("  for m = 3, 5, 7, 9, 11, ...")
print("  This is simply the integer sequence starting at 0.")
print()
print("  CONJECTURE: β_m^{orb} grows linearly with m because")
print("  the orbit complex gains one new independent cycle per")
print("  2 additional QR elements in the connection set.")
