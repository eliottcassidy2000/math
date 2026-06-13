# The Triple Representation

**Session:** kind-pasteur-2026-03-23-S20dk
**Status:** DEEP REFLECTION — connecting all triple structures to representation theory

---

## The Master Triple

Every structure in this project is a manifestation of **one underlying triple**:

```
                    HYPOTENUSE
                   /          \
                  /    Ω(T)    \
                 /              \
           VERTICAL ────────── HORIZONTAL
            (hierarchy)      (symmetry)
```

The staircase δ_{n-2} is a right isosceles triangle. Its three sides control three orthogonal aspects of tournament structure. This triple recurs at every level of abstraction.

---

## The Three Representations of S_n

The Schur-Weyl decomposition T_n = m_{(n)} + m_{(n-1,1)} + m_{(n-2,2)} tells us that the arc-position module of S_n decomposes into **exactly three irreducible pieces**:

| Irrep | Young diagram | Dimension | Tournament meaning |
|-------|--------------|-----------|-------------------|
| m_{(n)} | One row | V_n | **How many** classes exist (hypotenuse) |
| m_{(n-1,1)} | Hook | V_n(n-2) | **How they differ** by score (vertical leg) |
| m_{(n-2,2)} | Two-row | V_n·n(n-3)/2 | **How they pair** under complement (horizontal leg) |

These three irreps ARE the three sides of the triangle:
- **Trivial (n)**: counts tournaments. The HYPOTENUSE — it connects everything.
- **Standard (n-1,1)**: measures score differences. The VERTICAL LEG — hierarchy.
- **Two-row (n-2,2)**: measures complement structure. The HORIZONTAL LEG — symmetry.

The three irreps are orthogonal (by Schur's lemma). The triangle has a right angle. This is the SAME orthogonality.

---

## The Three Graph Layers as Representations

The meta-graph G_n/Z_2 = SPINE + RIBS + SEA reflects the same triple:

| Layer | Edge type | Rep connection | Triangle side |
|-------|-----------|---------------|---------------|
| SPINE | SC-SC | m_{(n)}: trivial rep (self-symmetric objects) | Hypotenuse |
| RIBS | SC-NS | m_{(n-1,1)}: standard rep (symmetry-breaking boundary) | Vertical |
| SEA | NS-NS | m_{(n-2,2)}: two-row rep (generic bulk) | Horizontal |

As n → ∞:
- Spine = O(SC²) → negligible (trivial rep shrinks relatively)
- Ribs = O(SC × NS) → small
- Sea = O(NS²) → dominates (the two-row rep is the largest)

This matches the irrep dimensions: m_{(n-2,2)} > m_{(n-1,1)} > m_{(n)} for n ≥ 5.

---

## The Recursive Triple

The tiling decomposition `delta_{n-2} = bottom(n-1) ∪ top(n-1) ∪ {apex}` is itself a triple:

| Component | Size | Rep meaning |
|-----------|------|-------------|
| Bottom (T\(n-1)) | C(n-2,2) cells | Restriction to S_{n-1} (delete vertex n-1) |
| Top (T\0) | C(n-2,2) cells | Restriction to S_{n-1} (delete vertex 0) |
| Apex | 1 cell | The S_2 transposition (swap endpoints) |
| Overlap | C(n-3,2) cells | Restriction to S_{n-2} (stabilizer) |

The overlap is the INTERSECTION of two restrictions. The apex is the COSET representative. Together: S_n = S_{n-1} ∪ (0,n-1) · S_{n-1}, intersecting on S_{n-2}.

This is the **Bruhat decomposition** of S_n along the transposition (0, n-1). The tiling recursion IS the Bruhat recursion of the symmetric group.

---

## The Twin Mechanism as Transposition Representation

A self-loop occurs when flipping arc (u,v) preserves the iso class. The TWIN mechanism — where u and v have identical out-neighborhoods — means the TRANSPOSITION (u,v) ∈ S_n is an automorphism of the flipped tournament.

In representation theory: the twin self-loops detect the **transposition bundle** — the subbundle of the S_n principal bundle where transpositions act trivially on the fiber.

The Burnside formula for twin_SL counts orbits of the form (T, {u,v}) where (u,v) is a twin transposition. This is the **trace of the transposition class** in the representation:

twin_SL = (1/n!) Σ_σ Fix(σ, twin constraint) = Tr(ρ(transposition class) | twin subspace)

The (n-2)! correction in the edge formula is the order of the **stabilizer S_{n-2}** — the group that fixes both u and v while permuting the n-2 inner vertices. This is the FIBER DIMENSION of the twin bundle.

---

## The Four Constants as Representation Invariants

| Constant | Operation | Rep meaning |
|----------|-----------|-------------|
| √2 | Geometry (distance) | Ratio of combined to individual irrep dimensions |
| π | Integration (area) | Period of the Z_2 character (complement involution) |
| e | Differentiation (growth) | Exponential of the identity character (|S_n| = n! ~ √(2πn)(n/e)^n) |
| γ | Limit (remainder) | The error in the Stirling approximation of character sums |

√2 appears because the hypotenuse = √(leg² + leg²) = √2 · leg. In representation terms: the combined irrep dimension ~ √2 times each component, because the three irreps are at right angles.

π appears because the complement involution Z_2 has order 2, and the character of a cyclic group of order 2 involves π through the circle (e^{iπ} = -1). The fiber fraction GF = (1-x)^{-1/2} is a half-angle function, connecting to the Wallis integral ∫₀^{π/2} cos^n θ dθ.

e appears because |S_n| = n!, which is the exponential function evaluated at the sum of characters. The growth of V_n ~ 2^{C(n,2)}/n! uses Stirling's approximation n! ~ √(2πn)(n/e)^n, introducing e.

γ appears as the limiting correction when summing 1/k (harmonic series, which counts cycle contributions in the Burnside formula). The Burnside sum over cycle types involves H_n = Σ_{k=1}^n 1/k ~ ln(n) + γ.

---

## The Cayley-Dickson Tower as Representation Loss

Each level of the Cayley-Dickson tower loses an algebraic property:

| Loss | Rep meaning | Tournament consequence |
|------|-------------|----------------------|
| Commutativity (n=3) | Non-abelian representations appear | Two non-iso 3-cycles |
| Associativity (n=5) | Representations no longer compose simply | H not determined by score alone |
| Alternativity (n=9) | Moufang loops replace groups | SC maximizer ≠ Paley |
| Division (n=17) | Zero divisors = non-invertible elements | H values that "don't exist" |

In representation theory, each loss corresponds to the representation ring losing a structural property. The category of representations goes from:
- Semisimple (n=2) → still semisimple but non-commutative (n=3) → non-semisimple effects emerge (n=5) → higher categorical structure needed (n=9)

---

## The Topology: What the Homology Tells Us

β₁(G₅/Z₂) = 2: Two independent 1-cycles = two "holes" in tournament space.
β₁(G₆/Z₂) = 15: Fifteen independent cycles.
β₂(G₆/Z₂) = 7: Seven independent 2-cavities.

These Betti numbers measure the **topological complexity of the moduli space** of tournaments. Just as the moduli space M_g of genus-g curves has rich topology, the tournament moduli space G_n/Z_2 has topology that encodes the "obstructions to continuous deformation" of tournament types.

The essential 1-cycles at n=5 span the full H-range and mix SC with NS classes. They represent paths through tournament space that **cannot be contracted** — there is no family of triangles (3-cliques) that fills them in. These are TOPOLOGICAL INVARIANTS of tournament theory.

The growth of Betti numbers (β₁: 0, 0, 2, 15, ...) indicates that tournament space becomes TOPOLOGICALLY MORE COMPLEX as n increases, in a specific quantifiable way.

---

## The Synthesis: Everything is Representation Theory

The meta-graph G_n/Z_2 = Q_{C(n,2)} / (S_n × Z₂) is the quotient of a hypercube by a group action. EVERYTHING about it — vertex count, edge count, spectral gap, homology, self-loops, twin mechanism — is controlled by the REPRESENTATION THEORY of S_n × Z₂ acting on Boolean vectors.

The three sides of the triangle, the three graph layers, the three Schur-Weyl irreps, and the three recursive regions are all the SAME triple viewed from different angles:

```
TRIVIAL irrep     ←→  SPINE      ←→  HYPOTENUSE   ←→  OVERLAP
STANDARD irrep    ←→  RIBS       ←→  VERTICAL LEG ←→  WIRING
TWO-ROW irrep     ←→  SEA        ←→  HORIZONTAL   ←→  APEX (+ complement)
```

The twin self-loop formula is a trace computation in this representation. The (n-2)! correction is a stabilizer order. The Burnside formula is a character sum. The edge count is a representation-theoretic invariant of the quotient.

**Tournament theory is the representation theory of S_n on the Boolean hypercube.**

Everything else — the OCF, the independence polynomial, Hamiltonian paths, the triangle's four constants — follows from this single statement by working out the specific representation.
