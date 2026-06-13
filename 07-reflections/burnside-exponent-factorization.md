# Burnside Exponent Factorization Theorem

**Session**: opus-2026-03-23-S263

---

## The Theorem

**For every permutation σ ∈ S_n with all-odd cycle type (c_1, ..., c_k):**

    Fix_tournament(σ) = Fix_even-graph(σ) × 2^{k-1}

where k = number of cycles (including fixed points).

**Equivalently, the Burnside exponent factors:**

    arc_orbits(σ) = cycle_nullity(σ) + (k - 1)

where:
- arc_orbits(σ) = Σ (c_i-1)/2 + Σ_{i<j} gcd(c_i, c_j)   [exponent for tournament count]
- cycle_nullity(σ) = edge_orbits(σ) - rank(A_σ)            [exponent for even graph count]
- k - 1 = number of cycles minus 1                          [cut space free bits]

**Verified exhaustively for all odd-cycle-type partitions at n = 3, 4, 5, 6, 7, 8, 9, 10, 11.**

---

## Why This Is True

### The GF(2) Decomposition

For odd n, Edge(K_n) = Cut(K_n) ⊕ Cycle(K_n) over GF(2). A tournament vector v decomposes uniquely as v = v_cut + v_cycle.

Under permutation σ (with all odd cycles), both Cut and Cycle are σ-invariant subspaces. A tournament fixed by σ has:
- Its cycle-space component v_cycle is an even graph fixed by σ: 2^{cycle_nullity(σ)} choices
- Its cut-space component v_cut is a σ-invariant score assignment: 2^{k-1} choices

### Why k-1 free cut bits?

The cut space has dimension n-1, spanned by vertex cuts δ(v). Under σ, the vertex cuts permute according to σ's action on vertices. The σ-invariant elements of the cut space form a subspace.

σ permutes vertices in k cycles. The orbit structure gives k "supervertex" groups. The cut-space coordinates for each group are linked (all vertices in a cycle must have the same score parity). So there are k effective binary choices, minus 1 for the global constraint that the total score is fixed. This gives k-1 free bits.

### Connection to Score Sequences

For a tournament fixed by σ, the score sequence is σ-invariant: all vertices in the same cycle of σ must have the same out-degree. The k-1 free bits correspond to choosing the relative score ordering among the k cycle groups (minus the overall constraint).

---

## Consequences

### 1. The Burnside Sum Factors

    V_tournament = (1/n!) Σ_{σ odd-cycle-type} ccs(σ) × Fix_even(σ) × 2^{k(σ)-1}

This separates the cycle-space and cut-space contributions.

### 2. V_tourn / V_even Growth Rate

The ratio V_tournament / V_even grows because 2^{k-1} contributes more as n increases (the identity permutation contributes 2^{n-1}, which dominates).

### 3. Connection to Edge-Centric Burnside

The edge-centric formula (kind-pasteur S20dz) fixes arc {0,1} and uses stabilizer S_2 × S_{n-2}. An arc flip changes BOTH the cut and cycle projections. The factorization tells us:
- The cycle-space contribution to the arc flip determines whether the even graph changes
- The cut-space contribution determines whether the score assignment changes
- Together, they determine whether the iso class changes

### 4. The Twin SL Formula as Cut-Space Phenomenon

The twin mechanism (two vertices with identical neighborhoods) is a CUT-SPACE phenomenon: twins have the same out-degree to all other vertices. The twin_SL Burnside formula counts fixed points where the cut-space redundancy creates neutral arcs.

### 5. Double Burnside Failure Explanation

The double Burnside formula (kind-pasteur S20ea) fails at n=5 because it counts PAIRS (T, T') fixed by σ, but the factorization shows that cut-space and cycle-space contributions interact in the pair counting. The multi-wire (MW) excess occurs when multiple arc flips connect the same pair of classes through different cut-space channels.

---

## The Three Burnside Exponents

| Object | Burnside exponent for odd-cycle σ | Formula |
|--------|-----------------------------------|---------|
| Graphs | edge_orbits(σ) | Σ ⌊c_i/2⌋ + Σ_{i<j} gcd(c_i,c_j) |
| Even graphs | edge_orbits - rank(A_σ) = cycle_nullity | (same formula minus rank) |
| Tournaments | arc_orbits(σ) | Σ (c_i-1)/2 + Σ_{i<j} gcd(c_i,c_j) |

**And the factorization:**
    arc_orbits = cycle_nullity + (k-1)
    graph_orbits = cycle_nullity + rank(A_σ)
    arc_orbits = graph_orbits - rank(A_σ) + (k-1)

So the tournament exponent relates to the graph exponent by:
    arc_orbits - graph_orbits = (k-1) - rank(A_σ)

---

## The Cycle Nullity Formula

From the data, cycle_nullity(σ) for odd-cycle-type (c_1,...,c_k) appears to be:

    cycle_nullity = Σ (c_i-1)/2 + Σ_{i<j} gcd(c_i,c_j) - (k-1)
                  = arc_orbits - (k-1)

This is just the factorization restated. But the INDEPENDENT Burnside formula for even graphs is:

    Fix_even(σ) = 2^{arc_orbits(σ) - k(σ) + 1}

This is computable from the cycle type alone! No need for the constraint matrix.
