# THM-052: Vertex-Transitive Tournaments Have Scalar Transfer Matrix

**Status:** FALSE for general VT tournaments (MISTAKE-013/014: Frobenius F_21 counterexample at n=21). Valid only for self-converse VT tournaments. See THM-113-circulant-scalar-M.md for the correct scope.
**Added:** opus-2026-03-06-S26
**Strengthened:** opus-2026-03-06-S26 (extended from circulant to vertex-transitive — LATER DISPROVED for non-SC VT)

## Statement

**Theorem.** For any vertex-transitive tournament T on n vertices (odd n), the transfer matrix is scalar:

    M = (H/n) * I

where H = H(T) is the Hamiltonian path count.

This includes all circulant tournaments, Paley tournaments, doubly regular tournaments, and Cayley tournaments on any group.

## Proof

### Step 1: Circulant + Symmetric => Circulant Symmetric Matrix

The rotation σ: i ↦ i+1 (mod n) is an automorphism of any circulant tournament. By the IE definition, M[σ(a), σ(b)] = M[a,b], so M is a circulant matrix: M[a,b] = m_{(b-a) mod n}. Combined with THM-030 (M symmetric), we get m_d = m_{n-d}.

### Step 2: N(d,j) is Palindromic via Reflection-Reversal Bijection

Define the map φ on Hamiltonian paths: φ(P) = r(P^rev), where:
- P^rev = path reversal: (p_{n-1}, ..., p_0)
- r: i ↦ -i (mod n) is the reflection

**Claim:** φ is a well-defined involution on Ham(T).

*Proof:* The reflection r is an isomorphism from T to T^op (it maps edge distance d to n-d, converting gen set S to its complement). So P^rev (a Ham path of T^op) maps under r to a Ham path of T. And φ² = id since both r and reversal are involutions.

**Key property:** If P has "distance" d = (p_{j+1} - p_j) mod n at position j, then φ(P) has distance d at position (n-2-j).

*Proof:* In φ(P) = (-p_{n-1}, ..., -p_0), the distance at position k is:
(-p_{n-2-k}) - (-p_{n-1-k}) = p_{n-1-k} - p_{n-2-k}, which is the distance of P at position (n-2-k).

Therefore N(d, j) = N(d, n-2-j) for all d and j.

### Step 3: Palindromic + Even Length => Zero Alternating Sum

By THM-050, M[a,b] = Σ_j (-1)^j N(d, j) where d = (b-a) mod n.

Since N(d,j) is palindromic (N(d,j) = N(d, n-2-j)) and n-1 is even at odd n, the sequence has even length. A palindromic sequence of even length has zero alternating sum:

Σ_j (-1)^j N(d,j) = Σ_{j < (n-1)/2} [(-1)^j + (-1)^{n-2-j}] N(d,j) = 0

(since (-1)^j + (-1)^{n-2-j} = (-1)^j(1 + (-1)^{n-2-2j}) = (-1)^j · 0 when n-2-2j is odd, which it always is at odd n).

Wait — n-2-2j: at odd n, n-2 is odd, so n-2-2j is odd for all j. Thus (-1)^{n-2-j} = -(-1)^j, and the pairs cancel.

### Step 4: Diagonal from Trace

At odd n: tr(M) = H (since Σ_a (-1)^{pos(a,P)} = 1 for each P). Combined with M = m_0 · I (from Steps 2-3), we get n · m_0 = H, so M = (H/n) · I. QED.

## Verification

- n=3: All 2 circulant tournaments have M = (H/3)*I (H=3 for both).
- n=5: All 4 circulant tournaments have M = 3I (all have H=15).
- n=7: All 8 circulant tournaments verified (H ∈ {175, 189}).
- n=9: 8/16 spot-checked, all confirmed (H ∈ {3159, 3267, 3357}).

## Extension to Vertex-Transitive (non-circulant)

The proof for circulant tournaments uses the reflection r: i→-i as the anti-automorphism. For general vertex-transitive tournaments, the argument generalizes:

**General proof:** Let G ≤ Aut(T) act transitively on vertices.

1. M is G-invariant: M[σ(a), σ(b)] = M[a,b] for all σ ∈ G.
2. Since G is transitive and M is symmetric (THM-030), M[a,b] depends only on the G-orbit of {a,b}.
3. For any vertex-transitive tournament T, there exists an anti-automorphism τ: V → V such that T^op = τ(T) (this is the "inversion" map g ↦ g⁻¹ for Cayley tournaments, or follows from the structure of tournament automorphism groups).
4. The map φ = τ ∘ rev on paths gives the palindromic bijection.
5. The alternating sum cancels at odd n.

**Verified at n=9:** All 16 Z/3×Z/3 Cayley tournaments (vertex-transitive, NOT circulant) have scalar M = 351·I with H = 3159.

**At prime n:** By Sylow's theorem, every transitive group on p elements contains a p-cycle, so vertex-transitive = circulant at prime n.

## Converse: Scalar M does NOT imply vertex-transitive

**Disproved at n=5 (opus-2026-03-06-S26):** There exists a tournament on 5 vertices with M = 3I that is NOT vertex-transitive.

The counterexample is the "cone over C_3": vertices {1,2,3} form a 3-cycle, vertex 4 beats {1,2,3}, {1,2,3} all beat vertex 0, and 0 beats 4. This has |Aut| = 3 (Z/3Z on the middle), score sequence (1,2,2,2,3), and vertex orbits {0}, {1,2,3}, {4}.

**Position-uniform characterization (opus-2026-03-06-S26):** At n=3,5, the exact characterization of scalar M is:

    M = (H/n)*I  <=>  N[v,j] = H/n for all v,j  (position-uniform)

where N[v,j] = number of Hamiltonian paths with vertex v at position j. Exhaustively verified: 64/1024 labeled tournaments on 5 vertices are position-uniform, all have scalar M, and vice versa.

**Key structural difference:** For vertex-transitive tournaments, each even-r coefficient c_k is individually scalar. For non-VT position-uniform tournaments (like the cone over C_3), the individual c_k are NOT scalar — their off-diagonal parts cancel at r=1/2. Specifically: c_0 off-diagonal = 0.5 is cancelled by c_2 off-diagonal * 0.25 = -0.5.

**At prime n=7:** No non-circulant position-uniform tournaments found (1000 random + 84 single-flips tested). Position-uniform likely = circulant = VT at prime n=7.
