# The King Vertex as a Window into Tournament Structure

**Session:** opus-2026-05-27-S1

## What is the King Theorem Really Saying?

The classical king theorem says: in any tournament, the vertex Q with maximum outdegree can reach every other vertex in ≤ 2 steps. This session revealed that this is the tip of an iceberg connecting score sequences, strong connectivity, the tiling model, and the OCF.

The key insight: **the number 2 in the king theorem and the number 2 in H = I(Ω,2) are the same 2.** Both arise from the binary structure of tournament arcs. The king theorem says: the squared adjacency matrix (2-step reachability) has no zero row at Q. The OCF specializes the independence polynomial at x=2. Both are evaluations at the "binary threshold" of tournament arc structure.

## The Court-Rivals Triangle

Every rival b ∈ N⁻(Q) (who beats Q) must be beaten by some court member a ∈ N⁺(Q) (who Q beats). This creates the fundamental triangle Q → a → b → Q. This is a directed 3-cycle through Q, and by Claim A it contributes 2μ(C) to H(T) − H(T−Q). Since μ(C) ≥ 1, we get the lower bound H(T)−H(T−Q) ≥ 2|rivals|.

This triangle (Q, court, rival) IS the "everything is the triangle" principle at the vertex level. The staircase triangle's three sides correspond exactly to: Q (source column / vertical leg), P (sink row / horizontal leg), and the apex (the corner connecting them). Every instance of the king theorem is an instance of the triangle principle.

## Strong Connectivity = Triangle Permeability Everywhere

In the tiling model, strong connectivity means: for EVERY cut k ∈ {1,...,n-1}, at least one tile goes "upward" (against the base path flow). This is equivalent to saying: the triangle structure is present at every level of the staircase, not just locally.

Non-strongly-connected tournaments = tilings with a "blocked cut" — a level where all flow goes downward, creating a one-way membrane. The tournament has a "feudal hierarchy": a dominant class that can only send information down, never receive from below.

The apex tile (n-1, 0) is the special tile that, when upward, instantly creates the "master cycle" n-1 → n-2 → ... → 0 → n-1. This Hamiltonian cycle is the ultimate expression of democracy: a cycle that includes everyone.

## The Q-P Axis as the Democracy Axis

The principal line of G_n/Z₂ runs from the transitive tournament (H=1, Q = absolute source, P = absolute sink, maximum hierarchy) to the regular tournament (H = max, Q = P = every vertex, zero hierarchy). This session shows this is the **Q-P democracy axis**:

- Transitive: Q has degree n-1 (absolute dominator), P has degree 0 (absolute submissive). gap = n-1.
- Regular: every vertex has the same degree. gap = 0.
- Everything between: partial hierarchy, partial democracy.

The formula H = 7 − 2×gap at n=4 makes this precise: H linearly tracks the democracy level. At n=5+, the relationship is nonlinear but remains the primary driver of H.

## The SC-Strict-Bound Dichotomy

One of the cleanest findings: for n ≥ 5, the bound H(T)−H(T−Q) = 2|rivals| (tight) implies the tournament is NOT strongly connected. The proof reveals WHY: if there's only one odd cycle through Q, the rival b is beaten by only one court member, and the other court members form a "transitive chain" that forces a vertex with outdegree 0 — a structural defect (absolute sink or source) that destroys strong connectivity.

This is the tournament's version of "you can't have it both ways." A tournament can be tightly minimal in its cycle structure OR strongly connected, but not both (for n ≥ 5). SC forces redundancy — extra triangles, extra cycles — that strictly exceeds the minimum required by the King theorem.

## The #Kings Rigidity at n=5

The discovery that #kings=4 → H=13 (exactly, for all 120 such tournaments at n=5) is striking. It means the count of king vertices is a remarkably informative invariant — more informative than the score sequence alone in some regimes. #kings measures how "democratic" the tournament is in terms of reach, which connects back to the Q-P axis.

At n=5: #kings = 5,4,3,1 (2 is impossible) corresponds to H taking values 15, 13, ?, 1-5. The impossible H=7 is the "missing" value in the H=1,3,5,9,11,13,15 sequence. #kings=4 precisely fills the gap between #kings=5 (H=15) and #kings=3 (H≤11).

## For Future Sessions

The most important open questions arising from this work:

1. **Prove or disprove:** min SC excess = 2(n-4) for all n ≥ 4. (HYP-1740)
2. **Identify** the non-SC tiling count sequence 1,3,14,121,1995,64648 in OEIS.
3. **Characterize** all tournaments where #kings = n (universal king property) for even n (where regular tournaments don't exist).
4. **Understand** why #kings = n-1 gives H=13 exactly at n=5 (and find the analog at n=7).
5. **Connect** the SC Cut Theorem to the Tutte polynomial / reliability polynomial of the tiling model.
6. **Prove** (rigorously, for all n) that tight King bound ↔ non-SC for n ≥ 5.
