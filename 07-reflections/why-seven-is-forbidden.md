# Why Seven Is Forbidden: The One Mechanism

*opus-2026-04-04-S10*

## The Question

7 and 21 are forbidden as H values (Hamiltonian path counts). They are also forbidden as tiling counts per isomorphism class. Why does the SAME number appear as forbidden in both senses? What is the single underlying mechanism?

## What Each Number Represents

Before asking why 7 is forbidden, ask: what does H = 7 *mean*?

H(T) = I(Ω(T), 2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...

Each term has a meaning:
- **1** = the empty independent set. Always present. The "ground state." Represents: the EXISTENCE of the tournament (there is always at least one Hamiltonian path).
- **2α₁** = the single-cycle contribution. α₁ odd cycles, each contributing weight 2. Represents: the CURVATURE of the tournament (odd cycles are the curvature elements of a directed graph).
- **4α₂** = the disjoint-pair contribution. Represents: INDEPENDENCE — two cycles living in non-overlapping parts of the tournament.
- **8α₃** = the triple-packing. Represents: DEPTH — three layers of independent structure coexisting.

**H = 7 = 1 + 2·3 + 4·0** means: three cycles, all conflicting, no independence. Pure curvature, zero depth.

## The Obstruction: Why Pure Curvature at Level 3 Is Impossible

The conflict graph Ω has α₁ = 3 vertices (cycles) forming a complete graph K₃ (all pairs share a tournament vertex).

But K₃ in Ω with no further vertices requires three odd cycles that are:
1. Pairwise vertex-sharing (adjacent in Ω)
2. The ONLY odd cycles (no other cycles exist)

Condition 1 forces the three cycles to share a COMMON vertex (by the pigeonhole argument: three 3-cycles on ≤5 vertices, pairwise sharing, must have a common vertex). This creates a "bouquet" — three petals meeting at one point.

But a bouquet of three 3-cycles at a common vertex spans exactly 5 vertices of the tournament. The tournament is COMPLETE on these 5 vertices (every pair has an arc). This completeness forces a 5-cycle through the same 5 vertices — an ADDITIONAL cycle that violates condition 2.

**The tournament's completeness is the cause.** In a general directed graph, you could have three cycles sharing a vertex without a 5-cycle. But in a TOURNAMENT, every pair of vertices has an arc, and this density forces the 5-cycle.

## What "Completeness Forces Additional Structure" Really Means

This is the deep point. A tournament is a COMPLETE directed graph — it has maximum density. Every pair of vertices interacts. There are no "empty" relationships.

When you try to create a specific structure (three conflicting cycles), the completeness means that the remaining arcs — the ones NOT used by the three cycles — are DETERMINED by the tournament property. And these determined arcs inevitably create new cycles.

**7 is forbidden because complete graphs cannot support "pure curvature" at the three-cycle level.** The completeness forces contamination — additional structure that pushes H above 7.

This is not a coincidence. It's a fundamental feature of complete graphs: they have too many arcs to allow isolated substructures. Every substructure bleeds into the whole.

## The Number 3

Why does the obstruction appear at EXACTLY 3 cycles?

- **1 cycle**: one 3-cycle on 3 vertices. The remaining n−3 vertices are unconstrained. No forced additional cycles. H = 3 is achievable.
- **2 cycles**: two 3-cycles sharing a vertex, on 5 vertices. The remaining arcs among these 5 vertices CAN be arranged to avoid a 5-cycle (e.g., if the non-cycle arcs form a transitive pattern). H = 5 is achievable.
- **3 cycles**: three 3-cycles sharing a common vertex, on 5 vertices. NOW the remaining arcs are fully determined (there are C(5,2) − 3·2 = 10 − 6 = 4 remaining arcs among 5 vertices, and they must form a consistent tournament). The completeness forces a 5-cycle. H = 7 is NOT achievable.

**3 is the critical count because at 3 cycles the remaining degrees of freedom drop to zero.** The system becomes RIGID, and the forced structure appears.

This is analogous to Ramsey theory: R(3,3) = 6 says that any 2-coloring of K₆ contains a monochromatic K₃. Here: any tournament on 5 vertices with three "bouquet" 3-cycles at a common vertex must contain a 5-cycle. The forcing threshold is 3 cycles on 5 vertices.

## From H-Forbidden to Tiling-Count-Forbidden

Now: why is tiling count = 7 also forbidden?

tc = H / |Aut|. For tc = 7:
- |Aut| = 1, H = 7: forbidden (the topological obstruction)
- |Aut| = 3, H = 21: forbidden (the composed obstruction)
- |Aut| = 5, H = 35: H exists, but NO tournament with H=35 has |Aut|=5
- |Aut| = 7, H = 49: H exists, but NO tournament with H=49 has |Aut|=7

Why can't H=35 coexist with |Aut|=5?

|Aut| = 5 means a 5-fold cyclic automorphism. On 7 vertices, a 5-cycle fixes 2 vertices and rotates 5. The 5 rotated vertices must form a tournament invariant under the rotation — a REGULAR tournament on 5 vertices.

A regular tournament on 5 vertices has H = 15 and α₁ = 6 (six 3-cycles). This is a HIGHLY SYMMETRIC, HIGH-CURVATURE tournament. To get the full 7-vertex tournament to have H = 35, the two fixed vertices must contribute precisely the right additional paths.

But H = 35 = 5·7 means the full tournament has 17 additional cycles (α₁) or specific independence structure. The 5-fold symmetry constrains the cycle structure so rigidly that H = 35 is not achievable with this symmetry. The symmetry either forces too many cycles (H > 35) or too few (H < 35).

**The same completeness that forbids H = 7 also prevents the symmetry combinations needed for tc = 7.** The tournament's density forces too much coherence for H to land on the exact multiples of 7 with the exact symmetry divisor.

## The One Mechanism

**The fundamental mechanism is: the completeness of tournaments creates topological pressure that prevents certain exact integer targets from being hit.**

This pressure manifests as:
1. **H-forbidden values** (7, 21): the conflict graph cannot be configured to produce I(Ω, 2) = target
2. **Tiling-count-forbidden values** (7, 21): the (H, |Aut|) combination needed to divide down to the target is blocked by the same pressure

The number 7 = 1 + 2·3 is the FIRST target that requires "pure curvature" (3 maximally conflicting cycles, no depth), and tournament completeness makes pure curvature at level 3 impossible.

21 = 3·7 is the SECOND target, living at the same "z = 0 floor" in the cuboid, requiring a higher-order version of the same pure-curvature obstruction.

## What This Tells Us About Tournaments

Tournaments are graphs where EVERY relationship is determined. There are no abstentions, no neutral positions. Every pair of vertices has a definite winner.

This totality — this completeness — is what makes tournaments special among directed graphs. And it is EXACTLY this completeness that creates the forbidden values. In a general directed graph, H = 7 is achievable (you can engineer the arc structure to produce exactly the right cycle packing). But in a tournament, you can't — because the arcs you don't explicitly choose are chosen for you, and they create structure you didn't intend.

**The forbidden values are the cost of completeness.** They are the integers that cannot be expressed as the independence polynomial of any conflict graph arising from a complete directed graph.

This is why 7 appears in both the H and tiling-count prohibitions: both are downstream consequences of the same upstream constraint (completeness), transmitted through the three-functor composition (tiling → tournament → conflict graph → integer).

## The Broader Principle

If this analysis is correct, then:

1. **Any complete combinatorial structure should have forbidden output values.** The completeness forces coherence that prevents certain exact targets.

2. **The forbidden values should be small numbers** — they are the first targets where the system becomes rigid enough for the forcing to manifest.

3. **The forbidden values should relate to Ramsey-type thresholds** — they appear at the count where "pure X without Y" first becomes impossible.

4. **42 = 2 · 3 · 7** is the product of:
   - 2 = the binary nature of tournament arcs (each arc has exactly 2 states)
   - 3 = the threshold for cycle-conflict rigidity (Ramsey threshold)
   - 7 = the first forbidden value (1 + 2·3 = the value at the rigidity threshold)
   
   It is the number that encodes ALL THREE structural constraints simultaneously.
