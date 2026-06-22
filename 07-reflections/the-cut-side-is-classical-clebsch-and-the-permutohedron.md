# The Cut Side Is Classical: the Clebsch Graph (Folded Cube) and the Permutohedron

*mac-mini-2026-06-22-S40. The owner offered three classical graphs (Unital, Clebsch, Truncated
Octahedral) as inspiration. Two of them turn out to BE the cut side of the project's Cut⊕Cycle frame
(S38) — verified. The even graph was the cycle side; these are the cut side. Companion to
[[the-even-graph-is-the-tournaments-cycle-half]].*

## Two verified identifications

**Clebsch graph = the cut-space Cayley graph of K₅.** The cut space of K_n over GF(2) is (ℤ/2)^{n−1},
spanned by the n vertex-stars δ_v with Σ_v δ_v = 0. For n=5 that is (ℤ/2)⁴ with generators
{δ₁,δ₂,δ₃,δ₄, δ₁+δ₂+δ₃+δ₄} — *exactly* the folded-5-cube generators. VERIFIED: the resulting Cayley
graph is SRG(16,5,0,2) — the Clebsch graph (16 vertices, 5-regular, λ=0, μ=2). Generally, the
cut-space Cayley graph of K_n is the **folded n-cube** (2^{n−1} vertices, n-regular); Clebsch is the n=5
member. This is the *bipartition* resolution of the cut side (each edge = flip one vertex across the cut).

**Truncated octahedron = the permutohedron of S₄ = the transitive tournaments at n=4.** VERIFIED: S₄
under adjacent transpositions is the 24-vertex 3-regular truncated octahedral graph. The transitive
tournaments are the pure-cut tournaments (H=1, zero cycle part); their adjacency under "swap two
consecutive in the order" is the permutohedron. This is the *ordering* resolution of the cut side (finer
than bipartitions: a full linear order, i.e. the score sequence realized).

## The cut side has two resolutions, both classical

| resolution | object at small n | what it records |
|---|---|---|
| ordering (fine) | **permutohedron** S_{n} = truncated octahedron (n=4) | full linear order / score sequence |
| bipartition (coarse) | **folded n-cube** = Clebsch (n=5) | the cut (which side of a 2-coloring) |

Both are quotients/Cayley graphs of the *cut* summand of E(K_n) = Cut(n−1) ⊕ Cycle(C(n−1,2)). The cycle
summand gave the even graphs / the conflict graph Ω / H = I(Ω,2) (S38). So the Cut⊕Cycle frame is now
populated by classical objects on **both** sides:

> **Cut side:** permutohedra (orderings) → folded cubes / Clebsch (bipartitions) — the score/hierarchy.
> **Cycle side:** even graphs / Ω / H = I(Ω,2) — the odd cycles, the apex-7, the LRC obstruction.

That the cut side is an SRG / a space-filling permutohedron, while the cycle side carries the
arithmetic apex-7 obstruction, is the geometric expression of "score is cut-side, H is cycle-side."

## The LRC connection: the winding walk lives on the permutohedron

For 13 runners the time axis [0,1) is cut by the resonance times t = k/(s_i−s_j) into chambers; within
each chamber the runners' cyclic order is fixed, and crossing a resonance is an adjacent transposition.
**So the winding tournament T(t) (HYP-2605) traces a closed walk on the permutohedron of S₁₃ — the cut
side is its state space**, and its odd cycles / H (the cycle side) record how twisted that order is.

Honest caveat (discipline): the LRC *safe* condition is NOT this subdivision. Safety is
observer-relative — runner s out of (−1/14,1/14) — so it cuts [0,1) at t = k/s_i ± 1/(14 s_i), a
*different* (finer, metric) subdivision than the pairwise-order resonances t = k/(s_i−s_j) that define
the permutohedron walk. The winding tournament records the order; the safe set records distance to the
observer. They are related (the order constrains, but does not determine, who is in the zone), and the
clean question is whether the cycle-side invariant (H / odd cycles of the winding tournament) controls
the observer-safe set. So this is a state-space identification for the *winding tournament*, placing it
on the cut side — a reframing of HYP-2605 in S38 coordinates, not yet a bridge to the safe condition.

## Open: the Unital (the third graph)

The unital 2-(q³+1, q+1, 1) (q=2: AG(2,3), 9 points; q=3: 2-(28,4,1)) is a *perfect pairwise covering
design*. The LRC counterexample is a covering system (unsafe residues cover ℤ/D). A unital-type rigidity
(each pair covered exactly once) would be a candidate non-realizable structure for the covering route
(HYP-2876/2878) — but I have not found a concrete map. Logged as a lead, not a result.

Verified core (S40): Clebsch = cut-space Cayley(K₅); truncated octahedron = permutohedron(S₄). Related:
[[the-even-graph-is-the-tournaments-cycle-half]], HYP-2605 (winding tournament), HYP-2888 (LRC coverage).
