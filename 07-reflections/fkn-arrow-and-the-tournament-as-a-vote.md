# FKN, Arrow, and the tournament as a vote: the tiling cube is the social-choice cube

**Source:** mac-mini-2026-06-15-S3. Dispatch: consider Friedgut–Kalai–Naor deeply; the
tiling default is transitive, perturbations are arc flips; the recursive A+B+C−D−E−F+G
decomposition. Canon: THM-512, T823, HYP-2534..2536, OPEN-Q-102. Builds on THM-511
(converse-parity = FKN level hierarchy), THM-163/478 (Walsh band limit), THM-002 (OCF).

## The recognition

A tournament IS a vote. The cube {±1}^{C(m,2)} of pairwise comparisons among m
candidates — Kalai's setting for the Fourier proof of Arrow's theorem — is *literally*
the tournament cube. A voter's strict ranking is a transitive tournament; pairwise-
majority aggregation outputs one arc per pair = the **majority tournament**; the society
is **rational** iff that tournament is **transitive** (acyclic), and suffers a
**Condorcet paradox** iff it contains a **directed 3-cycle** — the minimal odd cycle,
the seed of the whole OCF. So the project's central dichotomy (transitive/boring vs
cyclic/structured) is the social-choice dichotomy (rational vs paradoxical), and the
odd cycles the project counts are Condorcet cycles.

FKN is the theorem that makes this quantitative. A Boolean function with ≥1−ε of its
Fourier mass on level 1 is O(ε)-close to a dictator. Arrow's classical theorem is the
ε=0 case (always-rational ⟹ dictatorship); Kalai/Mossel upgrade it: **few Condorcet
cycles ⟹ near-dictator ⟹ (project side) near-transitive**, near the H=1 corner. The
"FKN-defect" of THM-511 — the fraction of an invariant's Fourier mass above level 1 — is
exactly "how Condorcet-cyclic, vs rational/ranking, the tournament is," and it grows
with n (0.25→0.51→0.68): **large tournaments stop being rankings**, the precise sense in
which generic pairwise data is paradoxical.

## Three things that snap into place

1. **The transitive tournament is Arrow's dictator and the ground state at once.**
   c3=0, H=1, score (0,…,n−1), the all-aligned tiling, |Aut|=1, self-complementary —
   THM-511's converse-symmetric reference where all cyclic (even-level) mass vanishes
   and only the level-1 score line survives. It is "the rational outcome," and Arrow says
   the only aggregation rules that always reach it are dictatorships.

2. **c3 is the Guilbaud quadratic.** Var(c3) = 3·C(m,3)/16 (verified m≤7), entirely at
   Walsh level 2 — each triple's cyclicity has Fourier mass only on its three arc-pairs,
   no degree-3 term (its two orientations are full reversals). The probability that
   3-candidate majority is rational tends to Guilbaud's number ≈ 0.9123 =
   ¾ + ¾·(2/π)arcsin(1/3) = (3/2π)arccos(−1/3) — an arcsine, a level-2 object. The project already knew
   c3 ⊂ level 2 (THM-163); naming it the Condorcet/Guilbaud quadratic is the bridge.
   And P(random m-tournament transitive) = m!/2^{C(m,2)} → 0 fast: rationality is
   exponentially rare; the ground state is a measure-zero island.

3. **The Möbius sieve is the Walsh dual.** The dispatch's A+B+C−D−E−F+G is the
   inclusion–exclusion Σ_{∅≠A}(−1)^{|A|+1}φ(T−A) = φ(T) over the n−1 non-anchor vertices
   (codim-c count C(n−1,c), sign (−1)^{c+1}, sieving to 1). "A tournament is comprised of
   its sub-tournaments" is the Möbius function of the Boolean lattice — and it is *dual*
   to the Walsh-Fourier expansion (a cube function is the signed sum of its sub-structure
   contributions, organized by degree). The anchor vertex is the chosen base-path /
   symmetry-breaking reference — the same choice that turns the score into a genuine
   level-1 (dictator) signal (THM-511). Reconstruction-by-sieve and expansion-by-degree
   are the two faces of one cube.

## The two perturbation scales, FKN-read

- **Flip one arc (wiggly, Hamming 1)** = the smallest ranking edit: shifts two scores by
  ±1, may create one 3-cycle. From the transitive ground state, the m single-flips are
  m pairwise-non-isomorphic "dictator atoms" (e=1 → m distinct iso classes), the level-1
  layer. The decisive single arcs (highest iso-influence) are the *local* (skip-2)
  comparisons — the near-neighbor preferences — not the global source→sink apex (the
  apex is the most *silent* flip, highest self-loop). Locally inconsistent pairs do more
  to change the tournament's identity than the one long-range verdict.
- **Flip a whole sub-tiling (a triangle of arcs, or the complement)** = a coordinated
  even-level move: a 3-cycle indicator is invariant under reversing its own triangle,
  which is *why* cyclicity has no odd-level (dictator) component. Condorcet content is
  intrinsically level-≥2; it cannot be read off any single arc.

## Where it points (OPEN-Q-102)

- The OCF `H = I(Ω,2)` and Kalai's `P_rational = ¾ − ¾·Stab_{−1/3}[f]` are two spectral
  encodings of odd-cycle/Condorcet content. Is the OCF specialization (p1→1, p_odd→2,
  p_even→0) a noise-stability evaluation at a specific ρ, mirroring ρ=−1/3 for the
  3-cycle? If so, the forbidden H-values become *forbidden Condorcet-cyclicity spectra*.
- KKL gives a free corollary: every balanced tournament invariant has a backward arc of
  influence Ω(log m / m) — a "decisive comparison" existence theorem. Friedgut's junta
  theorem: low-total-influence invariants depend on few arcs. Which project invariants
  are juntas?
- The multi-vertex Möbius sieve for H (not just c3): truncates at depth D=2⌊(n−1)/2⌋
  (the band limit) — a finite, well-posed identity to write down (HYP-2534).

The dispatch asked for the pattern structure of tournaments under perturbation. The
pattern is a vote: the boring ground state is consensus/rationality, the structure is
disagreement/Condorcet cyclicity, FKN says low-complexity outcomes are near-dictatorial,
and the project's whole odd-cycle apparatus is the spectroscopy of social paradox.
