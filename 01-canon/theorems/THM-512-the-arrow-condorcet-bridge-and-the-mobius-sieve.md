---
id: THM-512
title: The Arrow–Condorcet bridge (the tiling cube IS the social-choice cube; transitive = rational/dictator, 3-cycle = Condorcet paradox = the OCF obstruction; c3 = the Guilbaud quadratic) + the Möbius vertex-deletion sieve (A+B+C−D−E−F+G)
status: PROVED/VERIFIED + ADVERSARIALLY CHECKED (two external sign typos fixed: Kalai P_rational=¾−¾Stab_{-1/3}; Guilbaud=¾+¾(2/π)arcsin(1/3)≈0.9123; n=4 sieve overshoot characterized; all internal identities reconfirmed over ALL tournaments n≤6). The Arrow/Condorcet identifications are exact (Kalai 2002 Adv.Appl.Math 29:412; Mossel 2009/Keller 2010; FKN 2002 external; the project-side identities Var(c3)=3C(m,3)/16, P_transitive=m!/2^{C(m,2)}, c3 ⊂ Walsh-level-2 verified m≤7); the Möbius sieve Σ(-1)^{|A|+1}φ(T-A)=φ(T) verified n≤6 (n=4 edge case noted). Builds on THM-511.
source: mac-mini-2026-06-15-S3
depends_on:
  - THM-511   # converse-parity = Fourier-level-parity = ranking(level-1)/cyclic(even) split = FKN dictator hierarchy
  - THM-002   # H = I(Ω,2), OCF as the conflict-graph independence polynomial
  - THM-163   # c3 ⊂ Walsh level 2 (Var(c3)=3C(n,3)/16); H band-limited to D=2⌊(n-1)/2⌋
related:
  - HYP-2532  # H's even-level Walsh weights = OCF α_k strata
  - HYP-2533  # FKN-defect of H monotone in n
  - HYP-2534  # the Möbius sieve / Walsh-top duality
  - OPEN-Q-102
---

# THM-512 — the Arrow–Condorcet bridge and the Möbius sieve

The dispatch ("consider FKN deeply; the tiling default is transitive, perturbations are
arc flips") lands on the **home of FKN**: Kalai's Fourier proof of Arrow's theorem,
which lives on *exactly* the tournament cube. THM-511 already identified the tiling
model's converse-parity with the FKN level-1(score)/even(cycle) dictator hierarchy; this
theorem completes the bridge to social choice and adds the Möbius reconstruction the
dispatch describes.

## A. The Arrow–Condorcet identification (the deepest reframe)

Kalai's setting for Arrow/Condorcet is the cube **{±1}^{C(m,2)}** of pairwise
comparisons among m candidates — which IS the (full) tournament cube. The dictionary is
a literal identity, not an analogy:

| social choice (Kalai/Arrow/Condorcet) | tournament / project |
|---|---|
| the pairwise-comparison cube {±1}^{C(m,2)} | the arc cube / tiling cube Q_{C(n-1,2)} |
| a **rational** (transitive) outcome | the **transitive tournament** (c3=0, H=1, the ground state) |
| a **Condorcet paradox** | a **directed 3-cycle** = the minimal odd cycle = the OCF obstruction |
| Arrow's **dictator** | the H=1 transitive corner (one decisive coordinate; THM-511 level-1) |
| "irrationality" / cyclicity | the OCF content: H = I(Ω,2), the odd-cycle conflict graph |
| Kalai's P_rational = ¾ − ¾·Stab_{−1/3}[f] | a noise-stability sibling of the OCF |

**Robust Arrow = FKN = the project's transitive/cyclic split.** Kalai (via FKN): if a
neutral aggregation rule produces a non-transitive (Condorcet-cyclic) outcome with
probability ≤ ε, it is O(ε)-close to a dictator ±x_i; Mossel (2009, arXiv:0903.2574) makes this uniform in m, n; Keller
(arXiv:1003.3956, tight version, δ=Cε³) sharpens it. Translated: **a tournament with few 3-cycles is near-transitive (near the
H=1 corner)** — the FKN-defect of THM-511 (level-≥2 mass / variance, 0.25→0.51→0.68 at
n=4,5,6) is precisely "how Condorcet-cyclic vs rational." Arrow's ε=0 ("always rational
⟹ dictator") is the statement that c3≡0 forces transitivity.

**c3 is the Guilbaud quadratic.** Over uniform tournaments: E[c3]=C(m,3)/4,
**Var(c3)=3·C(m,3)/16** (verified m=3..7), sitting *entirely at Walsh level 2*
(THM-163/THM-511; each triple's cyclicity has Fourier mass only on its 3 arc-pairs, no
degree-3 term). This is the project face of Guilbaud's theorem: the probability that
3-candidate majority is rational tends to **Guilbaud's number ≈ 0.9123** =
¾ + ¾·(2/π)arcsin(1/3) = (3/2π)arccos(−1/3) = ¾ − ¾·Stab_{−1/3}(Maj_∞), an
arcsine/level-2 quantity. P(a random m-tournament is
transitive) = m!/2^{C(m,2)} → 0 fast (0.75, 0.375, 0.117, 0.022, 0.0024 for m=3..7):
**Condorcet cyclicity is generic; transitivity (rationality) is exponentially rare** —
the social-choice reading of "almost all tournaments are far from the ground state."

## B. The Möbius vertex-deletion sieve (the dispatch's A+B+C−D−E−F+G)

For an additive sub-structure invariant φ (e.g. c3 = #directed triangles), anchoring one
vertex and sieving over the other n−1:

> **Σ_{∅≠A ⊆ [n−1]} (−1)^{|A|+1} φ(T−A) = φ(T)** for n ≥ 5 (verified over ALL
> tournaments at n=5,6). At n=4 the sieve OVERSHOOTS by exactly the indicator
> [the all-deletable triple {1,2,3} is cyclic]: sieve = c3(T) + [deletable-triple
> cyclic] ∈ {c3, c3+1} (residual +1 on exactly 16/64 tournaments) — the lone triangle
> using all n−1=3 deletable vertices, which the c≥4 codim terms cannot reach.

The codim-c terms number **C(n−1,c)** with sign **(−1)^{c+1}**, and
Σ_c (−1)^{c+1} C(n−1,c) = 1 — the inclusion–exclusion / Möbius function of the Boolean
lattice B_{n−1} (one anchor + n−1 sieved vertices). At n=4 this is exactly
+3·(size-3) − 3·(size-2) + 1·(size-1) = **A+B+C−D−E−F+G**. So "a tournament is comprised
of its sub-tournaments" is literally the Möbius reconstruction, **dual to the
Walsh-Fourier expansion** (the cube function = Σ over subsets of its Fourier/Möbius
coefficients): the anchor = the chosen base-path / symmetry-breaking reference, and the
alternating sieve extracts the same content the Walsh transform organizes by degree.
(Verified n≤6 for φ=c3; the H version is the open multi-vertex sieve, HYP-2534, whose
truncation depth is the band limit D=2⌊(n−1)/2⌋ — no even-Walsh character survives past
it, THM-163/478.)

## C. The FKN dictator layer (ground-state perturbation structure)

In the fixed-base-path tiling cube {±1}^m, m=C(n−1,2), the all-aligned point is the
unique transitive ground state. At energy e=1 (one flipped arc) the number of distinct
iso classes equals **m = C(n−1,2)** (n=4,5,6 → 3,6,10): the m single-flips are pairwise
non-isomorphic — the **m level-1 "dictator atoms"** (organized in transpose pairs). The
iso-class diversity by backward-arc energy at n=6 is **1,10,25,40,43,41,36,31,23,6,1**
(peaks mid-band, mirroring the FKN-defect funnel). Per-tile influence is nearly uniform
(~0.93–0.98 at n=6) with NO dictator tile — the iso-class map is a near-even junta of
short-range arcs; the *decisive* single arcs (highest iso-influence) are the **local
(skip-2)** comparisons, not the global source→sink apex (which is highest self-loop).

## Provenance / what is new

THM-511 (kind-pasteur) established the converse-parity = FKN-level hierarchy and the
ranking/cyclic split; HYP-2532/2533 the α_k-strata and FKN-defect monotonicity. NEW
here: (i) the explicit **Arrow/Condorcet/Kalai social-choice bridge** (the tournament
cube IS the social-choice cube; 3-cycle = Condorcet paradox = OCF obstruction; c3 = the
Guilbaud level-2 quadratic; robust Arrow = FKN = transitive/cyclic split) — not
previously in canon; (ii) the **Möbius vertex-deletion sieve** formalizing the
dispatch's A+B+C−D−E−F+G as the Walsh-dual inclusion–exclusion; (iii) the **e=1 = m
dictator-atom** layer and the backward-arc energy sequences.
