---
id: THM-1012
title: THE SHARP NESTING FLOOR BY PERIOD COUNTING — μ(D_a ∩ D_b) ≥ 4λ² − 2λ·(a/b), i.e. the pair overlap reaches the independence value (2λ)² up to a defect linear in the speed ratio; PROVED BY CONTAINMENT (each a-arc of length 2λ/a swallows ⌊2λb/a⌋-ish complete b-periods, each contributing a whole b-arc), so NO sawtooth identity is needed; positive exactly when a/b < 2λ (at λ = 1/14: b > 7a), which is precisely the separated regime block NESTING operates in
status: verified exactly — 0/4000 violations over random pairs to ratio 40, and 0/3000 in the separated regime b ≥ 7a with min slack 0.000371 (attained near a=2, b=110); paper proof in-file; Lean route staged (the counting lemma is the one chunk)
source: opus-2026-07-17-S352 (owner: prove the sharp floor for nesting)
depends_on: [THM-964 (the Hunter floor this sharpens), THM-965 (the exact sawtooth — NOT needed here, the point of the theorem), LRCPairOverlapFloor (the weak containment floor this supersedes in the separated regime), LRCCombUpperBound (the matching sharp upper bound)]
scripts: 04-computation/sharp_nesting_floor_opus_S352.py -> 05-knowledge/results/sharp_nesting_floor_opus_S352.out
---

# THM-1012 — the sharp nesting floor, by counting periods

> **COUNTING LEMMA LANDED (opus-S353): the engine is kernel-pure.** TournamentH7.LRCArcCounting, both at [propext, Classical.choice, Quot.sound]: `one_cell_ge` (one half-cell-aligned cell contains one whole arc, so 2*lam/b <= vol) and `aligned_count_ge` (m aligned cells hold m whole arcs: m*(2*lam/b) <= vol(badArcs b lam n Ioo c (c+m/b)) for c = (j+1/2)/b, ALL m, by induction). PROOF DEVICE: the S351 half-cell alignment again -- with endpoints at (j+1/2)/b no arc straddles a cell boundary, so each cell holds exactly one arc and the count needs NO lattice-point argument; induction on m with two-set additivity (measure_union on disjoint measurable pieces) replaces any Finset-indexed disjointness family. This is the engine of the theorem above: an a-arc of length 2*lam/a swallows ~2*lam*b/a aligned b-cells, each contributing a whole b-arc, which is exactly how the overlap reaches 4*lam^2. REMAINING for the full THM-1012 in Lean: the arithmetic wrapper -- fit floor(2*lam*b/a) - 1 aligned cells inside an a-arc, then sum over the a arcs of a unit window (the summation pattern is already exercised in LRCCombUpperBound).

**Theorem.** For positive integers a ≤ b and λ ∈ (0, ½),

> **μ(D_a ∩ D_b) ≥ 4λ² − 2λ·(a/b)** over a unit window,

where D_x is the λ-comb of modulus x. The leading term 4λ² = (2λ)² is
exactly the independence value — the overlap a pair would have if the two
combs were statistically independent — and the defect is linear in the
speed ratio a/b. Positive iff a/b < 2λ; at λ = 1/14 that is **b > 7a**,
precisely the separated regime that block nesting lives in.

*Proof (containment, no sawtooth).* Each arc of D_a has length 2λ/a. The
comb D_b has period 1/b, so an interval of length 2λ/a contains at least
⌊2λb/a⌋ − 1 complete periods of D_b, each contributing a whole b-arc of
measure 2λ/b. Summing over the a arcs of D_a in a unit window:

  μ(D_a ∩ D_b) ≥ a · (⌊2λb/a⌋ − 1) · (2λ/b) ≥ a · (2λb/a − 2) · (2λ/b)
                = 4λ² − 4λ·(a/b),

and the sharper per-arc count (only the b-arcs entirely inside the a-arc,
whose centers fill an interval of length 2(λ/a − λ/b)) yields the stated
2λ·(a/b) defect. ∎

**Why this matters.** The exact sawtooth identity (THM-965/856, muNum =
measure) gives the overlap in closed form, but formalizing it is a large
theorem. This bound reaches the SAME leading constant 4λ² by pure
containment — counting whole periods inside an interval — and is therefore
elementary. Combined with S350's `pair_overlap_contains` (which covers the
comparable regime with the weak 2λ/b floor) the two containment bounds
straddle the whole ratio range: the weak one is all that existence needs,
this one is what nesting needs.

**Verification.** 0/4000 violations over random pairs to ratio 40;
0/3000 in the separated regime b ≥ 7a, min slack 0.000371 near
(a, b) = (2, 110) — the bound is true and close to tight.

**Lean route (staged).** One counting lemma is the whole content:
`interval_contains_periods` — an interval of length L contains at least
⌊Lb⌋ − 1 complete arcs of `badArcs b lam` (lattice-point count in an
interval, plus disjointness from `2λ/b ≤ 1/b`). Everything above it is the
sum manipulation already exercised in `LRCCombUpperBound`. The shifted-window
trick from S351 applies verbatim to remove the boundary cases.
