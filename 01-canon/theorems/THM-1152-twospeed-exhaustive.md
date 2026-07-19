---
id: THM-1152
title: THE TWO-SPEED SUBSTITUTION CASE IS EXHAUSTIVE — a density bound plus LRC(13) makes the double-substitution search FINITE, and it confirms no third tight family. Covering an interval of length ℓ by two combs forces ℓ(1−4λ) ≤ 2λ(1/r + 1/s), so the SMALLER new speed obeys r ≤ 2/(5·ℓ_max(E_{i,j})) — bounds 8…63 over the 78 pairs of {1,…,13}. The residue E′ = E_{i,j} \ D_r is never empty (0/2113 verified), because an empty residue would mean twelve speeds cover everything at radius 1/14, contradicting LRC(13); so the single-speed swap bound then caps the second speed at s ≤ 2λ/ℓ_max(E′) ≤ 62. Exhaustive enumeration of the resulting 6702 combinations yields 90 tight double substitutions but only TWO distinct families — {1,…,13} and {1,…,11,13,24} — with 90 = 78 identities + 12 routings through 12→24
status: PROVED for one- and two-speed substitutions on {1,…,13}. The density bound is formalised (`two_speed_density_bound`, kernel-pure); the enumeration is exhaustive because both speed bounds are proved, and the non-emptiness of E′ rests on LRC(13), a citation hypothesis for this project. NOT proved: three-or-more substitutions, or families far from {1,…,13} — so "exactly two tight families overall" remains open
source: opus-2026-07-17-S379 (owner: prove the two-speed case exhaustive too)
depends_on: [THM-1125 (the single-speed swap bound and essential-region criterion), THM-1120 (the tight locus this closes further), FragmentationLemma.fragmentation (the window bound supplying the density estimate), LRC(13) citation]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCEssentialRegion.lean
scripts: 04-computation/twospeed_bound_opus_S379.py -> 05-knowledge/results/twospeed_bound_opus_S379.out
---

# THM-1152 — closing the double-substitution neighbourhood

THM-1125 made SINGLE substitutions exhaustive. Two-speed is genuinely
harder: a component of the joint essential region

> E_{i,j} = I \ ⋃_{k ≠ i,j} D_k

may be covered by arcs of **both** new speeds, so the single-arc argument
fails outright. What replaces it is a density bound and a termination
argument.

## Step 1 — the density bound caps the smaller speed

The fragmentation window lemma gives μ(D_w ∩ I) ≤ 2λ|I| + 2λ/w. Covering a
component of length ℓ by two combs therefore forces

> ℓ ≤ (2λℓ + 2λ/r) + (2λℓ + 2λ/s),  i.e.  **ℓ(1 − 4λ) ≤ 2λ(1/r + 1/s)**

At λ = 1/14 this is ℓ ≤ (1/5)(1/r + 1/s). With r ≤ s it yields
**r ≤ 2/(5·ℓ_max(E_{i,j}))**. Over the 78 pairs of {1,…,13} the bounds run
**8 … 63**, worst at (i,j) = (6,10) with ℓ_max = 0.006279.

Note this bounds only the *smaller* speed — as s → ∞ the constraint
degenerates to ℓ ≤ 1/(5r), which is exactly right, since a very fine comb
has density 2λ and cannot help cover a specific interval.

## Step 2 — LRC(13) makes the recursion terminate

After fixing r, the residue E′ = E_{i,j} \ D_r must satisfy E′ ⊆ D_s. If
E′ were **empty**, the twelve speeds {1,…,13}∖{i,j} together with r would
cover [0,1] entirely at radius 1/14 — impossible, since LRC(13) guarantees
twelve speeds always leave a point of gap ≥ 1/13 > 1/14. Verified: **0
empty residues in 2113 (pair, r) combinations**, as predicted.

## Step 3 — the single-speed bound caps the second

E′ nonempty ⟹ **s ≤ 2λ/ℓ_max(E′)** by THM-1125. Maximum over all cases:
**62**. Both speeds bounded, so the search is finite.

## The exhaustive result

**6702** (i, j, r, s) combinations enumerated; **90** tight double
substitutions; **2** distinct families:

> **{1,…,13}** and **{1,…,11, 13, 24}**

and the count decomposes exactly: 90 = 78 identity pairs + 12 pairs routing
through 12→24, which is the arithmetic one expects and a useful check on
the enumeration.

## Status

**Proved:** on {1,…,13}, no substitution of one or two speeds produces a
tight family other than the two known. That upgrades the corresponding half
of THM-1120 from search to theorem.

**Not proved:** three-or-more substitutions, and families far from
{1,…,13}. So "exactly two tight families overall" remains open, and the
hill-climb evidence for it is still weak for the reason given in THM-1120 —
the tight locus has measure zero in family space.

## Lean

`two_speed_density_bound` in `LRCEssentialRegion.lean`: from
`Set.Icc x (x+L) ⊆ badArcs r lam ∪ badArcs s lam` it derives
`L * (1 - 4*lam) ≤ 2*lam*(1/r + 1/s)`, via measure subadditivity and two
applications of `fragmentation`. The module now carries **nine** theorems,
zero sorries, all depending only on `[propext, Classical.choice,
Quot.sound]`.
