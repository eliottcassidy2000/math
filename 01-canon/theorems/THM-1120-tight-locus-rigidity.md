---
id: THM-1120
title: THE TIGHT LOCUS OF LRC(14) IS EXTREMELY RIGID — exactly TWO families up to dilation survive every search: {1,…,13} and {1,…,11,13,24}, and they form a CLOSED PAIR under the involution 12 ↔ 24 (from {1,…,11,13,24} the only non-trivial tightness-preserving substitution is 24→12, straight back). Searches covering single substitutions to 120, two-speed multiplications, one-multiple-plus-one-free to 30, and 30 hill-climbs found no third. The mechanism is exact: {1,…,13} with i→r stays tight iff E_i ⊆ D_r where E_i = [0,1] \ ⋃_{j≠i} D_j is the ESSENTIAL REGION of speed i, and this criterion predicts the swappable speed to be exactly {12}, matching the search with no exceptions — while |E₆| = 0.0082 < |E₁₂| = 0.0122 shows size alone is NOT the criterion, position inside the finer comb is
status: the essential-region criterion is exact and its prediction matches the substitution search perfectly (predicts [12], search found [12]). "Only two tight families" is a SEARCH result over structured neighbourhoods, not a theorem — exhaustive search over 13-subsets is infeasible, and the hill-climb null result is weak evidence since the tight locus has measure zero in family space
source: opus-2026-07-17-S377 (owner: find more tight families)
depends_on: [THM-1115 (which found the second family), THM-1105 (the AP tightness classification this complements), THM-1035/1040 (the classical sieve at q=14, which puts every lonely point at p/14)]
scripts: 04-computation/tight_families_opus_S377.py, tight_chain_opus_S377.py, essential_region_opus_S377.py -> 05-knowledge/results/
---

# THM-1120 — the extremal objects, mapped

A family is **tight** when its uncovered measure is exactly 0: the danger
arcs cover [0,1] up to finitely many points, so the extremal gap is
exactly 1/14 and attained only on a null set. These are the extremal
objects of LRC(14). THM-1115 found a second one; this file maps the locus.

## The search

| search | scope | tight families found |
|---|---|---|
| single substitution i→r | r ≤ 120, all i | **only 12→24** |
| single multiplication i→i·m | m ≤ 12, all i | **only 12→24** |
| two-speed multiplication | m, m′ ∈ {2,3,4,5} | **none** |
| one multiple + one free | r ≤ 30 | nothing new |
| hill-climb to uncovered = 0 | 30 restarts × 300 steps | **none** |

Up to dilation, exactly **two** families survive:

> **{1, 2, …, 13}** (classical) and **{1,…,11, 13, 24}** (THM-1115)

and they form a **closed pair**: from {1,…,11,13,24} the only non-trivial
tightness-preserving substitution is 24→12, returning to the classical
family. The operation is dilation-equivariant — 2·{1,…,13} admits 24→48,
which is just 2·{1,…,11,13,24}.

## The mechanism, and why 12 is unique

Define the **essential region** of speed i in a family V:

> **E_i = [0,1] \ ⋃_{j≠i} D_j** — the part of the circle only speed i covers.

Then V with i replaced by r stays tight **iff E_i ⊆ D_r**. For {1,…,13}:

| i | 6 | 12 | 10 | 4 | 2 | 13 | 1 | 7 |
|---|---|---|---|---|---|---|---|---|
| \|E_i\| | 0.00816 | 0.01216 | 0.02410 | 0.02423 | 0.03022 | 0.03410 | 0.07143 | 0.08390 |
| E_i ⊆ D_2i | no | **YES** | no | no | no | no | no | no |

The criterion predicts the swappable speed to be exactly **{12}**, and the
independent substitution search found exactly **{12}**. No exceptions.

**Size is not the criterion.** Speed 6 has the *smallest* essential region
(0.00816, below 12's 0.01216) and is still not swappable. What matters is
whether E_i sits **inside the finer comb** D_{2i}: D_{2i} halves the arc
width at the old centres j/i while adding arcs at the midpoints, so the
swap survives only when the essential region happens to lie in the
surviving half. That is a positional condition, not a measure one.

## Status

The criterion is exact and verified. **"Exactly two tight families" is a
search result, not a theorem** — exhaustive search over 13-subsets is
infeasible, and the hill-climb null result carries little weight because
the tight locus has measure zero in family space, so a descent method
should not be expected to land on it. What the searches do establish is
that the locus is **thin and rigid**: structured neighbourhoods of the
classical family contain exactly one other point.
