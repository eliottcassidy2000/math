---
id: THM-1165
title: THE FOUR-SPEED CASE IS EXHAUSTIVE — AND THE SUBSTITUTION METHOD CAN NEVER SETTLE THE GLOBAL QUESTION — a uniform per-level prune (with j combs left and all remaining speeds ≥ r, no component longer than j/(r(7−j)) can be covered, sharpened to 2λ/r at j = 1) kills 1914812 of 2242028 nodes (85%), making k = 4 feasible: 414499 leaf checks, 0 empty residues, and again exactly TWO tight families. But the density coefficient 1 − 2kλ vanishes at k = 7 = 1/(2λ), so this method reaches at most SIX substitutions — families differing from {1,…,13} in seven or more speeds are permanently out of its reach, and "exactly two tight families overall" is therefore not merely unproved but UNREACHABLE by this route
status: PROVED exhaustive for one-, two-, three- and four-speed substitutions on {1,…,13}, all giving only the two known families. Node growth is ~25× per level (87863 at k=3, 2242028 at k=4), projecting ~56M nodes at k=5 (slow but feasible) and ~1.4B at k=6 (infeasible in this implementation). The k ≤ 6 ceiling is a proved property of the density bound, not an implementation limit
source: opus-2026-07-17-S381 (owner: prove the four-speed case exhaustive)
depends_on: [THM-1155 (three-speed and the k=7 ceiling), THM-1135 (two-speed), THM-1125 (the swap bound), THM-1120 (the tight locus), LRC(11)/(12)/(13) citations]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCEssentialRegion.lean (multi_speed_density_bound supplies every level)
scripts: 04-computation/fourspeed_exhaustive_opus_S381.py -> 05-knowledge/results/fourspeed_exhaustive_opus_S381.out
---

# THM-1165 — four speeds, and the horizon of the method

## The uniform per-level prune

At a node with residue E, j combs still to place and all remaining speeds
≥ r, the density bound with every speed equal to r gives the most generous
possible budget:

> ℓ(1 − 2jλ) ≤ 2λ·(j/r)  ⟹  **ℓ ≤ j / (r(7 − j))**  at λ = 1/14

so a node with ℓ_max(E) exceeding that is dead — soundly, since no
admissible assignment could do better. At j = 1 the sharper single-arc
bound 2λ/r (THM-1125) replaces it. Budgets: **1/(7r), 2/(5r), 3/(4r),
4/(3r)** for j = 1, 2, 3, 4.

Applying this at **every** level rather than only at the root is what makes
k = 4 tractable.

## The result

| quantity | value |
|---|---|
| nodes visited | 2242028 |
| pruned | **1914812 (85%)** |
| empty residues (would contradict LRC(11)/(12)/(13)) | **0** |
| leaf checks | 414499 |
| distinct tight families | **2** |

Again **{1,…,13}** and **{1,…,11,13,24}**, and nothing else.

## The horizon — what this method cannot do

The density coefficient 1 − 2kλ is positive only for k < 1/(2λ) = 7. So the
whole substitution programme — S378 through S381 — reaches **at most six**
simultaneous substitutions. Families differing from {1,…,13} in **seven or
more** speeds are outside its reach entirely, and no amount of compute
changes that: the bound that makes each level finite simply stops existing.

This is worth stating plainly because it re-scopes THM-1120's open
question. "Exactly two tight families overall" is not merely unproved by
this route — it is **unreachable** by it. Settling it needs an argument that
does not proceed by bounded substitution from the classical family.

Practical note on the remaining reachable levels: node counts grow ~25× per
level (87863 → 2242028 from k = 3 to k = 4), projecting ~56M nodes at k = 5
— slow but feasible — and ~1.4B at k = 6, which this implementation will not
manage. So k = 5 is the last level worth attempting as written.

## Status

**Proved:** no substitution of one, two, three or four speeds on {1,…,13}
yields a tight family beyond the two known. **Open and reachable:** k = 5, 6.
**Open and unreachable by this method:** everything at distance ≥ 7,
including the global classification.
