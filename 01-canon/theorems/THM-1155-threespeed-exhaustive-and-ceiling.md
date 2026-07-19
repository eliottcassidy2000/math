---
id: THM-1155
title: THE THREE-SPEED CASE IS EXHAUSTIVE, AND THE METHOD'S CEILING IS EXACTLY k = 7 = 1/(2λ) — the general k-comb density bound L(1 − 2kλ) ≤ 2λ·Σ(1/wᵢ) is now formalised for an arbitrary finite speed set, and it is informative precisely while 1 − 2kλ > 0, i.e. k < 1/(2λ): at λ = 1/14 the method reaches SIX simultaneous substitutions and degenerates exactly at SEVEN. At k = 3 it gives r ≤ 3/(4·ℓ_max(E_{i,j,k})) (bounds 14…102 over the 286 triples), then the S379 two-comb bound caps s and the S378 swap bound caps t; termination is LRC(12) and LRC(13) respectively (11 and 12 speeds cannot cover at radius 1/14), verified with 0 empty residues throughout. Exhaustive enumeration of 71890 combinations over 87863 branches yields the SAME TWO families — {1,…,13} and {1,…,11,13,24}
status: PROVED for one-, two- and three-speed substitutions on {1,…,13}. The general bound `multi_speed_density_bound` is formalised and kernel-pure; the enumeration is exhaustive because every speed bound in the chain is proved, with the non-emptiness of the residues resting on LRC(≤13) as a citation hypothesis. NOT proved: four-or-more substitutions (though the method survives to six), or families far from {1,…,13}
source: opus-2026-07-17-S380 (owner: prove the three-speed case exhaustive)
depends_on: [THM-1135 (two-speed), THM-1125 (the swap bound and essential-region criterion), THM-1120 (the tight locus), FragmentationLemma.fragmentation, LRC(12)/LRC(13) citations]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCEssentialRegion.lean
scripts: 04-computation/threespeed_bound_opus_S380.py, threespeed_exhaustive_opus_S380.py -> 05-knowledge/results/threespeed_bound_opus_S380.out
---

# THM-1155 — three speeds, and where the method stops

## The general bound and its ceiling

Each comb meets an interval of length L in measure at most 2λL + 2λ/w, so
covering by k combs forces

> **L(1 − 2kλ) ≤ 2λ · Σᵢ (1/wᵢ)**

which is informative only while **1 − 2kλ > 0**, i.e. **k < 1/(2λ)**. At
λ = 1/14:

| k | 1 2 3 4 5 6 | 7 |
|---|---|---|
| 1 − 2kλ | 6/7, 5/7, 4/7, 3/7, 2/7, 1/7 | **0** |

So the density method reaches **six** simultaneous substitutions and dies
**exactly at seven** — the same 7 = 1/(2λ) that governs S₁ = 13/7, the
forbidden window |W_q| = 1 ⟺ q ≤ 14, and every union bound this program has
tried. It is the same constant appearing yet again, now as the arity limit
of a covering argument.

## The three-speed chain

1. **k = 3 density bound:** ℓ ≤ (1/4)(1/r + 1/s + 1/t), so with r ≤ s ≤ t,
   **r ≤ 3/(4·ℓ_max(E_{i,j,k}))** — bounds **14…102** over the 286 triples,
   worst at (4,5,6) with ℓ_max = 0.007305.
2. **Residue E′ = E_{i,j,k} \ D_r** must be covered by two combs: the S379
   bound caps s. Termination by **LRC(12)** — eleven speeds cannot cover at
   radius 1/14, since they leave a gap ≥ 1/12 > 1/14.
3. **Residue E″ = E′ \ D_s** must be covered by one: the S378 swap bound
   caps t. Termination by **LRC(13)**.

A hard prune makes this feasible: with r ≤ s ≤ t the best two-comb budget
for E′ is at s = t = r, giving ℓ ≤ 2/(5r); if ℓ_max(E′) exceeds that, the
branch dies immediately. It killed **3496** (triple, r) pairs.

## The exhaustive result

| quantity | value |
|---|---|
| branches explored | 87863 |
| pruned by the 2/(5r) test | 3496 |
| empty E′ (would contradict LRC(12)) | **0** |
| empty E″ (would contradict LRC(13)) | **0** |
| (i,j,k,r,s,t) fully checked | **71890** |
| distinct tight families | **2** |

The two are **{1,…,13}** and **{1,…,11,13,24}** — no third. The zero empty
residues are a live check on the citation argument, not decoration: any
empty one would have falsified the LRC(≤13) reasoning that makes the
recursion terminate.

## Status

**Proved:** on {1,…,13}, no substitution of one, two or three speeds
produces a tight family other than the two known. **Not proved:** four or
more substitutions — though the density method survives to six, so k = 4, 5,
6 are reachable by the same chain at higher enumeration cost — nor families
far from {1,…,13}. "Exactly two tight families overall" remains open.

## Lean

`multi_speed_density_bound` (arbitrary `Finset ℕ` of speeds) joins the
module, proved by measure subadditivity over the finset plus a per-comb
application of `fragmentation`. Ten theorems, zero sorries, all depending
only on `[propext, Classical.choice, Quot.sound]`.
