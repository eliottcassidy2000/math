---
id: THM-1125
title: THE ESSENTIAL-REGION CRITERION, PROVED AND FORMALISED — AND THE SWAP BOUND THAT MAKES THM-1120 A THEOREM — replacing speed i by r leaves uncovered EXACTLY E_i \ D_r where E_i = I \ ⋃_{j≠i} D_j, so the swap covers iff E_i ⊆ D_r (an IFF at the level of sets, both directions formalised). The proof is Set.diff_diff. The separation lemma (consecutive arcs strictly apart, gap exactly (1−2λ)/w) then yields the SWAP BOUND: a closed interval inside badArcs r has length ≤ 2λ/r, so a component of E_i of length ℓ forces r ≤ 2λ/ℓ. Computing ℓ_max(E_i) over {1,…,13} gives bounds 4…52, all below the r ≤ 120 that S377 searched — so that search was EXHAUSTIVE and "12→24 is the only non-trivial single substitution" is now proved, not merely observed. Eight theorems land kernel-pure in Lean, zero sorries
status: PROVED and formalised. LRCEssentialRegion.lean builds green in the project tree (root import added); 8 theorems, 0 sorries, all depending only on [propext, Classical.choice, Quot.sound]. Cross-checks: criterion vs direct computation on 623 swaps = 0 mismatches; swap bounds computed exactly per speed. NOTE the exhaustiveness upgrade covers SINGLE substitutions on {1,…,13} only — the two-speed and global searches of THM-1120 remain searches
source: opus-2026-07-17-S378 (owner: prove the essential-region criterion; extended)
depends_on: [THM-1120 (found the criterion empirically; its single-substitution claim is upgraded here), THM-1115 (the second tight family), FragmentationLemma.lean (the badArcs definition)]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCEssentialRegion.lean
scripts: 04-computation/essential_criterion_check_opus_S378.py, swap_bound_opus_S378.py -> 05-knowledge/results/
---

# THM-1125 — the criterion, proved, and what it then proves

## Statement and proof

For an ambient interval I, family S, speed i ∈ S, replacement r, define the
**essential region** E_i = I \ ⋃_{j≠i} D_j — the part of I no *other* speed
covers. Then

> **I \ ⋃_{v ∈ (S∖{i})∪{r}} D_v  =  E_i \ D_r**

so the swapped family covers I **iff** E_i ⊆ D_r, and is tight iff
μ(E_i \ D_r) = 0.

*Proof.* With A = ⋃_{j≠i} D_j and B = D_r, `Set.diff_diff` gives
I \ (A ∪ B) = (I \ A) \ B = E_i \ B. ∎

The criterion is **De Morgan** — no depth. Its value is *reuse*: E_i is
computed once per speed, then tested against every candidate r. And it is an
**iff at the level of sets**, so a failed containment *proves* non-tightness.

## The swap bound — turning search into proof

Consecutive arcs of one modulus are strictly separated, with gap exactly
**(1 − 2λ)/w > 0**. Hence a closed interval contained in badArcs r lies in a
single arc, so its length is at most the arc length 2λ/r. Applied to the
longest component ℓ_max of E_i:

> **r ≤ 2λ / ℓ_max(E_i)** — an explicit finite bound on admissible replacements.

Computed on {1,…,13} (2λ = 1/7):

| i | 1 | 2 | 6 | 12 | 13 |
|---|---|---|---|---|---|
| ℓ_max(E_i) | 0.03571 | 0.01511 | 0.00271 | 0.00408 | 0.00595 |
| bound on r | 4 | 9 | **52** | 35 | 24 |
| admissible r | {1} | {2} | {6} | **{12, 24}** | {13} |

The largest bound is **r ≤ 52**. S377 searched r ≤ 120, which exceeds every
bound — so **that search was exhaustive**, and

> **on {1,…,13}, 12→24 is the ONLY non-trivial single substitution preserving
> tightness** — now a theorem, not a search result.

Scope: this upgrades the *single-substitution* claim. THM-1120's two-speed
and hill-climb searches remain searches, so "exactly two tight families
overall" is still unproved.

## Lean

`LRCEssentialRegion.lean`, in the build tree with a root import:

| theorem | content |
|---|---|
| `uncovered_swap_eq` | the set identity — the criterion |
| `tight_swap_iff` | measure-zero form |
| `swap_covers_of_subset` / `subset_of_swap_covers` | both directions of the iff |
| `badArcs_consecutive_separated` | arcs strictly apart when 2λ < 1 |
| `badArcs_gap_length` | gap is exactly (1 − 2λ)/w |
| `endpoint_not_mem_badArcs` | an arc endpoint escapes the whole union |
| `Icc_length_le_of_subset_badArcs` | **the swap bound**: interval ⊆ badArcs r ⟹ length ≤ 2λ/r |

Eight theorems, zero sorries, all depending only on `[propext,
Classical.choice, Quot.sound]`.

## Cross-check

Criterion vs direct computation of uncovered measure over all 623 swaps
(i ∈ {1,…,13}, r ≤ 60): **0 mismatches**, the tight swaps being the
identities plus the single non-trivial 12→24.
