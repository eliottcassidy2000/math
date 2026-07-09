---
source: opus-2026-07-09-S175
status: DISCHARGED klein-S203's hlink (the free-gap extraction) fully in Lean -- sorry-free, kernel-pure
  [propext, Classical.choice, Quot.sound]. The mergeSort ARGMAX (foldl_max_mem + mem_zipWith_sub_tail +
  exists_gap_decomp) + BOTH freeness branches (interior via kps-S101, the WRAPPING gap directly) assemble
  into hlink_extract: HasGoodPeriod E Vmax => a real interval (a,a+g), g>1/7, free of every integer
  tooth-translate. mreach_ge_of_goodPeriod_of_hembed then shows the good-period -> Mreach>=1/14 chain
  needs ONLY hembed (THM-527 Part A ruler embedding, the shared blocker). Completes the "mechanical
  assembly" kps-S104 flagged as opus-active; converges with kps's parallel LRCHlinkList toolkit.
tags:
  - lrc14
  - good-period
  - hlink
  - lean
  - mergesort-argmax
  - formalization
---

# hlink discharged: the mergeSort argmax and the wrapping gap

**opus-2026-07-09-S175.** Owner: keep working the mergeSort argmax and wrapping-gap case for `hlink`.
Done — `TournamentH7.LRCHlinkExtract`, sorry-free and kernel-pure.

## What was left, and what this closes

klein-S203's `LRCGoodPeriodReach.mreach_ge_of_goodPeriod` reduced the good-period leg to two links:
`hlink` (a good period yields a free residue gap `> 1/7`) and `hembed` (the ruler embedding, THM-527
Part A).  kps-S101 discharged the NON-wrapping freeness core; the remaining `hlink` content — the owner's
two named pieces — was the `mergeSort` ARGMAX extraction and the WRAPPING gap.  Both now proved.

### 1. The mergeSort argmax (the hard, novel part)

`maxCircGap E Vmax j = foldl max 0 (zipWith (·-·) cyc cyc.tail)`, `cyc = ps ++ [p0+Vmax]`, `ps` the
sorted residues.  Three reusable lemmas turn a POSITIVE `maxCircGap` into an adjacent gap pair:

- `foldl_max_mem` : `0 < foldl max 0 L → foldl max 0 L ∈ L` (induction, `max_choice`).
- `mem_zipWith_sub_tail` : `d ∈ zipWith (·-·) c c.tail → ∃ l₁ l₂ x y, c = l₁ ++ x::y::l₂ ∧ y−x = d`
  (an ADJACENT-pair decomposition — dodges the `List.get` index pain by structural recursion).
- `exists_gap_decomp` : unfolds `maxCircGap`'s `match`, applies the two above ⟹ `cyc = l₁ ++ x::y::l₂`
  with `y − x = maxCircGap`.

### 2. Both freeness branches

From `cyc` sorted (`pairwise_mergeSort'` + `pairwise_append`), the adjacent pair `x,y` has no residue
strictly between, and `y ∈ ps` (INTERIOR) or `y = p0+Vmax` (WRAPPING):

- **Interior** (`y < Vmax`): the interval `(x/Vmax, y/Vmax) ⊆ [0,1]` is tooth-free ⟹ kps-S101
  `free_translate_of_free_subInterval` gives the no-integer-translate conclusion.
- **Wrapping** (`y = p0+Vmax`): the interval `(ps.last/Vmax, p0/Vmax + 1)` overshoots `1`.  Directly:
  every residue `r` has `p0 ≤ r ≤ ps.last`, so `r/Vmax + n ∈ (a, a+g)` forces both `n ≥ 1` (from
  `r ≤ ps.last`) and `n ≤ 0` (from `p0 ≤ r`) — contradiction.  This is the piece kps-S101's
  non-wrapping lemma could not reach; the max/min-residue bound is what closes it.

`hlink_extract : (0 < Vmax) → HasGoodPeriod E Vmax → ∃ j a g, 1/7 < g ∧ ∀ c ∈ teeth, ∀ n, c+n ∉ (a,a+g)`
— exactly klein's `hlink`.

### 3. hlink discharged in the chain

`mreach_ge_of_goodPeriod_of_hembed` feeds `hlink_extract` into klein's chain: the good-period →
`Mreach ≥ 1/14` reduction now needs ONLY `hembed` (the ruler embedding, shared with the density route).
One of the good-period leg's two links is closed.

## Convergence + the paper

kps-S103/S104 built a PARALLEL toolkit (`LRCHlinkList`: `mem_zipWith_sub_adjacency`,
`sorted_adjacency_sep`, `tooth_not_in_gap` unifying both branches) and explicitly left "the mechanical
assembly (unfold `maxCircGap` match + teeth↔residue perm + internal/wrap dispatch)" as **opus-active** —
which is exactly what `hlink_extract` completes end-to-end.  Two independent routes to the same discharge
(kps's unified `tooth_not_in_gap` vs my two-branch dispatch); my file adds the full argmax→interval→chain
assembly.  The owner's paper (arXiv:2604.21187, *Doubly Saturated Ramsey Graphs*) is a case study in the
same SAT + LLM + **Lean formalization** paradigm this session enacted — decidable/finite structure driven
to a machine-checked proof.

## Ledger

- DISCHARGED klein-S203 `hlink`: `LRCHlinkExtract` (sorry-free, kernel-pure) — `foldl_max_mem`,
  `mem_zipWith_sub_tail`, `exists_gap_decomp` (mergeSort argmax) + `hlink_extract` (interior + WRAPPING
  freeness) + `mreach_ge_of_goodPeriod_of_hembed` (chain needs only `hembed`).  Wired into the aggregate.
- Converges with kps-S103/S104 `LRCHlinkList` (parallel toolkit; I completed the assembly they flagged
  opus-active).  Remaining good-period-leg blocker: `hembed` (THM-527 Part A ruler embedding).
- -> klein-S203 (`mreach_ge_of_goodPeriod`), kps-S101 (`free_translate_of_free_subInterval`),
  kps-S103/S104 (`LRCHlinkList`), `LRCGoodPeriodMaxgap` (`maxCircGap`/`HasGoodPeriod`).
