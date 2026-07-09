---
source: klein-2026-07-09-S203
status: The good-period leg of LRC(14)'s covering case is now ASSEMBLED end-to-end in Lean (dispatch +
  concrete predicate + reach bridge), and the exercise PINS the true remaining blocker: THM-527 Part A (the
  Vmax-ruler embedding / "positive witness density ⟹ M(S) ≥ 1/14") is the SHARED open node for BOTH the
  good-period route AND the density-floor route (the skeleton's `hpartA`, opaque `witnessG2`). Finishing the
  proof = formalizing Part A, NOT the good-period branches (LEM-012 proved, LEM-013 verified, dispatch done)
  and NOT the finite free-gap link. New: LRCGoodPeriodReach.lean (sorry-free, kernel-pure) proves
  HasGoodPeriod ⟹ Mreach ≥ 1/14 modulo exactly two named links — `hlink` (free-gap extraction, finite) and
  `hembed` (= Part A). Reduces the good-period leg to one shared analytic node.
tags: [lrc14, good-period, lean, assembly, thm-527-part-a, ruler-embedding, endgame]
---

# The good-period leg is assembled — Part A (the ruler embedding) is the shared blocker

**klein-2026-07-09-S203.** Owner: keep working the good-period dichotomy assembly, aim to finish the proof.
I traced the good-period leg's Lean assembly to the metal and wired it end-to-end. The result reframes the
endgame: the good-period dichotomy is essentially DONE, and what's left to *finish the proof* is a single
node shared with the other route.

## The good-period leg, fully mapped

The covering case's good-period leg is a chain of Lean pieces, now connected:

1. **The dichotomy** (kps-S99 `LRCGoodPeriodDispatch.hasGoodPeriod_of_dichotomy`): for `k ≥ 7`, the near-AP
   branch (`k−6 ≤ L`, LEM-012) and dissociated branch (`L ≤ k−7`, LEM-013) TILE all longest-AP `L`, so a good
   period exists. The two branches enter as hypotheses (the cited math lemmas); the tiling is proven.
2. **The concrete predicate** (klein `LRCGoodPeriodMaxgap`): `HasGoodPeriod E Vmax := ∃ j, 7·maxCircGap > Vmax`,
   decidable, with `native_decide` certs on the hard 7-structured clusters where the arc-count route dies.
3. **The reach core** (kps-S31 `LRCGapReach`): a free real gap `> 1/7` gives a fast phase with `nearInt`
   clearance `> 1/14` from every tooth — sorry-free.
4. **The reach bridge** (klein-S203 `LRCGoodPeriodReach`, NEW, sorry-free, kernel-pure): `mreach_ge_of_goodPeriod`
   proves `HasGoodPeriod E Vmax ⟹ 1/14 ≤ Mreach v`, composing (2)→(3)→`Mreach_ge_of_witness`, **modulo exactly
   two named links**:
   - `hlink` — the **free-gap extraction**: a good period's residue gap IS a free real interval `> 1/7`.
     FINITE / combinatorial (sortedness of the residues + the max adjacent gap is free). Laborious in Lean
     (the list-sort API has shifted — `List.Sorted` moved — so `maxCircGap`'s free-gap property is ~200 lines
     of `mergeSort`/`foldl`/wrap plumbing) but requires no analysis. A clean next Lean task.
   - `hembed` — the **ruler embedding**: a fast phase clearing all teeth is realized at a real time `τ` with
     `minReach v τ ≥ 1/14`. This is **THM-527 Part A**.

## The reframe: Part A is the SHARED blocker

The skeleton (`LRCFourteenSkeleton`) reduces LRC(14) to `hfloor` (density floor) + `hpartA`
(`0 < witnessG2 ⟹ 1/14 ≤ Mreach`), with `witnessG2` an **opaque** placeholder and `hpartA` an open node whose
"formalization needs the `Vmax`-ruler period structure, the fast-phase gap" (its own comment). That `hpartA` is
**THM-527 Part A — the same ruler embedding my good-period `hembed` needs.** So:

> **Both the good-period route and the density-floor route bottleneck on ONE node: THM-527 Part A (the
> `Vmax`-ruler embedding of the fast-phase gap into a real witness time).** Everything else — the dichotomy,
> the branches, the density-floor census (`hfloor`, machine-checked `k=8..13`), the gap→clearance geometry
> (GapReach), the witness→Mreach sup bound — is proven or cited sorry-free.

So "finishing the proof" is NOT about the good-period branches (LEM-012 is proved; LEM-013 is verified and its
regime is exact-check / density-floor covered per the klein-S201 2×2), NOR the finite free-gap link. It is the
**shared Part-A ruler embedding** — the slow-fast realization / equidistribution `ρ_K → ρ*` that kps-S31's
GapReach explicitly defers. That is the one substantive analytic node between the current sorry-free surface and
a complete Lean proof.

## What this session delivered

- `LRCGoodPeriodReach.lean` (sorry-free, `[propext, Classical.choice, Quot.sound]`): `teeth`,
  `reachMargin_of_residueGap` (residue gap → `nearInt` margin via GapReach), and `mreach_ge_of_goodPeriod`
  (`HasGoodPeriod ⟹ Mreach ≥ 1/14` modulo `hlink` + `hembed`). The first file connecting the good-period
  predicate through to the concrete `Mreach`.
- `foldl_max_mem` verified (kernel-pure, `[propext]`) — the seed helper for the `hlink` extraction, for whoever
  takes the free-gap plumbing.
- The endgame reframe above: **Part A is the shared bottleneck.** Every agent working the covering case should
  aim there; the two routes converge on it.

Files: `04-computation/lean/TournamentH7/TournamentH7/LRCGoodPeriodReach.lean`. Builds on kps-S99 dispatch,
kps-S31 GapReach, klein LRCGoodPeriodMaxgap, mac-mini/kps LRCMreachConcrete. Connects the klein-S201 2×2
(good-period `V≥Q+1` ∪ density-floor/exact-check) to the Lean assembly.
