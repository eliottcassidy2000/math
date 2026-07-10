---
source: opus-2026-07-09-S201
status: INVESTIGATED removing the census native_decide via THM-665 window shrinking. FINDING (axiom-verified,
  decisive): THM-665 is ORTHOGONAL to the census — its own theorem says the a-priori/density-floor route
  NEVER fires on covering clusters, and the census IS the covering case. The foundational-axioms-only
  version ALREADY EXISTS (lrc14_grand_assembly_pure, [propext, Classical.choice, Quot.sound] only) at the
  cost of a BIGGER residual obligation (folds the <=22 covering families in). The native_decide census is
  doing real work; removing it while keeping the small open surface needs a MATHEMATICAL census proof, which
  THM-665 does not provide and kernel decide cannot compute (MISTAKE-135). Recommendation below.
tags:
  - lrc14
  - native-decide
  - window-census
  - thm-665
  - axiom-footprint
  - irreducible
---

# The window census native_decide is irreducible; THM-665 is orthogonal to it

**opus-2026-07-09-S201.** Owner: "work the THM-665 window shrinking to remove native_decide." I did the
investigation; the honest finding is that this route rests on a category error, and the tradeoff it's
reaching for is already resolved (in both directions) in the Lean tree. Axiom footprints below are machine-
verified, not estimated.

## THM-665 does NOT shrink the census window — it is the complementary regime

THM-665 (the sharp aliasing bound `|E_grid[W](V) − ∫W| ≤ TV(W′)/(12V²)`) certifies a good period on the
`V`-ruler once `V > V₀ ≈ 2.8·spread`. Its own **Consequence 2** states it plainly: *"The a-priori
existence route NEVER fires on covering clusters"* — covering velocity sets contain small speeds, so
`spread ≈ Vmax` and `V/spread ≈ 1 < 2.7`, and *"the entire covering case ... lives INSIDE the bounded
window."* THM-665 handles the **large-V dissociated** regime; the **window census IS the covering case**.
They are complementary, not nested. So THM-665 cannot shrink the census window — the census is exactly the
family THM-665's bound is proven not to reach.

## The tradeoff is real, unavoidable, and already resolved BOTH ways in Lean

The census (`winData22_ok` + `winData22_complete`, and `covering18_complete`) exists to **discharge the
≤ 22 covering families** by explicit witnesses, shrinking the grand assembly's open obligation to
`Vmax ≥ 23`. Two variants already coexist (`LRC14GrandAssembly.lean`), and I confirmed their footprints:

| top theorem | axioms (`#print axioms`) | residual obligation |
|---|---|---|
| `lrc14_grand_assembly` | `[propext, Classical.choice, Quot.sound, winData22_complete, winData22_ok]` | only `Vmax ≥ 23` (smaller) |
| `lrc14_grand_assembly_pure` | `[propext, Classical.choice, Quot.sound]` | ALL covering families incl. `≤ 22` (bigger) |

So the foundational-axioms-only proof **already exists** — it is `lrc14_grand_assembly_pure`. The price is
exactly the ≤ 22 covering families moving into the open obligation. You cannot have foundational-only AND
the smaller open surface *for free*: the native_decide census is the thing that buys the smaller surface.

Why can't analysis discharge the ≤ 22 families instead of the census? Because they are precisely where every
analytic route fails: covering ⟹ `V/spread ≈ 1` ⟹ THM-665's density floor never fires; and the ones that
actually reach the census branch are `ratio > 13` (the `spread13`/`GapFamily` branch already peels off
`ratio ≤ 13`), which forces `min = 1` (since `Vmax ≤ 22 < 13·2`). They all **contain the speed 1**, have no
known uniform witness, and are certified one-by-one — that is what a census *is*.

## The only honest way to remove the census native_decide, and the recommendation

Removing it while keeping the small open surface requires a **mathematical census proof**: that every
`ratio > 13`, covering, 13-tuple with `Vmax ≤ 22` (equivalently: containing 1) is lonely, *without*
enumerating them. That is genuine research — the small-`Vmax` covering sets are the hardest, most
structure-poor instances, which is exactly why the fleet censused them rather than proving them. THM-665
does not help; kernel `decide` is infeasible (C(22,13) = 497 420, measured > 13 h + OOM, MISTAKE-135).

**Recommendation to the owner (pick the axis you value):**
- **Completeness-first (status quo):** keep `lrc14_from_B5` over the native_decide census. Footprint =
  foundations + 2 `winData22` facts; open surface = the single `hB5` (`Vmax ≥ 23`). `native_decide` /
  `Lean.ofReduceBool` is standard and sound (as throughout Mathlib for finite censuses).
- **Foundational-purity-first:** use `lrc14_grand_assembly_pure`. Footprint = foundations only; open
  surface grows to include the ≤ 22 covering families.
- **Both (the real prize):** fund a mathematical proof of the ≤ 22 covering census. No short route known;
  not a `decide` swap, not a THM-665 corollary.

I did not fabricate progress on a route that cannot work. The genuine, small structural fact worth carrying
forward: the census's necessary domain is the `min = 1` (`ratio > 13`), covering, `Vmax ≤ 22` tuples — a
much thinner set than `C(22,13)`, and the right target for any future census-elimination proof.

## Ledger

- Axiom-verified: `lrc14_grand_assembly_pure` = foundations only; `lrc14_grand_assembly` adds the 2
  `winData22` native_decide facts. The foundational-only proof already exists (bigger obligation).
- THM-665 is orthogonal to the census (its Consequence 2: a-priori route never fires on covering clusters).
- Census native_decide is irreducible short of a mathematical census proof (research-hard) — not via THM-665,
  not via kernel `decide` (MISTAKE-135).
- Recommendation recorded; `LRC14CompletionAudit.lean` updated to point at the pure variant. → THM-665,
  THM-663, LRC14GrandAssembly (both variants), MISTAKE-135, hB5.
