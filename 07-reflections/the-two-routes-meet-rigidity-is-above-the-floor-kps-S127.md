# The two routes meet — the primitive rigidity sits above opus's floor

*kind-pasteur-2026-07-10-S127. Owner: "wire the primitive collapse into opus's dissociated peel." This
note records where the two open routes to LRC(14) meet, and what the wiring reveals about their relative
strength.*

---

## Two agents, two strategies, one surface

opus and I have been closing LRC(14) from opposite directions onto the same object — the **primitive
residual family** (`tupleGcd v = 1`, covering, ratio `> 13`, …):

- **opus (S206–S209)** peels. `lrc14_grand_assembly_primitive` removes the dilates; `LRCDissociatedAssembly`
  hands the `d = 2,3` near-dilate μ-minimizers to THM-678 (`MultiDetunedDispatch`), leaving the analytic
  floor obligation on the *dissociated* remainder, where μ decorrelates toward `(6/7)^13`. The target
  shrinks to `SafeMeasureFloorPrimitive` — *"every primitive residual family has a positive-measure safe
  set"* — the weakest floor in the corpus.
- **I (S126–S127)** rigidify. `TightRigidity` says `μ = 0 ⟹` the family is a dilated AP; the difference-
  primitive collapse (`dilated_primitive_eq_range`) sharpens that on the primitive class to `μ = 0 ⟹ image
  |v| = {1,…,13}` *exactly*.

The wiring is one observation: **opus's clause 2 (`tupleGcd v = 1`) is literally my `Primitive v`**
(`primitive_of_tupleGcd_one`). So the collapse applies verbatim on opus's surface.

## What the collapse does on the floor

On the primitive residual, the collapse pins the entire tight locus to the single set `{1,…,13}`. But every
residual family is `GapFamily` (ratio `> 13`), and `{1,…,13}` has ratio exactly `13`
(`gapFamily_image_ne_range`, the AP-case of `not_dilated_of_gapFamily`). So **no primitive residual family
is tight** — and

```
PrimitiveTightRigidity  →  SafeMeasureFloorPrimitive       (safeMeasureFloorPrimitive_of_primitiveTightRigidity)
                        →  ResidualObligationPrimitive
                        →  ResidualObligationDissoc          (a fortiori — dissociation clause unused)
                        →  LRC14Statement                    (opus's dissociated assembly).
```

## The reveal: rigidity is strictly above the floor

The wiring makes the *relative strength* of the two routes precise, and it is not symmetric. My rigidity
discharges the **whole** primitive residual — dissociated or not — so the THM-678 dispatch is **redundant on
this route**: `lrc14_of_primitiveTightRigidity` closes LRC(14) with **no `hMD`** at all. The tight rigidity
is *sufficient overkill*. opus's `SafeMeasureFloorPrimitive` asks only for `μ > 0`; my
`PrimitiveTightRigidity` asks for the far stronger `μ = 0 ⟹ {1,…,13}` and then throws most of it away (it
only ever uses the `{1,…,13}`-vs-`GapFamily` contradiction).

That is the honest moral for the fleet, and a correction to any temptation to chase my own rigidity:
**`SafeMeasureFloorPrimitive` is the minimal analytic target.** Proving `μ > 0` on the primitive residual
is *weaker* than the extremal uniqueness and closes the same theorem. Prove the floor, not the rigidity.
opus's peeling strategy — shrink the surface until a plain positive-measure floor suffices — is aimed at
exactly the right, minimal object; my rigidity is a ceiling above it, useful for *locating* the difficulty
(`≥ LRC-hard`) but not the thing to spend proof effort on.

## The honest boundary

`PrimitiveTightRigidity` is the S127 open conjecture; everything here is a **reduction between open
hypotheses**, kernel-pure `[propext, Classical.choice, Quot.sound]`, not a proof of any floor. What is now a
machine-checked theorem is the *ordering*: `PrimitiveTightRigidity ⇒ SafeMeasureFloorPrimitive ⇒ opus's
dissociated obligation ⇒ LRC(14)`, with the collapse (`{1,…,13}`) load-bearing at the top and opus's floor
the weakest link — the one worth attacking.

*Files: `LRCDissociatedRigidity.lean` (the wiring), building on `LRCTightRigidity.lean` (mine,
`dilated_primitive_eq_range`) and opus's `LRCDissociatedAssembly.lean` / `LRCResidualMeasureFloorPrimitive.lean`.
Continues [[the-difference-primitive-case-collapses-the-conclusion-not-the-wall-kps-S127]].*
