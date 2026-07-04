# The two-killer ladder: the first multi-killer certified slice of CoveringFarLonely

*klein-2026-07-04-S130 (HYP-4091). Owner: push toward closing the proofs. After mapping the
active LRC(14) Lean route to its single remaining obligation, I formalized a new infinite certified
class of that obligation — the first drop-two / two-killer family — extending kps's single-swap
hexad into the multi-killer stratum. Kernel-pure, corpus-green.*

## Where the proof actually stands (the map)

The active top surface `lrc14_of_covering_far_of_window` (mac-mini) reduces LRC(14) to exactly two
pieces:
- `hwindow` — every family with all `|vᵢ| ≤ 20` is lonely. **DONE** (`hwindow20_closed`, a finite
  `native_decide`).
- **`CoveringFarLonely 20`** — every *covering* family with a *far* entry (`|vᵢ| > 20` for some `i`)
  is lonely. **The single remaining analytic obligation.**

`CoveringFarLonely` is exactly the covering-min statement (`M ≥ 1/14`, in fact `≥ 14/183`) my whole
track has mapped. It is discharged by certifying *infinite families* that tile the covering-far
domain: the AP block families `{V,…,V−12}` (mac-mini), the deep well, and kps's **single-swap
hexad** (drop one `j ∈ {8,…,13}`, add one killer `lcm(j,14)·k` — `drop8/9/10/11/12/13_lonely` in
`LRCOneSwapLadders.lean`). The hexad handles families with **one** killer. The next stratum —
families whose killer role is **split across two runners** — was uncertified.

## The result: the two-killer ladder

> **`twoKiller K := {1,…,11, 14, 156·K}` is lonely at `t = 13K/(156K+1)`, margin
> `M = 13K/(156K+1) > 1/14`, for every `K ≥ 1`.** (`LRCTwoKillerLadder.lean`, kernel-pure
> `[propext, Classical.choice, Quot.sound]`, no `native_decide`.)

This drops **both** `12` and `13` and adds **two** killers: `14` (covers `q=14`) and
`156K = 12·13·K` (covers `q=12` and `q=13`). It is covering for every `K` and has the far entry
`156K > 20`, so it is a `CoveringFarLonely 20` instance the hexad does not reach — the first
**multi-killer** certified class.

The certificate is the same engine as the hexad (kps's `residue_key`/`lattice_dist_ge`): one
rational witness `p/Q = 13K/(156K+1)`, floor `κ = 13K`, and a residue table linear in `K`:
```
  v = 1,…,11 : v·(13K) = 0·Q + 13Kv          (κ = 13K ≤ 13Kv ≤ 143K ≤ Q−κ)
  v = 14     : 14·(13K) = 1·Q + (26K−1)        (26K−1 ∈ [13K, 143K+1])
  v = 156K   : 156K·(13K) = (13K−1)·Q + (143K+1)   (143K+1 = Q−κ, the top binding)
```
with `14κ = 182K ≥ 156K+1 = Q ⟺ 26K ≥ 1 ⟺ K ≥ 1`. The **binding pair is runner 1 (bottom) and
killer 156K (top)** — precisely the 2-point equioscillation THM-618 identifies, now with the base
`{1,…,11}` (optimum `1/12`) and a *split* killer set. So the two-killer families sit on the same
"base-optimum minus killer-offset" geometry as the single-killer ladder; splitting the killer just
shifts the base from `{1,…,12}` to `{1,…,11}` and the value from `14k/(182k+1)` to `13k/(156k+1)`.

## Why this is proof progress

`CoveringFarLonely` is discharged family-class by family-class. Before: block families + deep well +
single-swap hexad. Now: **+ the two-killer (drop-{12,13}) ladder** — a genuinely new infinite slice
in the multi-killer stratum, closing one more of the classes an eventual full classification needs.
The value `13/157` was one of the low rungs of the covering-min spectrum (klein-S128); it is now
Lean-certified.

## Honest scope

One new certified class, not the whole obligation. `CoveringFarLonely` in full still needs the
general covering-min proof — the residual being the `m=2` folding (opus) and the remaining
multi-killer classes (drop-3, other split patterns). But the technique is now visibly uniform: every
covering-far family is (empirically, klein-S128/S129 + opus-S71) a base-optimum-minus-killer-offset
2-point equioscillation with a residue-table certificate, and each stratum falls to one `residue_key`
lemma. The two-killer ladder is the proof-of-concept that the multi-killer strata are as
formalizable as the single-killer hexad.

## Links

- Lean: `04-computation/lean/TournamentH7/TournamentH7/LRCTwoKillerLadder.lean`
  (`twoKiller_lonely`, `twoKiller156_lonely`; kernel-pure). HYP-4091.
- Search: `04-computation/lrc14_regimeC_witness_klein_S130.py` (probe of the far-cut route's
  regime C — the *other* route's residual).
- Builds on: kps `LRCOneSwapLadders` (single-swap hexad) + `LRCResidueLiar` (`residue_key`,
  `lattice_dist_ge`); mac-mini `LRC14CoveringFarSurface` (`CoveringFarLonely`, block families) +
  THM-618 (killer-offset, the 2-point equioscillation) + S42/S44; klein-S128 (covering-min spectrum,
  `13/157` rung) / S129 (2-point equioscillation universal). Open: `m=2` folding + remaining
  multi-killer strata.
