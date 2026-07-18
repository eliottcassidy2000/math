# THM-1000 — The far-element covering bound: a tight family has no runaway top speed (death-star-2026-07-17-S56)

**Status:**
- **Lemma B (PROVED):** for a tight family `V` and any `w ∈ V`, every connected component of the good
  set of `V∖{w}` (level `1/n`) has length `≤ 2/(nw)`.
- **Theorem T (PROVED):** every tight family `V` (`M(V)=1/n`, `|V|=n−1 ≥ 2`) satisfies
  `Vmax ≤ (n−1)·v₂`, where `v₂` is the second-largest speed. Uses only LRC(n−1) (settled for n≤14) +
  a Lipschitz/covering argument. Excludes the deep-well far family from tightness.
- **Consequence for the S56 residual:** this is genuine progress on "primitive covering tight ⟹
  small/comparable" — it kills the *runaway far element*. The **absolute** bound on `Vmax` (the full
  spread) is not closed; it reduces to bounding `v₂` recursively, i.e. to the tight-locus classification.

Source HYP-7305; extends THM-999. Script `04-computation/lrc_vmax_ratio_bound_deathstar_S56.py` (+`.out`).
Relates to THM-721 (u-escape), THM-733 (`{1..11,a,b}` tight ⟹ AP/GW), THM-755 (capped envelope),
THM-758 (far-count) — this is the clean general form, no case analysis. WLOG positive speeds.

---

## Lemma B — the covering width bound

**Statement.** Let `V` be tight (`M(V)=1/n`) and `w ∈ V`, `V' = V∖{w}`. Let
`G' = {t : min_{v∈V'} ‖vt‖ > 1/n}` be the good set of `V'`. Then every connected component of `G'`
has length `≤ 2/(nw)`.

**Proof.** The safe set of `V` is `{t : min_{v∈V}‖vt‖ > 1/n} = G' ∩ {‖wt‖ > 1/n}`. Tightness makes it
measure-zero, so `G' ⊆ {‖wt‖ ≤ 1/n}` up to measure zero. Now
`{‖wt‖ ≤ 1/n} = ⋃_k [k/w − 1/(nw), k/w + 1/(nw)]` — arcs of width `2/(nw)` centered at `k/w`, separated
by **gaps** of width `1/w − 2/(nw) = (n−2)/(nw) > 0` (for `n > 2`). A connected component `I` of `G'`
(an open interval) lies inside `⋃(arcs)` up to measure zero; if `I` met two different arcs it would
contain one of the positive-measure gaps between them, contradicting the inclusion. So `I` lies in a
single arc, giving `|I| ≤ 2/(nw)`. ∎

Verified: for AP the 12 good-components of `{1..12}` are all `≤ 2/(14·13)=1/91`; for GW the 4
good-components of `{1..11,13}` are all `≤ 2/(14·24)=1/168`. For the **non-tight** deep well
`{1..12,182}`, the covering **fails** — `{1..12}` has a good-component of width `1/168 > 2/(14·182) =
1/1274`, i.e. 182's arcs are too narrow to cover it — which is exactly why the deep well has a
positive-measure safe set and `M = 14/183 > 1/14`.

---

## Theorem T — the top-ratio bound `Vmax ≤ (n−1)·v₂`

**Statement.** For a tight family `V` with `n−1 ≥ 2` speeds, ordering `v₁ ≤ ⋯ ≤ v_{n−1} = Vmax`,
```
Vmax ≤ (n−1)·v_{n−2}.
```

**Proof.** Put `w = Vmax`, `V' = V∖{w} = {v₁,…,v_{n−2}}` (`n−2` speeds).
1. **Lower bound on a good-component of `V'`.** By LRC(n−1) (settled for `n ≤ 14`), `M(V') ≥ 1/(n−1)`.
   Let `t₀` attain it: `min_{v∈V'}‖vt₀‖ = M(V') ≥ 1/(n−1)`. The map `t ↦ min_{v∈V'}‖vt‖` is Lipschitz
   with constant `max(V') = v_{n−2}` (each `‖vt‖` has slope `≤ v`). Hence it stays `> 1/n` while
   `|t − t₀| < (1/(n−1) − 1/n)/v_{n−2} = 1/(n(n−1)v_{n−2})`. So the component of `G'` containing `t₀`
   has length `≥ 2/(n(n−1)v_{n−2})`.
2. **Upper bound from Lemma B.** That same component has length `≤ 2/(nw)`.
3. Combining: `2/(n(n−1)v_{n−2}) ≤ 2/(nw)`, i.e. `w ≤ (n−1)v_{n−2}`. ∎

**Verified:** AP `13 ≤ 13·12`, GW `24 ≤ 13·13`; the deep well `182 > 13·12 = 156` **violates** the
bound, re-deriving that it is not tight. The bound is loose (true top-ratios are `13/12` and `24/13`)
but general and exact in form; sharpening the constant needs a better good-component lower bound.

---

## What this settles, and the residual

**Settled:** a tight family has **no runaway far element** — the largest speed is within a factor
`n−1` of the second-largest. This is the clean general statement behind THM-733/755/758: it rules out
the deep-well/ladder family `{1..12,182m}` (and any single dominant outlier) from being *exactly* tight,
using only settled LRC(n−1) and elementary covering — no enumeration.

**Not settled (the residual = the classification).** Theorem T bounds `Vmax/v₂`, not `Vmax` absolutely.
Iterating (bound `v₂` by `v₃`, …) would need the peeled sub-family `V'` to be itself *tight* at level
`1/(n−1)`; this holds for the AP (`{1..12}` is LRC(13)-tight, so the AP is bounded by iterated
peeling) but **not** in general (GW's `{1..11,13}` is not tight), and the two-speed covering that would
bound `v_{n−2}` without sub-tightness only gives a vacuous inequality. So an **absolute** `Vmax` bound
over primitive tight families is equivalent to the full spread being bounded, i.e. to the tight-locus
classification `{AP, GW}` — the LRC-hard core (THM-995 VII).

**Net.** "Bound `Vmax` over primitive tight families" now splits cleanly: the *top ratio* is bounded
unconditionally (Theorem T, killing far elements); the *residual* is the bounded-spread/classification
core. Combined with THM-999 (R is a finite check once `Vmax` is bounded), a full absolute bound would
close R uniformly — so the three theorems THM-997/999/1000 reduce the S56 census residual to exactly
the one hard object the whole program is stuck on, with the far-element half now removed.

→ THM-995 (VII), THM-996, THM-997, THM-999, THM-721/733/755/758, HYP-7305.
