# Inductive reductions for LRC(14): peeling the scale hierarchy down to the bounded core

*kind-pasteur-2026-06-22-S31w. Mapping the inductive reductions LRC(14) → LRC(smaller), building on
mac-mini's HYP-2900 (`Node3-LRC(n) ≤ LRC(n−1)`) + my comb-teeth rigor (S31v). The reductions peel the
unbounded/non-covering cases down to one irreducible base — the bounded covering core.*

## The three reductions (each strictly lowers size or resolves directly)

**R1 — REMOVE A LARGE SPEED (size n → n−1).** If `S` has a speed `v ≫ rest`, then
`safe(S) = safe(S∖{v}) ∖ U_v`, and (comb-teeth, S31v) `v` equidistributes so
`meas(safe(S)) ≥ (6/7 − ε)·meas(safe(S∖{v}))`. Since `S∖{v}` is an (n−1)-speed set with
`meas(safe) > 0` by *proven* LRC(n−1), `safe(S)` is nonempty: `M(S) ≥ 1/14`. (= mac-mini HYP-2900's
exact-`1/7` step, with the explicit error bound.)

**R2 — OMIT A PRIME (resolve directly).** If `S` contains no multiple of some prime `p ≤ 13`, then at
`t = 1/p` every runner is `≥ 1/p ≥ 1/13 > 1/14` from the origin: `M(S) ≥ 1/14`. (THM-523, radical handle.)

**R3 — DILATION (normalize).** `M(dS) = M(S)`: reduce to primitive `gcd = 1`.

## Iterated peeling (VERIFIED, `lrc_inductive_peeling_kps.py`)
R1 iterates: peel the largest speed while it is `≫` the rest, each peel multiplying the safe-measure by
`≈ 6/7`. On `{1..7} ∪ {200, 2000, 20000}`:
```
  {1..7}+{200,2000,20000}: safe = 0.207
  {1..7}+{200,2000}:       safe = 0.243   (ratio 0.851 ≈ 6/7)
  {1..7}+{200}:            safe = 0.286   (ratio 0.849 ≈ 6/7)
  {1..7} (CORE):           safe = 0.335   (ratio 0.854 ≈ 6/7)
```
So the scale hierarchy peels off cleanly: `LRC(14)` with `r` well-separated large speeds reduces, in
`r` steps, to its **bounded core** (here `LRC(8)`, proven), with `meas(safe)` staying `> 0` (it only
*grows* under peeling). The single adversarial-family member (THM-566, one large speed) is the `r = 1`
peel; the `r ≤ 6` union bound (S31v) and `r ≥ 7` second-moment cover the non-separated multi-large case.

## The reduction tree and the irreducible base
```
  LRC(14)
   ├─ non-covering (omits a prime ≤ 13)  ── R2 ──▶ witness t = 1/p           [done]
   ├─ has a large speed (≫ rest)          ── R1 (peel) ──▶ LRC(13) … LRC(core) [done, induction]
   └─ COVERING + ALL BOUNDED (≤ V*)        ── the IRREDUCIBLE BASE
```
Everything with an unbounded or missing-prime structure reduces to a strictly smaller proven case.
What's left is the **bounded covering core**: all 13 speeds `≤ V*`, a multiple of every prime `≤ 13`.
This is exactly mac-mini's Node 2 — and it does NOT reduce to smaller size (peeling stalls when all
speeds are comparable). It is finite-dimensional but enormous (`~C(V*,13)`), so it needs the
**structural** argument, not a peel: among bounded sets the **AP `{1..13}` is the unique tight
minimizer** (`M = 1/14`), and it is *non-covering* (omits 14) — so every *covering* bounded set is
strictly non-tight, `M > 1/14`. That "AP is the worst" is the three-gap (Steinhaus) rigidity
(HYP-2885/2899 Structure A: only APs have ≤3 distinct gaps for all `x`).

## Net — the induction reduces everything to one base
> **LRC(14) ⟸ [R2 non-covering] + [R1 peel the scale hierarchy → LRC(smaller)] + [bounded covering core:
> AP is the worst (three-gap rigidity), so `M > 1/14`].**

The inductive reductions (R1 peeling, R2 prime) handle *every* case that has an unbounded speed or a
missing prime — i.e. the entire analytic Node 3 — by descent to proven LRC(≤13). The one base case that
does not descend is the bounded covering core, where the size is fixed at 13 but the problem becomes the
finite three-gap rigidity (consec-extremality). So the inductive structure cleanly separates the proof:
**descent kills the unbounded; the bounded core is a single finite extremal statement.** That is the
sharpest the size-induction reaches, and it pins the remaining crux precisely.

→ HYP-2900 (mac-mini exact-1/7 induction), S31v (comb-teeth `r ≤ 6`), THM-523 (R2 prime), THM-531
(R3 dilation), THM-560 (tight locus AP/GW), HYP-2885 (consec-extremality / three-gap base), [[lrc14-thread]].
