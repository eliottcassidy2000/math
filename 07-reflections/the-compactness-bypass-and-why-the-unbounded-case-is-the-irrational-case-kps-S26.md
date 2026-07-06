# The compactness bypass, why it fails cheaply, and how the unbounded case IS the irrational (three-gap) case

*kind-pasteur-2026-07-06-S26 — an honest attempt to bypass the unbounded-height
half of (G) by compactness, the precise reason it fails, and the positive bridge
it leaves: unbounded height ↔ irrational limit directions ↔ mac-mini's three-gap
quantisation (HYP-4412).*

## The setup: (G) lives on a compact space

`M(S) = max_t min_i ‖v_i t‖` is **scaling-invariant** (`M(cv) = M(v)`), so it is a
function on the **compact** projective space of directions `P¹¹(ℝ)` (a family `v`
is a rational point; its height is the size of the primitive representative). And
`M` is **continuous** there (`max` over compact `t` of a jointly-continuous
function). So the gap set `G = {d : M(d) ∈ (1/13, 2/25)}` is **open**.

The tempting bypass: an open set contains a rational point of *bounded* height
(Dirichlet), so if any gap direction existed — even an irrational, unbounded-
height limit — a **finite** gap member would exist near it, and a finite check
settles (G). The unbounded case would be bypassed by compactness.

## Why it fails cheaply: `M` is not equicontinuous

The bypass needs a *uniform* Lipschitz bound, and there isn't one. Under
normalisation `v ↦ v/|v|`, the loneliness parameter rescales: the optimal `t` for
the family becomes `s = |v|·t`, which ranges up to `|v| ~ height`. So

> `M` is Lipschitz on directions with constant `L(d) ~ height(d)/13`, **not**
> uniform — it **grows with height**.

Measured (`lrc_compactness_bypass_kps_S26`): the effective constant over
height-≤60 families is already `≈ 38`, not the naive `≤ 1/2`. Feeding it through
Dirichlet, a gap member whose `M` is `δ` from the nearest edge is only guaranteed
at height `~ (1/δ)¹¹`; for the deepest interior point (the mediant `3/38`,
`δ = 1/950`) that is `~10²⁰`. A finite bound in principle, astronomical in
practice, and no better than mac-mini's `q ≤ 2·max` lever. **The compactness
bypass does not give a usable finite check.**

And the reason it fails is exactly the crux: `M` **oscillates at frequency
`~height`** near the tight locus. That non-equicontinuity IS the density floor —
the same obstruction that makes MISTAKE-110 (unbounded witness modulus) and the
razor-thin no-low-order-Bonferroni cancellation (mac-mini HYP-4402). A caution
for the density-floor programme (mac-mini HYP-4452): any "height threshold" above
which a bounded check suffices is *non-uniform* — it must scale with the depth
into the gap, because `M`'s resolution near a direction is `1/height`, not a
constant.

## What the failure leaves — the positive bridge

Compactness still buys one real thing. If finite gap members existed at unbounded
height, their directions (compactness) have a convergent subsequence to a limit
direction `d*` with `M(d*) ∈ [1/13, 2/25]`. A generic accumulation `d*` is
**irrational** — a `{kα}`-type direction. So:

> the **unbounded-height** case of (G) IS the **irrational-direction** (loose-
> branch) case, reached as the limit of the rational gap members.

And the irrational directions are exactly where the **three-gap theorem**
(mac-mini HYP-4412) is *clean*: an irrational `d*` with small loneliness has a
witness whose `{v_i t*}` is a genuine `{kα}` orbit, with `≤ 3` gap lengths (Sós),
forcing `M(d*)` onto a continued-fraction rung — and the CF rungs are my
Stern–Brocot denominators `13a+25b` (S25), which **skip** the open interval. So
the two ends match: the *rational* gap values are Stern–Brocot descendants of
`(1/13, 2/25)` (S25); their *irrational accumulation points* are the CF limits,
governed by three-gap. The unbounded-rational case and the irrational-limit case
are **one object seen from the two sides of the rationals**, and mac-mini's
three-gap quantisation is the tool that applies on the irrational side.

This is the bypass, correctly placed: you do not check unbounded height directly;
you push it to the compact space's **irrational limit points**, where three-gap /
CF rigidity replaces the (hopeless) height induction. The residual is then the
**local isolation of two compact orbits** — the AP orbit (`M = 1/13`) and the
`{1,…,11,24}` orbit (`M = 2/25`) — each a `THM-615-Lemma-3`-flavoured
Lipschitz/derivative statement *at fixed compact directions*, where the constant
`L(d)` is finite (no longer competing with unbounded height). The density floor
is precisely that local isolation.

## The honest ledger

- **Bypass by compactness alone: NO** — `M` is continuous but not equicontinuous
  (`L ~ height`), so the finite bound is astronomical and useless.
- **What survives: a bridge** — unbounded height ⇒ irrational limit direction ⇒
  three-gap/CF regime (mac-mini HYP-4412), where the Stern–Brocot denominators
  (S25) skip the gap. The unbounded case is *transported* to the irrational side,
  not solved by height induction.
- **Residual: the two edge-orbit isolations** (AP at `1/13`, `{1..11,24}` at
  `2/25`), local and compact — the density floor, now with a *finite* local
  Lipschitz constant. This is where the Riesz/Selberg analysis (HYP-4452) does its
  work, and where it is genuinely needed.

## Pointers

- `lrc_compactness_bypass_kps_S26.py` / `.out` — the Lipschitz-constant estimate
  (`L ≈ 38`, growing with height) and the Dirichlet height calibration.
- kps S25 (Stern–Brocot denominators `13a+25b`, the rational gap values),
  HYP-4417 (residue split), S24 `gap_candidate_has_multiple`.
- mac-mini HYP-4412 (three-gap/CF quantisation — the irrational-side tool),
  HYP-4452 (Riesz-product density floor), HYP-4432 (the `q ≤ 2·max` lever).
- THM-615 Lemma 3 (the local Lipschitz+IVT non-extremity, = the edge-isolation
  shape).
