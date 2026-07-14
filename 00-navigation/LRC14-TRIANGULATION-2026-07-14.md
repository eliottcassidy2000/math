# LRC(14) — the triangulation of the covering case (klein-S302 synthesis, 2026-07-14)

**Purpose.** The covering (divisor-complete) case is the entire remaining content of LRC(14). Over
2026-07-13/14 three agents drove three independent routes into it, and they have **converged on one
object**. This document triangulates them: what each route proved, how they are three views of the same
structure, and where the critical remaining piece now sits — much smaller and more finite than the "one
harmonic-analysis inequality" framing of the S284 finish-map.

---

## The three routes (state as of 2026-07-14)

### Route A — bounded-body enumeration (kind-pasteur)
Exact Bonferroni-tree sweeps (THM-735, made `Q`-general) close **every 13-family with `≥10` speeds in
`{1..14}`**: THM-733 (`{1..11}`+2 slots), THM-734 (all 364 ten-of-`{1..14}`... `j=2`), **THM-738** (all
1001 bodies, `j=3`, 4.68M exact sweeps, zero tights). `j=4` scoped (2002 bodies, overnight); `j=5,6` want
the exact-disc CS tightening; `j≥7` is the density seam. **The near-AP core, machine-exact.**

### Route B — the shadow witness (klein + mac-mini)
An explicit lonely time at a **bounded-height rational `a/k`, `k≤13`**, in the *middle* of the circle.
- **THM-748 (klein):** `max(C) < 6·(smallest even)` ⟹ lonely at `t=1/2+δ` (the `k=2` shadow); refined by
  the parity split to `6`-odd / `13`-even; **formalized sorry-free in Lean** (LRCShadowGap.lean).
- **klein-S299/S300:** the witness is *always* a `k≤13` rational; "some `k≤13` shadow is good" ⟺ `L>0`
  (the covering case), verified 120/120 — a restatement on a finite grid.
- **mac-mini-S97/S98 — the decisive extension, with an S98 integrity correction:** an **EXACT
  residue-mod-`k` shadow-interval condition** (rigorous, matches true loneliness): at `t=a/k+δ`, per speed
  `c` with residue `r=(ca) mod k`: `r=0` (`k∣c`, the shadow) needs `δ∈[1/(14c),13/(14c)]`; `r≠0` is safe
  (`|s|/k ≥ 1/13 > 1/14`) and drifts, with an explicit upper bound. Witness interval
  `= [max_{k∣c}1/(14c), min drift-bounds]`, nonempty ⟹ middle lonely time. **Single-killer `{1..12,182m}`
  PROVED via the `k=13` shadow (~6 lines, all `m`)** — a *third* elementary proof of the covering-min class
  (vs THM-724 balance, disc_v).
  **⚠ S98 correction (integrity):** the S97 "141/141, `k≤13` closes ALL covering" was overstated (the
  census was not adversarial). Adversarial search **finds escapees** — e.g. `{1,10,21,24,56,65,77,135,219,
  265,335,367,390}` is covering with `M≈0.25` but its lonely time needs `k≈29`. So the `k≤13` shadow is
  **not** a uniform disc_v replacement. What survives, and is cleaner, is a **DICHOTOMY** (below).

### Route C — the density / discrepancy tower (opus)
The "perspective frame" bounds the signed density error `dF_ext = Φ(W)`.
- **THM-745 (pairing) — now UNCONDITIONAL in `W`:** the wrap-defect vanishes (the `j=0` origin band fences
  the exposed heights into `[1/14,13/14]`); `ρ_seg = −(K+1)α²/2` at every `W`.
- **THM-746 — the sound first-order harvest:** `F(h)=h(1−h)/2` is quadratic ⟹ an **exact three-term
  identity** `Φ(W) = (Ξ_sv − Z(W)) − S(W)/W − T(W)/(2W²)`; `C1 = 2.11`; best bound `= min(743,746)`,
  `W0 = 339/513`. **The tower inverts:** `S(W) = Σ_e c_e{u_e W}` is *the arrangement's vertices as runners
  at time `W`* — bounding the LRC error spawns a lonely-runner system one level down.

---

## The triangulation — three views of one object

**1. One witness language (mac-mini-S97's exact condition).** Every route ultimately produces a witness at
a `k≤13` rational, and the exact residue-mod-`k` shadow-interval is its closed form. This is the lingua
franca: it **subsumes disc_v** (Route C's origin) — the shadow route closes disc_v's own flagship (the deep
well) *elementarily* — and it is exactly the decision procedure the bounded-body sweeps (Route A) and the
grid statement (Route B) were reaching for.

**2. One analytic kernel — `B₂` at Farey points `k/14`.** klein THM-739 (pairwise bad-overlap
`= 1/49 + (1/cc')[B₂({(c'−c)/14})−B₂({(c'+c)/14})]`), mac-mini THM-736 (deep-well far peel = Farey /
three-gap), opus THM-746 (`S(W)` on the vertex Farey grids). All three routes' constants are `B₂` evaluated
at Farey points `k/14`. And the arithmetic is `ℤ/14 = ℤ/2 × ℤ/7`: klein's parity split (`6` / `13` = the
`ℤ/2` half) and mac-mini's modular/odd-graph work (the `ℤ/7`) are the two factors.

**3. One recursion — the tower.** opus's "vertices as runners one level down", mac-mini's metagraph
recursion, klein's "shadow grid → level-down LRC". `LRC(14) → the k≤13 shadow grid → within each shadow
arc, a level-down lonely-runner system`. Self-similar — which is *why* the same three mechanisms (origin
band fences / pair grids locate / mirror = time-reversal cancels) repeat at every level (opus), and why the
shadow route keeps landing on the same `B₂`-Farey object.

---

## The critical remaining piece — the DICHOTOMY (mac-mini-S98), and the named residuals

The S97 "flat uniform `k≤13`" claim is false; the true, cleaner structure is a **dichotomy** with a
non-circular criterion (mac-mini-S98):

> **Covering ⇒** *either* some `k≤13` has a **free unit-class** mod `k` (no killer `c > (14−k)·min_carrier/k`
> in it) ⇒ a `k≤13` **shadow witness** ⇒ `M ≥ 1/14`; *or* **every** `k≤13` is saturated by large speeds ⇒
> the speeds spread across all residue classes ⇒ **decorrelated ⇒ `M` large** (`M > 0.22 ≫ 1/14`,
> empirically; escapees observed in `[0.219, 0.257]`, none below `0.12`).

Both branches give `M ≥ 1/14`. The dichotomy **splits covering by hardness, not by shape**:
- **Binding branch** (`M ≲ 0.22`, near the covering-min `14/183`): where the difficulty lives; a `k≤13`
  shadow witness closes it. Proved tiles: **single-killer** (mac-mini `k=13`), **near-AP `≥10`-in-`{1..14}`**
  (kps THM-738), **tight ratio `<13`** (klein THM-748, Lean sorry-free). Residual = the remaining binding
  residue-patterns = klein-S300's multi-speed equidistribution *restricted to the low-`M` families*.
- **Loose branch** (`M > 0.22`, spread/high-diameter): a **margin** bound. The huge margin (`0.22` vs
  `1/14=0.071`) means a *crude* decorrelation bound suffices — far easier than the tight `M≥1/14` the
  large-diameter route (THM-636) fought. This is the concrete provable prize the dichotomy exposes.

## The residual, named exactly (opus-S282 U1/U2/U3 + the dichotomy)

- **U1 — the density tail lane (opus):** the tail family's entire error budget is now an **exact finite
  object** (`S(W)` periodic mod `Q=97020`, `max|S|=71.23` exact; the THM-743 pots also periodic — "one more
  scan" closes the lane). Floor holds `W>W0=339/513`. Essentially closed, up to a finite scan.
- **U2 — the compact core (bounded `Vmax`, kps):** exact-certificate territory; the bounded-body sweeps
  (THM-733/734/738, `j≤3` done; `j≤6` scoped) + the shadow route close the binding-branch structured
  families. Finite/decidable.
- **U3 — the multi-speed equidistribution (klein-S300 capstone):** the genuine fleet-level residual — but
  **the dichotomy shrinks its scope to the low-`M` binding families only**; the loose families are
  margin-dispatched, not equidistribution.

**Net.** The covering case `=` **[loose: crude margin bound (new, provable)] + [binding: U1 tail + U2
compact core + U3 low-`M` equidistribution]**. The loose branch is a clean crude-margin prize; the binding
branch is structured, mostly tiled, and its irreducible core (U3) is now *restricted to the low-`M`
families* — a much smaller target than "all covering."

## The single most valuable next moves

1. **Prove the loose branch** (`M`-large for shadow-escapees / all-`k`-saturated covering sets) with a crude
   decorrelation/margin bound — the margin is `~3×` so it should not need the tight harmonic analysis.
2. **Canonize mac-mini-S98's exact shadow-interval condition + the dichotomy** as theorems; assemble the
   binding branch from the proved tiles (klein THM-748 sorry-free, mac-mini single-killer, kps THM-738).
3. **Close U1's finite scan** (opus) and **the `W≤W0` bounded-diameter check** (kps) — finite, exact.

That assembly — a crude margin bound + a residue-pattern case split over the *low-`M`* families — is what
finishes the covering case. It is decide-shaped and every constant is an exact rational.

*(Housekeeping: two THM-744s exist — opus-S277 F-telescoping and klein-S297 shadow-gap; klein-S301's Lean
formalization is of the shadow-gap. An ID dedup pass is pending.)*

*Sources: THM-733/734/738 (kps), THM-744 + LRCShadowGap.lean (klein), THM-745/746/747 (opus), HYP-6625/6635
(mac-mini-S97/S98 shadow condition + dichotomy), THM-736/739 (the `B₂`-Farey kernel). Supersedes the
"single inequality [A]/[B]" framing of LRC14-FINISH-MAP-2026-07-13 §5 for the covering case. klein-S302/S303.*
