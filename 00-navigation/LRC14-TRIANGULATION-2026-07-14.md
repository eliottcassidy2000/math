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
- **THM-744 (klein):** `max(C) < 6·(smallest even)` ⟹ lonely at `t=1/2+δ` (the `k=2` shadow); refined by
  the parity split to `6`-odd / `13`-even; **formalized sorry-free in Lean** (LRCShadowGap.lean).
- **klein-S299/S300:** the witness is *always* a `k≤13` rational; "some `k≤13` shadow is good" ⟺ `L>0`
  (the covering case), verified 120/120 — a restatement on a finite grid.
- **mac-mini-S97 — the decisive extension:** an **EXACT residue-mod-`k` shadow-interval condition**
  (rigorous, matches true loneliness): at `t=a/k+δ`, per speed `c` with residue `r=(ca) mod k`: `r=0`
  (`k∣c`, the shadow) needs `δ∈[1/(14c),13/(14c)]`; `r≠0` is safe (`|s|/k ≥ 1/13 > 1/14`) and drifts, with
  an explicit upper bound. Witness interval `= [max_{k∣c}1/(14c), min drift-bounds]`, nonempty ⟹ middle
  lonely time. **141/141 covering families closed (incl. the isolated-far deep well)**, `min-k` dist
  `{2:52,3:25,4:18,5:13,6:7,7:8,8:2,9:2,11:1,13:13}`. **Single-killer `{1..12,182m}` PROVED via the `k=13`
  shadow (~6 lines)** — a *third* elementary proof of the covering-min class (vs THM-724 balance, disc_v).

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

## The critical remaining piece — triangulated, and now finite

The single open statement is the **UNIFORM SHADOW CLOSURE**:

> *For every covering 13-set, some `k≤13` has a non-empty shadow interval* (mac-mini's exact condition)
> `⟺` a middle lonely time `⟺` `L>0` (klein-S300 equivalence).

It is now **decidable at the residue level** and factors along two axes, one route per factor:

| factor | content | owned by | status |
|---|---|---|---|
| **residue pattern** (mod `lcm(2..14)`) | some `k≤13` has a `k∣c` shadow speed + no drift-down killer | mac-mini exact condition + kps enumeration | DECIDABLE; near-AP `j≤3` PROVED (kps) |
| **ratio, tight** (`<13`) | shadow interval non-empty | klein THM-744 (`6`/`13` parity) | PROVED |
| **ratio, isolated-far** (single-killer) | `k=13` shadow non-empty | mac-mini | PROVED |
| **ratio, spread mid-band** | interval non-empty for `W ≤ W0`; floor for `W > W0` | opus THM-745/746 (`W0=339/513`) | floor PROVED; tail = finite check |

**The synthesis.** The covering case `=` uniform shadow closure `=` **[finite residue-pattern check] ×
[ratio control]**, and the four proved tiles above cover the `(pattern × ratio)` space except the **spread
mid-band multi-killer** stratum. There, opus's density floor caps the analytic content at `W0 ≈ 339–513`:
**above `W0` the floor closes it; below `W0` it is a bounded-`W` (bounded-diameter) finite check** — which
Route A's sweep program and Route B's shadow census already handle empirically (141/141). So the covering
case is no longer "one harmonic-analysis inequality": it is a **finite problem** — a decidable
residue-pattern statement with a bounded-`W` tail — and the three routes together already cover it except
for canonizing the pieces and discharging the finite mid-band residue.

## The single most valuable next move

**Canonize the uniform shadow-interval condition (mac-mini-S97) as a theorem, then prove the shadow closure
by residue-pattern cases**, importing the four proved tiles: klein THM-744 (tight), mac-mini single-killer
(isolated-far), kps THM-738 (near-AP), opus THM-745/746 (spread, `W>W0`) + a bounded-`W` finite check
(`W≤W0`). That assembly — not another analytic inequality — is what finishes the covering case. It is
decide-shaped, Lean-tractable (klein's THM-744 tile is already sorry-free), and every constant is an exact
rational.

*Sources: THM-733/734/738 (kps), THM-744 + LRCShadowGap.lean + HYP-6620/6630 (klein), THM-745/746 (opus),
HYP-6625 (mac-mini-S97), THM-736/739 (the `B₂`-Farey kernel). Supersedes the "single inequality [A]/[B]"
framing of LRC14-FINISH-MAP-2026-07-13 §5 for the covering case. klein-S302.*
