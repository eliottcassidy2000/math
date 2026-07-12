# Divisor-completeness is one integer — so the coprime-core size is a diameter artifact

*kind-pasteur-2026-07-11-S127 cont.49. Owner: "search for more combinatorial connections and simplifications."
I searched, found the fleet (mac-mini cont.50, opus-S244) had just reached the same place, and can sharpen it
from the combinatorial-cover side with one crisp fact: the divisor-complete requirement `{2,…,14}` is met by a
single integer.*

---

## The one-integer cover

`lcm(2,…,14) = 360360`. A single runner equal to `360360` is divisible by every `d ∈ {2,…,14}`, so
**one integer satisfies the entire divisor-complete requirement.** (A tidy two-integer cover: `2520 = 2³·3²·5·7`
covers `{2,…,10,12,14}`, `143 = 11·13` covers `{11,13}`.) This is sharper than the fleet's "a multiple of 840
covers all even conditions" (mac-mini cont.50): `840` misses `9, 11, 13`; `360360` misses nothing.

The consequence is immediate. A divisor-complete family needs only **one** non-coprime-to-`30030` runner
(`30030 = 2·3·5·7·11·13`); the other twelve are free to be coprime to `30030`. Concretely,

`A = {360360} ∪ {17,19,23,29,31,37,41,43,47,53,59,61}`

is a **primitive, spread (longest-run 1), divisor-complete** 13-family with **12** coprime-to-`30030` runners,
`M ≈ 0.347`. And `B = {2520,143} ∪ {17,…,59}` has 11, `M ≈ 0.26`. Both far above `1/14`.

## What this confirms and sharpens

mac-mini cont.50 (MISTAKE-139) had already corrected the "≤6 distinct lifts" mechanism and, converging with
opus-S243, noted that the coprime-to-`30030` core is `≤6` only for **bounded-diameter** DC (mean 2.0), and
that adversarial families reach `12`. This session confirms both from a clean combinatorial angle and pins the
*reason*: the `≤6` is not a property of divisor-completeness — it is **availability**. Coprime-to-`30030`
integers have density `∏_{p≤13}(1−1/p) ≈ 0.192`, so only a handful exist below a small `Vmax` (exactly
`{1,17,19,23}` below `25`, hence `max = 4` there; `5,6,7` at `Vmax ≤ 40,80,200`). Since DC costs only one
runner, the core size is capped only by how many coprime values fit under the diameter — it grows with `Vmax`
to the structural ceiling `12`. The `≤6` is a bounded-diameter/availability artifact, the same shape as the
band-clearing window being bounded only under bounded diameter (cont.47). Two "small at bounded diameter,
unbounded in general" invariants, one cause.

## The simplification worth keeping

The universal, diameter-free statement is not `≤6` but:

> **The coprime-to-`30030` core is always `≤ 12`** — at least one runner is non-coprime (it carries the DC
> requirement) — so the core is a `≤12`-speed family, and **`LRC(≤13)` gives it `reach ≥ 1/13 > 1/14`
> regardless of its size.**

This reframes opus-S244's reduction (`LRC(14) ⟺` the coprime core does not blanket the good set `G'`): the
core's *size* is not the lever — it is bounded by availability and always `LRC(13)`-protected as a standalone
family. What can obstruct is only the **alignment** of the core's danger arcs with `G'`, not the number of
core speeds. So the honest crux is "the `LRC(13)`-good set of the core meets `G'`," an intersection/alignment
statement, and the "≤6" is a red herring for it — the core is loose on its own at any size; the question is
joint. (For the large-diameter families like `A`, the intersection is wide — `M ≈ 0.35` — so they are trivially
lonely; the obstruction lives only in the bounded-diameter near-tight neighborhood, as everywhere else in this
endgame.)

## Honest scope

This is a confirmation-plus-sharpening, not a new result. mac-mini cont.50 owns the correction; opus-S244 owns
the fold-to-core synthesis. What this adds is the one-integer cover (crisp, and cleaner than the 840 version),
the explicit `12`-coprime witnesses, the availability-not-cover-cost explanation of the `≤6`, and the
`≤12`/`LRC(13)`-protected reframe that removes the size bound from the crux. The lesson recurs: a "small
count" that grows with a free parameter (here the diameter) is not a structural bound — check the generic and
the adversarial, name the real cause (availability), and keep only the diameter-free invariant (`≤12`).

*Files: lrc14_dc_is_cheap_coprime_core_kps_S127.py/out. HYP-6140. Confirms/sharpens mac-mini cont.50
(MISTAKE-139) + opus-S243/S244 (coprime core) from the combinatorial-cover side; parallels
[[the-clearing-window-is-unbounded-the-finite-check-is-only-bounded-diameter-kps-S127]] (both invariants are
bounded-diameter artifacts).*
