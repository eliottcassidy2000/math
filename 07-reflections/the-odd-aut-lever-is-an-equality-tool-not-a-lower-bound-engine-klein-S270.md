# The odd-|Aut| tie-resolution lever is an equality tool, not a lower-bound engine — and where the lower bound actually lives

*klein-2026-07-12-S270. Owner: explore the odd-|Aut| tie-resolution lever (S269) as a proof route.
I did, and the honest answer is a demotion: the lever proves why the tight AP is pinned at exactly
1/14 (the equality case), but it is **not** a lower-bound engine for LRC(14), and I can show the exact
place its "you can't beat 1/14" claim fails. The genuine lower-bound content lives in two other
objects — the lander-exclusion count and the global antipodal degree — both open.*

---

## The lever, stated exactly

At `t = 1/14` the observer (speed 0) and 13 runners form the 14-grid `{0,1,…,13}/14`, which carries
the order-2 antipodal symmetry `x ↦ x+1/2` (`k ↦ k+7`). A tournament has **odd** `|Aut|` — two arc
states, no order-2 automorphism — so the winding tournament can't be invariant under `x↦x+1/2`; the 7
antipodal pairs at difference 7 are **tied**, and every one of the `2⁷` tie-resolutions leaves a
runner at distance exactly `1/14`. mac-mini cont.56's twin says the same thing arithmetically: the AP
is pinned at `1/14` because it **contains its own distance-1 landers** `{1,13} mod 14` (`v ≡ ±a⁻¹`).

## What it actually proves — verified

I tested the reach exactly (`lrc14_oddaut_lever_reach_klein_S270.py`). The decisive check is whether
"base 14 caps at `1/14`" holds for *all* families:

| family | best margin at base 14 | |
|---|---|---|
| AP `{1..13}` (full grid) | **`1/14`** | pinned — landers `±1` present |
| `{2,3,…,12,16,17,18}` (avoids `±1 mod 14`) | **`1/7`** | **beats `1/14`** at base 14 |
| a covering family (has a mult of 14) | **`0`** | blocked — a runner at residue 0 |

So "you can't beat `1/14` at base 14" is **false** the moment the residues leave the full grid. The
lever pins exactly `1/14` **only for the full-grid family** — i.e. it is a statement about *why the AP
achieves the tight value*, an **equality/uniqueness** fact. It is not a universal lower bound: a
family that dodges `±1 mod 14` does *better* at base 14 (it's non-covering, dispatched anyway), and a
covering family does *worse* (blocked, forced to a coarser base). The order-only winding tournament
forgets the metric (mac-mini-S57), and `H` is a class *selector*, not a loneliness scalar — the lever
inherits exactly that limitation.

## Why it can't touch the actual crux

The remaining conjecture is the covering case, `covering ⟹ M ≥ 14/183`, and it is realized at base
`183` (odd): the deep well's `M = 14/183` appears **only** at base 183 — its best over all bases
`q < 183` is `13/170 < 14/183` (verified). At the odd base 183 there are **no antipodal ties** (no
`q/2`), so the winding tournament is a genuine non-degenerate circulant and the odd-`|Aut|`
tie-resolution has *nothing to resolve*. The lever is a feature of the **even** base 14; the crux
lives at an **odd** base where it is silent. This is structural, not incidental: the tie-resolution is
the tournament shadow of the `2` in `14 = 2·7`, and the covering crux is the `7`-adic / far-base part.

## Where the lower bound actually lives (two open objects)

The exploration is not empty — it *locates* the lower-bound content precisely, in two places the
lever is a shadow of:

1. **The lander-exclusion count (mac-mini cont.56).** `M > 1/14` ⟺ the family, at its best base,
   contains no runner within `q/14` of 0 — in particular omits the distance-1 lander `v ≡ ±a⁻¹`. The
   covering case is then: *every covering family, at every base `q < 183`, contains a runner within
   `q/14` of 0* (so it can only escape at `q ≥ 183 = Φ₆(14)`). This is a finite, arithmetic
   lander-exclusion statement — the rotation-orbit form of the non-compressed peel — and it is the
   **open crux**. The odd-`|Aut|` lever supplies none of it (it's a base-14 equality fact).

2. **The antipodal degree — and why it is the *wrong category* for the floor (PROVED).** The deeper
   form of the parity: the involution `ι` (`a↦−a`, `t↦1−t`, `T↦T^op` = complement, THM-584) makes
   danger zones `ι`-symmetric (`‖−x‖=‖x‖`), so on odd moduli loneliness comes in `ι`-pairs — the
   witness count is **even** (verified: the deep well's good multipliers at base 183 are the pair
   `{14,169}`, `ι`-closed), the same mod-2 fact as `|Aut(T)|` odd / `H` odd. It is tempting (klein-S55)
   to hope a **global `ι`-equivariant degree** (antipodal Čech class / metagraph "odd index" /
   Borsuk–Ulam degree) vanishes ⟺ LRC(14). But the repo has **proved this degree is witness-side, not
   floor-side**: THM-582 — a lonely tournament has vertex 0 as a **source**, which the complement `R`
   turns into a **sink**, so it is **not self-converse**, and the odd/palindromic index does not apply
   to it; THM-581 — the floor (existence) is **σ-EVEN**, carries no sign-isotypic content, and does
   **not need** the Borsuk–Ulam certificate (the `ι`-fixed points `t=0,1/2` are danger for every
   covering set, no `p mod 4` split). So the antipodal degree is a genuine invariant of the
   **witness/Rédei half**, and the parity is a *constraint, not a detector* (klein-S55: it gives the
   *shape* of every hole, not that there are none — `drop-11 @ D=41` was lonely at radius 2, a perfect
   `ι`-pair). The topological-degree route to the *inequality* is not merely open — it is the wrong
   category. (And the naive apex-7 homological grade is refuted outright: `b₁⁻(7) = 1772`, `7 ∤ 1772`.)

**Where the floor's lower bound actually lives.** Not the degree but its **σ-even arithmetic
surrogate**: the 2-adic parity descent `meas(lonely S) = ∏ρ_j · ∏ meas(lonely O_j)` (THM-580) —
"compute the odd degree arithmetically instead of certifying it topologically" — leaving the analytic
**covering-min `inf_{covering} M ≥ 14/183`** (HYP-2566, open), an even-category discrepancy/Ostrowski
problem, exactly the lander-exclusion count of (1). Both open pieces are `σ`-even and arithmetic; the
antipodal degree sits on the other (witness) side of the `σ`-grading.

## Verdict

**The odd-`|Aut|` tie-resolution lever is an equality/uniqueness tool: it proves the full-grid AP is
pinned at exactly `1/14` (and, with the two mod-14 tie-resolutions, that the tight locus is
`{AP, GW}`). It is not a lower-bound engine — "you can't beat `1/14`" holds only on the full grid, the
covering crux lives at an odd base where the lever is silent, and its topological deepening (the
antipodal degree) is *proved* to be the wrong category for the floor (THM-581/582: the floor is
σ-even; the degree is witness-side).** The lower bound is the σ-even arithmetic content — the
lander-exclusion count / analytic covering-min `≥ 14/183` (HYP-2566, open) — of which the parity is
only the mod-2 shadow. Honest use of the lever: a clean proof of the *equality* half of the dichotomy
(why `1/14`, why `{AP,GW}`) — not, as I over-suggested in S269, "the load-bearing statement for why
you can't beat 1/14." A σ-grading verdict: **equality/witness = σ-odd (tie-resolution, degree);
inequality/floor = σ-even (descent, discrepancy). The lever lives entirely on the σ-odd side.**

## Two adjacent hopes, closed

- **"Kill the pressure-DAG counterexample-object with the symmetry lever" — no.** The pressure /
  two-clock "counterexample-shaped object" (a strict pressure SCC) is a *dormant reformulation*
  (codex, S470–S514, absent from the live proof map). Its SCC **never materializes** — every audit
  returns a DAG (4284/4284), and the realization clause it needs (THM-380) empirically fails
  (S504: the deletion-relief gauge carries only 2854/4284 owner-protection edges). What a covering
  counterexample *does* force is the weaker *owner-protection* cycle (THM-379, an elementary
  indegree-≥1 pigeonhole), not a pressure SCC. And the parity lever has **never** been attached to it,
  for a good reason: the pressure SCC is sampled at endpoint/gap times, so it **forgets M** (the same
  failure mode) — the lever belongs on the winding tournament *at the optimum*, which carries M. So
  the tempting "the counterexample-object must carry a symmetry a tournament can't hold" is not a live
  route; it would bolt the lever onto the wrong tournament.

- **The lever's real bottleneck (the actionable takeaway): it reduces to *consec-maximizes*.** On the
  winding tournament where it does live, "every tie-resolution gives M ≥ 1/14 / covering ⟹ M ≥ 14/183"
  is *logically equivalent* to "the extremal optimum tournament is the regular circulant ⟺ the
  residues are a `{kα}` three-gap orbit ⟺ consec/AP maximizes the coverage functional" — the
  three-gap / consec-maximizes-p0 crux the proof map is already grinding on (THM-716/717, the POS/BUNCH
  leg). The honest use of the lever is therefore not a new engine but a *pointer*: it says the equality
  rigidity and the inequality crux are the **same** three-gap statement, and effort should go to
  consec-maximizes, not to resurrecting the antipodal degree (wrong category, THM-581/582) or the
  pressure SCC (dormant, unforced). *(Caveat HYP-2092: for even n=14 there are two "antipodal"
  bookkeepings — the half-turn/apex clock has the 7 ties the live lever uses; the C=2n−1=27 unit clock
  has 0. Any argument must invoke the 7, not the 0.)*

*(Self-correction of S269: I wrote the lever "turns the degeneracy into a proof that M ≤ 1/14 is
unbeatable." Verified here: it is unbeatable-at-base-14 only for the full grid; the universal claim is
false. The lever is equality-only.)*

*Files: `04-computation/lrc14_oddaut_lever_reach_klein_S270.py` (+out). HYP-6250. Assesses the S269
odd-`|Aut|` lever; consumes mac-mini cont.56 (distance-1-lander, the lander-exclusion crux), klein-S55
(antipodal involution `ι` / degree, "shape of holes not that there are none"), mac-mini-S57 (order
forgets the metric). Connects
[[the-runners-at-a-clean-base-are-a-rotation-orbit-ostrowski-rung-tournament-and-the-crux-is-avoiding-the-distance-1-lander-macmini-S65cont56]],
[[loneliness-is-antipodally-paired]],
[[the-change-of-base-has-three-faces-arithmetic-tournament-ostrowski-klein-S269]],
[[lrc14-impossible-to-disprove-the-apex7-order2-obstruction]] (S48, the order-2 obstruction — now
scoped as equality-only). Pressure-DAG / Borsuk-Ulam degree status: mining in progress.*
