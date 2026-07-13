# The covering-min residual is "core runner 1"; my structural slice is its formalized extremal

*kind-pasteur-2026-07-11-S127 cont.62. Owner: "work the covering-min lower bound over all covering families."
Reading the fleet's freshly-pinned residual (mac-mini-S74, today: the equidistribution route closes `|core| ≥ 2`,
the residual is `|core| = 1`), I found that its **hard sub-case is core runner 1** — and that this is exactly
where my eight sessions of structural work + the single-killer Lean formalization live. The residual tiles:
structural (mine) over the small-core-runner extremal, analytic (opus/mac-mini) over the rest.*

---

## The residual, precisely (mac-mini-S74 / opus-S259)

A covering family splits into the **core** = speeds coprime to `30030 = 2·3·5·7·11·13` (the primes `≤ 13`) and
the **non-core** = the rest (each divisible by a prime `≤ 13`). LRC(14) for `v` ⟺ the core does not cover the
good set `G'` (`coreCover < 1`). The equidistribution union bound **closes `|core| ≥ 2`** (`coreCover ≈
1 − (6/7)^{|core|} < 1`, opus-S259). The residual is **`|core| = 1`**: a single coprime runner is one arc, does
not equidistribute, and its density on `G'` reaches `0.92` — so `coreCover < 1` there is the whole open crux.

## The hard sub-case is core runner 1 (verified)

The `|core| = 1` case still has a free parameter: *which* coprime runner. It is not uniformly hard — **the
smallest core runner is the unique hard case.** Swapping the core runner `1 → c` in `{c} ∪ {2,…,12, 182}`:

| core runner `c` | `M` | vs `14/183` |
|---|---|---|
| **1** | **14/183 ≈ 0.0765** | = (the floor) |
| 17 | 13/92 ≈ 0.141 | 2× looser |
| 19, 23 | 13/92 ≈ 0.141 | 2× looser |
| 29 | 4/31 ≈ 0.129 | looser |

A larger core runner has more arcs, equidistributes better, and its density on `G'` drops toward `1/7` — easy.
Only the **smallest coprime runner, `1`, concentrates** enough to approach `coreCover = 1`. So the covering-min
crux is not "`|core| = 1`" broadly; it is specifically **`|core| = 1` with core runner `1`** — the least-
equidistributed, most-concentrated single arc.

## That sub-case is exactly my structural domain — and partly formalized

Every covering-min extremal is `|core| = 1` with core `= {1}`:

- deep well `{1..12, 182}` → `M = 14/183`; ladder `S_2 = {1..12, 364}` → `28/365`; `S_3` → `42/547`;
  multi-killer `{1..11, 13, 84}` → `7/89` — all core `= {1}`.

These are precisely the **interval-core** families I have been analyzing, and the results transfer directly:

- **Single-killer** `{1,…,12, 182c}` (`c ≥ 1`): `reach ≥ 14c/(182c+1) ≥ 14/183`, **PROVED and machine-checked
  kernel-pure** (cont.60/61, `LRCSingleKillerLadder.lean`). This is the core-runner-1 extremal ladder — the exact
  families where the single arc (runner 1) is most concentrated — closed in Lean.
- **Multi-killer** `{1,…,k, outliers}` (`k ≤ 11`): all `≥ 14/183` by exact enumeration (cont.58), the covering-min
  monotone in core length `M ≥ 1/(k+1)`, reducing (cont.59) to LRC(13)-escape + finite check anchored by
  `1/13 > 14/183`.

So the residual's hard core (small core runner) is where the **structural** picture is sharpest and the
extremal is **formalized**, while the equidistribution / single-arc-density analytic machinery (opus-S259/S260,
mac-mini-S74) governs the larger-core-runner and `|core| ≥ 2` bulk.

## The tiling — and the honest open piece

> **The covering-min lower bound tiles by core structure:** `|core| ≥ 2` and large core runners →
> equidistribution (opus/mac-mini, coreCover `< 1`); `|core| = 1` with core runner `1` → the concentrated single
> arc, whose extremal (deep well) and single-killer ladder are structurally pinned and Lean-formalized (kps).

What remains genuinely open is the **general core-runner-1 `|core| = 1` family** — one runner `1` plus twelve
*arbitrary* `13`-smooth non-core runners (not necessarily `{2,…,12}` + a killer). My interval-core results cover
the slice whose non-core body is `{2,…,12}` (which contains the extremal), and the bounded such families are
klein's ILP (`≤ 182`); the general smooth body still needs the analytic bound `density(D_1 ∩ G') < 1` (opus's
mollified Erdős–Turán, open). So the frames are complementary: **structure pins and formalizes the extremal and
its ladder; analysis must close the concentrated single arc over all smooth bodies.** The one number that runs
through both is `1/13 > 14/183` — why the single arc has room to miss `G'`, and why multi-killer can never bind.

*Files: lrc14_core1_residual_kps_S127.py (+.out). Connects mac-mini-S74 (|core|=1 residual), opus-S259/S260
(equidistribution / Erdős–Turán) to kps cont.58–61 (multi-killer, monotonicity, single-killer ladder formalized).
HYP-6230.*
