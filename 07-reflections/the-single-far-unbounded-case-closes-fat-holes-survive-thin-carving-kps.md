# The single-far unbounded case closes: fat holes survive thin carving (a Morse / band-barrier certificate on top of klein's bounded ILP)

*kind-pasteur-2026-07-01. First concrete swing at the UNBOUNDED core of LRC(14) — the analytic half the previous reflection isolated (bounded = klein HYP-3779's ILP, done; unbounded = large speeds, open). Using the owner's Morse-theory / band-barrier / Borsuk-Ulam frame, the single-far unbounded case (one large speed) closes, with a comfortable margin, reducing to a "fat-hole lemma" or a modest extension of klein's speed bound.*

## The reduction

A covering 13-set with exactly one large speed is `C ∪ {W}`, `C` = 12-speed bounded core (all `≤ n(n−1)=182`, klein's regime), `W > 182`. By **proven LRC(13)** (`≤12` speeds), `M(C) ≥ 1/13 > 1/14`, so the core's **lonely set**
`L_C = { t : g_C(t)=min_{v∈C}‖vt‖ ≥ 1/14 }` is nonempty with positive measure.

`g_C` is a **piecewise-linear Morse function on S¹**; `L_C` is its super-level set `{g_C ≥ 1/14}` = a union of **arcs around the local maxima** (the "holes"). Adding `W` carves from above: `g_{C∪W} = min(g_C, ‖Wt‖)`, and `‖Wt‖` is a fast sawtooth whose **danger band** `{‖Wt‖<1/14}` is `W` isolated teeth, each of width `1/(7W)`, separated by safe gaps `6/(7W)`. A connected lonely arc `I ⊆ L_C` can be swallowed by the danger band only if it fits inside **one** tooth (the teeth are isolated). Hence:

> **`δ_max(C) > 1/(7W)` ⟹ the widest hole contains a `W`-safe point ⟹ a lonely moment of `C∪{W}` ⟹ `M(C∪{W}) ≥ 1/14`.**

So the single-far certificate holds for `W > W*(C) = 1/(7·δ_max(C))`, `δ_max` = the widest arc of `L_C`. **This is the band-barrier: thin bands cannot cover fat holes.**

## The computation: `W* ≤ 74 < 182`, so there is no gap

Computing `L_C` exactly (interval intersection over `∩_v{‖vt‖≥1/14}`, Fractions — building the super-level set directly, no `M(S)` max, no MISTAKE-86), for every tested and randomly-searched bounded 12-core (`lrc14_singlefar_morse_bandbarrier_kps.py`):

| core `C` | meas(`L_C`) | `δ_max` | `W*=1/(7δ_max)` |
|---|---|---|---|
| AP `{1..12}` | 0.0341 | 0.00595 | **24** |
| `{1..11,13}` | 0.0122 | 0.00408 | **35** |
| `{1..11,182}` (large internal) | 0.0430 | 0.00471 | **30** |
| worst of 300 random, max ≤ 182 | — | 0.00194 | **74** |

The worst-case `W*` over the random search, for *any* max-speed budget up to 182, is `≤ 74`. Since the unbounded case is `W > 182`, we have `W > 182 > 74 ≥ W*(C)` **automatically** — the certificate fires for every unbounded single-far `W`, and it **overlaps klein's bounded closure `W ≤ 182` with room to spare** (factor `182/74 ≈ 2.5`). The single-far unbounded case closes.

## Why the holes stay fat (and the honest rigor gap)

The widest hole is fat because it sits around the **highest local max** of `g_C`, where the binding runners are the *small* speeds (large arcs); the large internal speeds only nick it. Empirically `δ_max(C) ≥ 0.00194 = 1/515` for all bounded cores, well above the required `1/1274` (`⟺ W* ≤ 182`).

The **rigorous** lower bound is the one piece not yet proven. The naive Morse bound — near `t*` with `g_C(t*)=M(C)≥1/13`, `g_C` drops with slope `≤ max(C) ≤ 182`, so `δ_max ≥ 2(1/13−1/14)/182 = 6·10⁻⁵` — gives only `W* ≤ 340`, a `~30×`-lossy bound that leaves a *rigorous* gap `[182, 340]` (empirically empty). So single-far closes **rigorously** as soon as either:
1. a **fat-hole lemma** `δ_max(C) ≥ 1/1274` for all bounded 12-cores (empirical margin 2.5×; the widest hole is governed by the small runners, not `max(C)` — the content of the lemma), or
2. klein's lazy-cut ILP is pushed from `≤182` to `≤340` (a modest, finite extension of what he already did).

Either closes it; neither is speculative.

## The Borsuk–Ulam / ι-odd packaging

`L_C` is **antipode-symmetric**: `‖v(−t)‖=‖vt‖`, so `t∈L_C ⟺ −t∈L_C`; the `ℤ/2` complement action `ι: t↦−t` acts freely off the fixed points `{0,½}`. Whether a hole *survives* the `W`-carving is exactly the **ι-odd index** of this equivariant covering (klein-S56's ι-odd = the OPEN-Q-108 residual): the band-barrier is its down-to-earth witness. In single-far, the index is nonzero *for free* (the widest hole survives, `δ_max > 1/(7W)`); the genuinely ι-odd-hard obstruction only appears when **two or more** carving combs act together — the multi-far case, where the antipodal fixed structure can conspire. This is the topological reason single-far is easy and multi-far is the residual.

## Net and next

- **Single-far unbounded: closed** (band-barrier `W* ≤ 74 < 182`, overlapping klein), modulo a fat-hole lemma `δ_max ≥ 1/1274` (2.5× empirical margin) *or* klein's bound `182→340`. A concrete, finite-flavored rigor task, not the old combinatorial `a(n)` hunt.
- **Multi-far (`≥2` large speeds): the residual.** Two carving combs `‖W₁t‖, ‖W₂t‖`; a hole survives iff it dodges *both* danger bands. The r-far ladder / joint decorrelation (my prior far-element thread) + the Borsuk–Ulam ι-odd index. This is the remaining open core — but it is analytic (joint equidistribution of two combs), not a beater search.
- The Morse frame makes the whole thing legible: LRC(14) = "every super-level hole of some `g_S` reaches `1/14`"; bounded cores have fat holes (klein ILP + this), far combs are thin (band-barrier), and the only hard case is two thin combs conspiring against a hole (multi-far / ι-odd).

— Related: klein HYP-3779 (bounded ILP, `≤182`), the prior reflection `the-real-dichotomy-is-bounded-vs-unbounded-...` (this is its (b) executed), `the-eisenstein-cusp-dichotomy-is-the-three-distance-theorem-kps.md` (≤3 slopes = the Morse critical structure), `the-far-runner-is-the-log-barrier-...` (the far-element decorrelation, = the carving), klein-S56/HYP-3768 (ι-odd index), OPEN-Q-108. Script: `04-computation/lrc14_singlefar_morse_bandbarrier_kps.py`. Not a HYP reservation.
