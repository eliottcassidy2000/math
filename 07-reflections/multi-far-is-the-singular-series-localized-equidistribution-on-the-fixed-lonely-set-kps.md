# Multi-far is the singular series, localized: the far combs equidistribute on the fixed lonely set, single-far closes, and the residual is exactly the signed resonance correction

*kind-pasteur-2026-07-01. The two-comb swing on the single-far result. The right object turned out not to be the widest hole but the whole (fragmented) lonely set of the bounded core, with the far runners equidistributing relative to it — and that recovers the singular series, cleanly separating the closed single-far case from the OPEN-Q-108 residual.*

## The correction to the widest-hole picture

Single-far closed via the band-barrier on the *widest hole*. The natural next swing — does a fat hole survive *two* combs? — exposed that the widest hole is the wrong object for the spread cores. The worst core for two-far, `{1,26,74,94,122,130,161,172,174,176,177}`, has a thin widest hole (`δ_max ≈ 0.002`) yet **survives every comb pair exactly** (min jointly-safe `= 0.114 > 0`). Why: its lonely set `L_C = {t : g_C(t) ≥ 1/14}` is **fragmented into 300 thin holes of total measure `0.167`**, and two thin combs cannot cover 300 scattered holes. The certificate lives in the *whole* `L_C`, not one hole.

## The mechanism: equidistribution on the fixed `L_C` = the singular series

`L_C` is fixed (it depends only on the bounded core `C`, klein's regime). A far comb `‖Wt‖` (`W` large) equidistributes, so it keeps a `6/7`-fraction of `L_C` **regardless of where the holes are** (verified):

> `meas(L_C ∩ safe(W)) → (6/7)·meas(L_C)` as `W→∞` (Erdős–Turán/Koksma, rate `O(#holes / W)`).

`r` far combs keep the **main term** `(6/7)^r · meas(L_C) > 0` for every `r` (independent limit), minus the corrections from **resonances** (relations among the combs). Verified, two-comb worst case:

| core | `meas(L_C)` | `(6/7)²·meas` (main) | worst 2-comb survival | resonance correction |
|---|---|---|---|---|
| AP `{1..11}` | 0.0563 | 0.0414 | 0.0389 | **0.0025** |
| spread (worst) | 0.1667 | 0.1225 | 0.1137 | **0.0088** |

This is **exactly the singular series**, localized to the far combs on `L_C`:

> `survival(r) = (6/7)^r · meas(L_C) − [signed resonance correction]`  ⟷  `L = (6/7)^13 + Σ_T (−7/6)^{|T|} R_T`.

When all 13 speeds are far, `meas(L_C)=1` and the main term is `(6/7)^13` — the known singular-series bulk. The far-element/Morse frame *recovers the singular series* as the equidistribution of the combs on the bounded core's lonely set.

## The clean separation this gives

- **Single-far (`r=1`): CLOSED.** One comb cannot self-resonate, so there is **no correction**: `survival = (6/7)·meas(L_C) > 0`, exactly. (This *is* the band-barrier, now seen as `r=1` equidistribution — `meas(L_C ∩ safe(W)) → (6/7)meas`, verified `≥0.856·meas` down to the thinnest cores.)
- **Multi-far (`r=2..6`): the residual, now pinned.** The main term `(6/7)^r·meas(L_C)` is positive; the open content is the **signed resonance correction**, and it is *small* — `0.008` against a main term `0.122` at `r=2`, a comfortable factor. Bounding it below the main term is **OPEN-Q-108** = the r-far Dedekind ladder = the Riesz-product "factor of 2" (my prior sessions). The equidistribution *rate* is the far-element decorrelation `Δ_W = O(1/W)` (last session).
- **`r ≥ 7`: THM-573** (level-7 sieve: `≥7` far ⟹ `M>1/14`).

So the **whole unbounded case** = `[r=1 closed]` + `[r=2..6: main term positive, residual = the signed resonance correction, the known crux — now localized to a FIXED bounded core's L_C and shown small]` + `[r≥7 sieve]`.

## The reframe (what to pay attention to)

The question stops being "does a beater exist" (a config search over the disordered a(n)) and becomes **"is the signed resonance correction smaller than the main term `(6/7)^r·meas(L_C)`, uniformly over bounded cores `C` and far combs?"** — an *analytic* bound on a signed sum over a *fixed* set, which is precisely where the mature tools live (Riesz product, decorrelation, Erdős–Turán). The Morse/band-barrier picture the owner proposed and the singular series are the same object: holes of `g_C` = the lonely set; far combs carve by equidistribution; the correction is the ι-odd resonance obstruction (Borsuk–Ulam), which vanishes for one comb (single-far) and is the small signed residual for `≥2` (multi-far). This is real hope: the residual is not combinatorial and not large — it is a small, signed, analytic correction on a bounded object.

## Honest status

- **Verified:** single-comb equidistribution `→(6/7)meas` (single-far closes cleanly); two-comb survival `> 0` with the correction small (`0.008` vs main `0.122`); the `(6/7)^r` main term positive for `r=2..6`.
- **Not proved:** that the signed resonance correction stays below the main term *uniformly* over all bounded cores and all far combs — that is OPEN-Q-108/HYP-3132, unchanged in difficulty but now precisely localized and connected to the analytic tools (not a beater search). The equidistribution *rate* (how large `W`) is the finite constant-chase.
- **Union bound is too lossy** (`lrc14_rfar_unionbound_closure`: fails `r=2..6` on spread cores because boundary teeth eat the margin); the *exact* equidistribution is what survives — signalling that the correction, not the measure, is the content.

— Related: the prior reflections `the-real-dichotomy-is-bounded-vs-unbounded-...` and `the-single-far-unbounded-case-closes-...` (this is their multi-far continuation), klein HYP-3779/3781 (bounded ILP + Steinhaus tail), `lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md` + `the-far-runner-is-the-log-barrier-...` (the signed correction / decorrelation rate), THM-573 (level-7 sieve), HYP-3132 (r=2..6 floor), THM-504 (the singular series), OPEN-Q-108. Scripts: `04-computation/lrc14_{twofar_twocomb_resonance,rfar_unionbound_closure,multifar_equidistribution_singularseries}_kps.py`. Not a HYP reservation.
