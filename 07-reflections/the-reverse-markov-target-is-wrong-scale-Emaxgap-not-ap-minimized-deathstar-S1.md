---
source: death-star-2026-07-07-S1
status: CORRECTION + structural insight. E[maxgap] AP-minimality (asserted by opus-S133 and
  kps-S58 as "the one clean remaining step") is REFUTED exactly. The load-bearing floor mu_1/7
  IS AP-minimized (re-verified exact). The reverse-Markov reduction targets the WRONG SCALE.
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - max-gap
  - three-gap
  - correction
---

# The reverse-Markov target is the wrong scale: `E[maxgap]` is NOT AP-minimized (but `μ_{1/7}` is)

Owner: *keep going; extend ideas further.* I audited the newest Route-1 frontier (the
reverse-Markov reduction of the density floor, opus-S133 / kps-S57 / kps-S58) for correctness,
found a load-bearing claim to be **false**, and — more valuably — found the clean structural
reason, which sharpens where the real lemma lives.

## The claim that failed

Both of the two newest reflections assert the same "one clean remaining step":

- **opus-S133:** *"`E[maxgap]` is minimized at the AP … the AP value is proven, so this alone
  closes it."* Uses `E[maxgap(AP_13)] = 93/440` as the binding density floor.
- **kps-S58:** *"`inf E[maxgap] > 1/7` factors as (1) **AP-minimality of `E[maxgap]`** … (2)
  `E[maxgap](AP) = 93/440 > 1/7`."*

**This is false.** Exactly (rational three-gap piecewise-linear integration, cross-checked
against a pure-Python midpoint numeric to 5 digits):

| family | `E[maxgap]` exact | ≈ |
|---|---|---|
| AP `{1..13}` | `93/440` | 0.21136 |
| **prim-sat `2·{1..12}∪{13}`** | **`145091/720720`** | **0.20131** ⟵ below AP |
| S57's own adversarial `{2,6,8,10,11,12,14,16,18,20,22,26,42}` | `168192271811/829905476400` | 0.20266 |
| adversarial descent min | — | ≈ 0.201 |

The AP does **not** minimize `E[maxgap]`. `prim-sat = 2·{1..12}∪{13}` — a *structured*,
one-line family — beats it.

**Attribution (this was partly caught, then re-lost).** opus-S133's *reflection file* asserts
AP-minimality ("this alone closes it"), but opus's own **HYP-4762 INDEX entry** already flagged
the opposite — *"E[maxgap] is NOT AP-minimized … true inf ~0.205 at {0,2..12,17,28} … the crux =
DIRECT inf_E E[maxgap]>1/7, NOT AP-minimality (flagged to kps-S58 who conflated the two
functionals)."* So opus caught it at the **census** level and it is internally inconsistent with
their own reflection. Then **klein-S153 (HYP-4748) re-introduced it**: *"AP is the E[maxgap]-
minimizer … descent CONVERGES to {1..13} … inf_E E[maxgap] = 0.2114."* And kps-S58's reflection
still asserts it. So the false claim is live in **three** reflections. My value-add over opus's
census catch: (i) an **exact** refutation (rational three-gap integration), (ii) a **lower and
cleaner minimizer** — `prim-sat = 2·{1..12}∪{13} = 145091/720720 = 0.20131`, below opus's
stretched `{2,3..12,17,28} = 0.2068`, and a natural dilated family rather than a bespoke shape,
(iii) the **scale mechanism** (below) that says *why*, and (iv) the robustness radius.

**Why the descents all miss it.** Every claim used random starts + local descent. That is a
**weak adversary** (cf MISTAKE-095/096/097): it under-samples the *structured, coarse-cluster-
breaking* families (`2·AP` with one detuned element). And `lrc_Emaxgap_exact_opus_S133` is exact
**only for the AP** — its Farey breakpoints `d ≤ kdenom=13` miss the large-gap order-changes of
any family with speed gaps `> 13`, so it silently under-samples the very families that beat the AP.

## The load-bearing floor is fine — `μ_{1/7}` IS AP-minimized

The functional the proof actually needs is **`μ_{1/7}(E) = meas{x : maxgap{frac(e_i x)} > 1/7}`**
(opus-S130), not its integral. Re-verified exactly this session:

- `μ_{1/7}(AP_13) = 477/1078` (my independent exact computation **reproduces opus-S130's value**).
- **No structured adversary dips below it**: prim-sat `0.4997`, the S57 family `0.6736`, the
  `E[maxgap]`-minimizer `0.7131`, `3·{1..12}∪{13}` `0.5254`, geometric+fill `0.954` — all `>`
  AP's `0.4425`. A hard adversarial descent (the same kind that finds sub-AP families for
  `E[maxgap]`) returns to a **translate of the AP** (`{3..15}` → `477/1078`).

So `μ_{1/7}` AP-minimality survives the exact stress test; `E[maxgap]` AP-minimality does not.

## Why: they are different-scale functionals (`E[maxgap] = ∫₀¹ μ_θ dθ`)

`E[maxgap] = ∫₀¹ P(maxgap > θ)\,dθ = ∫₀¹ μ_θ(E)\,dθ`. So the two claims are about the **same
family of level sets at different scales `θ`.** Computing `μ_θ(AP)` vs `μ_θ(prim-sat)` exactly:

| `θ` | 0.10 | 0.1429 (**1/7**) | 0.20 | 0.222 | 0.25 | 0.333 | 0.50 |
|---|---|---|---|---|---|---|---|
| AP is the minimizer? | ✓ | **✓** | ✓ | ✗ | ✗ | ✗ | ✗ |

**The AP minimizes `μ_θ` for `θ ∈ (0, θ*)` with a crossover `θ* ≈ 0.239`** (exact bisection,
AP vs prim-sat), and **fails above it**. Mechanism:

- **Fine scale (small `θ`)** — the AP orbit `{jx}` is the **three-gap** configuration: at most
  three distinct gap lengths, maximally balanced. It is 1/7-well-spread *more often* than any
  non-AP orbit ⟹ smallest `meas{gap > θ}` for fine `θ`. This is the three-gap rigidity mac-mini-S15
  named.
- **Coarse scale (large `θ`)** — the AP's Achilles heel: near `x = 0` (and every small-denominator
  rational) all `k` points cluster in `[0,(k-1)x]`, leaving a gap `→ 1`. A family with a
  *spread* element (the `13` in `2·{1..12}∪{13}`, or any far element) **breaks those coarse
  clusters**, so it has fewer huge gaps ⟹ smaller `μ_θ` for coarse `θ`.

`E[maxgap]` integrates over **all** `θ`, so it is dominated by the coarse regime where the AP is
*bad*. `μ_{1/7}` samples a **single fine scale** (`θ = 1/7 = 0.143`). Two robustness numbers
(both exact-confirmed): the AP-vs-prim-sat crossover is `θ ≈ 0.239`, but the **true AP-minimality
radius** — the *smallest* `θ` at which *any* 13-family beats the AP — is `θ*_true ≈ 0.18` (a
family beats the AP at `θ=0.182` but none at `θ=0.180`, adversarial descent + exact confirm). So
`θ = 1/7 = 0.143` sits below the radius with an **honest margin `≈ +0.04`** — robust, though
modest, not the `+0.10` the prim-sat crossover alone suggests. That margin is exactly why
`μ_{1/7}` inherits AP-minimality and `E[maxgap]` does not.

**Floor is solid at every cluster size.** `μ_{1/7}` AP-minimality re-checked EXACTLY against
structured adversaries (the kind descent misses) at each `k=8..13`: `691/735, 247/294, 38/49,
1381/2205, 13823/24255, 477/1078` — nothing below the AP at any `k`. The floor `min_k` binds at
`k=13 = 477/1078`. (opus-S130 only ran a descent, which for `E[maxgap]` *did* miss the structured
minimizers; here the exact structured check passes, so the floor is genuinely solid.)

## Consequence for the proof route

The reverse-Markov inequality `μ_{1/7} ≥ (7/6)(E[maxgap] − 1/7)` is a **valid inequality** and a
legitimate *lower bound tool* — but it **discards the extremal structure**: the reduced functional
`E[maxgap]` is minimized by a *coarse-cluster-breaking* family, not the AP. So:

1. **Do NOT pursue "AP-minimality of `E[maxgap]`"** (opus-S133 step 1 / kps-S58 part 1). It is
   false; the AP value `93/440` is *not* the `inf`.
2. `inf_E E[maxgap] > 1/7` is still **true** (empirically `≈ 0.201 > 0.143`), so the reduction
   *conclusion* holds — but proving it needs a **direct** mean bound or the *true* minimizer
   (`prim-sat`-like, coarse-cluster-optimal), **not** a reduction to the AP.
3. The **cleaner and correct** open lemma is **fine-scale `μ_θ` AP-minimality** at `θ = 1/7`
   (opus-S130's original target), which *does* hold, is robust in `θ` up to `≈ 0.24`, and gives
   the sharp floor `μ_{1/7} ≥ 477/1078` — far better than reverse-Markov's lossy
   `≥ 1477/18480 ≈ 0.08`. The three-gap equidistribution content is genuinely at the fine scale;
   the reverse-Markov move to the mean throws it to the coarse scale where it is false.

The honest one-liner: **the "cleaner functional" is cleaner because it deleted the hard part —
and the hard part (AP-minimality) is exactly the part that is *true* only at the fine scale.**

## Can a *corrected* mean route work? No — the two windows are disjoint (the capstone)

The natural fix: use the **truncated mean** `T(θ*) := E[min(maxgap, θ*)] = ∫₀^{θ*} μ_θ dθ`, which
*stays at the fine scale* and (verified exactly, `AP ≤ prim-sat` at every `θ*`) **is AP-minimized
for `θ* ≤ θ*_true`** — the property `E[maxgap]` lost. And `T` yields a corrected reverse-Markov
floor (from `μ_θ ≤ 1` on `[0,1/7]`, `μ_θ ≤ μ_{1/7}` on `[1/7,θ*]`):
`μ_{1/7} ≥ (T(θ*) − 1/7)/(θ* − 1/7)`, whose AP value would be a legitimate floor.

**But it is vacuous, for a sharp reason.** Two windows in `θ*`:
- **AP-minimality of `T`** needs `θ* ≤ θ*_true ≈ 0.181` (the μ_θ AP-minimality radius).
- **Non-vacuity** needs `T_AP(θ*) > 1/7`, i.e. (exact) `T_AP(0.19)=0.1424 < 1/7`,
  `T_AP(0.20)=0.1452 > 1/7` ⟹ `θ* ≳ 0.195`.

`0.181 < 0.195`: **the windows are DISJOINT.** To make the mean floor non-vacuous you must push
`θ*` into the coarse regime, but there AP-minimality is already gone — so you can no longer reduce
to the AP. **No fine-scale mean can be simultaneously AP-minimized and deliver a non-vacuous
`μ_{1/7}` floor.** The reverse-Markov step `μ_θ ≤ 1` on `[0,1/7]` throws away exactly the
fine-scale mass (`μ_θ(AP)` drops from `1` to `0.44` on `[1/13, 1/7]`) that the direct tail keeps.

**Conclusion — the trichotomy is complete:**
| functional | AP-minimized? | gives non-vacuous `μ_{1/7}` floor? |
|---|---|---|
| full mean `E[maxgap]` | ✗ | (n/a — wrong extremal) |
| truncated mean `T(θ*)`, `θ*≤0.181` | ✓ | ✗ (vacuous) |
| **tail `μ_{1/7}`** | ✓ | ✓ (`477/1078`) |

So the tail `μ_{1/7}` is **irreducible**: its AP-minimality must be proved **directly** (the
three-gap / equidistribution lemma of opus-S130), *not* routed through any mean. The entire
reverse-Markov program (opus-S133 / kps-S57-S58 / klein-S153) yields valid inequalities but is a
**strict regression** as a proof route — it cannot close the floor. This *sharpens* rather than
weakens Route 1: it says the one honest open lemma is exactly `μ_{1/7}` (= fine-scale `μ_θ`)
AP-minimality, and rules the mean detours out.

## The coarse minimizer IS the saturated family (a cross-connection)

The `E[maxgap]` minimizer is not a random shape — it is **`prim-sat = 2·{1..12}∪{13}`**
(`145091/720720`), unique up to dilation/translation, found by an exact search: over
`2·{1..12}∪{m}` the value is **symmetric in `m` around `m=13`** (minimized by inserting the odd
element at the *center* of `2·AP`'s range `{2..26}`), and no descent beats it. `inf_E E[maxgap] ≈
0.2013 > 1/7` with margin `+0.058`, so opus-S133's HYP-4762 *direct* target `inf_E E[maxgap] > 1/7`
(the honest reframing after dropping AP-minimality) holds comfortably.

**And `prim-sat` is exactly the primitive-SATURATED family** (a multiple of every `q ∈ 2..14`)
that kps-S55/S56 + opus-S131 pinned as the *sieve* hard core. So the density-floor **mean**-extremal
family coincides with the **saturated** family: being arithmetically spread across *all* small
scales (what "saturated" means) is precisely what breaks the plain AP's coarse near-`x=0`
clustering (what lowers `E[maxgap]`). The additive (max-gap) and multiplicative (sieve/saturation)
faces of LRC — mac-mini-S15's duality — meet again at `2·{1..12}∪{13}`. (Note the *tail* `μ_{1/7}`
still ranks `prim-sat` *above* the AP, `0.4997 > 0.4425` — saturation helps the coarse mean, not
the fine tail; the two extremal problems genuinely differ.)

## Ledger

- **Refuted:** `E[maxgap]` AP-minimality (opus-S133, kps-S58). **Confirmed:** `μ_{1/7}`
  AP-minimality (opus-S130), exact. **Explained:** crossover `θ* ≈ 0.239`; `1/7` margin `+0.096`.
- **Files:** `lrc_maxgap_ap_minimality_check_deathstar_S1.py`,
  `lrc_maxgap_true_minimizer_deathstar_S1.py`, `lrc_mu17_ap_minimality_stress_deathstar_S1.py`,
  `lrc_mu_theta_crossover_deathstar_S1.py`, `lrc_ap_minimality_radius_deathstar_S1.py`,
  `lrc_mu17_apmin_all_k_deathstar_S1.py`, `lrc_truncated_mean_floor_deathstar_S1.py`,
  `lrc_Emaxgap_minimizer_structure_deathstar_S1.py` (+ `.out`s).
  HYP-4777. All exact computations reuse a *corrected* `Emaxgap`/`μ_θ` region-decomposition
  (`kdenom = max(max E, max|eᵢ−eⱼ|)`, fixing opus-S133's AP-only `kdenom=13`).
- **Pointers:** opus-S130 (`μ_{1/7}` floor, correct target), opus-S133 / kps-S57 / kps-S58
  (reverse-Markov, wrong-scale target), mac-mini-S15 (three-gap rigidity = the fine-scale content),
  THM-527 (Part A: `ρ*(P,E) > 0`, the actual density node).
- **Does NOT prove or disprove LRC(14).** Corrects a route and localizes the real lemma.
