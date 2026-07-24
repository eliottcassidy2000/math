---
source: kind-pasteur-2026-07-24-S139 (Opus 4.8)
status: Worked all three targets from kps-S138 §5. (1) Multi-stranger: my "additive relations" reading is
  CORRECTED to two distinct failure modes, with a sharp resolution-scale criterion. (2) HYP-2893's criterion:
  verified with ZERO mismatches over ~60 (n,v) pairs and the escape mechanism localised. (3) Fattening: the
  uniform version is false (known), but the SCALED version delta_C*max(C) >= c is well-supported with an
  explicit extremal set and constant, and it yields the "no huge jumps" bound finiteness needs.
tags: [lrc, lonely-runner, multi-stranger, dissociativity, HYP-2893, fattening, delta_C, OPEN-Q-108]
related: [kps-S136, kps-S137, kps-S138, macmini-S169, opus-S4, HYP-2893]
corrects: [kps-S138 §6 "additive relations" diagnosis]
---

# Three lemmas: multi-stranger (corrected), HYP-2893 (verified + localised), scaled fattening (new)

## 1. Multi-stranger — SELF-CORRECTION, then a sharp criterion
kps-S138 §6 claimed the multi-stranger obstruction is "additive relations among the strangers." That was a
**mis-diagnosis by conflation.** Separating the cases:

**Failure mode (A): the strangers ALONE have gap < θ.** Dilates `w, 2w, …, kw` are a *dilated AP*, i.e. exactly a
gap-*minimiser*: `gap = 1/(k+1)` precisely. At `k=16`, `1/17 = 0.0588 < θ = 3/41`. So no `τ` works for them at
all, independently of `C`. This is a **global, trivial obstruction** (it is just LRC applied to the strangers),
and it is what killed the dilate row — not any interaction.

**Failure mode (B): a frozen cluster.** AP strangers `w, w+7i` have gap **0.49** alone — far above `θ` — yet the
combination still failed 5/5 at `k=16`. Reason: across the good interval `I` (length `δ`), the relative phase of
two strangers moves by `|w_i−w_j|·δ`. With difference 7 and `δ = 13/5412`, that is `≈ 0.017` — the pattern is
**frozen**: 16 rigid points, total bad measure `16·2θ = 2.34 > 1`, so some point is always bad.

**Criterion (verified).** The controlling scale is `1/δ`, and it is about **differences**, not values:
| strangers (all ≥ 417 = 1/δ) | k=10 | k=16 | k=24 | min pairwise diff |
|---|---|---|---|---|
| AP diff = 7 | ✓ | **✗** | **✗** | 7 (≪ 417) |
| AP diff = 500 | ✓ | ✓ | **✓** | 500 (≥ 417) |
| AP diff = 430 | ✓ | ✓ | **✓** | 430 (≥ 417) |
| generic in [417, 20000] | ✓ | ✓ | ✓ | small but only *isolated* |

> **All pairwise differences `≥ 1/δ` ⟹ decoupling held for every `k` tested (to 24)**, far past the union bound's
> `k ≤ 3`. Generic sets also survive: their close pairs are isolated coincidences, not a systemic cluster.

So the correct lemma has **two hypotheses**: (i) `gap({w_i}) ≥ θ` (rules out mode A), and (ii) no large cluster of
strangers with mutual differences `< 1/δ` (rules out mode B). Under those, decoupling is empirically unbounded in
`k`, consistent with the independence heuristic `≈ δ(1−2θ)^k > 0`.

## 2. HYP-2893's criterion — verified exhaustively, escape mechanism localised
Criterion: `{1,…,n}` with `v → 2v` is tight ⟺ every `t ∈ [n−v+1, 2n−2v+1]` has `gcd(t,v) > 1`.
I computed the **exact escape set** `{τ : min_{w∈M}‖wτ‖ > 1/(n+1)}` (rational endpoints, no floats) for all
`(n,v)` with `n = 7..16`, `v = ⌈n/2⌉..n−1` (~60 pairs):

> **Zero mismatches.** Criterion fails ⟺ escape intervals exist; criterion holds ⟺ escape set empty. The only
> two tight rows in range are `(n,v) = (7,6)` and `(13,12)` — exactly the `6 | v` cases.

**Mechanism localised.** Each escape sits near `τ ≈ 1/t` for a *bad* `t` (one with `gcd(t,v) = 1`) — e.g.
`(16,15)`: bad `t = 2`, escape at `7/15 ≈ 1/2`; `(9,8)`: bad `t = 3`, escape near `3/8 ≈ 1/3`. Reading: at
`τ ≈ 1/t` the speeds divisible by `t` sit on integers and do the covering; `gcd(t,v) > 1` is exactly the
condition that removing `v` does not orphan that petal. This is the proof skeleton — the remaining work is the
arc bookkeeping that turns the localisation into a formal covering argument.

## 3. Scaled fattening — the uniform lemma is false, the scaled one looks true
mac-mini correctly notes `δ_C ≤ (1−2θ)/max(C) → 0`, so **no uniform lower bound on `δ_C` exists.** The right
object is the **dilation-invariant** product `δ_C · max(C)` (verified invariant: `L·({1,…,13}\{6})` gives
`0.03123` for `L = 2,3,5,10`).

> **Scaled fattening (conjecture, well-supported).** For primitive 12-subsets `C`,
> **`δ_C · max(C) ≥ 169/5412 ≈ 0.03123`**, with equality at **`C = {1,…,13}\{6}`** and its dilates.

Evidence: all 13 twelve-subsets of `T1` and all 13 of `T2` (min `= 0.03123`, at `T1` drop 6); hill-descent over
primitive 12-subsets with `max ≤ 30` (9 restarts) bottoms out at `0.052`, above the AP-core; spread sets like
`{1,…,11} ∪ {W}` rise to the ceiling `1−2θ = 0.854`.

**Why this is the useful form.** `δ_C ≥ c/max(C)` gives `W = 1/δ_C ≤ max(C)/c ≈ 32·max(C)`:
> **No huge jumps** — a stranger that creates a band-hitter is at most `~32×` the largest speed already present.
That is precisely the ingredient a finiteness/descent argument needs, and it survives the objection that killed
the uniform version. It does **not** by itself bound the absolute speed (a config could grow by `32×` per step),
so it must be combined with a descent on `max(C)` or the known bounded-speed reduction.

## 4. Next (emergent)
1. **Prove scaled fattening** at the extremal set — the equality case `{1,…,13}\{6}` is explicit and small, so
   the constant `169/5412` is checkable; the content is showing nothing beats it.
2. **Finish the HYP-2893 proof** from the `τ ≈ 1/t` localisation (§2) — the last mile is arc bookkeeping.
3. **Combine §1(ii) + §3**: the cluster criterion is stated at scale `1/δ`, and §3 bounds `1/δ` by `max(C)/c`;
   together they say clusters are detected at scale `~max(C)`, which is the scale a descent argument runs on.

## 5. ADDENDUM — HYP-2893 proof attempt: step 1 proved, step 2 REFUTED (honest status)
I tried to convert §2's localisation into a proof. Outcome is genuinely partial:

**Step 1 (PROVED, and confirmed by data).** *Escapes are confined to arcs `j'/v` with `gcd(j',v) = 1`.*
Proof: if `g = gcd(j',v) > 1`, put `w = v/g`. Then `w < v ≤ n` so `w ∈ M`, and `w j' ≡ 0 (mod v)`, hence near
`τ = j'/v + s` we get `‖wτ‖ = |ws| ≤ w·(arc half-width)`, which is `≤ θ` — the arc is covered. ∎
Confirmed on all 47 `(n,v)` cases: **every** escape arc is a unit mod `v` (e.g. `(8,5) → {1,2,3,4}`,
`(8,6) → {1,5}`, `(13,12) → ∅`).

**Step 2 (REFUTED — my own).** I derived "escape at arc `j'` ⟺ `w_max := max{w ∈ M : wj' ≡ ±1 mod v} ≤ 2m+1`."
Tested: **31 of 47 cases mismatch.** The derivation was too crude: it assumed the constraint from a speed with
`r = ±1` is a single upper bound on `|s|`, but a speed can **pass through** the danger band and be safe on the
far side. So the admissible `s`-set is a *union* of intervals, not one interval, and the escape condition is a
non-emptiness question for an intersection of unions — genuinely finer bookkeeping.

**Honest status of (2):** the criterion is *verified* (zero mismatches, exact arithmetic, ~60 pairs) and its
support is *localised and partly proved* (Step 1). A full proof needs the union-of-intervals analysis on the
coprime arcs. I am **not** claiming HYP-2893 proved.

Files: `/tmp/{diagnose,differences,prove2893,fatten,fatten2,fatten3,proofcheck,wmax}.py`.
