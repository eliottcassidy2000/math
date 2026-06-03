---
id: THM-398
title: LRC reduces to C′ (multiple-of-n ⟹ loose); the dominance-dodge lemma
status: PROVED (the reduction and the dodge lemma); C′ itself remains CONJECTURE
source: opus-2026-06-03-S572
depends_on:
  - THM-369   # divisibility / 1-clock witness
  - LRC(n-1)  # proven in the literature for n-1 <= 13
related:
  - HYP-2102  # the reduction, discovered S571
  - HYP-2097  # the 64-class worry container
  - HYP-2095  # the lift lemma (paired/anchored split)
---

# THM-398 — LRC reduces to C′, and the dominance-dodge lemma

## 0. Setup and definitions

Throughout `n ≥ 2`; a **speed set** `S = {v_1,…,v_{n-1}}` is `n-1` distinct positive
integers (the moving runners; the observer is runner `0`). For `t ∈ ℝ/ℤ` write
`‖x‖ = dist(x, ℤ)`. Define

```
M(S) = max_{t} min_{v∈S} ‖v t‖,
```

the maximal gap (max-collar). **LRC(n)** is the statement `M(S) ≥ 1/n` for every such
`S`. It is a theorem for `n ≤ 13` (literature) and open at `n = 14`.

The **safe set at level `1/n`** is `G(S) = { t : ‖v t‖ > 1/n for all v ∈ S }`, an open
subset of the circle. Note `μ(G(S)) > 0 ⟺ M(S) > 1/n` (one direction: a point of `G`
has all `‖v t‖ > 1/n` so `M > 1/n`; conversely the optimum `t*` with
`min ‖v t*‖ = M > 1/n` has an open neighbourhood still in `G` by continuity). We call
`S` **loose** if `M(S) > 1/n` and **tight** if `M(S) = 1/n`.

> **C′ (conjecture).** If `n | v` for some `v ∈ S`, then `M(S) > 1/n` (`S` is loose).

C′ was observed in S564 as a *property of tight sets* ("tight ⟹ no multiple of `n`").
Here we prove its converse direction is the entire conjecture.

## 1. Lemma A (the 1-clock witness; restatement of THM-369)

**If no element of `S` is a multiple of `n`, then `M(S) ≥ 1/n`.**

*Proof.* Evaluate at `t = 1/n`. For each `v ∈ S`, write `r = v mod n`. Since `n ∤ v`,
`r ∈ {1,…,n-1}`, so `‖v/n‖ = ‖r/n‖ = min(r, n-r)/n ≥ 1/n`. Thus
`min_{v∈S} ‖v/n‖ ≥ 1/n`, whence `M(S) ≥ 1/n`. ∎

## 2. Theorem (the reduction): **C′ ⟹ LRC(n)**

*Proof.* Let `S` be any speed set.
- If no `v ∈ S` is a multiple of `n`: `M(S) ≥ 1/n` by Lemma A.
- If some `v ∈ S` is a multiple of `n`: `M(S) > 1/n ≥ 1/n` by C′.

Either way `M(S) ≥ 1/n`; this is LRC(n). ∎

**Remark (why this is a genuine reduction).** The hypothesis C′ constrains only the
*proper subclass* of speed sets that contain a multiple of `n` — precisely the sets
on which the elementary `1/n`-clock witness fails. All of LRC's difficulty is thereby
isolated into a single structural statement about the distinguished runner that
vanishes on the whole `n`-clock. (Verified: 0 tight-with-multiple over exhaustive
small boxes and large samples, `n = 4..14`, every multiplier size — HYP-2102.)

**Remark (the Vitali handoff — S551o).** This split *is* the "Vitali wall" of S551o:
LRC = (positive-measure bulk, settled by *measure*) ∪ (measure-zero core, settled by
*construction*), with the Vitali set marking the handoff. THM-398 *locates the
handoff*: it is the equation **`n | v`**. Configs with `n ∤ v` are handled by Lemma A
— the `t=1/n` construction, which is *measure-blind* and so reaches the measure-zero
core (the worry-set has no multiple of `n`, S564). Configs with `n | v` are pushed to
the *measure* side (C′). On that side the danger of `v=nw` is a **periodic,
bounded-eccentricity arc family** (a genuine Vitali cover), so the **Vitali covering
lemma / Lebesgue density theorem** is the natural tool — see §3–§4 and HYP-2104.

## 3. Lemma B (the dominance-dodge): a partial proof of C′ — **and more**

For a speed set `S` and `v ∈ S` write `V'(v) = max(S \ {v})`.

> **Lemma B.** Assume LRC(n−1). If some `v ∈ S` satisfies `v > (n-1)·V'(v)`, then
> `M(S) > 1/n` (`S` is loose). *(No divisibility hypothesis: this holds for any
> dominant runner, multiple of `n` or not.)*

*Proof.* Put `S' = S \ {v}`, with `|S'| = n-2` distinct positive integers, and
`V' = V'(v) = max S'`. By LRC(n−1) applied to `S'` (proven for `n-1 ≤ 13`),
`M(S') ≥ 1/(n-1)`, so there is `t₀` with `‖u t₀‖ ≥ 1/(n-1)` for every `u ∈ S'`.

*The S′-window.* For `|t - t₀| < δ` with
```
δ := 1 / ( n (n-1) V' ),
```
each `u ∈ S'` satisfies `‖u t‖ ≥ ‖u t₀‖ - u·|t-t₀| > 1/(n-1) - V'·δ
   = 1/(n-1) - 1/(n(n-1)) = 1/n`.
So on the open interval `I = (t₀ - δ, t₀ + δ)` every runner of `S'` is `> 1/n`.

*Dodging v.* The danger set of `v` is `D_v = { t : ‖v t‖ < 1/n }`, a union of open
arcs of radius `ρ := 1/(n v)` centred at the points `k/v`, separated by safe gaps.
An interval can lie inside `D_v` only if it lies inside a *single* arc (a gap is
`v`-safe, so an interval spanning a gap meets the `v`-safe set). Hence if
`|I| = 2δ > 2ρ`, i.e. **`δ > ρ`**, then `I ⊄ D_v`, so `I` contains a point with
`‖v t‖ > 1/n` as well — and an open neighbourhood of it, since `I ∩ (circle ∖ D̄_v)`
is open and nonempty.

*The threshold.* `δ > ρ ⟺ 1/(n(n-1)V') > 1/(nv) ⟺ v > (n-1)V'`, exactly the
hypothesis. On the resulting open set every runner of `S` exceeds `1/n`, so
`μ(G(S)) > 0` and `M(S) > 1/n`. ∎

**Corollary B1 (large multiples).** If `n | v` and `v > (n-1)·V'(v)` then `S` is
loose. In particular every multiple-of-`n` config whose multiple dominates the rest
is handled — and for `n = 14` this draws on the *proven* `LRC(13)`.

**Corollary B2 (one dominant runner ⟹ loose).** Any speed set with a single runner
exceeding `(n-1)×` all others is loose, regardless of divisibility. *(Verified:
1500/1500 loose at each of `n = 6,8,10,12,14`, dominant runners chosen with arbitrary
residues — `lrc_dodge_formalization_s572.py`.)* This is the **clean general fact**
that the divisibility-flavoured C′ partial is a special case of: the dodge is about
**dominance**, not about multiples of `n`.

## 4. The sharpened criterion and the residual dichotomy

Inspecting the proof, Lemma B only used that *some* component interval of `G(S')`
(the safe set of the `n-2` runners at level `1/n`) is longer than one danger arc of
`v`. This is weaker than `v > (n-1)V'`:

> **Criterion B′.** If `G(S\{v})` has a connected component of length `> 2ρ = 2/(nv)`,
> then `S` is loose. *(Proof identical: that long interval cannot sit inside one
> `v`-arc, so it meets `G(S)`.)*

**Vitali-covering iff (HYP-2104).** Because `G(S\{v})` is a finite union of intervals
and `D_v`'s gaps are open and nonempty, *an interval of `G(S\{v})` lies in `D_v` iff
it fits inside a single arc*. So Criterion B′ is the Vitali-covering direction, and a
multiple-of-`n` config is tight (measure-0) *only if* every component of `G(S\{v})` is
both short **and** arc-aligned. Quantified: B′ already proves looseness for **72%
(n=6) → 96.8% (n=14)** of multiple-of-`n` configs, and the all-short residual is
**never tight** (0 across n=6..14) — `lrc_vitali_covering_residual_s573.py`.

This yields the honest **dichotomy** for a multiple-of-`n` config `S` (`v = nw`):

1. **Long-interval case** — `G(S\{v})` has a component longer than `2/(n·nw) = 2/(n²w)`:
   `S` is loose by Criterion B′. *(Always holds when `v` dominates; Lemma B.)*
2. **All-short case** (the **residual**) — every component of `G(S\{v})` is `≤ 2/(n²w)`:
   the dodge of a single interval is not guaranteed. Here `S` is still loose
   (verified), but for the equidistribution reason
   ```
   μ(G(S)) = μ(G(S')) − μ(G(S') ∩ D_v),   μ(G(S') ∩ D_v) ≈ (2/n)·μ(G(S')) < μ(G(S')),
   ```
   i.e. the single arithmetic progression of thin arcs `{k/(nw)}` (period `1/(nw)`,
   total danger `2/n`) cannot *align to cover* the fixed union of intervals `G(S')`.
   Proving this non-covering is the remaining open core of C′ — a discrepancy /
   three-distance statement about one AP against a fixed open set.

**The residual is sharpest at `v = n` (`w = 1`):** the arcs sit on the `n`-clock
points `k/n` (radius `1/n²`), and one must find a point in some clock gap
`(k/n + 1/n², (k+1)/n − 1/n²)` where all `n-2` other runners exceed `1/n`. The gap
*midpoint* `(2k+1)/(2n)` (where `v=n` is maximally safe, at `1/2`) suffices only about
half the time (773/1499 at `n=6`), so the witness genuinely ranges over the gap
interior — confirming the residual is a true interval-search, not a fixed sub-clock.

## 4½. The cover→congruence translator and Lemma C (small-owner uncoverability)

Each endpoint of a component of `G(S')` carries an **owner**: a left endpoint is
`a = (k_a n + 1)/(n u_a)` (runner `u_a` just turned safe), a right endpoint is
`b = (k_b n − 1)/(n u_b)`. The condition "(a,b) lies in the `v=nw` arc centred at
`j/(nw)`" *translates*, on clearing denominators (`× n u w`), into the **endpoint-owner
congruences**
```
|w(k_a n + 1) − j·u_a| < u_a/n     and     |w(k_b n − 1) − j·u_b| < u_b/n.
```
For an owner `u < n` the right side is `< 1`, and the bracket is an integer, so it is
**forced to 0** — the endpoint must equal the arc centre: `a = j/(nw)` (resp. `b`).
This rigidity gives a new w-free theorem:

> **Lemma C (small-owner uncoverability).** If a component `(a,b)` of `G(S')` has
> *both* owners `u_a < n` and `u_b < n`, then `S = S' ∪ {nw}` is loose for **every**
> `w`. *Proof.* If `(a,b)` sat inside one arc centre `c=j/(nw)`, the small-owner
> congruences force `a=c` and `b=c`, hence `a=b`, contradicting `a<b`. So `(a,b)` is
> not inside any arc; being connected with the arcs separated by nonempty `v`-safe
> gaps, it contains a `v`-safe point — a point of `G(S)`. ∎

(The "cross-relation" `u_b(k_a n+1) = u_a(k_b n−1)` that two small owners would need is
exactly `a=b`: verified `0/490, 0/1330, 0/2594, 0/5298` for `n=6,8,10,12`.) The slack
in the congruence appears **only** for owners `u ≥ n` (window `±u/n ≥ ±1`), so an
endpoint can sit off-centre inside an arc *only* when its owner is large.

**Translator verified:** the congruence-window characterisation matches direct
tight/loose computation **100%** (`2500/2500` each `n=6..14`). **Coverage by the two
proved criteria** (Lemma C ∪ Criterion B′), multiple-of-`n` configs:

| n | small-owner-component (Lemma C) | Lemma C ∪ B′ | residual |
|---|---|---|---|
| 6 | 8.2% | 73.4% | 531 |
| 8 | 18.8% | 81.7% | 366 |
| 10 | 33.4% | 92.0% | 160 |
| 12 | 56.2% | 96.0% | 80 |
| **14** | **81.3%** | **99.0%** | **19** |

The proved coverage **grows toward the frontier**: at `n=14` the two criteria
discharge `99%` of multiple-of-14 configs. The residual (`~1%`) is exactly the configs
where every component is short **and** every component has a **large** (`≥ n`) binding
owner — the only regime where the congruence slack permits an off-centre fit.

## 4¾. The endpoint-cover criterion and circuit positivity (Lemma D)

The cover condition has an exact one-line form per component. Let `C_i = (a_i, b_i)` be
a component of `G(S')`, with midpoint `m_i = (a_i+b_i)/2` and length `ℓ_i = b_i-a_i`.

> **Lemma D (endpoint-cover criterion).**
> `C_i ⊆ D_v  ⟺  ∃ j∈ℤ: ‖v a_i − j‖ ≤ 1/n and ‖v b_i − j‖ ≤ 1/n  ⟺  ‖v m_i‖ ≤ 1/n − (v/2)ℓ_i.`
> *Proof.* `C_i ⊆` arc `(j/v − 1/(nv), j/v + 1/(nv))` ⟺ `v a_i, v b_i ∈ (j−1/n, j+1/n)`;
> a common integer `j` exists iff the interval `[v b_i − 1/n, v a_i + 1/n]` — midpoint
> `v m_i`, half-length `1/n − (v/2)ℓ_i` — contains an integer, i.e. `‖v m_i‖ ≤ 1/n −
> (v/2)ℓ_i`. ∎

Hence **`S` is tight ⟺ every component satisfies Lemma D**, and the **circuit-positivity
margin**
```
P(S) := max_i ( ‖v m_i‖ + (v/2)ℓ_i − 1/n )      satisfies   P(S) > 0  ⟺  S is loose.
```
**Verified 100%** (criterion vs direct tight/loose, `2492/2492` … `700/700`, n=6..14;
`lrc_endpoint_cover_circuit_positivity_s575.py`), and `P(S) > 0` for **every** sampled
multiple-of-`n` config.

**The circuit.** The `M` components' arc indices `j_i` wind once around the circle —
`Σ_i (j_{i+1}-j_i) = v` — so the per-component conditions are not independent: one
integer `v` must *simultaneously* bring every midpoint within `1/n` (phase) of an arc
centre. **Summing** Lemma D over the circuit gives a provable necessary condition:
```
S tight  ⟹  Σ_i ‖v m_i‖  ≤  M/n − (v/2)·μ(G(S'))  <  M/n,
```
i.e. the components' midpoints must have **average v-phase distance `< 1/n`**. The
actual average is **`≈ 0.245`** (near the generic `1/4`) for all n=6..14, versus the
requirement `1/n` (`= 0.071` at n=14): a positivity margin that **grows toward the
frontier**. So tightness would demand an anomalous simultaneous resonance of `v` with
*all* component midpoints — the circuit obstruction, the sharpened residual of C′.

## 5. Status ledger

| Statement | status |
|---|---|
| Lemma A (no-multiple ⟹ `M ≥ 1/n`) | **PROVED** (THM-369) |
| Reduction **C′ ⟹ LRC(n)** | **PROVED** (§2) |
| Lemma B / Cor B1 (dominant multiple dodge, uses LRC(n−1)) | **PROVED** (§3) |
| Cor B2 (one dominant runner ⟹ loose) | **PROVED** (§3), verified 1500×5 |
| Criterion B′ (long component ⟹ loose) | **PROVED** (§4) |
| Lemma C (both-small-owner component ⟹ loose, w-free) | **PROVED** (§4½), covers 81% at n=14 |
| Lemma C ∪ B′ | **PROVED**, covers **99%** of multiple-of-14 configs |
| Lemma D (exact endpoint-cover criterion + circuit positivity `P>0 ⟺ loose`) | **PROVED** (§4¾), verified 100% n=6..14 |
| Summed corollary (tight ⟹ avg midpoint v-phase `< 1/n`) | **PROVED**; actual avg `≈0.245 ≫ 1/n` |
| C′ residual = `P(S) > 0` always (simultaneous-resonance infeasibility) | **OPEN** (~1% at n=14; the circuit Diophantine obstruction) |

So LRC(14) now sits on a single open assertion: *the thin evenly-spaced danger arcs
of a multiple of 14 cannot cover the safe set of the other twelve runners* (the
all-short case), with the long-interval/dominant case fully proved from the literature's
LRC(13).

**Artifacts:** `04-computation/lrc_dodge_formalization_s572.py` (+`.out`),
`lrc_lift_Cprime_residual_s571.py`, `lrc_lift_lemma_measure_bound_s571.py`.
See reflection `07-reflections/lrc-formalizing-the-Cprime-reduction-and-dominance-dodge-s572.md`.
