---
id: THM-1002
title: The pair-sum denominator bound (M(A)=val/q with q | v_i+v_j, hence q ≤ 2·max A) and the resulting PROVED bounded case of the analytic stability gap — every 12-set with max(A) ≤ 18 satisfies M(A)=1/13 or M(A) ≥ 2/25
status: PROVED (elementary; the bounded case max(A)≤18 only). The GENERAL stability bound "non-AP ⟹ M ≥ 2/25" is CRUX (C) and remains OPEN — this file proves a sub-case and pins the exact obstruction to extending it.
source: klein-2026-07-17-S313
depends_on:
  - THM-999   # death-star-S56 Lemma A — PRIOR ART for §1 (tight case); §1 here is the general-M form
  - THM-666   # witnesses live on pair-sum rulers (the maximizer-location fact used here)
related:
  - HYP-7310  # klein-S313 n=12 tight census (AP {1..12} unique; the Farey gap reproduced)
  - HYP-6820  # the n=12 uniformity audit (sporadic branch)
  - THM-996   # death-star resonance confinement (tight ⟹ lonely only at resonant times)
external: LRC(13) SETTLED (gives M(A) ≥ 1/13 for every 12-set).
---

# THM-1002 — the pair-sum denominator bound, and the bounded case of the stability gap

Owner directive: *prove the analytic stability bound `non-AP ⟹ M ≥ 2/25`*. That statement is
**CRUX (C)** of the PROOF-MAP ("no integer 12-family has `M ∈ (1/13, 2/25)`") and is **open** — it is
strictly stronger than n=12 AP-uniqueness (HYP-6820), which is itself Tao's optimistic conjecture
(HYP-7310). This file proves the **bounded case** and isolates exactly what blocks the general one.

Throughout `A = {v_1<…<v_12}` is a set of 12 positive integers and `M(A) = max_t min_{v∈A} ‖vt‖`.

## 1. The pair-sum denominator bound (the engine)

> **PRIOR ART — attribution.** This lemma is **death-star's THM-999 Lemma A**
> (`death-star-2026-07-17-S56`, pushed before this file), stated there for **tight** families
> (`M = 1/n`): *every loneliness time has reduced denominator dividing `v₊+v₋ ≤ 2·Vmax`*, with the
> same opposite-sides active-pair proof. Credit for the lemma is theirs. The only thing added here is
> that **tightness is not needed**: death-star's step 1 invokes tightness to make maxima isolated, but
> `f = min_v ‖vt‖` is a min of pieces of slope `±v ≠ 0`, hence *never locally constant*, so every
> maximum is a breakpoint for any family at any value of `M`. That general-`M` form is what the gap
> application below requires, since the gap concerns **non-tight** families with `M` just above `1/13`.

**Lemma (general-`M` form of THM-999 Lemma A).** Write the maximizer as `t* = a/q` in lowest terms.
Then `q | (v_i + v_j)` for some `i ≤ j`, and consequently **`q ≤ 2·max(A)`**. Moreover
`M(A) = val/q` where `val = min_v |va|_q` (`|x|_q := min(x mod q, q − x mod q)`).

*Proof.* `f(t) = min_v ‖vt‖` is piecewise linear; each `‖vt‖` rises with slope `v` until `vt ≡ 1/2`,
then falls with slope `−v`. At a local maximum of the minimum, an active rising branch meets an
active falling branch (a single branch's own peak cannot be the min, since there `‖vt‖ = 1/2` is
maximal and 11 other speeds are smaller). Equating them, `v_i t − k = k' − v_j t`, so
`(v_i + v_j)t = k + k'` and `t = (k+k')/(v_i+v_j)`. In lowest terms `t = a/q` with `q | (v_i+v_j)`,
hence `q ≤ v_i + v_j ≤ 2·max(A)`. ∎

*(Consistency: AP `{1..12}` → `q = 13 = 1+12`, `val = 1`, `M = 1/13`. `{1..11,24}` → `q = 25 = 1+24`,
`val = 2`, `M = 2/25`. Both verified by the exact pair-sum evaluator.)*

## 2. The gap admits no small numerator

`M = val/q ∈ (1/13, 2/25)` ⟺ `12.5·val < q < 13·val`. For `val = 1` this is `q ∈ (12.5,13)` and for
`val = 2` it is `q ∈ (25,26)` — **both integer-free**. Hence any violation has `val ≥ 3`, and then
`q > 37.5`, i.e. **`q ≥ 38`**. (First admissible pairs: `(val,q) = (3,38), (4,51), (5,63), (5,64),
(6,77), (7,88), (7,89), (7,90), …`)

## 3. The proved bounded case

**Theorem.** If `max(A) ≤ 18` then `M(A) = 1/13` or `M(A) ≥ 2/25`.

*Proof.* By §1, `q ≤ 2·max(A) ≤ 36`. By §2 a gap violation needs `q ≥ 38 > 36`. ∎

Equivalently: **a gap violation requires `max(A) ≥ 19`.**

## 4. Why this does not extend (the exact obstruction)

The natural next step — rule out each admissible `(val,q)` by a *discrete* necessary condition — **fails**.
WLOG the witness is `a = 1` (lowest terms ⟹ `a` invertible; relabel `A` by `a`), so `R = A mod q` lies in
`{val,…,q−val}`, and `M ≤ val/q` forces every `k ∈ Z/q` to have some `c ∈ R` with `|ck|_q ≤ val`, i.e.
the sets `B_c = {k : |ck|_q ≤ val}` must **cover** `Z/q` using `≤ 12` of them. Computed: for every
admissible `(val,q)` with `q ≤ 149`, a greedy cover of size **7–10 ≤ 12** exists. So the residue-covering
condition is satisfiable at every in-gap denominator and rules out none of them.

**Therefore the obstruction to the general bound is integer realizability, not residue arithmetic** — which
is precisely the hard LRC content of CRUX (C).

## 5. Erratum (recorded because it was caught by testing)

An earlier attempt claimed the crossing forces speeds `≡ ±1 (mod q)`, which would have given
`max(A) ≥ q−1 ≥ 37`. That is **FALSE**: the crossing gives `v_i a ≡ +val` and `v_j a ≡ −val (mod q)`,
so `v_i ≡ val·a^{-1}`, which equals `±1` only when `val = 1`. An empirical check refuted it at
**332/400** random primitive 12-sets. The surviving, correct consequence is `q | (v_i+v_j)` (§1).

## 6. Evidence beyond the proof

Exact pair-sum evaluation (validated against a fine grid; deviations ≤ 2.9e−4, within grid resolution
`max(A)/NG`): **zero** gap violations among ~174,000 primitive 12-sets — 119,998 with `max ≤ 36` and
53,740 with `max ≤ 60`. Consistent with CRUX (C) far beyond the proved range, but not a proof.

*Files: `04-computation/lrc14_gap_reduction_klein_S313b.py` (+ `lrc14_gap_hunt/_structure/_reduction` .out).
The n=12 twin of death-star's THM-996 resonance confinement: extremal times are pinned to pair-sum rulers.*
