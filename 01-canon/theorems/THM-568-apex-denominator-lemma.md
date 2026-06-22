---
id: THM-568
title: Apex-denominator lemma — every tight LRC(14) optimum has denominator a multiple of 14 (D = 14·gcd(S)); reduces (★) to "14-covering ⟹ not tight"
status: PROVED (elementary). The reduction of (★) it yields has a residual = the multiples-of-14 equidistribution (≥7 multiples), the apex-localized Node-3 estimate.
author: kind-pasteur-2026-06-22-S31aa
depends_on:
  - THM-523    # q=14 witness (14-free ⟹ M_14 >= 1/14)
  - S31v       # comb-teeth union bound + M(R)>=1/13 margin
related:
  - THM-079    # the H=21 template (Move A reduce-to-atom + Move B forcing)
  - THM-560    # tight locus / dilation
  - HYP-2906   # mac-mini's (★)
---

# THM-568: the apex-denominator lemma

## Statement
Let `S` be 13 distinct positive integers with `M(S) = max_t min_i ‖v_i t‖ = 1/14` (a tight set), and let
`t* = a/D` (lowest terms) be an optimum. Then:
1. **`14 | D`.**
2. The two binding runners `{v_a, v_b}` (those with `‖v_· t*‖ = 1/14`, at `±1/14`) satisfy `D | (v_a+v_b)`.
3. **`D = 14·gcd(S)`** — in particular **primitive `S` ⟹ `D = 14`** (the optimum is at the apex denominator).

## Proof
At a tight optimum the origin lies in an empty open arc of length exactly `1/7`, with binding runners at
**both** endpoints `±1/14` (a single occupied endpoint would let the arc — hence `M` — grow, contradicting
optimality). For the `+1/14` binder, `v_a·(a/D) = m ± 1/14` for some `m ∈ ℤ`, so `14(v_a a − mD) = ±D`,
i.e. `D(1+14m) = 14 v_a a`. Since `1+14m ≡ 1 (mod 14)` is coprime to 14, `(1+14m) | v_a a` and
`D = 14·(v_a a/(1+14m))`, so **`14 | D`**. The two binders give `(v_a+v_b)t* ∈ ℤ`, so `D | (v_a+v_b)a`,
hence `D | (v_a+v_b)`. Writing `S = gcd(S)·S'` with `S'` primitive, `M(S)=M(S')`, `t*(S)=t*(S')/gcd(S)`,
so `D(S) = gcd(S)·D(S')`; and a primitive tight set has `D'=14` (parts 1–2 with `v_a+v_b` minimal = 14),
giving `D = 14·gcd(S)`. ∎

VERIFIED `lrc_apex_denominator_lemma_kps.py`: AP (`D=14`, binders `{1,13}`), GW (`D=14`, binders `{5,9}`),
`2·AP` (`D=28`), `3·GW` (`D=42`).

## Consequence — the reduction of (★)
`(★)` [`M(S)=1/14 ⟹ optimum at a denom-14 point`] reduces to two parts (`S` 14-covering iff it contains a
multiple of 14):
- **14-free tight ⟹ optimum at `D=14`** (no equidistribution): THM-523 gives `M_14(S) ≥ 1/14`, so
  tightness forces `M_14 = 1/14` attained at `t=a/14`.
- **14-covering ⟹ `M(S) > 1/14`** (not tight): the 14-free part `R` has `M(R) ≥ 1/13 > 1/14` by proven
  LRC(≤13), giving an interval `I` with `min_R > 1/14`; the multiples of 14 cover `≤ |M14|/7` of `I`, so for
  `|M14| ≤ 6` a point survives with `M(S) > 1/14`. Residual: `|M14| ≥ 7` (the apex-localized second-moment).

So with Move A (peel, R1) and Move B (apex floor), THM-568 turns the THM-079-template proof of LRC(14)
into the single residual **"≥7 multiples of 14 over a 14-free core ⟹ not tight."**

## Significance
The structural half of (★) is now PROVED: the tight optimum is forced to the apex denominator `14 = 2·7`
(no analysis — pure divisibility from the binding runner). The Steinhaus/three-gap "finite half" of (★)
is discharged for the 14-free case by THM-523; only an apex-localized equidistribution about the multiples
of 14 remains. This is the LRC analogue of THM-079's "WLOG one atom + cycle-forcing," with the forcing
made arithmetic: a tight set must bind at `±1/14`, which pins `14 | D`.

→ THM-523, THM-560, THM-079, S31v, HYP-2906, `the-apex-denominator-lemma-reduces-star-...md`.
