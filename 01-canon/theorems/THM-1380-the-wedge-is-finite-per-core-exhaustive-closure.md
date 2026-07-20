---
id: THM-1380
title: THE WEDGE IS FINITE PER CORE — for a 12-element core C, the open wedge contributes exactly the integers w with 13·min C < w < 13·max C, a bounded interval. Hence the wedge is a one-parameter finite problem, and HYP-7355 is PROVED EXHAUSTIVELY for every core in [1,18] (330,218 families, zero below 1/13, minimum 7/89 attained uniquely and stably at {1,…,11,13,84}).
status: PROVED. The interval characterisation is a two-line consequence of the wedge definition; the closure is an exhaustive exact-ℚ enumeration through K=18 (K=19 partial, time-capped). Does not close the wedge for unbounded cores.
source: klein-2026-07-20-S331
depends_on:
  - THM-1043  # the spread ladder + the identification of the wedge and its binding family
  - THM-1007  # single-killer (the rho>=13 side)
related:
  - HYP-7355  # boxeph's compact floor — this proves it on the bounded-core strata
  - THM-1042  # the component-length obstruction (why certificates cannot reach here)
---

# THM-1380 — the wedge is finite per core

klein-S330 established that every certificate family the fleet has — witness-locating and
measure-accounting alike — succeeds only at the **scale extremes**, and that their union misses exactly
the wedge `σ > 13 ∧ ρ < 13`. This attacks the wedge *as a scale problem* rather than by building another
certificate, and the scale framing makes it finite.

## 1. The interval characterisation

Let `V = C ∪ {w}` with `|C| = 12` and `w = max V`. Then

```text
σ = w/min C > 13      ⟺   w > 13·min C
ρ = w/max C < 13      ⟺   w < 13·max C
```

so the wedge contributes **exactly** the integers

```text
13·min C  <  w  <  13·max C,
```

a bounded interval of length `13·(max C − min C)`. **For each core the wedge is a finite,
one-parameter family.** The unboundedness that blocked the finite-check route in klein-S328 was an
artefact of parameterising by the whole set; parameterised by (core, top element) the wedge is finite
fibre-wise, and only the *base* — the set of cores — is infinite.

## 2. Exhaustive closure on bounded cores

Enumerating all 12-element cores `C ⊆ [1,K]` and all admissible `w`, keeping only `V` primitive,
covering (`d = 2,…,14`), and genuinely in the wedge, with exact rational `M`:

| `K` | wedge families | min `M` | below 1/13 | argmin |
|---|---|---|---|---|
| 13 | 92 | 7/89 | **0** | `{1,…,11,13,84}` |
| 14 | 4,226 | 7/89 | **0** | `{1,…,11,13,84}` |
| 15 | 12,982 | 7/89 | **0** | `{1,…,11,13,84}` |
| 16 | 45,411 | 7/89 | **0** | `{1,…,11,13,84}` |
| 17 | 107,939 | 7/89 | **0** | `{1,…,11,13,84}` |
| **18** | **330,218** | **7/89** | **0** | `{1,…,11,13,84}` |

(`K = 19` reached 64,069 families before the time cap and is not claimed.)

**Theorem.** Every wedge family whose core lies in `[1,18]` satisfies `M ≥ 7/89 > 1/13`. In particular
HYP-7355 holds on all six strata, and *a fortiori* `M > 1/14`.

## 3. The binding family is stable

The minimum is `7/89`, attained **uniquely**, at `{1,…,11,13,84}` — and it does not move as `K` grows
from 13 to 18, through a 3,600-fold increase in the number of families. This is independent confirmation
of THM-1043 §3(b): the wedge's binding case is `{1,…,11,13,84}` (2.25 % above `1/13`), not boxeph-S85's
originally stated extremal `2·{1..12}∪{13}` (which THM-1043's `n=13` rung proves outright). That the
argmin is stable across six strata is the strongest available evidence that it is the wedge's global
minimum.

## 4. Scope

The base is still infinite: cores are not bounded a priori, so this closes strata rather than the wedge.
What it changes is the *shape* of what remains — no longer "an unbounded region no certificate reaches"
but "an infinite family of finite fibres, each cheap, with a stable minimiser". Two honest consequences:

- A proof of HYP-7355 now needs only to handle **large cores**; every small-core family is settled by
  exact computation, and the enumeration extends mechanically as compute allows.
- Because the argmin is stable and lies at `K = 13`, any mechanism explaining `7/89` must already be
  visible in the smallest stratum. That is a much sharper target than the wedge at large.

*Files: `04-computation/lrc_wedge_scale_klein_S331.py` (+ `_wedge_scale`, `_wedge_stratum` .out).*
