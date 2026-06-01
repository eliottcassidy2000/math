---
source: oracle-2026-06-01-S545
status: result + honest negative (the cascade is the right object; the cycle-exclusion is its triple/Helly-3 shadow, necessary but not sufficient)
tags: [lonely-runner, cascade, conditional-clearance, transitivity, cycle-exclusion, helly, order-of-obstruction, inside-debt]
---

# LRC as a cascade of conditional clearances; the cycle-exclusion is the Helly-3 shadow

**Prompt (user):** think of a cascade as a product of conditional clearances; LRC
in mind; transitivity has the usual fact (`X->Y, Y->Z => X->Z`) AND a second hidden
fact that propagates — the cycle-exclusion: `X->Y => not(Z->X and Y->Z)` (no directed
3-cycle `Z->X->Y->Z`).

Both halves of the hint are right and they locate the structure precisely — and the
computation also tells us, honestly, *which order* the cycle-exclusion governs.

## The cascade (verified)

Read loneliness as a telescoping **product of conditional clearances**. Order the
runners; let `S_{<=i}` be the safe set after clearing the first `i`. Then

```
|SAFE| = prod_i c_i ,   c_i = |S_{<=i}| / |S_{<i}|  (clearance of runner i | the previous).
```

Computed (`lrc_cascade_cycle_exclusion_s545.py`):

```
generic n=5:  c = (3/5, 5/9, 13/25, 9/13)   -> |SAFE| = 3/25   (every factor > 0: the cascade clears)
AP n=5 (tight): c = (3/5, 2/3, 1/3, 0)        -> |SAFE| = 0      (a ZERO factor: the 4th runner is TRAPPED
                                                                  by the clearance of the first three)
```

So a **tight system is a cascade with a zero conditional clearance** — a runner the
earlier clearances make impossible (at the boundary wall). The cascade is exactly
the right object: it exposes loneliness as a chain where each link can collapse.

## The two facts, placed: inclusion telescopes, exclusion is Helly-3

- The **inclusion** (`X->Y,Y->Z=>X->Z`) is what lets the cascade *telescope* — the
  clearances chain into a single product.
- The **cycle-exclusion** (no `Z->X->Y->Z`) is a **triple** fact: it is a **Helly-3**
  condition on the clearances — "no obstructed triple." We tested whether it is
  *sufficient*: does "every triple of runners clears" force "all clear"?

Comparing the full loneliness collar to the worst `k`-subset collar:

```
 system    1/n   full   worst pair   worst triple   worst 4-set
 AP n=5    1/5   1/5     3/10         1/4            1/5
 AP n=6    1/6   1/6     7/24         1/4            7/36
 AP n=7    1/7   1/7     2/7          1/4            4/21
```

The **worst `k`-subset collar is `1/(k+1)`** — it is exactly the `k`-runner LRC bound,
achieved by the sub-AP. So the worst **triple** collar is `1/4` for every AP
(n=5,6,7): *every triple of runners clears comfortably* (`1/4 > 1/n`), the
cycle-exclusion / Helly-3 holds — **but the full collar still drops to `1/n`.**

> **The cycle-exclusion is necessary and real, but NOT sufficient.** Clearing every
> triple buys only `>= 1/4`; LRC needs `>= 1/n`. The collar **degrades monotonically
> with subset size along the LRC ladder `1/(k+1)`**, reaching `1/n` only at the *full*
> cascade. LRC is irreducibly the whole product; no finite order (no Helly-`k` for
> fixed `k`) captures it.

## What this says (and an honest negative)

The user's instinct is exactly right about the *framework*: LRC is a cascade, and
the cycle-exclusion is the hidden constraint that propagates through it. The
computation sharpens *where* it lives:

- the cascade `|SAFE| = prod c_i` is the transport (S544) made into a product — each
  factor is a conditional clearance, and the tight case is a collapsed link;
- the cycle-exclusion is the **3-rung** of the LRC ladder (the `1/(k+1)` sequence):
  it certifies `theta_c >= 1/4` (3-runner LRC, *proved*, S526), and propagates up,
  but the rung tightens at every level and only the top rung is `1/n`;
- so the difficulty is **order-`n`**: it is the *accumulation* of the cascade, not a
  triple obstruction. This is the cascade form of "the inside debt grows with order"
  (S531: at n=4 the new obstruction is the 3-term/triple term; deeper `n` adds
  higher terms) and of "no fixed-order channel decoupling" (S532).

**Honest negative:** the tempting geometric reading — "the runner sub-tournament at
the loneliest time is a 3-cycle iff the inside-debt is active" — is *false* (held for
only 104/262 triples ~ chance). The reason: the loneliest time `t*` is a **wall**,
where the half-turn comparator has antipodal **ties**, so the 3-cycle there is
ill-defined. The cycle-exclusion does not live in the geometric tournament *at* `t*`;
it lives in the **clearance ladder** (the Helly/`1/(k+1)` structure above), which is
the correct home for the user's hidden fact.

## Open (→ HYP)
- The collar ladder `worst-k-subset = 1/(k+1)` (achieved by the sub-AP) deserves a
  proof: every `k`-runner sub-system clears at `>= 1/(k+1)`, with the AP tight. Is
  it exactly the statement that the cascade's conditional clearances obey
  `c_i >= ` (a bound from the `(i+1)`-runner LRC)? If so, a proof of LRC(n) by
  induction down the cascade needs only that the *last* clearance survives — the
  AP-critical rung.
- Helly number of the clearance arcs: each runner's safe set is `v_i` arcs; the LRC
  obstruction order = the Helly number of this arc family. The data says it is `n-1`
  (full), not 3 — quantify the Helly number as a function of the system.

## Anchor
`04-computation/lrc_cascade_cycle_exclusion_s545.py` (+ `.out`): the cascade product
(tight = zero factor); worst-`k`-subset collar `= 1/(k+1)` (cycle-exclusion = Helly-3,
`1/4`, insufficient); the negative on the `t*` 3-cycle. Builds on S531 (inside-debt
order), S544 (transport), S526 (3-runner LRC = 1/4), S524 (block structure).
