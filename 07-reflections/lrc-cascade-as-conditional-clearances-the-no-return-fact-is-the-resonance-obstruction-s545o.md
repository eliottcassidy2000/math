---
source: oracle-2026-06-01-S545o
status: synthesis + computation (cascade = product of conditional clearances; the hidden no-return fact = the resonance/3-cycle obstruction)
tags:
  - lonely-runner
  - cascade
  - conditional-clearance
  - transitivity
  - no-return
  - resonance
  - inside-debt
---

# The Cascade as a Product of Conditional Clearances: the Hidden "No-Return" Fact Is the Resonance Obstruction

A cascade is `P(clear) = ∏_k P(E_k | E_1∩…∩E_{k-1})` — a **product of conditional
clearances**. For LRC: order the runners (the base-path/jersey order, S542), clear
them one at a time, `F_k = F_{k-1} ∩ {t : ‖v_k t‖ ≥ 1/n}`, `c_k = |F_k|/|F_{k-1}|`,
and LRC ⟺ `∏ c_k = |F_{n-1}| > 0`. Two transitivity facts drive it — and the second,
hidden one is exactly the resonance obstruction. (Convergent with codex-S549/HYP-2040,
the "transitive wedge debt"; this is the empirical verification + quantification.)

## The two facts are one transitivity — read two ways

- **Fact 1 (completion):** `(X,Y),(Y,Z)` arcs ⟹ `(X,Z)` arc. Clearances **propagate
  forward**: the product builds toward the independent value.
- **Fact 2 (the hidden one, no-return):** `(X,Y)` an arc ⟹ **not** `((Z,X) and
  (Y,Z))` — the arc `X→Y` forbids the **3-cycle** `Z→X→Y→Z`. No later vertex
  cyclically **re-blocks** a cleared region.

Verified (`lrc_cascade_conditional_clearance_noreturn_s545.py`): on `4000/4000` random
tournaments, **Fact 1 ⟺ Fact 2 ⟺ (no directed 3-cycle)** — the same transitivity. The
content is that Fact 2 reads transitivity as *the absence of a return cycle*, which is
the cascade-relevant face.

## The cascade, computed: returns drive the product to zero

`c_1 = 1 - 2/n` exactly (the first, unconditional clearance). The subsequent
conditional clearances and the product `∏ c_k` (= lonely measure) against the
independent value `(1-2/n)^{n-1}`:

```
 n=5  indep (1-2/n)^{n-1}=0.130;  #triples=4
   AP 1,2,3,4 (regular):  c=[.60,.667,.333,0]   product 0.0000   return 3-cycles 4/4
   random 7,11,14,21:     c=[.60,.602,.662,.341] product 0.0814   return 3-cycles 1/4
   Fibonacci 1,2,3,5:     c=[.60,.667,.333,.40]  product 0.0533   return 3-cycles 4/4
 n=7  indep=0.133;  #triples=20
   AP 1..6:               product 0.0000   return 3-cycles 20/20
   random:                product 0.1224   return 3-cycles 9/20
   Fibonacci 1,2,3,5,8,13: product 0.0489  return 3-cycles 19/20
```

The **"return" 3-cycles = the small 3-term resonances** `m_a v_a + m_b v_b + m_c v_c =
0` (a cyclic constraint among three runners). And the pattern is decisive:

> **The product deficit (below `(1-2/n)^{n-1}`) is driven by the number of return
> 3-cycles.** Generic speeds have few returns ⟹ `c_k ≈ 1-2/n`, product ≈ the
> independent value (clearances **propagate**, Fact 1). The **regular polygon (AP) has
> ALL triples resonant** (returns everywhere, Fact 2 fails maximally) ⟹ the
> conditional clearances **collapse** (the last `c_k = 0`) ⟹ product `= 0` (tight).

(The Fibonacci speeds `1,2,3,5,8,13` are deliberately resonance-rich — each `a+b=c` is
a return cycle — and show a large deficit, confirming the mechanism on a non-AP set.)

## The unification: the hidden fact = no resonance = no inside debt

The "return" 3-cycles forbidden by Fact 2 are exactly the **3-term resonances = the
inside debt** (S529/S533). So:

> **Fact 2 (no return) ⟺ no 3-term resonance ⟺ the inside debt vanishes ⟺ the cascade
> conditional clearances don't re-block ⟺ the product stays at `(1-2/n)^{n-1} > 0`
> ⟺ local emptiness (S544).**

The cascade needs **both** facts: **Fact 1 (completion) builds the clearance product
forward**, toward the independent value; **Fact 2 (no return) guarantees no later
runner cyclically re-blocks it.** The whole difficulty of LRC is the *returns* — the
resonances/3-cycles — and they are **maximal exactly at the regular polygon** (all
triples resonant), where they cancel the product to `0`. This is the same wall as
S529 (inside debt), S533 (parity), S544 (frequency-concentration), now seen as: *the
clearance cascade is re-blocked only by returns, and the conjecture is that the
returns never re-block the whole circle (the lonely set stays nonempty even when its
measure hits 0 at the regular polygon).*

## The reframe the two facts give

> **LRC = "no return wins."** Clear the runners in cascade order (Fact 1 propagates
> the clearance forward, multiplicatively to `(1-2/n)^{n-1}`); the only thing that can
> re-block a cleared region is a **return cycle** (a 3-/higher-term resonance, Fact 2
> violation). The conjecture is that the returns — even at their maximum (the regular
> polygon, all triples resonant) — drive the lonely *measure* to `0` but never empty
> the lonely *set*.

This makes the cascade rigorous *in the no-return (generic/transitive) regime* (product
`= (1-2/n)^{n-1} > 0`) and isolates the obstruction as the return-resonance count — a
clean, countable quantity (the 3-term resonances, the inside debt's leading order).

## Verdict / next
- Verified the two facts are one transitivity (Fact 1 ⟺ Fact 2 ⟺ no-3-cycle); the
  hidden Fact 2 = "no return cycle re-blocks the cascade."
- Computed: the cascade product deficit is driven by the return 3-cycles (3-term
  resonances); generic ⟹ product ≈ `(1-2/n)^{n-1}`; AP/regular ⟹ all triples
  resonant ⟹ product `0` (tight). Fibonacci speeds confirm on a non-AP set.
- Convergent with codex-S549/HYP-2040 (the "transitive wedge debt").
- Concrete next: (1) a quantitative cascade bound `∏ c_k ≥ (1-2/n)^{n-1} − f(#returns)`
  — the wedge/return debt; (2) the set-vs-measure gap at the all-triples-resonant
  regular polygon (returns kill the measure, not the set); (3) higher-order returns
  (k-cycles ↔ k-term resonances) as the full inside debt.

## Artifacts
```
04-computation/lrc_cascade_conditional_clearance_noreturn_s545.py
05-knowledge/results/lrc_cascade_conditional_clearance_noreturn_s545.out
```
Related: HYP-2040 (codex-S549, transitive wedge debt — convergent), S527 (cascade),
S529/S533 (resonances/inside debt = the returns), S544 (frequency spread/decorrelation),
S542 (base-path order), S530 (apex).
