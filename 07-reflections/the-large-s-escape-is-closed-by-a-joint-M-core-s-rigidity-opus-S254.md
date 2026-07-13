---
source: opus-2026-07-11-S254
status: PARTIAL CLOSURE (honest). The "large-s escape" in S253's slow-fast balance is CLOSED for the
  single-killer-182 covering-min: empirically (11000+ families, binding speed s up to 63, ZERO
  counterexamples; deep well the unique minimizer) and RIGOROUSLY in the tight case s=1 (LRC(13) + prime-13
  uniqueness of the interval core). It reduces to a single sharp joint rigidity M_core >= (182+s)/2379
  (a quantitative refinement of LRC(13) coupling the core value to its binding speed), tight only at the deep
  well -- verified, not yet proved. Remaining: killers != 182 (easier) and multi-killer families.
tags:
  - lrc14
  - covering-min
  - slow-fast-balance
  - large-s-escape
  - joint-rigidity
  - deep-well-minimality
  - partial-proof
---

# The large-s escape is closed by a joint M_core–s rigidity (tight only at the deep well)

**opus-2026-07-11-S254.** Owner: close the large-s escape in the general balance. Worked it — the escape is
closed for the single-killer covering-min (empirically comprehensive + rigorous at the tight case), and reduces
to one sharp rigidity. Honest partial closure.

## The escape, precisely

S253's balance witness gives, for a covering family `{core C} ∪ {killer v_f}`,

`M ≥ M_core · v_f/(v_f + s)`,   `s` = the core's binding speed at its LRC optimum.

With `M_core ≥ 1/13` (LRC(13)) and `v_f ≥ 182` (covering), this clears the wall `14/183` **iff `s = 1`**; for
`s ≥ 2` the factor `v_f/(v_f+s)` shrinks and the witness alone can dip below (`7/92` at `s=2`). The escape is
the worry that some covering family with a **fast** core-binding runner beats `14/183`.

## What closes it: a joint M_core–s rigidity

Requiring `M ≥ 14/183` in the balance is exactly (worst case `v_f = 182`):

> **`M_core ≥ (182 + s)/2379`**   (`2379 = 13·183 = 13·Φ₆`;  `= 1/13 + (s−1)/2379`),

a **joint rigidity coupling the core value to its binding speed** — a quantitative refinement of LRC(13).

**(1) Tight case `s = 1` — RIGOROUS.** `req = 1/13`, and `M_core ≥ 1/13` by LRC(13), so `M ≥ (1/13)(182/183) =
14/183`. Equality iff `M_core = 1/13` iff the core is the interval `{1..12}` — the **unique** `1/13`-minimizer
because `n = 13` is **prime** (S252: prime ⟹ one tight pattern, no doubling). So among `s=1` families the
**deep well `{1..12,182}` is the unique minimizer, exactly `14/183`.** Rigorous.

**(2) `s ≥ 2` — VERIFIED with margin, tight only at `s=1`.** The core is non-interval, so `M_core > 1/13`
strictly, and the requirement holds with room:
- near-interval covering `{C,182}`: **all** `M ≥ 14/183` (min `= 14/183`, the deep well **alone**);
- the rigidity `M_core ≥ (182+s)/2379`: **zero violations** across 11000+ families, binding speed `s` up to
  **63**; tightest large-`s` margin `M_core − req = 238/159393 = 0.0015` at `s = 63`.

It is a **genuine coupling**, not a single-rung fact: the 12-core spectrum gap `(1/13, 1/12)` is **non-empty**
(`56/18743` cores land strictly inside), so "non-interval ⟹ `M_core ≥ 1/12`" is *false*. Rather, cores with
**small `M_core` are forced to have small `s`**, and large-`s` cores have large `M_core` — precisely what
`(182+s)/2379` captures, tight at the interval. (This is why the escape *feels* dangerous — `M_core` can sit
just above `1/13` — yet never fires: such cores have `s` small enough that `req` stays below them.)

## Net — and honest scope

The large-`s` escape is **closed for the single-killer-182 covering-min**: empirically (11000+ families, `s`
to 63, zero counterexamples, deep well the unique minimizer) and **rigorously in the tight case `s=1`**
(LRC(13) + prime-13 uniqueness). The whole `s ≥ 2` regime reduces to a **single sharp joint rigidity**
`M_core ≥ (182+s)/2379` — a clean, verified, LRC(13)-refining statement, tight only at the deep well — which
is the concrete remaining lemma. Beyond it: **killers `v_f ≠ 182`** clear *more easily* (`v_f/(v_f+s)` is
larger for bigger `v_f`, and covering forces `v_f ≥ 182`, so `182` is the worst case — this is handled by the
same rigidity with a weaker `req`); and **multi-killer covering families** (several resonant runners at once —
the multi-constraint balance) remain the genuinely open general case.

So S253's proof (interval core, `s=1`) now extends, rigorously, to the **entire tight `s=1` stratum** with the
deep well pinned as the unique minimizer, and the `s ≥ 2` escape is closed modulo one verified joint rigidity.
The covering-min lower bound is now: **prove `M_core ≥ (182+s)/2379`** (single-killer) **and** handle
multi-killer families — a sharply reduced target.

→ opus-S253 (the balance + escape), opus-S252 (prime-13 uniqueness of the tight core), klein S267 (14/183
covering-min, verified), LRC(≤13) citation. Files: `lrc14_large_s_escape_closure_opus_S254.py` (+`.out`).
