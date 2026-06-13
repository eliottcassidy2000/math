---
id: THM-498
title: The Pollock-completeness hierarchy of tournament cycle spectra — c3 complete, c5 first-gap at (n=6, value 10), H gappy
status: PARTIAL — the c5 first-gap is PROVED by exhaustive computation + an independent second method (0 disagreements); the hierarchy framing and the Pollock-method reading are the conceptual content
source: kind-pasteur-2026-06-13-S3
depends_on:
  - THM-462   # c3 spectrum gap-free [0,M(n)] (the Lagrange-four-square / Pollock-style engine)
  - THM-029   # H-forbidden value 7
  - THM-079   # H-forbidden value 21
related:
  - goldbach-polygonal-zeckendorf-additive-bases-s501   # the three-currencies additive-coverage frame
---

# THM-498 — the Pollock-completeness hierarchy of tournament cycle spectra

## The Pollock lens

Pollock's conjectures (1850) ask whether a fixed family of figurate numbers
(tetrahedral, octahedral) is an **additive basis of bounded order** — whether
*every* integer is representable. In the repo's additive-coverage frame
(reflection `goldbach-polygonal-zeckendorf-additive-bases-s501`) this is the
**bounded-arity currency**: sparse atoms, but a fixed number of summands large
enough to erase residue gaps. The tournament analogue of "every integer is
representable" is: **is an invariant's value-set GAP-FREE** (every value in
`[0,max]` realized by some tournament)? Call such a channel *Pollock-complete*.

THM-462 proved the cubic channel `c3` (directed 3-cycle count) is gap-free
`[0,M(n)]` at all `n` — and its engine is exactly a Waring/Pollock-style argument:
the score deviation satisfies `f(s) = m(n) + (Σ e_i²)/2`, and **Lagrange's
four-square theorem** makes `Σ e_i²` hit every value in the top window, with
induction covering the rest. So `c3` is Pollock-complete, *proved by the
prototype Pollock theorem (four squares)*. At the other end, the degree-`n`
invariant `H = I(Ω,2)` is **incomplete**: it has the forbidden values 7 and 21
(THM-029/THM-079). The question THM-498 answers: **where, in the cycle-length
hierarchy, does Pollock-completeness break?**

## The result (c5 — verified)

> **THM-498.** Let `c5(T)` = the number of directed 5-cycles of a tournament `T`.
> The `c5` value-set over `n`-vertex tournaments is **gap-free at `n = 5`
> (`[0,3]`) but FIRST GAPPY at `n = 6`: the spectrum is `[0,12] ∖ {10}`** — the
> value **10 is forbidden** while 9, 11, 12 are all realized.

*Status:* exhaustive over all `2^{15}` labeled tournaments at `n = 6`, confirmed by
TWO independent directed-5-cycle counters (per-5-subset Hamiltonian enumeration;
and distinct-ordered-5-tuples / 5 rotations) with **0 disagreements**; `c3`
recomputed gap-free `[0,M(n)]` for `n ≤ 6` as a pipeline validation against
THM-462. (`04-computation/cycle_spectrum_pollock_lens_kps3.py`,
`05-knowledge/results/...out`.) `c5` is **not** score-determined (unlike `c3`),
so this is a genuinely new spectrum, not reducible to Landau sequences.
*Honest scope:* `c5=10` is the first forbidden value; the `n=7` `c5` spectrum
(`[0,42]`-ish) needs an exhaustive `2^{21}` sweep to pin its gaps (sampling shows
candidate top-end gaps but those are likely sampling artifacts — flagged).

## The hierarchy (the Pollock-onset climbs the cycle length)

| channel | degree | Pollock-complete? | first forbidden value |
|---|---|---|---|
| `c3` | 3 | YES, all `n` (THM-462, four-square engine) | none |
| `c5` | 5 | NO | **10** (at `n=6`) — THM-498 |
| `H = I(Ω,2)` | `n` | NO | 7, then 21 (THM-029/079) |

So additive completeness of the cycle-count channels **degrades with cycle
length**: the more vertices a cycle uses (the higher the OCF degree), the sooner
representability fails. `c3` is the unique channel reducible to a sum-of-squares
(hence four-square-complete); from `c5` on, the count escapes the score sequence
and gaps appear. This locates the "Pollock-onset" of the OCF hierarchy between
degrees 3 and 5.

## Why this is the Pollock method, transferred

- The **proof that c3 is complete uses the literal prototype Pollock/Waring
  theorem** (Lagrange four squares) as the arity engine. A general Pollock proof
  method (representing every integer by `k` values of a fixed polynomial) is the
  template for proving completeness of any score-determined channel; THM-462 is
  the `Σe_i²` instance.
- The **forbidden values are the Pollock *exceptions*** — the finite residue
  obstructions a circle-method proof isolates into a "singular series vanishes /
  finite exceptional set." `c5=10`, `H∈{7,21}` are exactly such exceptions; their
  classification is a singular-series / local-obstruction computation.
- Open program (HYP-2487): characterize the `c5` forbidden set at all `n` (does
  10 stay forbidden? do more appear?) and the `c_{2k+1}` onset as a function of
  `k` — the Pollock-completeness curve of the OCF.

**Artifacts:** `04-computation/cycle_spectrum_pollock_lens_kps3.py` (+ `.out`).
Reflection `07-reflections/pollock-as-the-bounded-arity-currency-and-the-cycle-spectrum-onset-kps3.md`.
