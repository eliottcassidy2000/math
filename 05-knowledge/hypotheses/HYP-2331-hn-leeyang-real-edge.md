---
id: HYP-2331
title: Hadwiger–Nelson as a Lee–Yang edge — χ(plane) = 1 + rightmost REAL chromatic zero of UD graphs; the integrality gap is the real-edge/complex-bulk lag
status: OPEN (reframe; NO elimination of {5,6,7}); real-edge=χ−1 + complex bulk VERIFIED on gadgets
source: claudebox-2026-06-06-S640
related:
  - chromatic-number-of-the-plane-narrowing-attempts-and-the-combinatorial-wall-s699  # the Vitali wall χ_f≤4.36<5
  - hn-moser-spindle-field-tower-fta-s687  # field tower / FTA root-locus (geometry)
  - HYP-2326  # FTA: coefficients↔roots; HYP-2316 glass transition (Lee-Yang)
  - HYP-2295  # chromatic polynomial = zero-T antiferromagnetic Potts
---

# HYP-2331 — the Lee–Yang face of Hadwiger–Nelson

A complement to S699 (Vitali wall) and S687 (FTA on the unit-distance *geometry*): the Lee–Yang lens
on the chromatic *polynomial's* zeros (the `q`-plane). **No elimination of {5,6,7}** — a reframe that
locates the problem precisely in the zero-locus language and ties it to the fleet's integrality gap.

## The reframe

`P(G,q)` is the zero-temperature antiferromagnetic Potts partition function; its zeros are the
**chromatic zeros = Lee–Yang zeros in the state-number `q`**. `χ(G) = ` smallest integer `q` with
`P(G,q)>0`, so:

> **`χ(plane) = 1 + (rightmost REAL chromatic-zero accumulation point of finite unit-distance graphs)`**
> — a real Lee–Yang edge in the `q`-plane. de Grey ⟹ the edge `≥ 4` (so `χ≥5`); `{5,6,7}` = whether
> the real edge sits at `4, 5,` or `6`.

**Verified on gadgets** (`hn_leeyang_s640.out`): the rightmost *real* chromatic zero is exactly
`χ−1` — `W6` (Eisenstein floor, χ=3) → real edge `2`; Moser spindle (χ=4) → real edge `3`. The
**complex** zeros form a *bulk extending further right* (W6: real part up to `2.81 > 2`).

## The integrality gap is the real-edge/complex-bulk lag = the Vitali wall

S699's rigorous result: `χ_f(plane) = 1/m₁ ≤ 4.36 < 5`, so no fractional/spectral/measure argument
reaches `χ≥5` — narrowing `{5,6,7}` is irreducibly combinatorial (the LRC Vitali wall / integrality
gap). In Lee–Yang language:

> The **complex** chromatic-zero bulk reaches the *analytic* value (≈ the fractional/measure bound
> `χ_f ≤ 4.36`); the **real** Lee–Yang edge is `χ−1 ≥ 4`. The integrality gap `χ − χ_f` is the
> **lag of the real edge behind the complex accumulation** — the real zeros sit to the *left* of where
> the bulk reaches. Analytic methods see the complex bulk (≤ 4.36); the integer `χ` is the *real*
> edge, which only a combinatorial gadget (de Grey) can push to 4, and to 5 (=`χ≥6`) if ever.

So the "combinatorial wall" (S699) is exactly: **you cannot move the real Lee–Yang edge by analytic
means; the complex bulk is pinned below 5.** This is the glass-transition picture (HYP-2316: real-axis
pinch of partition-function zeros) and the FTA coeff↔root duality (HYP-2326) applied to HN.

## Why this is honest and where it sits
- It does **not** eliminate 5/6/7 (that needs a 6-coloring, a 6-chromatic gadget, or a 5-coloring).
- It is a *clean restatement* of S699's wall in zero-locus terms, with the gadget zeros computed,
  unifying HN with the LRC Lee–Yang/glass picture (the real-axis-pinch = collapse/transition motif).
- It sharpens the target: **find a UD graph with a real chromatic zero in `[5,6)` (i.e. `P(G,5)=0`)
  ⟺ `χ≥6`** — the chromatic-root version of "kill 5".

## To do
1. Track the real chromatic-zero edge (`= χ−1`) and the complex-bulk right-edge across a growing UD
   family (Moser ring / field-tower gadgets, S687); does the complex edge climb toward `4.36` while
   the real edge stays integer-stepped? (would visualize the integrality gap as a zero-locus picture.)
2. Is the complex-bulk right-edge provably `≤ χ_f` (or related)? If so the Vitali wall becomes a
   Lee–Yang statement: complex bulk `≤ 4.36`, real edge `= χ−1`, gap = integrality.
3. The Loeschian/field-tower (S699/S687) in zero-locus terms: `W6` zeros are `2 + (5th roots of −1)`
   — the cyclotomic floor; the `√−11` (spindle, disc −11) zeros are the post-cyclotomic edge.
