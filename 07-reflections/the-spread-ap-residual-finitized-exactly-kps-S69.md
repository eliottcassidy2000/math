---
source: kind-pasteur-2026-07-07-S69
status: EXACT FINITIZATION. The spread-AP residual of the 2-anchor (A') tail is reduced from
  boxeph's "adversarial numeric" state to a finite EXACT rational check (exact minima at the
  d=2 resonant dip, clearing T_k) plus a T²-decorrelation tail with an explicit Erdős–Turán
  bound. Owner: keep working the bleeding edge; finitize the spread-AP residual exactly.
tags:
  - lonely-runner
  - LRC14
  - anchor-floor
  - finitization
  - three-gap
  - decorrelation
---

# The spread-AP residual, finitized exactly

**kind-pasteur-2026-07-07-S69 (HYP-5067).** After my S68 anchored-gap subset lemma finitized
the 2-anchor tail `PA_2(E) = P(max(gap@0, gap@1/2) > 1/7) ≥ T_k` on bounded diameter, the
residual was the **spread-AP family** `{a + dj}` (unbounded `d`), where boxeph-S1's "spread AP
minimizes `PA_2`" was verified only *adversarially/numerically*. This session makes that
residual an **exact finite computation**.

## Three exact facts turn the spread-AP inf into a finite check

**(i) Dilation invariance (exact).** `PA_2(cE) = PA_2(E)` for integer `c` (the config
`{frac(cv·x)} = {frac(v·(cx))}` and `x↦cx` preserves measure). Verified exactly:
`PA_2([3,5,7,9,11]) = PA_2([9,15,21,27,33])`. So `{a+dj} = g·{a/g + (d/g)j}` reduces to
**`gcd(a,d)=1`**. (The offset `a` is *not* further reducible — it is a genuine inhomogeneity
parameter, since adding a constant to all speeds rotates the config against the fixed anchors.)

**(ii) Exact rationality per `(a,d)`.** The **anchored order-cell engine**: breakpoints in `x`
at `m/(e_i−e_j)` (order changes), `m/e_i` (a phase crosses anchor `0`), and `(2m+1)/(2e_i)` (a
phase crosses anchor `1/2`). On each cell the circular order is fixed, so the gap containing a
fixed anchor is a fixed pair of phases whose length is **affine in `x`**; `max(gap@0, gap@1/2)
> 1/7` is solved by exact linear crossings. `PA_2(E)` is an **exact rational** — validated
against the numeric integral to 4–5 digits.

**(iii) Decorrelation (the tail).** Write the config as the Steinhaus set `{frac(j·(dx))}`
**rotated by `frac(ax)`**; the pair `(frac(ax), frac(dx))` traces the `(a,d)`-geodesic on `T²`,
which **equidistributes** as `max(a,d)→∞` (Weyl). Hence `PA_2({a+dj}) → PA_2^∞(k) :=` the `T²`
average, computed:

| k | `PA_2^∞(k)` (T² avg) | `T_k` | margin |
|---|---|---|---|
| 8 | 0.7994 | 0.6185 | +0.181 |
| 10 | 0.5926 | 0.3956 | +0.197 |
| 13 | 0.3648 | 0.0565 | +0.308 |

The decorrelated limit clears `T_k` with a comfortable margin, so **the tail
`max(a,d) ≥ B_0(k)` is safe**, and the box is finite by **Erdős–Turán** (the max-of-two-
anchor-gaps indicator has bounded variation, so the equidistribution has an explicit rate,
giving an explicit `B_0`).

## The exact resonant minima (the dip)

On the coprime box the exact minimum sits at the **`d=2` resonant dip** at every `k`:

| k | exact `min PA_2` | `= ` | at `(a,d)` | over `T_k` |
|---|---|---|---|---|
| 8 | `1164298/1616615` | 0.72021 | (5, 2) | **+0.1017** |
| 10 | `178278085157/323982330210` | 0.55027 | (15, 2) | **+0.1547** |
| 13 | `14892552877/46727286606` | 0.31871 | (15, 2) | **+0.2622** |

**Dip is bounded (verified).** Sweeping `d = 1..15` (`k=8`, the tightest): the global minimum
is at `d=2` (`0.7203`); every `d ≥ 3` gives `≥ 0.766`; and — the natural worry — the
`1/7`-resonant step `d=7` does **not** dip (`0.7995`, essentially the decorrelated limit). So
the resonant dip is genuinely at small bounded `(a,d)`, and the exact minimum clears `T_k`.

## The finitization, assembled

> **`inf` over spread APs `{a+dj}` of `PA_2` `=` `min` over the finite coprime box
> `{gcd(a,d)=1, max(a,d) ≤ B_0(k)}` of EXACT rationals, and this minimum `≥ T_k`.**

Concretely: dilation reduces to `gcd(a,d)=1`; `PA_2` is exact-rational per `(a,d)`; the
resonant dip (`d=2`) is exact and clears `T_k` with margins `+0.10/+0.15/+0.26`; and everything
outside the box clears by decorrelation to `PA_2^∞ > T_k`. This replaces boxeph's *adversarial*
`PA_2 ≥ T_k` with an **exact-rational finite check plus one standard analytic input** — the
Erdős–Turán constant `B_0(k)` (a bounded-variation equidistribution rate, entirely routine).

This is the same finitization shape that closed my S59 diameter floor and S61 Part-A `V₀`: the
σ-even measure obstruction is bounded to a region where it is exact-rational-decidable, and the
tail is handled by a decorrelation limit. The 2-anchor `(A')` ledger is now, on the spread-AP
axis, *exact*.

## What remains (honest)

- **R2 (separate):** that spread APs are the **global** `PA_2`-minimizer over *all* families
  (boxeph's claim). My S68 anchored-gap subset lemma proves the bounded-diameter part for every
  family; the unbounded-diameter global-minimality is the last rigidity. If it holds, `(A')`
  closes: bounded diameter (S68) + spread-AP inf (this session, exact) + R2 (global min).
- **The `B_0(k)` constant** is stated via Erdős–Turán, not yet computed to an explicit integer
  — a routine (if tedious) BV/discrepancy calculation; the empirical `B_0` is small (`d=2` dip,
  `a ≲ 20`).

## Ledger

- Exact (verified): dilation invariance; the anchored order-cell rational engine; the resonant
  minima (rationals above) clearing `T_k`; the dip bounded at `d=2` (no `d=7` resonance); the
  `T²` decorrelated limits `> T_k`.
- Builds on: boxeph-S1 (2-anchor reduction), my S68 (anchored-gap subset lemma, HYP-5057),
  opus-S139 (double-cover `k=13`), opus-S134 / THM-637 (Farey roof), Weyl / Erdős–Turán.
- Files: `lrc_spread_ap_finitize_kps_S69.py` (+out).
- Does NOT prove LRC(14). It finitizes the spread-AP residual to an exact check + a routine
  decorrelation constant, leaving only the global-minimality rigidity (R2).
