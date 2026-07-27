---
id: THM-2570
title: "The word-depth collision law and the owner-clock host array: r = nu_13(deepest in-word blocker), and the first MISTAKE-281-clean 13x7 live-branch host"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT (two exact companions, each with
  double-route verification, exact gates, and python/-O agreement);
  independent hostile audit REQUESTED.  SCOPE: canonical TYPED row
  THM-2309 (25) only (not an asserted scalar cover); the host array W
  is NOT claimed to be THM-2512's d_A, no transplant theorem, no
  physical current, no row removed; LRC(14) OPEN.
source: opus-2026-07-28 (batons 2 and 3 of the THM-2560 broadcast)
depends_on:
  - THM-2471 (collision index (23)-(25); colours (8)-(12); sidecar (44)-(47))
  - THM-2449 (owner-clock partition, Sec 4)
  - THM-2560 (the r = 3 Stage-1 anchor, reproduced exactly)
related:
  - HYP-9032 (prediction P2 CONFIRMED: the C_7 re-homing is the owner clock)
  - THM-2512 (the bridge test now LIVE at sigma = {b} via the r-law)
  - MISTAKE-281 (common-base discipline, obeyed)
scripts:
  - 04-computation/lrc14_r_index_census_opus_20260728.py
  - 04-computation/lrc14_base_only_bridge_opus_20260728.py
outputs:
  - 05-knowledge/results/lrc14_r_index_census_opus_20260728.out
  - 05-knowledge/results/lrc14_base_only_bridge_opus_20260728.out
---

# THM-2570 -- the word decides the collision depth; the owner clock hosts

Row: `THM-2309 (25)`, `w = (1,14,27,40,53,66,13,2197,742586)`, owner
`j = 1` (the only lawful owner).

## (A) The word-depth collision law [r-census, 9 lawful configurations]

The THM-2471 collision index `r = min{s >= 1 : I_s > 0}` was computed
exactly at every lawful `(sigma, K)`, `sigma in {{a},{b},{a,b}}`,
`K in {2,3,4}` (all returns positive; empty word unlawful, recorded as
diagnostic only):

```text
r({a},   K) = 3   for K = 2,3,4     I_3({a},2) = 9926558757352/109707098520974955
r({b},   K) = 5   for K = 2,3,4     I_5({b},2) = 48602521488933856/337437093630814766589
r({a,b}, K) = 5   for K = 2,3,4     I_5({a,b},2) = 3919820186838/303177981698845253
```

**`r` depends ONLY on the word, never on the packet clock, and equals
the 13-adic depth of the deepest blocker appearing in the word:**
`r = nu_13(c_2) = 3` for the `a`-word, `r = nu_13(c_3) = 5` for every
word containing `b`. Consequence: THM-2560 (A)'s Stage-1 STOP was a
`sigma = {a}` artifact -- at six lawful configurations the collision
sits at `d = 13^5` where THM-2471 (47) supplies the unique affine
invariant `theta = t - 2u`: **THM-2512 Section 5's u-slaved bridge
test is LIVE**, canonically at `(j=1, sigma={b}, K=2)` (largest first
service, `I_5 ~ 1.44e-4`).

## (B) The owner-clock host array [base-only bridge, r = 3 packet]

On the common physical circle, pairing the `r = 3` collision-colour
data (`U = P_d f`, `V = P_d e`, `d = 13^3`) against the native
owner-clock partition (`ell in F_7`, the seven translates of
`d(c_1 y - ell/7)`, THM-2449 Sec 4) INSIDE one integral yields the
exact `13 x 7` array `W(k, ell)` (integer table over
`DENC = 2^4 3^3 5 7^2 11 13^14 53`), with exact gates:

- `sum_ell W(k, ell) = J(k)` for every `k` (reproduces the THM-2471
  colours; `J(0) = I_3` matches THM-2560's triple-verified value);
- `sum_k W(k, ell) = 0` for every `ell` (cellwise first-collision
  disjointness);
- Hermitian symmetry `W(13-k, ell) = conj W(k, ell)`;
- `W` genuinely `(k, ell)`-mixed: rank `>= 2`, `1548/1638` nonzero
  `2x2` minors; every one of the 7 cells fires all 12 nonzero colours.
  The only 90 vanishing minors are forced exactly by the reflection
  symmetry `C_ell(s) = C_{(7-ell) mod 7}(-s)` (negation-symmetric
  `E_1, Q`), which makes columns `ell = 3, 4` identical.

The centred array `Wc(k, ell) = W(k, ell) - J(k)/7` satisfies ALL
THREE transplant inputs of HYP-9032: **(T1)** row-zero over `F_7` AND
column-zero over `F_13`; **(T2)** the `k`-independence dichotomy holds
structurally and is non-vacuous (`91/91` entries nonzero); **(T3)** an
exact universal denominator floor. This is the first
MISTAKE-281-clean common-base object on a live-branch row carrying
both a collision-root `F_13` colour index and a native `F_7` label --
and the label is the OWNER CLOCK, exactly HYP-9032's prediction P2.

## What follows and what does not

The transplant program has its candidate physical host and its live
u-slaved bridge window, on one row. Open burdens: (i) consuming the
THM-2506/2507/2508-style conclusions from a `Q(zeta_13)`-valued
(rather than integral) array; (ii) running the Stage-2 slaved
contraction at `sigma = {b}` (the `r = 5` window); (iii) everything
uniform over the 165 rows. NOT implied: any transplant theorem,
physical current, row exclusion, or LRC(14) progress beyond the cited
statements.
