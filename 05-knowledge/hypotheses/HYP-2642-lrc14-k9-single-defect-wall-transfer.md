---
id: HYP-2642
title: LRC(14) k=9 single-defect wall-transfer certificate
status: PARTIAL proof scaffold; exact certificate for the binding non-AP row
source: codex-2026-06-19-S29
depends_on:
  - HYP-2638
  - HYP-2637
  - HYP-2607
  - THM-534
related:
  - HYP-2640
  - HYP-2639
  - HYP-2635
  - OPEN-Q-108
---

# HYP-2642 - LRC(14) k=9 Single-Defect Wall Transfer

> Namespace note: this target was first pushed as `HYP-2641/T889`, but KPS S13
> independently claimed `HYP-2641` for the compatible far-element plateau
> recursion.  The wall-transfer target is renumbered here as `HYP-2642/T890`.

## Claim Being Tested

After the KPS S12 correction, the binding non-AP row is no longer a wide GAP
or a mod-27 stratum phenomenon.  It is the bounded k=9 single-end-defect row

```text
A = (0,1,2,3,4,5,6,7,8)
D = (0,1,2,3,4,5,6,7,9).
```

HYP-2638 already proves by exact finite enumeration that `D` is below the
cap.  This note asks for a proof-shaped explanation of the margin: replacing
runner `8` by `9` should force enough common-wall mass to move from higher
`g_9(N)` values to lower ones that the AP-to-defect loss exceeds the tiny
AP-to-cap slack.

## Computation

Script:

- `04-computation/lrc14_k9_single_defect_wall_transfer_codex_s29.py`
- output: `05-knowledge/results/lrc14_k9_single_defect_wall_transfer_codex_s29.out`

The script uses exact rational breakpoint integration.  It refines the unit
interval by all walls of both `A` and `D`, records the missed sector sets on
each common interval, and aggregates the transition

```text
N_A(x) -> N_D(x).
```

For k=9,

```text
g_9(t) = -(t-2)(t-3)(t-6)/36
       = [1, 5/18, 0, 0, 1/9, 1/6, 0].
```

## Exact Wall Certificate

The exact distributions are:

```text
t   p_A                         p_D                         p_D-p_A
0   2447/5880                   1181/2940                  -17/1176
1   653/2940                    23/90                       59/1764
2   27/196                      1139/8820                  -19/2205
3   23/245                      38/441                     -17/2205
4   13/196                      317/4410                    1/180
5   9/196                       5/126                      -11/1764
6   1/56                        1/63                       -1/504
```

Therefore

```text
L_y(A)       = 26083/52920      = 0.492876039
L_y(D)       = 38681/79380      = 0.487288990
L_y(A)-L_y(D)= 887/158760       = 0.005587050
cap_9-L_y(D) = 39541/5675670    = 0.006966755
cap_9-L_y(A) = 10441/7567560    = 0.001379705
```

On the common-wall refinement, the weighted transfers split as:

```text
weighted positive transfers = 9749/158760  = 0.061407155
weighted negative transfers = 2659/39690   = 0.066994205
negative - positive         = 887/158760   = 0.005587050
zero-weight transfer mass   = 4393/5880    = 0.747108844
```

So the exact AP-to-defect loss is a wall-transfer surplus.  A proof can aim to
pair every positive transfer against negative transfer mass and retain at
least `887/158760`; retaining only the smaller surplus
`cap_9-L_y(A)=10441/7567560` already proves `D < cap_9`.

The largest negative moves are the clean AP-cover losses

```text
() -> (6,)  mass=1/45
() -> (1,)  mass=31/1470
() -> (5,)  mass=1/56
```

each with `dg=-13/18`.  The largest compensating gains are smaller or have
less favorable weights.

## One-Gap Envelope

For one-gap near-APs

```text
E(L,s) = {0,...,L-1} union {L-1+s,...,L-1+s+8-L},  s>=2,
```

the S29 exact scan through `s<=30` finds:

1. for each tested `s`, the best primitive row is the endpoint defect `L=8`;
2. among those endpoint defects, the global best is `s=2`, namely `D`;
3. larger endpoint gaps oscillate by residue, so naive monotonicity in `s` is
   false, but none approach the `s=2` value.

This extends the earlier KPS `lrc14_single_defect_family.py` scan, which found
the same `s=2` endpoint maximizer for k=8,9,10 through `s<=60`.

## Endpoint Gap Asymptotic

The endpoint family has an exact independent-extra-runner limit.  Write

```text
F_s = (0,1,2,3,4,5,6,7,7+s).
```

For the base row `B=(0,1,2,3,4,5,6,7)`, the missed-count distribution is

```text
[481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49].
```

If the last runner is treated as equidistributed relative to the base row, the
limit is

```text
lim_{s->infty} L_y(F_s) = 20698/46305 = 0.446992765.
```

The binding defect sits well above this limit:

```text
L_y(F_2) - 20698/46305 = 22391/555660 = 0.040296224.
```

The scan through `s<=120` has empirical

```text
max |(7+s) * (L_y(F_s)-20698/46305)| = 2143/3087 = 0.694201490.
```

Therefore a clean discrepancy lemma of the form

```text
|L_y(F_s)-20698/46305| <= 1/(7+s)
```

would reduce endpoint minimal-gap dominance to the finite check `s<=17`.  A
stronger `2/(3(7+s))` bound would leave only `s<=9`.  This converts the
oscillating residue sequence into a standard rotation-discrepancy target.

## Proposed Lemma Split

To turn the finite near-AP check into a theorem-shaped argument, prove:

1. **Endpoint dominance.** For k=9 one-gap rows, `L_y(E(L,s))` is maximized at
   the endpoint `L=8` for every `s>=2`.
2. **Minimal-gap dominance.** For endpoint rows
   `F_s=(0,1,2,3,4,5,6,7,7+s)`, `L_y(F_s) <= L_y(F_2)` for every `s>=2`.
   Because the sequence in `s` oscillates, the proof should use wall-transfer
   pairing, residue-class envelopes, or the endpoint discrepancy limit above
   rather than monotone descent.
3. **AP-to-defect wall surplus.** On the common-wall refinement for `A` and
   `F_2`, pair positive transfers against negative transfers and keep surplus
   at least `10441/7567560`.  The exact surplus is larger:
   `887/158760`.

Together with HYP-2638's Freiman `3k-4` finite pocket, these lemmas would give
a structural explanation of the tight bounded non-AP case.  They do not by
themselves prove LRC(14); the wide-spread signed tail / large-doubling side
still belongs to HYP-2639, HYP-2636, and HYP-2633.

## Assumption Challenge

The tested tournament vertices were runners, gaps, fixed sections, wall
crossing events, residue classes, relation vectors, and proof obligations.
The chosen quotient is the common-wall transfer of missed-sector counts.  It
preserves the exact `L_y` difference and the cap margin while discarding
runner labels inside intervals.  This challenges the assumption, refuted again
by HYP-2640, that saturated raw relation rank should linearly predict the
correction.
