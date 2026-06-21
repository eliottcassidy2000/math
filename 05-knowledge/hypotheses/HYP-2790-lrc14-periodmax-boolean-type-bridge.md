---
id: HYP-2790
title: THM-563 bounded-base period maxima are not explained by scalar Boolean/type slack; endpoint and q0 coordinates survive
status: TESTED/PARTIALLY-TRUE; scalar bridge weak, Phi_low base-transfer refuted at k=8
source: codex-2026-06-21
depends_on:
  - THM-563
  - HYP-2788
  - HYP-2791
  - HYP-2752
  - HYP-2789
related:
  - HYP-2780
  - HYP-2781
  - HYP-2792
  - HYP-2770
  - HYP-2744
  - THM-534
  - OPEN-Q-108
---

# HYP-2790: Period-Max / Boolean-Type Bridge

## Claim Under Test

`THM-563` turns the single-far signed deviation into a finite endpoint-period
maximum.  `HYP-2788` says the wide near-cap branch routes to
single-perturbation bounded scaffolds.  The remaining finite object is therefore
a list of bounded bases `B` and the inequality

```text
periodmax(B) < 15 * (cap_k - Plat(B)).
```

This hypothesis tests whether the same Boolean/type slack discovered in
`HYP-2791` and the containment-cut refinements of `HYP-2752` explain why the
dangerous bounded bases satisfy the period-max inequality.

The first S76 coordinate overlay gives a correction: the three-term
`Phi_low` cut does **not** transfer naively to the one-far bounded-base ledger.
At base size `7` (final row `k=8`) some non-AP frontier bases have negative
`Phi_low-AP` gap.  The safer bridge is:

```text
period-max direct scan / THM-563
  -> AP and dilation orbit filter
  -> q0 cover-atom slack on bounded bases
  -> Phi_low only on final-row Boolean laws, or on a size-shifted k>=9 ledger
```

## Planned Finite Ledger

For each bounded base `B` with `0 in B`, `max(B)<=14`, and final row size
`k=|B|+1`, compute:

```text
Plat(B) = p0(B) + p1(B)/7,
margin(B) = cap_k - Plat(B),
periodmax(B) = max_w w*Delta_w(B),
ratio(B) = periodmax(B) / margin(B),
```

then compare `ratio(B)` with Boolean/type features of the missed-sector law of
`B`, especially:

```text
T1=(1,(1)), T2sep=(2,(1,1)), T2adj=(2,(2)),
Phi_low = 21*T1 + 57*T2sep + 2*T2adj.
```

The first goal is not to prove the inequality.  It is to determine whether the
top period-max pressure rows form a small family already visible to the
Boolean/type quotient, sorted-cell leak quotient, or affine gap-word quotient.

## Result: Scalar Bridge Refuted, Endpoint Coordinate Promoted

The exact scout
`04-computation/lrc14_periodmax_type_bridge_codex_20260621.py` uses the
`THM-563` endpoint identity but evaluates the period scan by integer numerators
on the common `7*period` grid, since the centering constants cancel in every
arc difference.  It uses the canonical cap table

```text
cap_8..13 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7, 1.
```

The checked high-plateau bounded bases gave:

```text
k=8: 36 checked, 0 failures, worst ratio 10.818767 at (0,2,4,6,8,10,12)
k=9: 32 checked, 0 failures, worst ratio 13.280470 at (0,2,4,6,8,10,12,14)
k=10: 30 checked, 0 failures, worst ratio 9.780630 at (0,2,4,5,6,8,10,12,14)
k=11: 22 checked, 0 failures, worst ratio 8.540750 at (0,2,4,5,6,7,8,10,12,14)
k=12: 15 checked, 0 failures, worst ratio 5.882825 at (0,2,3,4,6,7,8,9,10,12,14)
```

Total: `135` checked rows, `0` failures of
`periodmax(B)<15*(cap_k-Plat(B))`.  The global checked pressure is the k=9 AP
dilation, where Boolean/type and containment deficits are exactly zero.  This
is correct behavior for a dilation-invariant quotient, but it also means these
quotients cannot by themselves explain the period numerator.

For non-AP rows, the diagnostic correlations point away from the
Boolean/type scalar:

```text
k=8:  corr(ratio,arc_count)=+0.6556, corr(ratio,KPS3def)=-0.3512
k=9:  corr(ratio,arc_count)=+0.7331, corr(ratio,KPS3def)=-0.0576
k=10: corr(ratio,arc_count)=+0.6924, corr(ratio,KPS3def)=+0.4088
k=11: corr(ratio,arc_count)=+0.4971, corr(ratio,KPS3def)=+0.2700
k=12: corr(ratio,arc_count)=+0.7114, corr(ratio,KPS3def)=+0.4963
```

The signs and strengths are not stable enough for a finite Boolean/type slack
certificate.  The surviving coordinate is the endpoint-period numerator itself:
arc count, endpoint period, and the signed sawtooth/Dedekind sum in `HYP-2792`.

## Proof Consequence

This does not close the single-far bounded-base finite check.  It narrows the
right proof object.  A plausible proof should not try to force
`periodmax(B)` through the HYP-2791/HYP-2752 Boolean/type cut.  Instead, prove a
uniform signed endpoint-period estimate, likely a generalized Dedekind
reciprocity bound of the form

```text
periodmax(B) <= C(B)
```

where `C(B)` is controlled by the signed endpoint orbit, arc pairing, or
reciprocity data and then compare `C(B)` with `15*(cap_k-Plat(B))`.  The
Boolean/type quotient remains useful for bounded plateau/cap certification; it
is not the natural carrier for the single-far oscillatory numerator.

## S76 Coordinate Overlay

Script:

```text
04-computation/lrc14_periodmax_boolean_bridge_codex_s76.py
```

Stored output:

```text
05-knowledge/results/lrc14_periodmax_boolean_bridge_codex_s76.out
```

The script is deliberately a coordinate overlay.  Exact period maxima are
delegated to the incoming S6/S7 period-max scripts, especially
`lrc_periodmax_general_macmini_0621s6.py`,
`lrc_periodmax_dangerous_scan_macmini_0621s7.py`, and
`lrc_periodmax_skipped_audit_thread5.py`.  The overlay attaches `Plat(B)`,
canonical and strict margins, endpoint period `P=7*lcm(B)`, AP/dilation type,
`q0`, and `Phi_low` coordinates to the published frontier rows.

Findings:

```text
global non-AP frontier min Phi_low gap = -4153/3080
global non-AP frontier min q0 gap      = 71/5880
k=8  min Phi_low gap=-4153/3080; min q0 gap=1/40
k=9  min Phi_low gap=23/12;      min q0 gap=3/56
k=10 min Phi_low gap=10867/5880; min q0 gap=71/5880
```

Thus the naive statement

```text
non-AP bounded base => Phi_low(B) > Phi_low(AP_base)
```

is false at the `k=8` base level.  Concrete negative witnesses in the S6/S21
frontier include:

```text
B=(0,1,2,4,5,6,10):       Phi_low gap = -283/420
B=(0,2,3,5,6,8,11):       Phi_low gap = -4153/3080
B=(0,1,2,10,11,12,13):    Phi_low gap = -33499/60060
B=(0,1,3,5,9,11,13):      Phi_low gap = -7669/16380
```

All non-AP frontier rows in the overlay have positive `q0` gap, including the
S6 broad-scan worst checked row

```text
B=(0,4,6,8,10,12,14): Phi_low gap=493/294, q0 gap=5/49.
```

The updated proof interpretation is that `q0`/cover slack is a better
base-ledger coordinate, while HYP-2791's `Phi_low` cut remains a final-row
Boolean law or possibly a size-shifted `k>=9` subledger.  This is a useful
negative result: it prevents using the finite k=8 final-row certificate at the
wrong projection level.

The overlay also records the cap-normalization guardrail: the canonical final
LRC cap has `cap_10=55/91`, while older period-max scripts sometimes report a
strict comparison against the proved floor `4/7`.  The latter is harder and
should not be silently conflated with the final cap ledger.

## Why This Is Not A Duplicate

`HYP-2791` is a full k=8 bounded-row atom-run-type certificate.  `HYP-2752`
shows fixed-k low-level Mobius/containment cuts exist but are non-uniform and
not structurally free.  `HYP-2790` asks a different question: after the wide
branch has been reduced to single-far bounded scaffolds, can the finite
period-max pressure itself be ranked or certified by the Boolean/type slack?

If yes, this gives a route from:

```text
wide near-cap -> single-perturbation bounded scaffold
              -> THM-563 endpoint-period max
              -> Boolean/type finite ledger
```

instead of a separate ad hoc period scan for every bounded base.

## Assumption Challenge

The vertices for this test are bounded bases/scaffolds, not final rows or
runner arcs.  The quotient preserves missed-sector cyclic shape and endpoint
period pressure; it destroys speed ownership inside wall atoms and the full
Fourier phase distribution.  If the test fails, the likely missing coordinate
is the endpoint-period numerator itself, not another scalar energy statistic.

This is now exactly what the scout indicates.  Alternate vertices considered:
bounded bases, endpoint arcs, missing-sector type atoms, containment atoms,
sorted cell widths, and explanatory-lens tournament nodes.  The quotient
preserved enough to recognize AP/dilation extremals but destroyed the phase
ordering of endpoints under multiplication by `w`.
