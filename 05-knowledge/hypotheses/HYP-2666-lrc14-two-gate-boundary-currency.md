---
id: HYP-2666
title: LRC14 two-gate boundary currency - shell-1 gate plus refined p1 tax
status: OPEN; exact bounded-bank evidence, raw tax constant corrected by HYP-2667
source: codex-2026-06-20-S40
tangent: T908
depends_on:
  - HYP-2661
  - HYP-2662
  - HYP-2663
  - HYP-2664
  - HYP-2665
related:
  - HYP-2667
  - HYP-2648
  - HYP-2655
  - HYP-2658
  - HYP-2660
  - OPEN-Q-108
---

# HYP-2666 - LRC14 Two-Gate Boundary Currency

## Claim

The two active S39 HYP-2664 routes, refined by HYP-2665, are not separate
proof obligations.  They are two gauges of one boundary currency:

1. First apply the HYP-2661 shell-1 gate.  If a row damages the tower
   `{1,2,4,8}`, it should be discharged by tower-deletion or mouth-ownership
   rigidity before any far residual estimate is needed.
2. On the remaining shell-1-full or genuinely nonlocal rows, pay the far
   HYP-2662/HYP-2664 endpoint imbalance with the single-missed-sector mass
   `p1(E')`, using the HYP-2665/HYP-2667 correction that the raw `p1/3` and
   `3p1/8` constants are false; `2p1/5` or a packet split is the viable target.

The scalar target is therefore not a naked residual bound, but the ordered
certificate

```text
p0(E') + (1/7 + c) p1(E') <= cap_9,
cap_9 = 1979/4004,
```

after quotienting first by the shell-1 packet.

**S41 correction:** HYP-2667 shows that the cap bank can have room for `c=3/8`
while the actual far discrepancy still violates `Delta_w^+ <= 3*p1(E')/8`.
The worst full `B=13` row,
`E'=(0,1,2,4,6,7,8,10), w=12`, is shell-1-full and has
`Delta_w^+/p1=997/2562 > 3/8`.  Thus the two-gate proof order remains useful,
but the raw far-tax target must be `2p1/5` or a generic/dyadic packet split.

## Exact Bounded Evidence

`04-computation/lrc14_two_gate_boundary_currency_codex_s40.py` scans

```text
E' = {0} union A, |A|=7, A subset [1,B],
```

with primitive `E'`, then stratifies each row by the missing shell-1 packet
`{1,2,4,8} \ E'`.  The S40 version uses exact integer common-wall coordinates
for the missed-sector distribution, so the `[1,18]` bank is now practical.

Stored run:

```text
B=18
primitive rows=31788
c in {1/4, 13/51, 1/3, 3/8}
c=3/8 violations=0
```

Global maximum at `c=3/8`:

```text
E'=(0,1,2,3,4,5,6,7)
shell1_missing=(8,)
value=12449/27440
slack=159213/3923920
minimum allowed c=388929/718718 ~= 0.541143
```

Shell-1-full stratum:

```text
rows=364
max row=(0,1,2,3,4,5,6,8)
value=6389/15680
slack=194613/2242240
minimum allowed c=194226/323323 ~= 0.600718
```

The widened bank keeps the same qualitative shape as the earlier `[1,15]`
fallback but adds a stronger check: rows newly admitted by `[16,18]` can enter
the top thirty, especially all-damaged `3`-multiple packets, yet none beat the
AP8 boundary row and none threaten `c=3/8`.

At the older `c=1/3` checkpoint, the same max rows have slack
`448069/8828820` globally and `51871/504504` in the shell-1-full stratum.
HYP-2665 shows that `c=1/3` is not enough for the actual far discrepancy, and
HYP-2667 shows that `3/8` is not enough either.  This cap bank still has room
for `3/8`; that is a cap-slack fact, not a proof of the actual raw tax.

## Interpretation

The global top row is AP8 missing the high tower bit `8`.  That is exactly the
kind of row HYP-2661 says should be handled by shell-1/tower deletion before
the far-residual p1 tax is invoked.  After imposing shell-1-full, the margin
gets much larger:

```text
all rows slack at c=3/8:        159213/3923920 ~= 0.04058
shell-1-full slack at c=3/8:    194613/2242240 ~= 0.08679
```

So the proof order looks like:

```text
shell1_gate > p1_boundary_tax > missing_packet > cap_slack > raw_value
```

The corresponding tournament vertices are proof obligations, not runners:
`shell1_gate`, `p1_boundary_tax`, `missing_packet`, `cap_slack`, and
`raw_value`.  The observable is preservation of the cap inequality under the
switch that splits rows by shell-1 missing packet before scalarizing `p1`.

## Next Proof Target

Prove the ordered lemma:

```text
Either E' damages shell 1 and pays the HYP-2661 tower-deletion/mouth tax,
or Delta_w^+(E') <= 2*p1(E')/5
```

with the residual interpreted through the HYP-2662 endpoint-weight/Galois
phase decomposition and HYP-2665/HYP-2667's phase-packet cancellation warning.
A more precise alternative is to prove generic phase packets below `3p1/8` and
classify the dyadic-even frontier below `2p1/5`, while keeping the bounded cap
bank charged with its exact slack ledger.

## Guardrails

This is not a proof of LRC14.  It is a proof-order reduction supported by an
exact bounded bank.  The computation also warns against scalarizing too early:
some all-damaged `3`-multiple packets have high `p1`-tax values in the wider
bank, but their address is different.  They should route through the shell
gate or a `3`-ramified/fibered packet lemma before being compared as raw
numbers.  It also absorbs HYP-2665's correction: the phrase "p1 tax" does not
mean a universal `p1/3` endpoint bound.

## Artifacts

- `04-computation/lrc14_two_gate_boundary_currency_codex_s40.py`
- `05-knowledge/results/lrc14_two_gate_boundary_currency_codex_s40.out`
- `07-reflections/lrc14-two-gate-boundary-currency-codex-s40.md`
