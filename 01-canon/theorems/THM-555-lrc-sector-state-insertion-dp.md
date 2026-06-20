---
id: THM-555
title: LRC sector-state insertion is an exact one-sector deletion dynamic program
status: PROVED; verified by assertions in `lrc14_state_transfer_dp_codex_s58.py`
source: codex-2026-06-20-S58
depends_on:
  - THM-551
related:
  - THM-548
  - HYP-2675
  - HYP-2683
  - HYP-2684
  - HYP-2690
  - HYP-2691
---

# THM-555 - LRC Sector-State Insertion DP

## Statement

Work in the LRC(14) seven-sector model, and let the six inner sectors be
`1,...,6`.  For a finite speed row `P`, write

```text
M_P(t) = {inner sectors not occupied by any speed in P at time t}.
```

For a new speed `e`, let `s_e(t)` be the sector occupied by `e` at time `t`.
Then pointwise

```text
M_{P union {e}}(t) =
    M_P(t) \ {s_e(t)}    if s_e(t) is an inner sector in M_P(t),
    M_P(t)               otherwise.
```

Consequently adding one runner can remove at most one missed inner sector on
each wall atom.  The transfer kernel on missed-sector states is lower
triangular in inclusion:

```text
K_e(A,B) = meas{t: M_P(t)=A and M_{P union {e}}(t)=B}
```

can be nonzero only when `B=A` or `B=A\{s}` for one sector `s`.

In particular, for the all-covered sector mass

```text
p0(P) = meas{t: M_P(t)=empty},
```

the insertion increment is exactly

```text
p0(P union {e}) - p0(P)
  = sum_{s=1}^6 meas{t: M_P(t)={s} and s_e(t)=s}.
```

Equivalently, the only atoms that newly become all-covered are atoms where the
prefix missed exactly one inner sector and the new speed lands in that sector.

## Proof

At a fixed time `t`, the new speed occupies exactly one of the seven sectors.
If it occupies sector `0`, it does not affect the inner-sector missed set.  If
it occupies an inner sector already occupied by a prefix speed, it also does
not affect the missed set.  The only possible change is when `e` occupies an
inner sector that was previously missed; then that one sector is deleted from
`M_P(t)`.  No other missed sector can change because one speed has one sector
position at time `t`.

This proves the pointwise state formula and the lower-triangular support of
`K_e`.  The formula for `p0` follows by taking the measure of the transition
into the empty state and subtracting the atoms that were already empty before
the insertion.  Those newly empty atoms are exactly the atoms with
`M_P(t)={s}` and `s_e(t)=s` for some `s in {1,...,6}`.

## Computational Check

`04-computation/lrc14_state_transfer_dp_codex_s58.py` computes the transfer
kernel on the exact common wall refinement.  Its internal assertions verify
for every tested prefix and insertion that:

```text
M_after subset M_before,
|M_before|-|M_after| in {0,1},
sum closure_mass = p0(after)-p0(before).
```

The stored output
`05-knowledge/results/lrc14_state_transfer_dp_codex_s58.out` applies the
identity to AP, boundary, dyadic, true-wide, and AP-triple phase rows.
