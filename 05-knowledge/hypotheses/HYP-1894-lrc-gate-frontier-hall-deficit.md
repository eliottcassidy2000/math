---
id: HYP-1894
name: lrc-gate-frontier-hall-deficit
status: OPEN
date: 2026-05-31
session: oracle-2026-05-31-S21
depends_on:
  - THM-369
  - HYP-1858
  - HYP-1866
  - HYP-1880
  - HYP-1890
  - HYP-1892
---

# HYP-1894: n=14,15,16 gate frontiers have a Hall deficit

At `n=14,15,16`, once THM-369 forces an `n`-gate, the local endpoint row for
that gate and the coarse denominator rows `q <= n` appear to be incompatible
with a primitive `n-1`-column all-protected repair.

The finite certificate should be a Hall/dual-weight inequality combining:

- endpoint rows owned by the mandatory `n`-gate;
- coarse divisor rows from THM-369;
- primitivity/gcd rows that forbid imprimitive even-window repairs.

## Evidence

`lrc_n14_15_16_gate_frontier_s21.py` audits the exact local endpoint set cover
for the `n`-gate endpoints using candidates `1..2n-1`.

- `n=14`: minimum window cover size is `7`; there are `24` minimum covers. The
  `6` covers that already pay every coarse divisor row are all gcd-`2` even
  covers. Primitive minimum covers leave divisor debt.
- `n=15`: minimum window cover size is `6`; all `32` minimum covers are
  primitive, but every one leaves at least two coarse divisor rows unpaid.
- `n=16`: the minimum window cover is unique, `(8,18,22,26,30)` plus the
  `16`-gate; it has gcd `2` and still misses rows `7,12,14`.

Greedy deterministic completions that pay all missing coarse rows remain
positive-gap for all six lower/window starts across `n=14,15,16`.

## Prediction

A proof should assign nonnegative weights to gate-endpoint rows, coarse divisor
rows, and primitivity rows so that every speed column has weight at most `1`,
but the total required row weight exceeds `n-1` unless an endpoint remains
unprotected or a positive interval gap remains.

## Artifacts

- `04-computation/lrc_n14_15_16_gate_frontier_s21.py`
- `05-knowledge/results/lrc_n14_15_16_gate_frontier_s21.out`
- `07-reflections/lrc-n14-n15-n16-gate-frontier-s21.md`
