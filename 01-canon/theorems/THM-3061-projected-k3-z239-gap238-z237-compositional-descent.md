---
id: THM-3061
title: "Projected k3 z239 gap238 z237 compositional descent"
status: CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT
source: codex-lrc-z239-z237-composition-2026-08-01
depends_on:
  - THM-3052-projected-k3-z242-gap241-z240-compositional-descent
  - THM-3033-projected-k3-z246-to-z244-descent-and-z243-high-floor-addendum
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.py
output: 05-knowledge/results/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.out
script_sha256: 93a90c48ebed4bc782bca31f378ebb0d4f7ee19ef471a60220fcda3e8927e2fb
output_sha256: b30961378985a543a86a82681f8556353243703aca9bacf09bb6b4f61648274c
semantic_sha256: 2b825c5f12a59048b497b73dba234b79301d78e2db65325c63a580870eaccc88
hash_basis: LF-normalized bytes
---

# THM-3061 -- projected k3 z239 / gap238 / z237 compositional descent

**CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

This file remains outside the proved dependency graph until an independent
hostile audit and explicit status promotion. No navigation or ledger file is
changed by this candidate.

## Statement

In the lossless projected `k=3` body atlas inherited from THM-3052 and
THM-3033, the exact compositor closes

```text
all 4 occupied z_1=239 rows,
the zero-row atlas gap z_1=238,
all 44 occupied z_1=237 rows.                         (1)
```

The occupied-row decrement in `(1)` is `4+44=48`; the gap costs zero. Thus,
**if this candidate is promoted**, THM-3052's proved projected necessary-row
ledger and cap would update by the exact arithmetic

```text
374828-48=374780,
projected k=3 cap: z_1<=236.                          (2)
```

Equation `(2)` is an audit target, not a live ledger or navigation mutation
while the theorem has candidate status. This is not LRC(14), says nothing
about `k<=1` or the final rung, and closes only the named projected necessary
sector.

## 1. Pinned disjoint universe and the `238` gap

The compositor pins LF-normalized THM-3052 source/output/semantic hashes and
the THM-2941 projected-body atlas hash. It parses every one of the `6,060`
`row=` records by the structured row grammar. A separate literal `;z1=...;`
scan gives the same neighboring census:

```text
z_1=240: 52 rows = 49 below floor +  3 at/above floor;
z_1=239:  4 rows =  2 below floor +  2 at/above floor;
z_1=238:  0 rows;
z_1=237: 44 rows = 33 below floor + 11 at/above floor;
z_1=236:  1 row  =  1 below floor +  0 at/above floor. (3)
```

The positive `240/239/237/236` neighbors prevent a parser miss from
masquerading as an empty layer. The combined occupied bank has exactly

```text
48 rows = 35 wall + 13 order,                         (4)
```

with row-order digest
`76af1f045822eb75d507e98ecaae9eaf0b65c67ae6a1b8ce8432c9be40d05787`.
The `z_1=238` conclusion is exact absence from the pinned atlas universe, not
an inference from a zero-measure carrier.

## 2. Exact screen and first failed implication

The inherited crude and status gates give

```text
combined: 4018 states = 1823 crude + 1825 status + 370 residual;
order:      82 states =   71 crude +   11 status +   0 residual. (5)
```

Levelwise, `z_1=239` has only `13=1+12+0` states. At `z_1=237`, the
screen has

```text
4005 states = 1822 crude + 1813 status + 370 residual, (6)
```

and the `370` residuals occupy exactly `14` wall bodies. Therefore the first
failed implication is explicit: crude plus status closes `3,648` of `4,018`
states, not the remaining `370`.

## 3. Literal terminal closure

For each of the `14` residual bodies the THM-3052 terminal compositor rebuilds
the literal safe tables and proves a strictly positive duplicate-two-high gap.
The remaining split is exact:

```text
zero-high scalar cases: 341;
one-high translated cases: 377.                       (7)
```

Every one of the `377` one-high cases closes by strict complete-cell
cardinality. No max-gap fallback and no unit-phase fallback is used. The
minimum printed strict slack is `1`; the fourteen exact rational two-high gaps
and per-body case hashes are frozen in the companion transcript. Thus all
fourteen terminal bodies close, with zero failed case.

The direction is the inherited necessary one: every literal cover maps into
the relaxation, and the terminal cells retain actual safe residues. No
null-set or representative-root inference is made.

## 4. Exact evidence

The companion binds its `48` screen and `14` terminal certificates to cache
fingerprint

```text
1ab69b9536f8105c53f8b9df2499695b49307768c4089d9cf5fd5240bd7a6bd1.
```

Fresh disjoint normal and optimized runs produced `62/62` byte-identical
checkpoint envelopes and byte-identical transcripts:

```text
python 04-computation/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.py --processes 4 --checkpoint-dir <fresh-normal>
python -O 04-computation/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.py --processes 4 --checkpoint-dir <fresh-optimized> --output <optimized-output>
```

The stored output ends with `all_exact_controls=PASS`. Its semantic digest is
`2b825c5f12a59048b497b73dba234b79301d78e2db65325c63a580870eaccc88`.

## 5. Scope boundary

This theorem candidate is a compositional descent inside the already-defined
projected `k=3` necessary atlas. It does not upgrade the projected rows to
physical covers, does not change the `k<=1` debt, and does not prove LRC(14).
The exact consequence `(2)` becomes canonical only after independent hostile
audit and promotion.
