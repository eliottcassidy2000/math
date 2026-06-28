---
id: HYP-3425
title: LRC14 additive energy needs a sheet sidecar
status: SYNTHESIS / exact counterexample-and-repair scout; not an LRC14 proof
source: monad-explorer-2026-06-28 pass on HYP-3424's add/mult transfer rule, integrating HYP-3140 fiber PGF, HYP-3129 signed SPEC, HYP-2272 additive energy, and HYP-3422/HYP-3421 covering packets
tangent: T1386
related:
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3418
  - HYP-3415
  - HYP-3140
  - HYP-3129
  - HYP-2272
  - HYP-2129
  - HYP-2128
  - THM-414
  - OPEN-Q-108
script: 04-computation/lrc14_additive_energy_sheet_sidecar_codex_20260628.py
result: 05-knowledge/results/lrc14_additive_energy_sheet_sidecar_codex_20260628.out
reflection: 07-reflections/lrc14-additive-energy-sheet-sidecar-codex-20260628.md
---

# HYP-3425: LRC14 Additive Energy Needs A Sheet Sidecar

## Claim

HYP-3424's add/mult transfer rule is real, but only in packet form.

The old scalar summaries

```text
full additive energy E_+
R-only additive energy
odd-sector additive energy
```

do not by themselves control the covering-floor bias

```text
Rprime - 1 = E[N_R | Q-lonely] / E[N_R] - 1.
```

Exact covering packets collide on those scalar energies while having different
`Rprime` values, different signs of `delta = E[N_R | Q] - E[N_R]`, and
different sheet behavior.  The first legal repair is therefore an
energy-plus-sheet packet, not another scalar.

The right reading of HYP-2272 after HYP-3418/HYP-3424 is:

```text
additive energy
-> possible floor sidecar only after keeping a sheet-profile coordinate
   such as q_zero_mass, q_sheet_range, or another HYP-3140 packet field.
```

## Executable Readout

Script:

```text
04-computation/lrc14_additive_energy_sheet_sidecar_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_additive_energy_sheet_sidecar_codex_20260628.out
```

The scout combines:

```text
- HYP-3140 canonical Rprime rows,
- HYP-3422 covering_AP_with_84 and covering_AP_with_12_and_84,
- a small exact random covering sample (seed 2, max_speed 90, Qsize <= 3).
```

## Exact Counterexamples

### 1. Full additive energy collides

The curated bank contains:

```text
canonical_r3 and canonical_r5 with fullE = 389
```

but

```text
canonical_r3: Rprime = 0.967327..., delta = -0.025767...
canonical_r5: Rprime = 1.039567..., delta = +0.100299...
```

So a raw full-energy scalar does not even determine the sign of the sheet bias.

### 2. R-only additive energy collides

```text
canonical_r1_drop12 and covering_AP_with_84 both have RE = 246
```

but

```text
canonical_r1_drop12: Rprime = 7/6,       delta > 0, q_zero_mass = 26/35
covering_AP_with_84: Rprime = 0.513954..., delta < 0, q_zero_mass = 83/105
```

This is the cleanest rejection of "R-energy alone controls the floor."

### 3. Odd-sector energy collides

```text
canonical_r1_drop12 and covering_AP_with_84 both have oddE = 47
```

with the same sign split as above; similarly the `oddE = 29` class contains
three rows with visibly different `Rprime`.

So HYP-3424's "odd data becomes phase debt" is not just rhetoric: odd-sector
energy alone is too coarse to carry the floor.

## Minimal Repairs On The Curated Bank

On the exact curated bank, no single scalar among the tested energy summaries
survives.  But several two-coordinate packets do.  The scout finds sign-
separating pairs such as:

```text
(RE, q_zero_mass)
(oddE, q_zero_mass)
(fullE, q_range_hi)
(RE, sumv2_even)
```

and full row-separating pairs such as:

```text
(fullE, RE)
(fullE, q_range_hi)
(RE, sumv2_even)
```

Interpretation: the energy coordinate becomes useful only after a sheet-sidecar
or another explicit floor-facing packet coordinate is retained.

## Small Random Covering Sample

The small exact random bank is only a scout, not a theorem, but it supports
the negative half:

```text
corr(RE, delta)    = +0.628
corr(oddE, delta)  = +0.134
corr(evenE, delta) = -0.047
```

Thus even on fresh primitive covering rows, the raw odd/even energy scalars are
weak.  `RE` carries some signal, but the exact collision pairs show it is not
safe without a sidecar.

## Candidate Packet Theorem

For every primitive covering packet routed through HYP-3424's additive branch,
the proof should keep a packet of the form

```text
(additive energy coordinate, sheet-profile sidecar)
```

where the sheet-profile sidecar is one of:

```text
q_zero_mass,
q_sheet_range,
Qsize / far-sheet depth,
or another HYP-3140-equivalent fiber field.
```

Then the packet should feed one of:

```text
low-SPEC penalty,
signed-SPEC/fiber-PGF inequality,
odd phase-cover debt,
or named terminal debt.
```

The forbidden move is:

```text
energy scalar -> terminal floor claim
```

without a sheet sidecar.

## Assumption Challenge

Alternate scalar summaries tested by the scout:

```text
fullE, RE, oddE, evenE, shell count, Qsize, q_zero_mass, q_range_hi, sumv2_even
```

The useful object is not a single winner from that list.  The useful object is
the packet that keeps one energy coordinate and one sheet coordinate together.
