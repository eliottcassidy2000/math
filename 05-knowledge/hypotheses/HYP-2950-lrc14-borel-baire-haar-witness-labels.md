---
id: HYP-2950
title: LRC14 witness sets need Borel code, Baire core, and Haar mass
status: COMPUTATIONAL SCOUT / proof-interface label proposal; not a proof
source: codex-2026-06-24
tags: [lrc14, borel, baire, haar, measure, category, boundary-witness]
related:
  - HYP-2951
  - HYP-2952
  - HYP-2949
  - HYP-2947
  - HYP-2943
  - HYP-2248
  - HYP-2247
results:
  - 04-computation/lrc14_borel_baire_haar_path_codex.py
  - 05-knowledge/results/lrc14_borel_baire_haar_path_codex.out
external:
  - https://en.wikipedia.org/wiki/Haar_measure
  - https://en.wikipedia.org/wiki/Baire_measure
---

# HYP-2950: Borel-Baire-Haar witness labels for LRC14

Haar's theorem gives the invariant volume coordinate on locally compact
groups: up to scale, there is a translation-invariant regular Borel measure.
For compact models this is the right way to count torus direction sets without
choosing an arbitrary coordinate chart.

But the LRC14 proof carrier should not collapse to Haar mass alone.

The finite scout in
[lrc14_borel_baire_haar_path_codex.py](/home/bigo/math/04-computation/lrc14_borel_baire_haar_path_codex.py:1)
models a compact direction group by `C_80`.  The output is stored in
[lrc14_borel_baire_haar_path_codex.out](/home/bigo/math/05-knowledge/results/lrc14_borel_baire_haar_path_codex.out:1).

The key toy witness:

```text
contiguous Borel arc       size=20 mass=1/4 robust=18 boundary=2
20-point Baire dust        size=20 mass=1/4 robust=0  boundary=20
```

These two sets have the same Haar mass in the finite quotient, but they have
opposite Baire behavior.  One contains a robust interval core; the other is
pure boundary dust.

## Claim

An LRC14 witness label should be a triple:

```text
(Borel code, Baire core/boundary split, Haar mass)
```

The intended readings are:

```text
Borel code:
  finite Boolean formula for the safe/danger event, such as torus inequalities,
  Farey branch labels, AP/Goddyn-Wong marks, and C27 owner constraints.

Baire core/boundary split:
  which parts of the event are robust under small perturbation, and which
  survivors exist only on equality walls or tangencies.

Haar mass:
  invariant volume of the event under the ambient compact direction/action
  group.
```

This order matters.  Haar mass can certify that an open danger region covers
almost all directions, but it cannot by itself decide a tight equality witness
living on the boundary.

## LRC14 Transfer

Earlier AP/Goddyn-Wong work suggests the `n=14` branch may cover open danger
arcs up to measure zero while leaving equality-time survivors.  HYP-2950 says
those survivors should not be discarded as measure artifacts.  They should be
classified as Baire-boundary witnesses with explicit Borel addresses.

The proof target becomes:

```text
Every candidate LRC14 residual either
  (a) has positive Haar danger mass,
  (b) loses its robust Baire core under the AP/GW/C27 labels,
  (c) or is forced onto a named Borel boundary wall where the exact equality
      calculation can be discharged.
```

This also sharpens HYP-2248.  The finite invariant-selector lesson was that a
selector needs an address tax before it can be meaningful.  Here the address
tax has three coordinates: Borel formula, Baire boundary status, and invariant
Haar mass.

## Guardrail

Do not use:

```text
Haar mass zero => irrelevant.
```

Use:

```text
Haar mass zero + named Borel boundary + nonempty Baire boundary witness
  => exact equality case, not noise.
```

This is the same reason a tournament proof cannot forget sink labels merely
because the exceptional sink class is small.
