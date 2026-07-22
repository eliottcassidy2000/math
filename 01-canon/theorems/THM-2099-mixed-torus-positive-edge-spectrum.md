---
id: THM-2099
title: "The mixed-torus positive-edge spectrum and rank-eight tree gate"
status: >
  RESERVED / IN PROGRESS. The proved seed is that a transverse pair whose
  joint danger band meets the guard complement is governed by the primitive
  three-character relation among the guard and the pair. The zero layer is
  THM-2097 pair rigidity. This file is reserved to quantify the positive
  layers and test whether their exact threshold graph forces a spanning-tree
  weight above THM-2098's rank-eight budget 5/49. No positive-edge floor,
  tree bound, rank-eight exclusion, or LRC(14) conclusion is claimed yet.
source: codex-2026-07-22-LRC-mixed-torus-positive-spectrum
depends_on:
  - THM-2097
  - THM-2098
related:
  - THM-1221
  - THM-2096
---

# THM-2099 -- reserved positive-edge spectrum

## Exact target

For transverse terminal characters `f_i,f_j` and guard `g`, determine or
sharply bound

```text
w_ij=measure({||f_i||<1/14, ||f_j||<1/14, ||g||>1/7})
```

from the primitive integer relation among `(g,f_i,f_j)`, including the
finite-kernel case. Build the nested graphs obtained by thresholding these
weights, and minimize their maximum spanning-tree weight over eight distinct
positive specializations.

The desired strict conclusion would be

```text
tau>5/49,
```

which would contradict THM-2098 for an eight-band transverse cover. Equality
and counterexample families must be classified before any promotion.

## Current proof seed

When `g=a f_i+b f_j` with both coefficients nonzero, rescale the terminal
box to `[-1,1]^2`. The conditional outside-guard fraction is the exact area
fraction

```text
area{(x,y) in [-1,1]^2:|a x+b y|>2}/4.
```

The `|a|+|b|=2` layer is zero (THM-2097). The first positive relation
`{|a|,|b|}={1,2}` has fraction `1/8`, hence weight `1/392`. These are proved
only for the stated integral-factor case; the finite-kernel and dependent-
terminal cases remain to be incorporated.

## Assumption challenge and carrier

The challenged assumption is that all positive pair intersections can be
collapsed to one bit. Their actual relation height and weight are the likely
sidecars. Candidate vertices are terminal characters, relation types,
positive edges, and spanning-tree obligations. A tournament orientation by
weight is only a scheduler: it forgets the threshold colors and the graphic-
matroid maximum. The faithful provisional carrier is the edge-weighted
complete graph decorated by primitive three-character relation data.
