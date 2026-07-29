---
id: THM-2898
title: "Unique maximum gate six-seed parity closure"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.  A scratch
  exact computation closes the unique THM-2896 root with K=25 through six
  marked paid branches, 287 H4-pair residuals, and scalar cascades.  No
  result is claimed until the verifier is self-contained, hash-pinned,
  replayed, and independently audited.
source: root-2026-07-29
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
---

# THM-2898 -- unique maximum gate six-seed parity closure

**RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.**

The candidate root is

```text
E=(1,8,10,11,12,13,14),
```

the unique THM-2896 root with least gate size `K=25`.  A scratch
closed-state breadth-first search finds a six-paid-edge marked path

```text
(23,27,19,46,18,17).
```

All six paid carriers satisfy the THM-2895 singleton-cap gate.  Their H4
cores have sizes

```text
(13,13,10,11,7,5),
```

giving `287` pair residuals.  The direct top-three certificate closes
`286`; the `B2+B1` certificate closes all `287`; no recursive row remains.
Scalar cascades then activate all `25` gate apices.

This is only a reservation of the exact statement and intended mechanism.
The scratch program still imports unpinned scratch search modules and does
not yet provide a canonical ledger.  Until those defects are repaired and
audited, the theorem remains outside the proved graph and does not reduce
the official `3428`-root residual.
