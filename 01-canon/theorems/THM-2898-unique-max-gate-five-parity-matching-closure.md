---
id: THM-2898
title: "Unique maximum gate five-parity matching closure"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.  A scratch
  exact computation closes the unique THM-2896 root with K=25 through five
  parity branches, 277 H4-pair residuals, one matching repair, and scalar
  cascades.  No
  result is claimed until the verifier is self-contained, hash-pinned,
  replayed, and independently audited.
source: root-2026-07-29
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
---

# THM-2898 -- unique maximum gate five-parity matching closure

**RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.**

The candidate root is

```text
E=(1,8,10,11,12,13,14),
```

the unique THM-2896 root with least gate size `K=25`.  The scalar
closed-state search first isolates the six-apex hostile core

```text
(23,27,19,46,18,17).
```

THM-2895 parity certifies the first five on their actual prefixes.  Their
H4 cores have sizes

```text
(13,13,10,11,7),
```

giving `277` pair residuals.  THM-2897's matching repair certifies apex `17`.
At its live carrier the heavy-edge cutoff is `1844`; the heavy graph has
only three edges.  Its sole disjoint pair lies below the required threshold
by `69186919/39666627000`, so `q_5+M_(2,2)<h`.  Scalar cascades then activate
all `25` gate apices.

This is only a reservation of the strengthened exact statement and intended
mechanism.  Until the combined matching/parity verifier is frozen and
independently audited, the theorem remains outside the proved graph and does
not reduce the official `3428`-root residual.
