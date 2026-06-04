---
id: THM-408
title: Moser layered slabs have explicit unit Hamiltonian spines
status: PROVED
source: codex-2026-06-04-S628
depends_on:
  - HYP-2203
related:
  - HYP-2202
  - HYP-2201
  - HYP-2188
---

# THM-408 - Moser Layered Slabs Have Unit Spines

## Statement

Work in the rank-4 Moser coordinate carrier used by S614/S617/S626.  Its unit
shell contains the following vectors:

```text
(-1,1,0,0), (0,-1,0,0), (0,0,-1,0),
(0,0,1,0),  (0,0,0,-1), (0,0,0,1)
```

For every integer `m >= 0`, define the `plus` layer word

```text
W_a^+ =
(a,1,1,-1), (a,1,0,-1), (a,1,0,0), (a,1,1,0),
(a,0,1,0),  (a,0,0,0),  (a,0,0,-1), (a,0,1,-1)
```

and the cap

```text
C^+ =
(0,1,1,-1), (0,1,0,-1), (0,1,0,0),
(0,1,1,0),  (0,0,1,0),  (0,0,0,0).
```

Let `P_m^+` be the set listed by

```text
W_m^+, W_{m-1}^+, ..., W_1^+, C^+.
```

Then the unit-distance graph on `P_m^+` has a Hamiltonian path consisting only
of unit edges.  In particular `|P_m^+| = 8m+6`; the S626 `n=14` and `n=22`
Moser witnesses are `P_1^+` and `P_2^+`.

There is also a shifted `minus` family.  Define

```text
W_a^- =
(a,1,0,-1),  (a,1,-1,-1), (a,1,-1,0), (a,1,0,0),
(a,0,0,0),   (a,0,-1,0),  (a,0,-1,-1), (a,0,0,-1)

C^- =
(0,1,0,-1), (0,1,-1,-1), (0,1,-1,0), (0,1,0,0), (0,0,0,0).
```

The set `P_m^-` listed by `W_m^-,...,W_1^-,C^-` also has a Hamiltonian unit
spine.  Here `|P_m^-|=8m+5`; the S626 exact `n=21` witness is `P_2^-`.

## Proof

The displayed order lists every vertex of `P_m^+` exactly once: the layers have
distinct first coordinate `a=1,...,m`, and the cap lies in the separate
`a=0` layer.  Thus it remains only to check that each consecutive difference is
a Moser unit vector.

Inside a full `plus` layer the consecutive differences are

```text
(0,0,-1,0), (0,0,0,1), (0,0,1,0), (0,-1,0,0),
(0,0,-1,0), (0,0,0,-1), (0,0,1,0).
```

All are in the Moser unit shell.  The bridge from the last vertex of `W_a^+`
to the first vertex of `W_{a-1}^+` is

```text
(a-1,1,1,-1) - (a,0,1,-1) = (-1,1,0,0),
```

also a Moser unit.  The bridge from `W_1^+` to `C^+` is the same vector.  The
cap `C^+` uses only the first five internal step types from the full layer.
Therefore the displayed order is a Hamiltonian unit path on `P_m^+`.

The `minus` proof is identical after translating the `(c,d)` square.  Its
internal differences are the same list up to sign/order:

```text
(0,0,-1,0), (0,0,0,1), (0,0,1,0), (0,-1,0,0),
(0,0,-1,0), (0,0,0,-1), (0,0,1,0),
```

and the bridge is again `(-1,1,0,0)`.  The cap `C^-` is a prefix of the same
layer word plus the final origin attachment by `(0,-1,0,0)`.  Hence `P_m^-`
also has a Hamiltonian unit path.  This proves both families.  QED.

## Edge Counts

The spine theorem does not need edge maximality, but the sampled slabs satisfy
the counts used in S626:

```text
P_1^+ : n=14, edges=33
P_2^- : n=21, edges=57
P_2^+ : n=22, edges=60
```

The certifier checks the formulas

```text
E(P_m^+) = 27m+6 for m>=1, with E(P_0^+)=8,
E(P_m^-) = 27m+3 for m>=1, with E(P_0^-)=6.
```

## Abstract Reading

The theorem says that the S626 unit spine is a **traceable section** over a
layer quotient.  Each 8-vertex fiber is a Gray-code slab, and the bridge vector
`(-1,1,0,0)` is the gluing morphism between adjacent fibers.  This is why a
fixed point-order tournament can lie: the traceable section exists upstream of
the point-order quotient, and the quotient can forget the section.

This places unit spines next to three existing repo abstractions:

- **ear/slab gluing:** a Hamiltonian path persists when traceable blocks have
  endpoint-compatible unit bridges;
- **certificate sheaves:** the spine is a local section that glues across
  carrier fibers, echoing the LRC n=14 certificate-sheaf language;
- **Tournament Analysis:** the useful tournament vertices are proof routes,
  slabs, ears, or traceability obligations, not only raw points.

## Verification

`04-computation/unit_distance_spine_ladder_s628.py` verifies the step alphabet,
samples both families through `m=8`, checks the edge-count formulas above, and
records the proof-route tournament.

**Artifacts:** `04-computation/unit_distance_spine_ladder_s628.py`,
`05-knowledge/results/unit_distance_spine_ladder_s628.out`,
`07-reflections/unit-spines-as-traceable-sections-s628.md`.  Builds on
HYP-2203, HYP-2202, HYP-2201, and the Moser carrier from HYP-2188.
