# LRC14 Labelled Packet Theorem Gauntlet

- Created: 2026-06-24T10:15:00Z
- Coordinator: codex
- Cycle: manual-user-request
- External input: arXiv:2606.22636

## The Breakthrough Shape

The gauntlet plus boundary-moment bridge now wants a labelled packet theorem:

```text
Every primitive LRC14 row is in F0-F5,
or it is an F6 non-migrating boundary-moment kernel,
unless a new F7 Johnson sector exists.
```

This is not yet a proof.  It is a cleaner classification target for all
possible counterexamples.

## Families

```text
F0 qdiv witness discharge
F1 AP/GW boundary atoms
F2 unit-petal/two-block discharge
F3 K33 state-lift packet
F4 unlabelled q14 positive front
F5 covering positive / boundary-moment strictness
F6 covering zero-open non-migration kernel
F7 new Johnson-harmonic packet sector
```

An actual counterexample `M(S)<1/14` cannot lie in F0-F5.  It has to be F6,
unless F7 reveals a missing non-scalar packet sector.

## Exact Seeds

S151 audits the named seeds exactly:

```text
AP, GW 12->24            -> F1
12->36, P10+K33          -> F3
P10, P13, P10+GW         -> F2
12->26                   -> F0
12->96                   -> F4
12->84                   -> F5
```

S150's enlarged AP-neighborhood gauntlet remains the key local support:

```text
covered qdiv>=14 rows:         0
non-AP/GW boundary-only rows:  0
```

## Fixed-Margin Import

The new arXiv paper on binary fixed-margin swap chains contributes a proof
architecture:

```text
fixed-margin binary relation
  -> local swaps / two-row heat-bath comparison
  -> three-row reduction
  -> scalar count sector + Johnson harmonic sectors
```

For LRC14:

```text
scalar count sector     = qdiv, exact M/Farey node, Haar mass
Johnson harmonic sector = boundary owners, C27/unital labels, K33 debt
```

So the next theorem is a packet-sector theorem, not a scalar estimate:

```text
Every fixed-margin source-spectrum packet component reaches F0-F5 by swaps,
unless a three-row Johnson sector is a genuine F7 falsifier.
```

Then the boundary-moment bridge only has to close F6.

## Files

```text
04-computation/lrc14_labelled_packet_theorem_gauntlet_codex_s151.py
05-knowledge/results/lrc14_labelled_packet_theorem_gauntlet_codex_s151.out
05-knowledge/hypotheses/HYP-2956-lrc14-labelled-packet-classification-theorem.md
07-reflections/lrc14-labelled-packet-classification-theorem-codex-s151.md
```
