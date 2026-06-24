---
id: HYP-2969
title: LRC14 boundary-moment packet ledger
status: PROOF-INTERFACE / finite boundary-moment proxy ledger; not a proof
source: codex-2026-06-24-S154-ledger
related:
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2950
  - HYP-2947
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_boundary_moment_packet_ledger_codex_s154.py
  - 05-knowledge/results/lrc14_boundary_moment_packet_ledger_codex_s154.out
external:
  - https://arxiv.org/abs/2606.22636
---

# HYP-2969: LRC14 Boundary-Moment Packet Ledger

HYP-2969 turns the phrase `COVERING-MOMENT` into an executable packet ledger.
It complements the HYP-2964 moon-core proof skeleton, the HYP-2965 boundary-gap
bridge, the HYP-2966 NORK pinch-template atlas, the HYP-2967 apex-aperture
comb gate, and the HYP-2968 few-apex lift-packet bridge by auditing the
covering boundary-moment residual after those local positive-open certificates.
It does not prove the HYP-2954 boundary-moment
bridge, but it adds a concrete
finite proxy for the missing map

```text
adaptive exact-period packet boundary
  -> missed-depth sector vector
  -> L_y / gK8-style moment readout.
```

For a labelled packet row `S` and a chosen exact-period denominator `D`, scan
unit packets `a in (Z/DZ)^*` and record:

```text
covered / boundary / strict threshold state at t=a/D,
six-sector hit support of {s*a/D : s in S},
missed-depth histogram q_0,...,q_6,
L_y = 10*q_0 + q_3 + 10*q_6.
```

This is intentionally a scout ledger.  It preserves the packet theorem's
labels but still destroys continuous torus geometry between the sampled
denominators.  Therefore a single all-covered denominator chart is not a
counterexample certificate.

## Computation

Script:

```text
04-computation/lrc14_boundary_moment_packet_ledger_codex_s154.py
```

Stored output:

```text
05-knowledge/results/lrc14_boundary_moment_packet_ledger_codex_s154.out
```

Default run is a curated theorem-facing bank, not the full HYP-2963 audit.  It
contains named adversaries, AP/GW, unit petals, S138 splices, direct-q decoys,
covering lcm tails, a K33 larger marker, and fattening/comb rows.

## Main Readout

The default run audited `35` source packets and emitted `29` moment ledgers:

```text
below-threshold packets = 0
zero-open packets       = 2
dangerous moment-kernel rows = 0
```

Kernel flags:

```text
AP/GW-zero-open-equality    2
named-K33-state-lift        3
positive-Haar-open         10
q-witness-discharge        10
unit-petal-named            4
```

Route-level moment statistics:

```text
Q-WITNESS                n=10, all-covered@D=0
BOUNDARY-AP-GW           n= 2, all-covered@D=0
BOUNDARY-PETAL-SPORADIC  n= 4, all-covered@D=0
K33-STATE-LIFT           n= 3, all-covered@D=0
COVERING-MOMENT          n=10, all-covered@D=7
```

The last line is the useful warning.  Several covering rows are all-covered at
the selected exact-period chart while still having positive Haar-open mass.
Thus the final bridge cannot be:

```text
one denominator chart all-covered -> obstruction.
```

It must be a labelled multi-chart or feasible-region theorem:

```text
all relevant exact-period charts + retained packet labels
  -> positive boundary-moment/gK8 image
  or named K33/TournamentStateLift debt.
```

## Relation To arXiv:2606.22636

Fu-Qin-Wang prove a spectral gap for binary fixed-margin swap chains by keeping
fixed margins, comparing swap moves to two-row heat-bath moves, reducing to
three rows, and splitting the scalar count sector from Johnson harmonic
sectors.  HYP-2969 uses this only as proof architecture:

```text
fixed margins        -> labelled LRC packet signatures
swap fiber           -> packet-preserving mutations
scalar sector        -> qdiv / exact M / Haar-open status
Johnson sectors      -> C27, K33, source-spectrum, boundary-moment
three-row reduction  -> triples of proof obligations, not runner triples
```

The lesson is conservative: scalar count evidence is valid only after the
non-scalar packet sectors have been conditioned on.

## Tournament Analysis

Vertices are proof layers / packet sectors, not runners:

```text
qdiv_gate
exact_M_Farey
Haar_Baire_front
C27_K33_labels
fixed_margin_packet
boundary_moment_ledger
TournamentStateLift
raw_scalar_or_residue
```

Pairwise observable:

```text
which layer preserves strict-counterexample status and owner labels before
scalarization.
```

Fingerprint:

```text
score_hist={0:1,1:1,3:2,4:1,5:1,6:2}
directed_3cycles=4
Hamiltonian_paths=17
SCCs=(6,1,1)
```

Unlike the transitive fixed-margin classifier, the boundary-moment ledger sits
inside a nontrivial SCC.  That is the right warning: qdiv, exact M, Haar,
C27/K33, fixed-margin packet, and boundary-moment data must be transported
together until the final proof chooses a discharge.

## Theorem Target

The sharpened target is:

```text
Every primitive LRC14 residual emits a fixed-margin labelled packet.  If the
packet is strict-bad, then qdiv>14 and it lies in the covering
boundary-moment fiber.  That fiber either:
  (1) has positive Haar/moment image in the true multi-chart B_D map,
  (2) carries a named K33 state-lift debt discharged by HYP-2908/THM-572,
  (3) or exposes a new Johnson-harmonic sector.
```

The default HYP-2969 ledger finds no row of type (3).  The remaining work is
to replace the finite sector proxy by the true gK8/L_y feasible-region map.
