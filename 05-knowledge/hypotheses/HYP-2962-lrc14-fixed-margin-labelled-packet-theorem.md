---
id: HYP-2962
title: LRC14 fixed-margin labelled packet theorem
status: CLASSIFICATION THEOREM TARGET / fixed-margin family scaffold, not a proof
source: codex-2026-06-24-S150
related:
  - HYP-2961
  - HYP-2960
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2947
  - HYP-2942
  - HYP-2937
  - HYP-2908
  - THM-523
  - THM-566
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.py
  - 05-knowledge/results/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.out
external:
  - https://arxiv.org/abs/2606.22636
  - https://arxiv.org/html/2606.22636
---

# HYP-2962: LRC14 Fixed-Margin Labelled Packet Theorem

## Claim

The breakthrough shape should be a **labelled packet theorem**.

HYP-2961 classifies possible strict LRC14 counterexamples into live families
and sporadic buckets.  HYP-2962 adds the missing grammar for the word
"family":

```text
family = fixed-margin labelled packet class
sporadic = bounded singleton, after family parameters are removed
```

The packet is not a raw row, raw residue set, or raw tournament class.  It is
an incidence object retaining exactly the proof data that the current LRC14
stack says cannot be dropped.

## External Proof-Pattern Import

This pass used arXiv:2606.22636, "Spectral Gap for the Binary Fixed-Margin
Swap Chain", by Fu, Qin, and Wang, submitted June 21, 2026, only as a proof
architecture analogy.  The paper studies binary matrices with prescribed row
and column sums, compares the swap chain with a two-row heat-bath chain,
reduces to a three-row model, and separates a scalar count sector from Johnson
harmonic non-scalar sectors:

- arXiv abstract: <https://arxiv.org/abs/2606.22636>
- HTML paper: <https://arxiv.org/html/2606.22636>

No theorem from that paper is imported into LRC.  The useful pattern is:

```text
fixed margins              -> packet families
swap moves                 -> packet-preserving mutations
two-row heat bath          -> conditional resampling of two packet rows
three-row reduction        -> local triples of proof obligations
scalar count sector        -> qdiv/q-cover margins
Johnson harmonic sectors   -> owner/source/C27/K33/moment labels
```

## Labelled Packet

For a primitive 13-row `S`, define a fixed-margin labelled packet `P(S)` with
rows/features:

```text
Q_d(S) for d=2..14:
  whether S contains a multiple of d.

U_u^0(S), U_u^1(S) for u in {1,3,5,9,11,13}:
  denominator-14 apex zero contacts and boundary contacts.

C27(S):
  C27 shell state counts and marked hole/double transfer.

O(S):
  strict Haar-open safe mass and component count.

F(S):
  exact M/Farey branch and source-spectrum binding labels.

K(S):
  C27/unital atom labels and K33/state-lift flag.

Y(S):
  boundary-moment / gK8 / L_y image, when the row is in qdiv>14.
```

The **fixed margins** are the row sums and typed totals in this packet.  A
packet-preserving mutation is any row mutation preserving these margins while
possibly changing raw runner names or exact representatives.  This is the LRC
analogue of staying inside a prescribed-margin binary-matrix fiber.

## Theorem Target

The theorem target is:

```text
Every primitive 13-row S maps to a fixed-margin labelled packet P(S).
If S is a strict LRC14 counterexample, then P(S) lies in exactly one live
HYP-2961 family:

  L1 apex-multiple residual,
  L2 genuine-wide zero-moment,
  L3 bounded covering core,
  L4 K33 zero-open state-lift,
  L5 unnamed source-kernel.

All other fixed-margin packet classes are discharged by the q-witness,
Haar-open migration, AP/GW boundary equality, unit-petal/GW strip discharge,
positive K33 state-lift labels, or existing wide/decorrelation machinery.
```

In this form a future proof can attack whole packet families rather than rows
one by one.  The desired "boundary-moment bridge" becomes the map from the
non-scalar labelled packet sectors into the `Y(S)` moment image.

## Computation

`04-computation/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.py`
implements the small representative classifier.  It imports existing exact
repo machinery:

```text
S138 exact M/Farey and C27 transfer,
S145 measurable packet rank,
S147 Haar/Baire strict-open mass.
```

It then emits a fixed-margin signature made from:

```text
qdiv bucket,
q-cover vector for d=2..14,
U14 zero profile,
U14 boundary profile,
C27 shell state counts,
packet atom keys.
```

The stored output is
`05-knowledge/results/lrc14_fixed_margin_labelled_packet_classifier_codex_s150.out`.

Representative run:

```text
rows audited                         = 18
fixed-margin classes                 = 16
shared family signatures             = 2
singleton sporadic signatures        = 14
strict counterexamples found         = 0
qdiv>14 representatives              = 5
qdiv>14 zero-open representatives    = 0
```

Bucket counts:

```text
D0 direct q-witness family                 4
D2 positive Haar-open row                  2
D2 positive covering Haar-open family      5
D3 positive unit-petal/GW strip            3
D4 positive K33/state-lift labelled row    2
S0 equality sporadic AP/GW                 2
```

Two shared fixed-margin signatures appear:

```text
direct-q apex decoy family:
  12->28 and 12->56

covering lcm-tail family:
  12->84 and 12->168
```

The second is the important one.  It shows how `qdiv>14` covering rows can
belong to an honest packet family while still being discharged by positive
Haar-open mass in this bank.

## Counterexample Readout

The run found no actual strict counterexample and no zero-open representative
in the `qdiv>14` covering branch.  In HYP-2961 language:

```text
L1 apex-multiple residual:
  not seriously tested; needs many 14-multiples.

L2 genuine-wide zero-moment:
  moment image not computed in this script.

L3 bounded covering core:
  represented by covering repairs and lcm tails; all positive-open.

L4 K33 zero-open state-lift:
  no zero-open representative found.

L5 unnamed source-kernel:
  no zero-open non-K33 representative found.
```

Thus HYP-2962 does not close LRC14.  It narrows the next computation:

```text
search for qdiv>14, zero strict-open, zero/nonpositive Y(S), non-K33,
non-unit-petal, non-AP/GW fixed-margin packet classes.
```

## Tournament Analysis

Assumption challenge: tournament vertices need not be runners, residues, arcs,
or one tournament isomorphism class.

Candidate vertex sets considered:

```text
runners,
gaps,
divisor clocks,
fixed circle sections,
section boundaries,
wall-crossing events,
residues,
cover arcs,
Fourier/moment modes,
matroid-like circuits,
K33/state-lift obligations,
fixed-margin packet families.
```

Chosen vertices:

```text
qdiv/q-cover count sector
Haar open-vs-boundary sector
source-spectrum pullback
C27/unital owner sector
K33 state-lift sector
gK8/L_y moment bridge
fixed-margin mutation family
unnamed source-kernel bucket
raw residue/tournament shadow
```

Pairwise observable:

```text
strict predicate retention,
owner-label retention,
source/kernel visibility,
named proof exit,
finite-family pressure,
anti-scalar guard.
```

Fingerprint from the script:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed 3-cycles = 0
SCCs = nine singleton SCCs
Hamiltonian paths = 1
```

Hamiltonian path:

```text
qdiv/q-cover count sector
> Haar open-vs-boundary sector
> source-spectrum pullback
> fixed-margin mutation family
> C27/unital owner sector
> K33 state-lift sector
> gK8/L_y moment bridge
> unnamed source-kernel bucket
> raw residue/tournament shadow
```

What this quotient preserves:

```text
strict-counterexample status,
open-front discharge,
boundary-owner labels,
source-spectrum ownership,
state-lift address,
and whether a live family parameter remains.
```

What it destroys:

```text
raw row identity,
fine wall-crossing order,
and runner names after conversion to proof packets.
```

## Next Work

1. Add the `Y(S)` moment image to the classifier, using the gK8/L_y scripts.
2. Build an L1-specific bank with many 14-multiples and fixed residual `R`.
3. Add local packet swaps inside each fixed-margin class and test whether
   positive-open discharge is invariant under those swaps.
4. Try to prove a three-feature reduction:

```text
bad fixed-margin packet
  -> one of three interacting rows/features:
       q-cover count,
       boundary/source owner,
       moment/K33 output.
```

That is the exact place where the arXiv fixed-margin proof pattern is most
suggestive.

