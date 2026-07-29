---
id: THM-2886
title: "Stepped-origin provenance transport on the V4 horn"
status: >
  PROVED + VERIFIED-EXACT.  The THM-2884 diagonal H canonically selects
  the stepped-origin q3 Q0 subcopy and the empty zero-origin component on
  every one of the 449 horn sheets.  Transporting those two marked
  pieces, rather than reapplying H locally, gives the unique selector
  path q3:(notE3,notE3), q11:(notE3,E3),
  q7:(E3,notE3).  The empty/full origin occupancies remain (0,1),
  J3(a)+8U=J11(a), and J11(a)+9U=J7(a+1).  The fixed-source signed raw
  endpoint current is therefore the same nonzero value at q3,q11,q7 in
  both exact fields for all 26 samples and on all 449 sheets.  This
  closes selector parity and marked sheetwise coefficient transport.
  It does not yet descend to the canonical unmarked current: on the
  wrapped edge the Prony pair transforms by diag(omega^4,1), so U+V has
  a nonzero scalar-descent defect on all 26 field/section states.  No row
  exclusion or LRC(14) proof follows.
source: root/lrc-stepped-origin-provenance-2026-07-29
depends_on:
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2882-event-twisted-all-q-coefficient-carry-lift
  - THM-2884-macro-semantic-diagonal-horn-carrier-and-origin-even-boundary
related:
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2878-endpoint-factor-exit-carry-transducer-and-flat-vertex-boundary
  - THM-2880-q0-q3-one-fibre-selector-provenance-obstruction
script: 04-computation/lrc14_stepped_origin_v4_provenance_transport_thm2886.py
output: 05-knowledge/results/lrc14_stepped_origin_v4_provenance_transport_thm2886.out
script_sha256: 1f2cbffb8151c0c74bb22beff58e0ace5715eba3d4a9afc59481e7ab5e0d6dc9
output_sha256: 60ba517f1b8fa92fdaeba48b68d10da97f23f00a9fb1526de9e39f1311e49929
hash_basis: LF-normalized bytes
---

# THM-2886 -- stepped-origin provenance transport on the V4 horn

**PROVED + VERIFIED-EXACT.**

THM-2884 found the physical diagonal

```text
H=ker(e_E3 XOR u_7 XOR v_8) ~= V4                         (1)
```

and showed that local evaluation of its projector is origin-odd only at
`q3`; it cancels at `q11` and `q7`.  The next move is not to evaluate the
same projector again.  It is to retain the **image provenance** of the
unique stepped-origin `q3` piece under the two literal target
translations.

This produces an exact positive result.  The selector parity, fixed
source, semantic word, endpoint atom, endpoint origin, and raw coefficient
all transport together through the entire `q3 -> q11 -> q7` horn.  The
remaining failure occurs one level later: the two Prony channels do not
descend through their scalar recombination.

## 1. The local selector baseline

Order the two marked endpoint origins as

```text
(zero, stepped)=(o0,o1),                                 (2)
```

and encode `E3=1`, `not-E3=0`.  The target truth supports are

```text
o0: {0,3,11},             o1: {0,11}.                    (3)
```

Both origins retain the fixed source in `E3`.

An exhaustive census covers all four selector pairs at every residue and
all `32` object-independent componentwise Boolean recharts, including
origin swap.  It proves:

```text
target-only selector support             all 13 q,
same-E3-grade source/target support       {0,3,11},
same-grade q7 support                     empty.          (4)
```

The unique bijective rechart carrying the positive `q11` target selector
to the positive `q7` selector is global truth complementation.  It sends
the source truth pair `(1,1)` to `(0,0)`.  Thus no object-independent
Boolean operation simultaneously transports target, selector, and fixed
source within the old `E3/not-E3` grading.

Equation `(4)` remains the sharp baseline.  The result below uses the
finer diagonal grading `(1)`, not a same-`E3`-grade current.

## 2. The canonical q3 seed

On all `449` THM-2835 labels, the relevant diagonal states are

```text
fixed source QB : 101,
q3 Q0           : o0=100, o1=000,
q11 QA          : 110 at both origins,
q7 QAB          : 011 at both origins.                  (5)
```

Hence the local `H` memberships at `(q3,q11,q7)` are

```text
q3 :(0,1),
q11:(1,1),
q7 :(1,1).                                               (6)
```

Local signed evaluation therefore gives `(-1,0,0)`.  But `(6)` contains
more information than that signed vector: it canonically identifies a
full stepped-origin `q3` atom and an empty zero-origin component.

The `H` rule at `q3` selects

```text
s(3)=(not-E3,not-E3),                                    (7)
```

whose occupancy is `(0,1)`.  No endpoint origin is chosen by hand:
`(7)` is one common truth label, and the physical truth asymmetry in
`(3)` makes only the stepped component live.

## 3. Unique empty/full-preserving provenance transport

Transport each of the two marked pieces in `(7)`, including the empty
one.  Requiring empty to remain empty and full to remain full uniquely
forces

```text
s(3) =(not-E3,not-E3),
s(11)=(not-E3,E3),
s(7) =(E3,not-E3).                                      (8)
```

Their selector parities are

```text
(0,1,1)=1+delta_3                                       (9)
```

on the `q3/q11/q7` seam.  More strongly, exact weighted-set translation
gives, at each origin `o`,

```text
tau_(8U)(J3 intersect s_o(3))
  =J11 intersect s_o(11),

tau_(9U)(J11 intersect s_o(11))
  =J7 intersect s_o(7).                                (10)
```

The origin occupancies in `(10)` are `(0,1)` at all three checkpoints.
Thus the same stepped subcopy is followed throughout; the zero-origin
component is not discarded and recreated, but transported as the empty
subobject.

This does not conflict with THM-2880.  No affine endpoint-address map
changes the marked origin pair.  The endpoint origins in `(2)` are fixed;
only the target atom is translated, and its image provenance is retained.

## 4. Exact physical and semantic transport

Before reduction modulo the physical period, every `a` in the `449`-set
obeys

```text
J3(a)+8U = J11(a),
J11(a)+9U = J7(a+1).                                    (11)
```

The companion checks both identities on all labels: `898` literal
cylinder equalities.  It also reconstructs

```text
Q0(q3,a) -> QA(q11,a) -> QAB(q7,a+1)                    (12)
```

on all `449/449` sheets, with

```text
min a=59306,       max a=311961,
a=156 mod169,      a+1=157 mod169.                      (13)
```

The source remains the same fixed `I,a,QB,E3` object.

On the control cell

```text
(s,t,clock)=(0,5,1),                                    (14)
```

the six-factor signatures of source, `q3`, and `q11` are all ones, while
`q7` has

```text
(0,1,1,1,1,1)                                          (15)
```

in order `(E3,clock,q1,q2,c2,c3)`.  Thus only the old macro grade changes.
THM-2884 proves that the same statement holds on its exact `20`-cell
bank after replacing `E3` by `H`.  Since the endpoint/word transport is
cell-independent, `(10)`--`(12)` combine with all

```text
20*449=8980                                             (16)
```

THM-2884 packets.

## 5. The nonzero transported raw current

Let `C_m(q,o)` be the exact right-endpoint value on the selected piece,
and let `P` be the inherited nonzero fixed-source coefficient.  For each
of the `26` lawful unit samples, `(8)` gives

```text
(C_m(q,o0),C_m(q,o1))=(0,C_m)
```

at `q=3,11,7`.  Target translation by a multiple of `U` has zero endpoint
phase at these frequencies.  Consequently

```text
P(C_m(q,o0)-C_m(q,o1))=-P C_m !=0
```

and the value is literally the same at all three checkpoints.             (17)

The companion proves `(17)` in both exact fields

```text
352341050142921841,
956354278959359281                                      (18)
```

for all `52` field/sample rows.  Multiplication by `449` remains nonzero
in both fields, so the unit-weighted sum of the `449` identical
origin-resolved endpoint copies is nonzero.

The last sentence is deliberately sheetwise.  It does not identify the
unit-weighted sum with the inherited canonical semantic Fourier current;
registering the marked `H` subcopy in that current interface remains a
separate descent obligation.

## 6. First failed implication: scalar charged descent

Prony splitting writes the raw endpoint value as

```text
S_r=U_r+V_r,                                             (19)
```

where frequency-section transport has

```text
U_(r+1)=omega^3 U_r,       V_(r+1)=V_r.                 (20)
```

On the wrapped `q11 -> q7` edge, ancestry carry changes the section by the
amount giving

```text
E=omega^4,
(U,V) -> (E U,V).                                       (21)
```

Both channels are nonzero.  Therefore

```text
S' - E S
 =(E U+V)-E(U+V)
 =(1-E)V !=0.                                           (22)
```

The companion verifies `(22)` for all `13` ancestry residues in both
fields: `26/26` field/section failures.  Thus the marked raw current in
`(17)` is not a scalar `C169`-charged local system after recombination.
The charged line `U` transports, and the trivial line `V` transports, but
their sum must be retained as a two-channel object.

## 7. Exact consequence and boundary

This theorem proves:

1. a canonical origin-independent `H` seed at `q3`;
2. the unique marked-origin, empty/full-preserving selector transport
   `(8)`;
3. literal target, ancestry, word, and fixed-source transport on all
   `449` sheets;
4. compatibility with the `20` THM-2884 physical cells;
5. the same nonzero signed raw endpoint current at `q3,q11,q7`;
6. the exact two-channel scalar-descent obstruction `(22)`.

It does **not** prove:

- that the marked path-space subcopy has already descended to the
  canonical unmarked THM-2334/2365 current;
- scalar transport of `U+V`;
- a row exclusion;
- LRC(14).

The next live object is not another scalar selector.  It is a typed
two-channel or Hermitian observable on the marked `H` path object, together
with an explicit descent from that object to the canonical current
interface.

## 8. Reproduction

Run

```text
python3 04-computation/lrc14_stepped_origin_v4_provenance_transport_thm2886.py
python3 -O 04-computation/lrc14_stepped_origin_v4_provenance_transport_thm2886.py
```

Both executions have no executable Python `assert` and byte-match

```text
05-knowledge/results/lrc14_stepped_origin_v4_provenance_transport_thm2886.out.
```
