---
id: HYP-2884
status: COMPUTATIONAL SIGNAL / proof-target
source: codex-2026-06-22-S102
tags: [lrc14, coimage, character-kernel, reciprocal-tail, packet-graph, wall-ledger, abel-summation, tournament-analysis]
related:
  - HYP-2883
  - HYP-2887
  - HYP-2632
  - HYP-2633
  - HYP-2634
  - HYP-2636
  - HYP-2617
  - HYP-2828
  - OPEN-Q-108
---

# HYP-2884: the lifted packet divergence defect is the LRC local-current lemma

HYP-2883 turns the HYP-2632 repeated-residue kernel into a finite signed graph
with exact local balance:

```text
loop(a) + sum_b edge(a,b) = 0
```

for `a in {0,2,3,4,5,6}`.  HYP-2633 warns that finite packet signs do not
automatically survive the actual reciprocal relation-lattice lift.  The right
bridge between those statements is a local divergence defect.

The conjectural LRC proof lemma should be:

```text
finite packet current
  -> integer-lift divergence defect
  -> finite low-height wall deletion
  -> additive-frequency Abel summation
  -> signed reciprocal-tail bound.
```

## Computation

Script:

- `04-computation/lrc14_packet_balance_lift_probe_codex_s102.py`
- output: `05-knowledge/results/lrc14_packet_balance_lift_probe_codex_s102.out`

The script fixes the four residue-`1` core

```text
CORE = (1,8,15,22)
```

and compares two integer-lift gauges for the finite graph:

```text
start_aligned: edges use residues in layer 0, loops use layers 0 and 1
raised_pair:   edges use residues in layer 1, loops use layers 1 and 2
```

For each lifted packet, it computes the actual `d=9` support-six reciprocal
tail through height `H=12`, using the HYP-2614/S12 residue-cusp machinery.
Then it computes the lifted divergence

```text
div_H(a) = lift_H(loop a) + sum_{b != a} lift_H(edge a,b).
```

## Findings

The finite packet graph still has exact integer-unit balance:

```text
v   loop  incident  balance
0     -4         4        0
2    -25        25        0
3    -18        18        0
4    -25        25        0
5    -18        18        0
6    -18        18        0
```

But the reciprocal lift has nonzero divergence.  At height `H=12`:

```text
start_aligned:
  max |div| = 0.00512112
  L1 div    = 0.0193444
  sum div   = 0.0193444

raised_pair:
  max |div| = 0.00191161
  L1 div    = 0.00610376
  sum div   = 0.00430821
```

Raising the repeated pair does not prove the lemma, but it cuts the measured
`H=12` local-divergence mass by roughly a factor of `3.17`.  This matches the
HYP-2634 warning: much of the finite/lift sign mismatch is an integer-height
placement phenomenon, not a failure of the mod-7 packet identity.

The start-aligned lifted supports have `LOW_HEIGHT=2` relation-wall counts in
the range `62..86`; the raised-pair supports reduce that range to `58..68`.
The largest start-aligned wall count is the loop at residue `0`:

```text
(1,7,8,14,15,22): 86 low-height wall directions.
```

This suggests that the finite packet current is not the final proof object by
itself.  The proof object is the pair:

```text
(finite current, integer-height wall filtration).
```

## Proof target

The next lemma should not assert that local balance is preserved under every
integer lift.  It should assert a controlled divergence after deleting the
finite low-height wall ledger:

```text
For each repeated-packet lift gauge G and height cutoff H, write

  div_H^G(a) = loop_H^G(a) + sum_b edge_H^G(a,b).

After removing coefficient-height <=2 wall directions and grouping the
remaining relation lattice by HYP-2636 additive-frequency shells,

  sum_a |AbelTail(div^G(a))|

is bounded by the available LRC(14) support-six margin.
```

The reason this is more promising than the scalar HYP-2632 ledger is that
vertex divergence is a local object.  It can be attacked with wall deletion,
summation by parts, and relation-depth/Freiman structure at the residue vertex
where the defect occurs.

## Winding-tournament interpretation

Incoming KPS S31i adds a useful runner-shadow interpretation.  In the winding
tournament model (HYP-2605), lonely phases have a scale-`1/7` local sink, while
non-lonely phases have no such sink and carry higher directed-`3`-cycle
content.  The AP probes report non-lonely phases with about `1.22..1.25x` the
mean cyclic-triangle content of lonely phases.

HYP-2884 is the support-six reciprocal-tail version of the same statement.
The finite mod-7 packet graph is locally balanced, but an integer lift creates
cyclic excess measured as divergence.  Thus:

```text
winding tournament no-sink / high C3 content
  <-> cyclic packet current not absorbed by a local sink
  <-> lifted reciprocal divergence after choosing integer representatives.
```

This does not prove LRC(14), but it gives a common observable for two formerly
separate pictures: winding-tournament cyclicity and support-six coimage tails.
The finite wall ledger should remove coherent low-height cyclic scars; the
remaining spread cyclicity should be small by additive-frequency Abel
summation.

## LRC route

The live LRC chain becomes:

```text
HYP-2617 finite coimage atlas
  -> HYP-2632 finite signed packet graph
  -> HYP-2883 local packet balance
  -> HYP-2884 lifted divergence defect
  -> HYP-2636 additive-frequency block transfer
  -> HYP-2828 relation-depth separator
  -> OPEN-Q-108 support-six tail bound.
```

For the genuine-wide branch, this gives a more surgical target than "bound the
wide support-six correction."  Bound the divergence defects of packet vertices.
High additive energy should concentrate the defect into finitely many
low-height walls; low additive energy should make the HYP-2636 channelwise
Abel/Cauchy bound small.

## Assumption challenge

Candidate vertices considered: runners, raw speed supports, residue packets,
finite loops, finite edges, integer lifts, low-height wall motifs, lifted
divergence defects, additive-frequency shells, relation-depth buckets, and
proof obligations.

Chosen quotient: lifted packet divergence.  It preserves the LRC predicate
needed by HYP-2633: signed reciprocal-tail cancellation after finite wall
deletion.  It destroys raw runner labels and raw finite packet scalar sums.

The challenged assumption is that the exact finite mod-7 balance should lift
unchanged to integer relation lattices.  The computation refutes that naive
claim but replaces it with a sharper proof object: divergence plus wall
filtration.

## HYP-2887 Follow-Up

HYP-2887 identifies the finite topology behind this divergence defect.  The
nonzero repeated-packet graph is the octahedron `L(K4)`: residues are the six
edges of a tetrahedron, and the affine-zero lane is the opposite-edge matching.
Thus an integer lift gauge is a height cochain on an octahedral carrier, not an
arbitrary packet orientation.

This sharpens the local-current lemma again.  After the height-2 wall ledger is
deleted, the remaining defect should be decomposed into:

```text
vertex divergence on L(K4)
+ triangular face curl on the 8 octahedron faces
+ HYP-2636 additive-frequency tails.
```

Since the octahedron is spherical, there is no independent harmonic
one-current once divergence and face curl are controlled.  This is the most
concrete current-realizability route found so far.
