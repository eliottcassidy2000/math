---
id: THM-3042
title: "Reflected spread-five closure and three-D cone"
status: >
  PROVED + VERIFIED-EXACT + SCOPE-AUDITED.  In the reflected k=1 sufficient
  residual family of THM-2941, every packet of spread D<=5 closes at every
  positive minimum level.  Its bodywise bank closes arbitrary positive levels
  on 2,421 of the 3,003 six-label bodies and leaves 582 bodies uncovered; it
  does not assert that those 582 fail other certificates.  Every packet with
  D>=6 and m>=3D also closes.  Thus the assembled certificate-failure locus
  is confined to those 582 bodies with D>=6 and 1<=m<3D.  This is not a
  physical-survivor classification and not LRC(14).
source: codex-reflected-chain-recover-2026-08-01
audit: >
  The promotion audit checked the reflected-family typing, the three-way
  logical intersection, all body and spread partitions, strict overlap/debt
  directions, the two exceptional chromatic graphs, the D=5 head/tail seam,
  the robust edge-nine body partition, the unique zero (3,4) phase channel,
  all three reverse ladders, and the remaining-wedge quantifiers.  The
  assembly pins LF-normalized sources and outputs for all ten terminal inputs
  and requires their ordinary/-O and exact-control tokens.  Fresh ordinary
  and optimized assembly replays are byte-identical to each other and to the
  stored transcript; the promoted sources compile and contain no
  optimization-stripped assertions.
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1226-gcd-period-projective-charge-obstruction
related:
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
script: 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
output: 05-knowledge/results/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.out
script_sha256: 2e242a99b8b95ab2b03f0b13f2347aaa6faf9f09e1b987c28b96ef5d25f5c906
output_sha256: 1a3788f5e224d3be247104b6963cf8fe9d71af45650a70798e53dc66174c942e
semantic_sha256: 9bc1da061c14d44ef1016a9b887097666fd50ea6af242ea8302822a6dbe664b1
hash_basis: LF-normalized bytes
---

# THM-3042 -- reflected spread-five closure and three-D cone

**PROVED + VERIFIED-EXACT + SCOPE-AUDITED.**

## Statement

Let

```text
E subset {1,...,14},   |E|=6,   L=14*lcm(E),
Z(E,q)={q_e L-e:e in E},
m=min_e q_e,           D=max_e q_e-min_e q_e,
```

where all six levels `q_e` are positive integers.  Restrict to the reflected
`k=1` residual family and its sufficient projection-mass certificate from
THM-2941.  Then:

1. every packet with `D<=5` closes for every `m>=1`;
2. the assembled bodywise bank closes every positive level assignment on
   `2,421` of the `3,003` bodies, with no spread restriction; and
3. every packet with `D>=6` and `m>=3D` closes on all `3,003` bodies.

Intersecting these statements confines possible **certificate failures** to

```text
582 bodies,   D>=6,   1<=m<3D.                         (1)
```

The word “bank” is essential: the other `582` bodies are not covered by this
bodywise bank, rather than proved to lack some other arbitrary-level
certificate.  This theorem neither identifies the
physical LRC survivor set nor proves that every physical configuration enters
this reflected family.  It does not close arbitrary `k<=1`, the final rung,
or LRC(14).

## 1. The graph carried by reflected pair overlap

For a body `E`, attach a graph `G_E` to its six labels.  An edge `a--b` is
good when an exact same-level pair-overlap floor exceeds the largest possible
six-clause singleton debt.  A repeated level on a good edge gives, by strict
one-edge Bonferroni,

```text
mu(union of the six reflected clauses)<6/7,             (2)
```

and therefore closes the THM-2941 projection-mass residual.

The exact signed-chart census checks all `45,045` labelled body edges.  Its
graph anatomy is

```text
3,001 bodies: G_E=K6;
E=(1,2,7,9,11,13): bad-edge path 7--9--11--13;
E=(2,4,7,9,11,13): bad matching {7--11,9--13}.          (3)
```

Both exceptional good graphs have chromatic number four.  Their five bad
edges are structural: an exact `593,992`-cell incidence census finds no
body-safe cell supporting the required pair overlap.  Thus this is not a
near-threshold numerical defect.

Every word using at most three level values has a repeated good edge.  This
proves all `D=0,1,2` sectors at once.  For `D=3`, only the two proper
four-colour exceptional rays survive the graph gate; an exact finite head and
the nearest-level analytic tail close both.  For `D=4`, the exceptional graph
leaves `432+312` proper words.  Exactly `72` on each exceptional body use only
four colours, so “spread four forces all five offsets” is false.  Retaining
the correct nearest-offset bound `Delta<=2`, the exact finite head and tail
close every such word.  This explains both the success and the precise limit
of the chromatic mechanism.

## 2. Spread five: transport, orientation, and the median cell

At spread five, robust good-edge theorems first remove every body with at
least fourteen robust edges, leaving `727` bodies.  The remaining word must
be a permutation of the six offsets `0,...,5`.

For labels `a,b` at levels `p,q`, the cross-determinant transport through
`(a,p)` has target slope

```text
A=p(qL-b)/(pL-a).                                      (4)
```

The reverse slope `A'` satisfies `AA'=pq`; hence at least one orientation has
the one-sided slope `>=1` required by the overlap bound.  This closes every
`D=5` packet with `m>=16`.

For `1<=m<=15`, the exact universe has

```text
727*15*6! = 7,851,600 assignments.                      (5)
```

Cross-determinant transport closes `7,835,524`.  The residual counts for
`m=1,...,15` are

```text
(12985,1511,773,427,235,53,43,23,8,9,5,2,0,0,2),       (6)
```

totalling `16,076`.  Order each body's safe cells and select its upper median.
On that single cell, one exact pair overlap exceeds the actual singleton debt
for every residual word; direct six-arc union reconstruction independently
checks `(2)`.  The selector is canonical under the fixed-point-free reflection
`j -> L-1-j`: lower and upper medians have equal singleton, pair, and union
masses.  Equations `(4)`--`(6)` therefore close all of `D=5`.

## 3. Bodywise arbitrary-level closure

The robust-edge graph, its threshold-graph complement, and low-phase clique
classification close increasingly sparse body classes without restricting
the positive levels.  The edge-ten block contributes `32` bodies and
`188,544` exact certificate rows.  Its projectively infinite complement is a
`K_(1,5)` star; a primitive five-leaf template controls the free center
because dilation strengthens the selected pair and decreases all leaf debt.

The next threshold block has `35` bodies and two creation codes, `DIIID` and
`IDIDI`.  Besides connected projective banks, it has star-`K5` and
star-`(K5-e)` cylinders.  The latter needs a new exact classification of `168`
labelled primitive words, or `14` projective orbits.  In all, `440,352`
positive-margin rows and `105` direct transport controls close that block.
The bodywise bank therefore has the exact coverage ledger

```text
covered through robust edge ten: 2,386 bodies;
covered by robust edge nine:         35 bodies;
not covered by this bank:           582 bodies;
total:                            3,003 bodies.         (7)
```

No implication in `(7)` runs backwards: the `582` bodies are outside this
assembled bodywise bank, not proved survivor bodies and not proved to lack
some different arbitrary-level certificate.

## 4. The three-D cone and its three reverse ladders

Suppose `D>=6` and `m>=3D`.  For a normalized pair, reduce its levels to
coprime `P<Q`.  Then `P>=3(Q-P)`.  The exact primitive phase fibre has the
sharp hierarchy

```text
(P,Q)=(3,4):  floor 0, uniquely;
(P,Q)=(4,5):  floor 1/70, uniquely off (3,4);
all other channels: floor at least 1/55.                (8)
```

The complete `PQ<110` bank and analytic correction tail prove `(8)`.  On the
broader `649`-body complement of the earlier edge-eleven theorem, the
ordinary nonzero-channel invoice at `D=6,m=18` closes `647` bodies.  Its
strict worst margin is

```text
1173003125846269/758117895348950710>0.                  (9)
```

The two remaining bodies use the stronger off-`(4,5)` floor and the already
proved central `(4,5)` ladder.  This deliberately retains a broader audit
universe than the `582` bodies left by `(7)`.

The zero channel can occur in `[m,m+D]` only at `m=3D`, with endpoint levels
`3D,4D`.  On each body's upper-median safe cell, retain the located primitive
phase instead of minimizing it away.  At `D=6` this closes `1,295` of the
`1,298` oriented body cases.  Located phase is scale-independent, whereas
transport loss and singleton debt decrease with `D`.

Exactly three reverse `(4D,3D)` assignments remain, on

```text
H1=(1,2,3,4,6,12),    cell 90,  labels (1,2);
H2=(1,3,4,6,8,12),    cell 174, labels (3,4);
H3=(2,4,6,8,9,12),    cell 540, labels (8,9).           (10)
```

Their exact located overlaps are

```text
30D(91D-2)/[(252D-1)(672D-1)],
14D(91D-2)/[(252D-1)(448D-1)],
   D(924D-17)/[(336D-1)(504D-1)].                       (11)
```

Each interval sweep consists of exactly `D` intersections in one arithmetic
ladder.  The three functions in `(11)` increase for `D>=1` and already exceed
their exact body debt at `D=6`.  Thus every packet with `D>=6,m>=3D` closes.
Intersecting this cone with `(7)` and the spread-five theorem gives precisely
the scoped wedge `(1)`.

## 5. Exact evidence and audit boundary

The assembly compositor pins ten terminal sources and stored transcripts:
the universal chromatic theorem, nearest-level tail, exceptional `D=3` and
`D=4` closures, `D=5` head and tail, arbitrary-level edge-ten and edge-nine
closures, and both the four-D control cone and stronger three-D cone.  It
verifies their exact conclusion tokens and rejects any change in their
LF-normalized hashes.  Retaining the weaker four-D theorem makes the sharpened
cone an audited extension rather than silently replacing its hostile controls.
The component chain additionally pins all earlier robust-edge and low-phase
inputs.

The frozen digests are

```text
assembly source: 2e242a99b8b95ab2b03f0b13f2347aaa6faf9f09e1b987c28b96ef5d25f5c906
assembly output: 1a3788f5e224d3be247104b6963cf8fe9d71af45650a70798e53dc66174c942e
semantic:        9bc1da061c14d44ef1016a9b887097666fd50ea6af242ea8302822a6dbe664b1
```

Canonical replay:

```text
python3 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
```

Fresh ordinary and optimized assembly outputs are byte-identical to the
stored output.  This audit establishes exactly the scoped conclusions above;
the residual wedge `(1)` is the next reflected sufficient-certificate target.
