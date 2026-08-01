---
id: THM-3032
title: "Reflected spread-five closure and four-D cone"
status: >
  PROVED + VERIFIED-EXACT + SCOPE-AUDITED.  In the reflected k=1 sufficient
  residual family of THM-2941, every packet of spread D<=5 closes at every
  positive minimum level; 2,386 of the 3,003 six-label bodies close for
  arbitrary positive levels; and every packet with D>=6 and m>=4D closes.
  Thus certificate failure is confined to 617 bodies with D>=6 and
  1<=m<4D.  This is not a physical-survivor classification and not LRC(14).
source: codex-reflected-chain-recover-2026-08-01
audit: >
  The promotion audit checked the reflected-family typing, the three-way
  logical intersection, all body and spread partitions, strict overlap/debt
  directions, the two exceptional chromatic graphs, the D=5 head/tail seam,
  the unique hostile (4,5) cone orientation, and the remaining-wedge
  quantifiers.  The assembly pins LF-normalized sources and outputs for all
  eight terminal inputs and requires their ordinary/-O and exact-control
  tokens.  Fresh ordinary and optimized assembly replays are byte-identical
  to each other and to the stored transcript; all 17 added Python sources
  compile and contain no optimization-stripped assertions.
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1226-gcd-period-projective-charge-obstruction
related:
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
script: 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
output: 05-knowledge/results/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.out
script_sha256: 1e52808545a4bdf747d338400170e006a82d5300bac95b93c3aa7c274fcca93b
output_sha256: 2274166ec43bcf517a318cae0e6b46be3e16a93453cfb3d73f30ff618768b0ab
semantic_sha256: 077c8a3df606eec6ad68ee4f1826b1cb69851a7f806a6f16566cd5ce39d7e4ca
hash_basis: LF-normalized bytes
---

# THM-3032 -- reflected spread-five closure and four-D cone

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
2. on exactly `2,386` of the `3,003` bodies, every positive level assignment
   closes, with no spread restriction; and
3. every packet with `D>=6` and `m>=4D` closes on all `3,003` bodies.

Intersecting these statements confines possible **certificate failures** to

```text
617 bodies,   D>=6,   1<=m<4D.                         (1)
```

The word “certificate” is essential.  This theorem neither identifies the
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
the positive levels.  The final edge-ten block contains `32` bodies and
`188,544` exact certificate rows.  Its only projectively infinite complement
is a `K_(1,5)` star; a primitive five-leaf template controls the free center
because dilation strengthens the selected pair and decreases all leaf debt.

The cumulative exact partition is

```text
arbitrary-level closed: 2,386 bodies;
remaining:                617 bodies;
total:                   3,003 bodies.                  (7)
```

No implication in `(7)` runs backwards: the `617` bodies are failures of the
assembled sufficient certificates, not proved survivor bodies.

## 4. The four-D cone and its unique hostile channel

Suppose `D>=6` and `m>=4D`.  For the closest normalized pair, reduce its
levels to coprime `P<Q`.  Then `P>=4(Q-P)`, and the exact phase fibre satisfies

```text
F_(P,Q)>=1/70.                                          (8)
```

Equality occurs only at `(P,Q)=(4,5)`; every other channel has the stronger
floor `1/55`.  The finite small-channel bank and the analytic correction tail
prove this uniformly.

For every body except

```text
H=(1,2,3,4,6,12),                                      (9)
```

the closest-pair label ratio is at most `1/168`.  At the boundary `D=6,
m=24`, the phase floor minus transport loss and singleton debt is already
strictly positive; both losses decrease as `D` or `m` increases.

On `H`, labels `1,2` close every non-`(4,5)` channel using the `1/55` floor.
The equality channel fits in the interval `[m,m+D]` only at `m=4D` with
levels `4D,5D`.  The forward orientation closes directly.  In the sole
reverse exception, a located certificate on safe cell `j=90` has exact pair
intersection

```text
I_D=36D(322D-1)/[(840D-1)(672D-2)]
   >23/1120>1/70,                                      (10)
```

while total singleton debt is at most `1/(168D-3)<1/70`.  This closes the
last cone packet and proves statement 3.

## 5. Exact evidence and audit boundary

The assembly compositor pins eight terminal sources and stored transcripts:
the universal chromatic theorem, nearest-level tail, exceptional `D=3` and
`D=4` closures, `D=5` head and tail, arbitrary-level edge-ten closure, and
the four-D cone.  It verifies their exact conclusion tokens and rejects any
change in their LF-normalized hashes.  The component chain additionally pins
all earlier robust-edge and low-phase inputs.

The frozen digests are

```text
assembly source: 1e52808545a4bdf747d338400170e006a82d5300bac95b93c3aa7c274fcca93b
assembly output: 2274166ec43bcf517a318cae0e6b46be3e16a93453cfb3d73f30ff618768b0ab
semantic:        077c8a3df606eec6ad68ee4f1826b1cb69851a7f806a6f16566cd5ce39d7e4ca
```

Canonical replay:

```text
python3 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_d0_d5_wedge_assembly_thm2941.py
```

Fresh ordinary and optimized assembly outputs are byte-identical to the
stored output.  This audit establishes exactly the scoped conclusions above;
the residual wedge `(1)` is the next reflected sufficient-certificate target.
