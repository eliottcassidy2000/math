---
id: THM-857
title: Scale-one Hamming-six exact component closure
status: PROVED + FINITE-EXACT — every proper scale-one Hamming-six packet is strictly loose except the doubled AP equality 2[12]; 580,919,164-node logical census with independent closed-danger replay
source: codex-2026-07-15-S10 H6 continuation
depends_on: [THM-815]
related: [THM-770, THM-820, THM-845, THM-856, THM-858, HYP-6820]
verification:
  - 04-computation/lrc13_scale_one_hamming_six_component_recursion_codex_S10.cpp
  - 05-knowledge/results/lrc13_scale_one_hamming_six_component_recursion_codex_S10.out
  - 04-computation/lrc13_hamming_six_closed_danger_union_replay_codex_S10.cpp
---

# THM-857 — the scale-one Hamming-six chart closes

For a six-set `R subset [12]`, put

```text
P=[12] minus R,
W=P union {r+13h_r:r in R},                    h_r>=1.   (1)
```

Thus `W` is a proper residue-preserving scale-one Hamming-six lift of
`[12]`.  For a finite positive speed set `Q`, write

```text
E(Q)={t in R/Z: ||qt||>1/13 for every q in Q}.           (2)
```

The open set `E(Q)` is represented as its exact sorted union of rational
components.  At every materialized node, the stored intervals are exactly the
connected components of the prefix strict-safe set, by induction under exact
intersection with the next speed's open safe teeth.  Thus `E(W)` is nonempty
exactly when `M(W)>1/13`, while `E(W)` is empty exactly when the closed danger
combs cover the circle.

## Theorem

Over all `C(12,6)=924` choices of `R` and all unbounded proper heights in
(1), the exact recursion below has one and only one covering terminal:

```text
R={1,3,5,7,9,11},
(r:r+13h_r)=(1:14,3:16,5:18,7:20,9:22,11:24),
W={2,4,...,24}=2[12].                                  (3)
```

Multiplication by two is surjective on `R/Z`, so
`M(2[12])=M([12])=1/13`; equivalently, the strict-safe set in (2) is empty
and the usual AP witness attains the boundary.  Every other packet (1),
primitive or nonprimitive, has

```text
M(W)>1/13.                                               (4)
```

In particular, the `923` primitive retained-core rows have zero covering
terminals.  The exceptional missing-odd row also has no primitive
mixed-height cover; its only cover is the all-even equality (3).  This closes
the entire proper **scale-one** Hamming-six chart.

## 1. Complete finite recursion

Numerically order the six proper lifts, with `x_0=0` for the empty prefix:

```text
x_1<...<x_6.                                             (5)
```

At a depth-`j` prefix let `L` be the length of its longest component of
`E(P union {x_1,...,x_j})`, and put `r=6-j`.  A depth-`j<6` prefix has
`6+j<=11` speeds.  The settled lower-runner theorem gives
`M>=1/(7+j)>=1/12>1/13`, hence a positive-length strict-safe component.
THM-815's sharp one-component danger-discrepancy inequality
then says that every tight completion obeys

```text
x_(j+1) <= floor(22r/[13(13-2r)L]).                     (6)
```

The denominator is positive for every `1<=r<=6`.  If a completed packet
covers, its remaining `r` danger combs cover every longest prefix component,
and all their speeds are at least `x_(j+1)`.  Therefore (6) applies to the
actual next lift, which is consequently enumerated.

For each unused label the
recursion enumerates every congruent proper lift above `x_j` and at most (6),
then sorts these candidates by speed.  Distinct labels give distinct residues
modulo thirteen, so every completed packet has one unique numerical order.
Induction on `j` proves that the tree contains every possible terminal with
empty strict-safe set; there is no height cutoff.

THM-815 recorded `194,040` first-prefix states by using the worst root cap
`468` uniformly.  Applying (6) to each deletion core's own longest component
reduces the exact first layer to

```text
83,881.                                                  (7)
```

## 2. Two theorem-safe shortcuts

The logical tree is large, but most children have an immediate exact
certificate.

### 2.1 A complete safe tooth

The open safe teeth of a candidate speed `u` are

```text
B_k(u)=((13k+1)/(13u),(13k+12)/(13u)),       0<=k<u.    (8)
```

Suppose at parent depth `j>=2` one component contains a whole `B_k(u)`.
After adjoining `u`, the child retains an interval of length `11/(13u)`.
There are `r=5-j<=3` later combs.  Applying (6) to that interval bounds the
next speed by

```text
2ru/(13-2r) <= 6u/7 < u,                                (9)
```

contradicting (5).  At `j=5`, the retained tooth directly proves that the
terminal is loose.  Endpoint equality is allowed in the containment test
because both safe intervals are open.

### 2.2 Streaming the child cap

Otherwise intersect the parent components with the safe teeth of `u` in
endpoint order.  Let `y_min` be the least proper congruent lift above `u`
among all unused labels, and let `r=5-j>0`.  As soon as an emerging child
piece `J` satisfies

```text
22r/[13(13-2r)|J|] < y_min,                              (10)
```

the child cannot have an outgoing tight edge.  Indeed its actual longest
component has length at least `|J|`, so its cap is no larger than the left
side of (10).  At the last insertion, the first nonempty child piece certifies
looseness; if no piece occurs, the strict-safe set is empty and the covering
terminal is retained.  All comparisons in (8)--(10) are exact integer cross
products.

For either nonterminal shortcut, the certified child's exact next cap is
below its least legal future lift, so the unshortened recursion has zero
outgoing edges there.  At depth six the shortcut instead certifies a loose
terminal.  Recording the child and its incoming edge therefore preserves the
complete logical node and edge census; only the certified child is not
materialized further.

## 3. Frozen exact census

The complete logical tree is:

| depth | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| nodes | 924 | 83,881 | 8,906,315 | 559,202,706 | 12,671,505 | 53,812 | 21 |
| dead / no possible tight continuation | 0 | 0 | 0 | 555,565,824 | 12,638,291 | 53,792 | 0 |
| complete-tooth certificates | 0 | 0 | 0 | 495,797,163 | 0 | 0 | 0 |
| streaming-cap certificates | 0 | 0 | 0 | 59,768,661 | 12,638,291 | 53,792 | 20 |

Within this finite cap-admissible logical tree there are exactly

```text
logical nodes              580,919,164
candidate edges            580,918,240
covering terminals                   1
loose terminals                     20.                 (11)
```

The materialized nonterminal prefixes, excluding shortcut-certified
children, have the following exact extrema:

| depth | prefixes | minimum longest component | multiplicity | maximum next cap |
|---:|---:|---:|---:|---:|
| 0 | 924 | `31/1430` | 1 | 468 |
| 1 | 83,881 | `11/6071` | 1 | 1,556 |
| 2 | 8,906,315 | `11/20215` | 1 | 2,488 |
| 3 | 3,636,882 | `2/5941` | 1 | 2,154 |
| 4 | 33,214 | `2251/2834949` | 1 | 473 |
| 5 | 20 | `1/325` | 1 | 50 |

The primary source asserts every number in (7), (11), both tables, and the
literal row (3).  Its rational subtraction first cancels the common
denominator factor, and all cap and containment decisions use signed
128-bit cross products.

## 4. Independent reconstruction

The companion replay shares no interval-intersection implementation with the
primary source.  At every expanded prefix it instead:

1. constructs every **closed** danger tooth from all prefix speeds;
2. sorts and merges touching closed intervals;
3. takes the open complementary gaps;
4. independently subtracts the next closed comb; and
5. checks the resulting gaps endpoint-for-endpoint against a reconstruction
   from the full child prefix.

For every complete-tooth or streaming-cap shortcut it freezes the parent-gap
index, tooth index, exact witness interval, remaining-comb count, least future
speed, and cap.  Expanded nodes, terminal verdicts, and shortcut witnesses
feed a canonical SHA-256 stream per root; the `924` root hashes are combined
in root-index order.  This makes the certificate sensitive to the exact
labelled tree without storing its hundreds of millions of logical nodes as
text.

The independent replay and frozen source/output hashes are recorded after the
reproduction commands below.  [Certificate manifest pending final all-root
replay; no mathematical claim in this theorem depends on an unrecorded hash.]

## 5. Tournament analysis and challenged vertices

At each deletion root, take the six missing labels as vertices.  The raw
gauge orders vertices by their least proper lifts `r+13`.  The
core-conditioned gauge orders them by the exact five-comb cap left after
adjoining that least lift; ties use increasing label.  Both gauges are scalar
rankings and hence give the same transitive fingerprint on all `924` roots:

```text
score histogram                 {0:1,1:1,2:1,3:1,4:1,5:1}
directed triangles                                      0
SCC-size histogram                                  {1:6}
Hamiltonian paths                                        1
conditioned ties                                       406
raw/conditioned edge flips                           4,500
flip histogram  {0:47,1:61,2:93,3:115,4:146,5:108,6:110,
                 7:78,8:49,9:46,10:33,11:18,12:16,13:4}. (12)
```

The switch moves thousands of edges but cannot distinguish the exceptional
root containing the unique cover from the `923` primitive retained-core roots,
each of which has zero covers.  The pairwise planning order is not the cover
predicate.

Alternate vertices explain why earlier tournament threads were still useful.
THM-815(C.2) uses oriented AP cusps as owner vertices and provider-start ratios
as directed edge weights.  A pair may carry both arrows or neither, so the
object is not a tournament; its min-plus cycle products close fourteen roots.
THM-856 uses remaining combs as vertices but retains reduced ratio, common
scale, the mod-thirteen pair sawtooth, and endpoint discrepancy on every edge;
its max-plus spanning-tree character crosses the seven-comb first-moment wall.
The exact H6 recursion subsumes the former scale-one partial closure, while
the latter remains the relevant radius-seven lead.

The predicate-preserving carrier here is

```text
(literal strict-safe components,
 remaining labelled residue progressions,
 last speed,
 exact shortcut witness).                              (13)
```

Using components, AP cusps, boundaries, or proof obligations as vertices can
expose a useful cycle or tree algebra.  Quotienting (13) to bare runners,
root caps, antipodal-pair counts, or a binary tournament destroys exact
component coverage.

## Scope guardrail

This theorem closes proper residue-preserving Hamming-six lifts at common
scale one.  It does not transport the result across arbitrary AP scale or
arbitrary common sheets, classify the finite Hamming-five bank left by
THM-858, cross the seven-comb wall, or prove global `n=12` sporadic-branch
emptiness.

Reproduce the primary census with

```bash
c++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc13_scale_one_hamming_six_component_recursion_codex_S10.cpp \
  -o /tmp/thm857-primary
/tmp/thm857-primary > /tmp/thm857-primary.out
cmp /tmp/thm857-primary.out \
  05-knowledge/results/lrc13_scale_one_hamming_six_component_recursion_codex_S10.out
```
