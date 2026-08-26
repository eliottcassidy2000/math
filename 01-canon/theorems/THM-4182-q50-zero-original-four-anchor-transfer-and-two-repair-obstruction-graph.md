---
id: THM-4182
title: "q=50 zero-original four-anchor transfer and two-repair obstruction graph"
status: >
  PROVED RELATIVE TO THM-4150/4179 + VERIFIED-EXACT PYTHON POOL-WALL AND
  INDEPENDENT C++20 JOINT-WALL AUDITS; LRC(14) OPEN. Exactly 32 primitive
  divisor-complete four-anchors avoid the original three anchors. Twenty-four
  close at deletion depth six; eight cross sharply from tau=6 to tau=7 at
  depth seven. Their 27 labelled blockers collapse to five bodies on a
  C4-plus-leaf selector graph and have a sharp common repair-bank size two.
  The zero-original slice contains 1,138,494 distinct transferred q=50 cores.
source: codex-lrc-anchor-incidence-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
script: 04-computation/lrc14_q50_zero_original_four_anchor_transfer.py
output: 05-knowledge/results/lrc14_q50_zero_original_four_anchor_transfer.out
independent_audit_script: 04-computation/lrc14_q50_zero_original_four_anchor_transfer_independent.cpp
independent_audit_output: 05-knowledge/results/lrc14_q50_zero_original_four_anchor_transfer_independent.out
script_sha256: 1e1ec137b27b8555b2c9759a3e283746e09f445a8bd6bdba25ccd94187b05958
output_sha256: 2b1cbefeb13893840b0216b860f8925afadbdd39324d15a2e5dd34e643770d05
independent_audit_script_sha256: b6bd4ef87ff0fda05bf020dde7b1a92df3513e76e6f31fb714a9b1f6312c7eec
independent_audit_output_sha256: df68de24eb4b5f1a4510739c80123d921d7e2fe7f7a1ce4fe6c967bac6f12b06
independent_geometry_dependency_sha256: 47c0d8f1107d35b770d63348a21fcac04ffc08dc4a611b46b221f0be0ac5aa3b
hash_basis: raw LF bytes
primary_audit: >
  PASS. Frozen THM-4178 pool-wall failure atoms are subset-zeta summed over
  every one of the C(30,6)+C(30,7)=2,629,575 global deletion masks. A new
  Python-big-int incidence recursion searches every four-anchor row through
  cover budget six, enumerates all 27 blockers, verifies the two-repair bank
  and all eight depth-seven covers, and counts the body union by 3,230,304
  anchor presentations. Normal, optimized, and hash-seeded outputs byte-match.
independent_audit: >
  PASS. A warning-clean signed-integer C++20 path independently inserts q=50
  into the joint wall arrangement, reconstructs both global repair layers,
  searches covers with a distinct frequency/packing recursion, independently
  enumerates all blockers, and counts the union by scanning all C(27,10)
  zero-original bodies rather than presentations. O2, O0, and UBSan streams
  byte-match. The paths share only the fixed pool, threshold, and elementary
  blocker-descent theorem; their wall geometries, arithmetic representations,
  cover searches, and body-census directions differ.
---

# THM-4182 -- zero-original four-anchor transfer and obstruction graph

**PROVED RELATIVE TO THM-4150/4179 + VERIFIED-EXACT; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the fixed pool, original anchor, and newcomer

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143},                 U=P\A_0,             q=50.       (1)
```

Call `A in binom(U,4)` a **zero-original four-anchor** when

```text
gcd(A)=1
and for every 2<=m<=14 some a in A satisfies m|a.                 (2)
```

There are exactly `32` such anchors. They are the following disjoint
parametric rows, with `x` ranging over the displayed set:

```text
{63,168,x,286},      x in {10,20,30,40,60,80,170,190,240,290};
{63,84,x,286},       x in {40,80,240};
{63,252,x,286},      x in {40,80,240};
{15,252,x,286},      x in {40,80,240};
{40,252,x,286},      x in {85,95,145,193};
{80,252,x,286},      x in {85,95,145,193};
{240,252,x,286},     x in {85,95,145,193};
{42,63,240,286}.                                               (3)
```

In particular, every one contains `286`. For such an `A`, put

```text
V_A=P\A,
E_d^A={R in binom(V_A,d):mu(G_((P union {50})\R))>=4/63}.        (4)
```

> **Four-anchor transfer theorem.** For 24 of the anchors in `(3)`,
>
> ```text
> tau(E_6^A)>6.                                                  (5)
> ```
>
> The remaining eight anchors are
>
> ```text
> L_x={80,x,252,286},      R_x={x,240,252,286},
> x in X={85,95,145,193}.                                      (6)
> ```
>
> For every one of these rows,
>
> ```text
> tau(E_6^A)=6,                  tau(E_7^A)=7.                  (7)
> ```
>
> Consequently, for every anchor `A` in `(3)`, every
> `K in binom(V_A,6)`, every positive integer `c`, and every two distinct
> positive odd integers `a,b`, there is a point `y in R/Z` such that, for
>
> ```text
> H_(A,K)=A union {50} union K,
> ```
>
> one has
>
> ```text
> min_(v in 2cH_(A,K) union {a,b}) ||vy|| >= 1/14.              (8)
> ```

Restricting `K` to `binom(U\A,6)` gives bodies containing none of the three
original anchors. The `32 binom(23,6)=3,230,304` anchor presentations collapse
to exactly

```text
1,138,494 distinct zero-original q=50 core bodies.              (9)
```

The complete anchor-presentation multiplicity histogram is

```text
1:262761, 2:332137, 3:204900, 4:178178, 5:59560,
6:66399, 7:13188, 8:11229, 9:6149, 10:2666,
11:658, 12:493, 13:140, 15:36.                                 (10)
```

There is no multiplicity-14 body. Equation `(9)` is the union of bodies
containing at least one anchor in `(3)`, not a classification of all
zero-original divisor-complete bodies.

## 2. Inheritance and the body-level coordinate

The closest mechanism is
[THM-4179](THM-4179-q50-seventh-deletion-primitive-anchor-completion.md):
one-step blocker descent says that, for a monotone deletion statistic on an
`n`-vertex ground set, every `k`-cover of `E_(d+1)` is a `k`-cover of `E_d`
when `d+k<n`. Here

```text
n=|V_A|=26,       d=6,       k=6,       d+k=12<26.             (11)
```

Thus every cover of `E_7^A` through size six must already occur in the
complete blocker list of `E_6^A`. This lets the depth-seven proof inspect only
the surviving blockers, rather than repeat a blind transversal search.

The canonical hostiles are the eight rows `(6)`. The corrected near miss is
again that a shallow repair cover is a body obstruction. It is not: the
27 labelled depth-six blockers below are only presentations of five distinct
body states, and each state has many depth-seven repairs. The least-used
sidecar is therefore the map

```text
(anchor A, blocker K) -> body S=A union K -> repairs disjoint from S.       (12)
```

It preserves exactly the coordinate on which repair existence depends and
forgets the redundant anchor presentation only after its role in the repair
restriction has been checked.

The source-to-target contract is

```text
source:       q=50 labelled d=6,7 repair layers on all thirty pool labels
target:       every four-anchor/six-choice body and odd-tail transfer
map:          anchor restriction -> exact blocker -> body collapse
              -> d=7 missed repair -> THM-4150
preserved:    q, body labels, threshold 4/63, divisor pins, content
destroyed:    component order, canonical safe phase, selected repair after use
sidecar:      anchor presentation until blocker enumeration, then body state
positive:     24 rows have tau(E6)>6; eight have tau(E7)=7
hostile:      the eight rows have 27 exact d6 blockers
decisive test: all 2,629,575 global d=6,7 masks and two exact cover paths. (13)
```

## 3. Exact anchor and repair-layer classifications

Since `143` is excluded and `286` is the only remaining pool label divisible
by `13`, condition `(2)` forces `286`. Exhausting the remaining
`binom(26,3)=2,600` triples attached to it, or equivalently all
`binom(27,4)=17,550` four-subsets of `U`, against the thirteen divisibility
predicates and the gcd gives exactly `(3)`. Both audits perform this literal
classification.

The primary audit retains the THM-4178 pool-wall arrangement

```text
L=18,241,159,416,480,       7,134 walls,       7,133 cells.     (14)
```

It integrates the `q=50` safe comb inside each cell and subset-zeta sums the
complete labelled failure atoms. Globally, before restricting any anchor, the
two repair layers are

| depth `d` | candidates | repairs | threshold equalities | minimum edge delta | maximum nonedge delta |
|---:|---:|---:|---:|---:|---:|
|6|`593,775`|`85,324`|`0`|`348,508,440`|`-2,639,430,360`|
|7|`2,035,800`|`821,737`|`0`|`348,508,440`|`-402,312,960`|

Here the integer comparison is the THM-4178 normalization

```text
9N(R)-8qL >= 0.                                               (15)
```

Filtering the global edges by `R intersection A=empty`, and then exhaustively
searching every row through cover budget six, gives exactly the `24/8` split
in `(5)--(6)`. The full artifact prints all 32 row edge counts and search
levels. Every cover found in the eight residual rows has size exactly six;
there is no smaller cover.

## 4. All 27 depth-six blockers

The blocker list admits a compact exact description. Put

```text
X={85,95,145,193},       Q_L={88,168,240},       Q_R={88,168}. (16)
```

For the left rows `L_x={80,x,252,286}`:

```text
x=85:
  K=Q_L union (X\{85});                                      (17)

x in X\{85}:
  K=Q_L union (X\{85,x}) union {t},
  t in {8,85,170}.                                           (18)
```

These are `1+3*3=10` blocker presentations. For the right rows
`R_x={x,240,252,286}`:

```text
x=85:
  K=Q_R union (X\{85}) union {t},
  t in {16,80};                                              (19)

x in X\{85}:
  K=Q_R union (X\{85,x}) union T,
  T in {{8,80},{16,85},{16,170},{80,85},{80,170}}.           (20)
```

These are `2+3*5=17` more. Equations `(17)--(20)` therefore give all
`27`, not selected witnesses. Direct substitution reproduces the exact
per-row blocker counts

```text
L_85:1;   L_95,L_145,L_193:3 each;
R_85:2;   R_95,R_145,R_193:5 each.                            (21)
```

The primary first-uncovered-edge recursion and the independent C++
frequency/packing decision followed by a separately written exhaustive
enumerator agree on every labelled set in `(17)--(20)`.

## 5. The five-state obstruction graph and sharp repair bank

Now pass from presentations to bodies. Define the common eight-label spine

```text
D={88,95,145,168,193,240,252,286}                           (22)
```

and the selector graph

```text
Gamma={{8,80},{16,85},{16,170},{80,85},{80,170}}.           (23)
```

The last four edges form `K_(2,2)` on
`{16,80} x {85,170}` and `{8,80}` is a leaf edge. The 27 presentations in
Section 4 collapse **exactly** to the five bodies `D union e`, `e in Gamma`:

| selector edge `e` | anchor presentations | depth-seven repairs disjoint from `D union e` |
|:---|---:|---:|
|`{8,80}`|`6`|`243`|
|`{16,85}`|`4`|`439`|
|`{16,170}`|`3`|`736`|
|`{80,85}`|`8`|`278`|
|`{80,170}`|`6`|`514`|

Thus the blocker multiplicities sum to `6+4+3+8+6=27`, while the repair
counts are body invariants independent of which anchor presentation led to
the body.

Index the five table rows by `0,...,4`. The complete nonempty incidence
patterns of a depth-seven repair against these bodies are

| bodies repaired | number of repairs |
|:---|---:|
|`{1,2,3,4}`|`27`|
|`{0,2,4}`|`5`|
|`{1,2}`|`284`|
|`{1,3}`|`11`|
|`{2,4}`|`91`|
|`{3,4}`|`172`|
|`{0}`|`238`|
|`{1}`|`117`|
|`{2}`|`329`|
|`{3}`|`68`|
|`{4}`|`219`|

There is no `{0,1,2,3,4}` row: no single depth-seven repair is disjoint from
all five obstruction bodies. The minimum simultaneous repair bank therefore
has size at least two. It has size exactly two, witnessed by

```text
R_*={8,10,15,20,143,176,290},
R_0={10,15,20,85,120,143,176}.                            (24)
```

The first repairs bodies `{1,2,3,4}`, while the second repairs `{0,2,4}`.
Their exact safe masses and strict margins above `4/63` are

```text
mu(G_((P union {50})\R_*))
  =2912475704927/45602898541200,
margin=1894841703/5066988726800;

mu(G_((P union {50})\R_0))
  =2566056877/40320865200,
margin=18005831/120962595600.                            (25)
```

This is the mechanism behind the depth-seven crossing. Every size-six cover
of an `E_7^A` row would, by `(11)`, be one of the depth-six blockers. Its body
is one of the five rows above, so one member of `(24)` is an active
depth-seven repair disjoint from that blocker. This contradiction proves

```text
tau(E_7^A)>6                                             (26)
```

for every residual anchor. The following exact seven-covers prove the reverse
bounds and hence `(7)`:

| anchor | one cover of `E_7^A` | `|E_7^A|` |
|:---|:---|---:|
|`{80,85,252,286}`|`{8,88,95,145,168,193,240}`|`202,280`|
|`{80,95,252,286}`|`{8,85,88,145,168,193,240}`|`175,101`|
|`{80,145,252,286}`|`{8,85,88,95,168,193,240}`|`192,626`|
|`{80,193,252,286}`|`{8,85,88,95,145,168,240}`|`178,835`|
|`{85,240,252,286}`|`{8,80,88,95,145,168,193}`|`194,302`|
|`{95,240,252,286}`|`{8,80,85,88,145,168,193}`|`166,779`|
|`{145,240,252,286}`|`{8,80,85,88,95,168,193}`|`180,334`|
|`{193,240,252,286}`|`{8,80,85,88,95,145,168}`|`168,585`|

Depth seven is sharp for these eight rows because `(21)` proves
`tau(E_6^A)=6`.

## 6. Haar transfer and body union

Fix `A` and `K` as in Section 1. In a row satisfying `(5)`, `K` misses some
six-repair. In a residual row, `(7)` says it misses some seven-repair. In
either case there is an `R subset V_A\K` with

```text
mu(G_((P union {50})\R))>=4/63.                         (27)
```

Because `R` misses both `A` and `K`,

```text
H_(A,K) subset (P union {50})\R,
G_((P union {50})\R) subset G_(H_(A,K)).               (28)
```

Thus `mu(G_(H_(A,K)))>=4/63`. Every body is primitive and divisor-complete
because its anchor is, and `50 notin P` makes its eleven labels distinct.
The circle endomorphism `z -> cz` preserves Haar measure of the safe-set
pullback. THM-4150 then proves `(8)` for every content and every distinct
positive odd tail pair; doubling separates body labels from the odd tails.

For the zero-original slice, each anchor has exactly `23` eligible remaining
labels and therefore `binom(23,6)=100,947` presentations. The primary audit
inserts all `3,230,304` presentation masks into a multiplicity map. The
independent audit instead scans all `binom(27,10)=8,436,285` zero-original
ten-subsets and counts which of the 32 anchors they contain. Both give
`(9)--(10)` exactly.

The new union is disjoint by original-anchor count from the fixed-q=50 slices
of THM-4174, THM-4175, and THM-4179: those contain respectively three, two,
and exactly one member of `A_0`, while this one contains zero. Consequently
these four declared slices contain exactly

```text
888,030+6,660,225+1,071,961+1,138,494
  =9,758,710                                               (29)
```

distinct scale-one q=50 cores. This is an overlap-safe count for these four
named mechanisms, not a global novelty claim against every canonical family.

## 7. Independent audits and replay

The primary Python path hash-pins the THM-4178 pool-wall dependency. It builds
both global repair layers once, uses Python integers as edge-incidence bitsets
to enumerate covers through budget six, and counts bodies from anchor
presentations.

The independent C++20 source includes the frozen THM-4179 joint-geometry
implementation, whose declared dependency hash is recorded in the
frontmatter. It inserts the newcomer walls into

```text
D=91,205,797,082,400,       7,214 walls,       7,213 atoms,      (30)
```

with labelled-atom FNV-1a-64 control `0ccd305ae47c79ea`. It then reconstructs
the new global layers by signed-integer subset sums, decides covers with a
frequency-ordered matching-pruned recursion, separately enumerates every
surviving blocker by first uncovered edges, and performs the reverse-direction
`binom(27,10)` body census. The two paths agree on all anchor labels, global
and row edge counts, zero equalities, all 27 blockers, five body states,
eleven witness-incidence patterns, eight exact seven-covers, and the complete
body histogram.

Primary replay:

```bash
python3 -B 04-computation/lrc14_q50_zero_original_four_anchor_transfer.py
python3 -B -O 04-computation/lrc14_q50_zero_original_four_anchor_transfer.py
PYTHONHASHSEED=4182 python3 -B \
  04-computation/lrc14_q50_zero_original_four_anchor_transfer.py
```

Independent replay:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_zero_original_four_anchor_transfer_independent.cpp \
  -o /tmp/lrc4182-o2
/tmp/lrc4182-o2

g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_zero_original_four_anchor_transfer_independent.cpp \
  -o /tmp/lrc4182-o0
/tmp/lrc4182-o0

g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined \
  04-computation/lrc14_q50_zero_original_four_anchor_transfer_independent.cpp \
  -o /tmp/lrc4182-ubsan
/tmp/lrc4182-ubsan
```

All three Python streams and all three C++ streams byte-match their respective
frozen outputs.

## 8. Scope and frontier

This theorem is `q=50`-only and fixed-pool-only. It treats primitive
divisor-complete four-anchors, doubled bodies, and distinct positive odd tails.
The transfer conclusion applies to every `K in binom(P\A,6)`, but the exact
new union `(9)` deliberately counts only `K subset U\A`, so “zero-original”
is a body predicate rather than an anchor-only slogan.

It does **not** say every zero-original divisor-complete body contains one of
the 32 anchors, prove an all-`q` four-anchor theorem, provide physical entry,
classify arbitrary bodies, handle mixed/even tails, prove necessity of the
Haar threshold, or prove LRC(14). **QED.**
