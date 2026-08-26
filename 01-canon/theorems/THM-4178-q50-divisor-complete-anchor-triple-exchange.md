---
id: THM-4178
title: "q=50 primitive-anchor exchange profile and genuine two-anchor Haar transfer"
status: >
  VERIFIED-EXACT q=50 SIX-ANCHOR d=5,6 TRANSVERSAL CLASSIFICATION +
  INDEPENDENT C++20 JOINT-WALL AUDIT; PROVED RELATIVE TO THM-4150 FOR ONE
  GENUINE TWO-ANCHOR 480,700-CORE ODD-TAIL FAMILY; THM-4175 TOMOGRAPHY
  SPECIALIZED; LRC(14) OPEN.
source: codex-open-frontiers-anchor-exchange-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
related:
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4172-multideletion-support-tomography-and-same-parity-johnson-holonomy
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
script: 04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py
output: 05-knowledge/results/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.out
independent_audit_script: 04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178_independent_audit.out
script_sha256: ff9748f07cbc1bab20832860ff5f52bca59573262725a4e69c5b42ab2d078251
output_sha256: b5411413f063c3eb9ad22ff2b69e7795c2b744c03de54c2871573f51c7439a13
independent_audit_script_sha256: 96483771eb2cc38b89e765ccc6d1d8e03fc2e983a279ecf94099c20d32fae806
independent_audit_output_sha256: 1f082c7e1d24ef267ec9adf1383c9e6e1c235631a581f02409ba5299fc657edd
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction pool-wall geometry retains all 30 labelled failure
  masks and integrates the q=50 safe comb by signed-integer prefix
  differences. Exhaustive subset-zeta sums and exact cover searches through
  budget eight reproduce all twelve d=5,6 rows. Normal, optimized, and
  hash-seeded Python streams byte-match.
independent_audit: >
  PASS. A warning-clean C++20 path instead builds the joint P-union-{50}
  arrangement on denominator 91,205,797,082,400, with 7,214 walls and 7,213
  atoms, sums their lengths directly, and independently searches every cover
  budget zero through eight. O0 and O2 outputs byte-match. All edge counts,
  equality counts, and exact transversal numbers agree; minimum covers need
  not be unique. The labelled-atom FNV-1a-64 control is 0ccd305ae47c79ea.
---

# THM-4178 -- q=50 primitive-anchor exchange profile and genuine two-anchor Haar transfer

**VERIFIED-EXACT q=50 CLASSIFICATION + PROVED RELATIVE TRANSFER +
INDEPENDENTLY AUDITED; THM-4175 TOMOGRAPHY SPECIALIZED; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the THM-4156 pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                      (1)
```

For a finite positive set `S`, write

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.              (2)
```

For `y in R/Z`, define the complete labelled failure set of the pool by

```text
F_P(y)={p in P:||py||<1/14}.                          (3)
```

For a positive integer `q notin P` and every `F subset P`, define its
`q`-safe failure-atom weight

```text
w_q(F)=mu({y:F_P(y)=F and ||qy||>=1/14}).             (4)
```

For `R subset P`, put

```text
M_q(R)=mu(G_((P union {q})\R)).                       (5)
```

> **THM-4175 specialization.** For every `q notin P` and every
> `R subset P`,
>
> ```text
> M_q(R)=sum_(F subset R) w_q(F).                     (6)
> ```
>
> Consequently the complete deletion deck and the complete labelled atom
> vector determine one another by Boolean-lattice inversion:
>
> ```text
> w_q(F)=sum_(E subset F)(-1)^(|F|-|E|)M_q(E).        (7)
> ```

Call a three-set `A subset P` a **primitive divisor-complete anchor** when

```text
gcd(A)=1
and for every 2<=m<=14 some a in A satisfies m|a.      (8)
```

There are exactly six such anchors:

```text
Acal={
 {40,143,252}, {80,143,252}, {120,126,143},
 {120,143,252}, {126,143,240}, {143,240,252}
}.                                                     (9)
```

For `A in Acal`, set `O_A=P\A`. For `d>=0`, define the `q=50`
deletion-repair hypergraph

```text
E_d^A={R in binom(O_A,d):M_50(R)>=4/63}.              (10)
```

The exact `d=5,6` ledger is:

| anchor `A` | `|E_5^A|` | `tau(E_5^A)` | one minimum cover of `E_5^A` | `|E_6^A|` | `tau(E_6^A)` | one minimum cover of `E_6^A` |
|:---|---:|---:|:---|---:|---:|:---|
| `{40,143,252}` | 1,962 | 5 | `{95,168,193,240,290}` | 44,562 | 8 | `{16,85,88,95,145,168,193,240}` |
| `{80,143,252}` | 1,398 | 5 | `{85,95,145,168,193}` | 33,636 | 7 | `{16,88,95,145,168,193,240}` |
| `{120,126,143}` | 2,063 | 5 | `{95,168,193,240,290}` | 46,261 | 8 | `{16,85,88,95,168,193,240,290}` |
| `{120,143,252}` | 1,664 | 5 | `{95,168,193,240,290}` | 37,637 | 8 | `{8,80,85,88,95,145,168,193}` |
| `{126,143,240}` | 1,558 | 4 | `{95,168,193,290}` | 38,915 | 8 | `{8,80,85,88,95,145,168,193}` |
| `{143,240,252}` | 1,234 | 4 | `{95,168,193,290}` | 31,234 | 7 | `{16,85,88,95,145,168,193}` |

Every threshold-equality count in the table is zero. Thus the five
nonoriginal anchors split exactly as

```text
A_+={{40,143,252},{120,143,252},{126,143,240}},
A_blk={{80,143,252},{143,240,252}},                    (11)
```

where

```text
A in A_+  iff tau(E_6^A)=8,
A in A_blk iff tau(E_6^A)=7,                          (12)
```

for `A in Acal minus {{120,126,143}}`.

For every `A in A_+`, every `K in binom(O_A,7)`, every positive integer
`c`, and every two distinct positive odd integers `a,b`, put

```text
H_(A,K)=A union {50} union K.                          (13)
```

Then there exists `x in R/Z` such that

```text
min_(v in 2cH_(A,K) union {a,b})||vx||>=1/14.          (14)
```

For each fixed `A in A_+`, this gives exactly

```text
binom(27,7)=888,030                                    (15)
```

`K`-labelled core bodies. Equation `(15)` is a within-anchor count. Bodies
from different anchor rows may overlap, so the three counts must not be added
without a separate overlap census.

Relative to the original anchor `A_0={120,126,143}`, the rows
`{120,143,252}` and `{126,143,240}` omit one original anchor. For either row,
the `888,030` choices split as

```text
binom(26,7)=657,800  one-anchor bodies from THM-4175,
binom(26,6)=230,230  all-anchor bodies from THM-4174. (15a)
```

Their transfer conclusions are therefore subsumed, although their uniform
anchor-preserving edge and transversal ledgers are new exact controls. The
remaining positive row

```text
A_2={40,143,252}                                     (15b)
```

exchanges both `120` and `126`. Restricting to

```text
K in binom(O_(A_2)\{120,126},7)
```

gives exactly

```text
binom(25,7)=480,700                                  (15c)
```

bodies containing exactly one member of `A_0`. They are therefore disjoint
from both the all-anchor THM-4174 family and the one-anchor-exchange THM-4175
family. The other choices in this row split as

```text
2 binom(25,6)=354,200  one-anchor overlaps,
  binom(25,5)= 53,130  all-anchor overlaps.           (15d)
```

Among the two-anchor replacement triples, `A_2` is positive while
`{80,143,252}` and `{143,240,252}` are exact seven-cover certificate blocks.
THM-4175's projected total-depth hypergraph and the present anchor-preserving
uniform-deletion hypergraph have different vertex sets and edge arities; this
overlap ledger compares target bodies, not repair decks.

## 2. THM-4175 failure-atom specialization

Equation `(6)` is the `q`-safe specialization of THM-4175's finite-measure
deletion tomography. The short direct verification is retained to fix the
boundary convention and the exact labelled sidecar used by both audits.

For every `y`, membership in the deletion safe set in `(5)` is equivalent to

```text
||qy||>=1/14
and ||py||>=1/14 for every p in P\R.                  (16)
```

By definition `(3)`, the second condition says exactly

```text
F_P(y) subset R.                                      (17)
```

The labelled atoms `{y:F_P(y)=F}`, `F subset P`, partition the circle.
Intersecting this partition with `G_{50}` or, more generally, `G_q`, and
summing the measures of the atoms satisfying `(17)` proves `(6)`. This uses
the strict failure convention in `(3)`, so boundary points with clearance
exactly `1/14` are safe; in any case the rational walls have Haar measure
zero. Formula `(7)` is ordinary Mobius inversion on the Boolean subset
lattice. **QED for `(6)--(7)`.**

Two consequences are type-critical:

1. If `|R|<=6`, only atom labels `F` with `|F|<=6` can contribute to `(6)`.
2. Aggregating `w_q(F)` only by `|F|` destroys the labelled overlaps that
   determine `tau(E_d^A)`. The full mask, or an equivalent incidence sidecar,
   is essential; THM-4172's symmetric-layer warning applies literally.

## 3. Exact anchor classification

The classification `(9)` is a complete finite decision over the explicit
universe

```text
binom(P,3),                    |binom(P,3)|=4,060.      (18)
```

For each triple the primary audit evaluates the fourteen integer predicates in `(8)`
and its exact integer gcd. Exactly six triples pass. Equivalently, after the
required `13`-divisibility coordinate is imposed, the only surviving partner
pairs are

```text
{40,252}, {80,252}, {120,126}, {120,252},
{126,240}, {240,252},                                  (19)
```

each attached to `143`. Direct divisibility verifies that every triple in
`(9)` covers each modulus `2,...,14`, and direct gcd calculation gives one.
The exhaustive rejection of the other `4,054` triples proves completeness.

The nearby type hostile

```text
{120,126,286}                                          (20)
```

is divisor-complete but has gcd two. It shows why the anchor predicate cannot
be replaced by the Haar/transversal computation alone.

## 4. Exact q=50 computation

The fixed pool has the THM-4174 wall lattice

```text
L=lcm{14p:p in P}=18,241,159,416,480,
7,134 walls and 7,133 open cells.                      (21)
```

On each cell the complete 30-label mask `F_P` is constant. Unlike the
THM-4174 original-anchor specialization, the primary audit retains every mask,
including all `3,177` cells on which at least one of `120,126,143` fails.
There are `2,950` distinct labelled masks, `2,605` of which have positive
`q=50` safe mass. On denominator `14*50*L`, the total `q`-safe numerator is

```text
10,944,695,649,888,000,                                (22)
```

and the part on original-anchor-failing cells has numerator

```text
4,072,752,399,148,720.                                 (23)
```

Thus the old-anchor-failing atoms are materially present, not merely stored
as zero-weight labels.

For a wall tick `t`, the q-safe prefix numerator used by the primary audit is

```text
J_q(t)=12kL+
  0,                 if 14r<=L,
  14r-L,             if L<14r<13L,
  12L,               if 13L<=14r,

where qt=kL+r and 0<=r<L.                              (24)
```

Cellwise differences of `(24)` give the exact integer atom weights. The
subset-zeta sum `(6)` is then evaluated once for every global deletion mask:

```text
d=5: binom(30,5)=142,506 masks, 3,017 repairs,
d=6: binom(30,6)=593,775 masks, 85,324 repairs.         (25)
```

There are no threshold equalities in either global layer. Restricting masks
to `R subset O_A` produces the table in Section 1.

For each restricted hypergraph, an exhaustive recursion searches cover
budgets `0,...,7`. A greedy packing of pairwise disjoint uncovered edges is a
valid lower-bound prune. If no seven-cover exists, budget eight is searched;
every returned witness is checked against every active edge. Hence each
reported value `4,5,7,8` is an exact transversal number, not merely an upper
bound. The original-anchor controls reproduce THM-4174 exactly:

```text
d=5: 2,063 edges and tau=5,
d=6: 46,261 edges and tau=8.                           (26)
```

Normal and optimized Python streams byte-match the frozen output.

## 5. Haar transfer for exactly the positive rows

Fix `A in A_+` and `K in binom(O_A,7)`. By `(12)`,

```text
tau(E_6^A)=8>7.                                       (27)
```

Therefore `K` is not a transversal of `E_6^A`: some repair
`R in E_6^A` is disjoint from `K`. Since `R subset O_A`, it is also disjoint
from `A`, and therefore

```text
H_(A,K) subset (P union {50})\R.                      (28)
```

Safe-set inclusion reverses under inclusion of speed sets, so `(10)` and
`(28)` give

```text
G_((P union {50})\R) subset G_(H_(A,K)),
mu(G_(H_(A,K)))>=4/63.                                (29)
```

For positive `c`, the surjective circle endomorphism `m_c(y)=cy` preserves
Haar measure and satisfies

```text
G_(cH)=m_c^(-1)(G_H),             mu(G_(cH))=mu(G_H). (30)
```

THM-4150 applied to `(29)--(30)` proves `(14)` for every distinct positive
odd tail pair. Every core is primitive because `gcd(A)=1`, and it is
divisor-complete because `A` already supplies a multiple of each modulus
`2,...,14`. Also `50 notin P`, so each core has eleven distinct labels.
After doubling, body labels are even and cannot collide with the odd tails.
This proves `(14)--(15)` for the three rows in `A_+`. **QED relative to
THM-4150.**

For `A_2`, the optional ground set has 27 labels, of which `120,126` are
exactly the two omitted members of `A_0`. Excluding them leaves 25 choices,
so `(15c)` is `binom(25,7)`. Choosing exactly one or both omitted anchors gives
the two counts in `(15d)`. Every body in the `(15c)` slice contains `143` and
neither `120` nor `126`; hence it has exactly one original anchor. THM-4174
bodies have all three original anchors and THM-4175 bodies have exactly two.
This proves the overlap partition and novelty assertion following `(15c)`.

For `A in A_blk`, the displayed seven-set in the table meets every six-repair.
Thus the implication from `(27)` to `(28)` is unavailable for that body.
This is an exact failure of the six-deletion transversal certificate only.
It does **not** prove that the corresponding body is unsafe, refute the Haar
criterion, or supply an LRC counterexample.

## 6. Source-to-target contract and failure anatomy

```text
source:       full q=50-labelled failure-atom measure vector on P
target:       each re-anchored d=5,6 repair hypergraph, then odd-tail rows
map:          subset-zeta sum -> 4/63 threshold -> anchor restriction
              -> exact transversal -> missed repair -> THM-4150
preserved:    q=50, all pool labels, failure overlaps, deletion arity,
              exact Haar threshold, anchor divisor pins and content
destroyed:    interval/component order, a canonical safe phase, the selected
              repair after existential transfer, body-specific extra phases
sidecar:      labelled anchor triple plus full repair-incidence deck
positive:     the three rows A_+ have tau(E_6^A)=8; only the A_2 slice in
              (15c) is new beyond THM-4174/4175
blocks:       both A_blk rows have explicit seven-covers; every d=5 row has
              tau<=5; {120,126,286} is divisor-complete but nonprimitive
decisive test: all 736,281 global d=5,6 masks plus exact covers through 8.
                                                               (31)
```

The mechanism is not edge abundance. For example, the hostile anchor
`{80,143,252}` has `33,636` six-repairs but a seven-cover, while the positive
anchor `{120,143,252}` has `37,637` and needs eight. The exact coordinate is
labelled transversal dispersion, not an arity or edge-count histogram.

## 7. Audit architecture, replay, and scope

The two exact paths share no implementation; both use the elementary
failure-atom identity and declared pool. The primary audit keeps the 7,134
pool walls and integrates the
`q=50` comb exactly inside each pool cell by the prefix formula `(24)`. The
independent C++20 audit inserts the `q=50` walls into the arrangement, works
directly on

```text
D=lcm({14p:p in P} union {700})=91,205,797,082,400,
7,214 joint walls and 7,213 open atoms,                  (32)
```

and sums joint-atom lengths. It then rebuilds each re-anchored hypergraph and
searches cover budgets `0,...,8` independently. Its labelled-atom
FNV-1a-64 control is `0ccd305ae47c79ea`. All twelve edge counts, zero-equality
counts, and exact transversal numbers agree. Several displayed minimum covers
differ because minimum hitting sets are nonunique; every witness is checked
against every active edge.

Primary replay:

```bash
python3 -B 04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py
python3 -B -O 04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py
PYTHONHASHSEED=4178 python3 -B \
  04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178.py
```

Independent replay:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178_independent_audit.cpp \
  -o /tmp/lrc4178-o2
/tmp/lrc4178-o2
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_divisor_complete_anchor_triple_exchange_thm4178_independent_audit.cpp \
  -o /tmp/lrc4178-o0
/tmp/lrc4178-o0
```

Both Python modes and both C++ optimization levels byte-match their respective
frozen outputs. The C++ path forces binary stdout on Windows, the Python path
forces LF, and repository attributes pin raw artifact bytes.

This theorem specializes and directly rechecks THM-4175's elementary
transform, then proves a **q=50-only** finite-exact anchor classification. It
establishes no all-`q` anchor exchange, eventual
discrepancy tail, or statement about anchors outside the fixed pool. It treats
only deletion arities five and six and distinct positive odd tails after
doubling. It does not provide physical entry, classify arbitrary
divisor-complete bodies, handle mixed/even tails, prove the two hostile rows
unsafe, or prove LRC(14).
