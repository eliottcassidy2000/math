---
id: THM-4185
title: "q=50 complete zero-original ten-body anchor hierarchy and repair transfer"
status: >
  PROVED RELATIVE TO THM-4150/4179/4182 + VERIFIED-EXACT PYTHON POOL-WALL
  AND INDEPENDENT C++20 JOINT-WALL AUDITS; LRC(14) OPEN. Every one of the
  1,491,665 primitive divisor-complete ten-subsets of the 27-label
  zero-original ground set contains a certified anchor of size four, five,
  or six. Complete blocker descent gives depth at most seven for the
  size-five rows and depth six for the size-six rows, closing the entire
  declared zero-original ten-body slice at q=50.
source: codex-lrc-anchor-incidence-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
  - THM-4182-q50-zero-original-four-anchor-transfer-and-two-repair-obstruction-graph
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
script: 04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy.py
output: 05-knowledge/results/lrc14_q50_zero_original_minimal_anchor_hierarchy.out
independent_audit_script: 04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy_independent.cpp
independent_audit_output: 05-knowledge/results/lrc14_q50_zero_original_minimal_anchor_hierarchy_independent.out
script_sha256: 23cc6d018178bfec4adb25e5149b69e6401dd1cb8c9c552d039e4e35fb46df13
output_sha256: f2ed373f70c073aa576c30aa97094eebed8687100784a925c8d1bbd2879196c9
independent_audit_script_sha256: 0c26df20b1a1c047af4d673ee6649b3f918df628ce73c6d87793b932f441140b
independent_audit_output_sha256: 0d231f3cd589ded4b42f6d42da036f29fec92f68439f409cdf491dafd43e491c
primary_dependency_sha256: 1e1ec137b27b8555b2c9759a3e283746e09f445a8bd6bdba25ccd94187b05958
independent_geometry_dependency_sha256: 47c0d8f1107d35b770d63348a21fcac04ffc08dc4a611b46b221f0be0ac5aa3b
hash_basis: raw LF bytes
primary_audit: >
  PASS. A hash-pinned THM-4182 pool-wall path exhausts all size-one through
  size-six candidate anchors, all three global deletion layers, every
  labelled cover through budgets five and four by literal combination
  enumeration, and all C(27,10) zero-original bodies. Normal, optimized, and
  hostile-hash-seed streams byte-match. The load-bearing complete size-five
  ledger has 21,506 covers, not the 8,075 first-hit states from the corrected
  exploratory near miss.
independent_audit: >
  PASS. A warning-clean signed-integer C++20 path rebuilds q=50 in the joint
  wall arrangement, scans every cover candidate directly against active
  edges, and classifies each of the C(27,10) bodies by its smallest contained
  size-four/five/six anchor rather than taking an anchor-presentation union.
  O2, O0, and UBSan streams byte-match. The two paths agree on both row
  ledgers, all residual presentations, witness counts, and the body partition.
---

# THM-4185 -- complete zero-original ten-body anchor hierarchy and repair transfer

**PROVED RELATIVE TO THM-4150/4179/4182 + VERIFIED-EXACT; LRC(14) REMAINS
OPEN.**

## 1. Exact ten-body statement

Retain the fixed pool, original anchor, zero-original ground set, and newcomer

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143},                 U=P\A_0,             q=50.       (1)
```

Call a nonempty set `S subset U` **good** when

```text
gcd(S)=1
and for every 2<=m<=14 some s in S satisfies m|s.                  (2)
```

For `1<=j<=6`, let `M_j` be the good `j`-subsets of `U` containing no
proper good subset. Exact enumeration gives

| `j` | `|M_j|` | SHA-256 of ordered little-endian pool masks |
|---:|---:|:---|
|1|`0`|`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`|
|2|`0`|same empty hash|
|3|`0`|same empty hash|
|4|`32`|`359881969de938d942ddf69501c0d66fee783fe12910e348b915a23c82e3b560`|
|5|`297`|`fdde1bd87afcb906da126ee5a8cd43c7f3eb0cabc1a67f844c150fb3f7cc96d6`|
|6|`24`|`34564111d2efd367ca6a23f02eef4099760d055b08e5cd9e4e0b30f981ce07d5`|

The 32 members of `M_4` are exactly the four-anchors of THM-4182. Every
member of `M_4 union M_5` contains `286`. The size-six layer has the complete
product classification

```text
M_6={{42,63,132,286,e,f}:
     e in {8,16,88,176},
     f in {10,20,30,170,190,290}}.                              (3)
```

Now define the declared ten-body universe

```text
G_10={B in binom(U,10):B is good}.                              (4)
```

> **Complete zero-original ten-body hierarchy.** One has
>
> ```text
> |G_10|=1,491,665.                                             (5)
> ```
>
> Every `B in G_10` contains a member of `M_4 union M_5 union M_6`.
> Partitioning by the least such anchor size gives exactly
>
> ```text
> size 4: 1,138,494,
> size 5:   348,712,
> size 6:     4,459.                                           (6)
> ```

Thus the cumulative unions after the size-four, size-five, and size-six
layers have cardinalities

```text
1,138,494,             1,487,206,             1,491,665.       (7)
```

The scope is deliberately exact: `(3)--(7)` classify every good ten-body in
`U`. The two promoted audits enumerate candidate anchors only through size
six and make no claim here about globally inclusion-minimal good subsets of
`U` having size seven or more.

For `A in M_j`, put `V_A=P\A` and

```text
E_d^A={R in binom(V_A,d):
       mu(G_((P union {50})\R))>=4/63}.                         (8)
```

> **Anchor-size repair staircase.** Uniformly over the declared anchors,
>
> ```text
> A in M_4  implies tau(E_7^A)>6,
> A in M_5  implies tau(E_7^A)>5,
> A in M_6  implies tau(E_6^A)>4.                              (9)
> ```

The first row is THM-4182 plus one-step blocker descent. The second and third
rows are the new exact repair-incidence theorem proved below.

Consequently, for every `B in G_10`, every positive integer `c`, and every
two distinct positive odd integers `a,b`, there is a point `y in R/Z` such
that

```text
min_(v in 2c(B union {50}) union {a,b}) ||vy|| >= 1/14.         (10)
```

This closes the entire primitive divisor-complete zero-original ten-body
slice of the fixed `q=50` pool. It is not a theorem about arbitrary integer
bodies and does not prove `LRC(14)`.

## 2. Inheritance, mechanism, and exact map

The closest mechanism is
[THM-4179](THM-4179-q50-seventh-deletion-primitive-anchor-completion.md):
for monotone repair layers on an `n`-vertex ground set, every size-at-most-`k`
cover of `E_(d+1)` is a cover of `E_d` whenever `d+k<n`. The canonical
hostiles inherited from THM-4182 are its eight depth-six four-anchor rows.
The corrected near miss is that a first-hit search ledger is a complete cover
ledger; it is not. The least-used sidecar is the exact labelled
anchor--cover--body--repair incidence, retained across both deletion steps.

The new source-to-target contract is

```text
source:       every good B in binom(U,10) and global q=50 E5,E6,E7 layers
target:       a threshold repair for B union {50}, then THM-4150 odd-tail transfer
map:          B -> least contained anchor size -> (A,K=B\A)
              -> complete shallow cover ledger -> missed deeper repair
preserved:    all body labels, q, gcd, divisor pins, threshold 4/63, content
destroyed:    anchor presentation after union, component address, chosen phase
sidecar:      anchor mask, every cover through budget, residual body, repairs
positive:     41 size-five and 22 size-six rows have no shallow cover
hostile:      15 size-five E6 covers and four size-six E5 covers
decisive test: all C(27,10) bodies and every labelled cover candidate.      (11)
```

This is an anchor-replacement closure rather than a repeat of a fixed-anchor
count. THM-4182 covers the bodies containing a four-anchor. Equation `(6)`
finds the exact `353,171` good holdouts and proves that size-five anchors
cover all but `4,459`; the size-six product layer covers the remainder.

## 3. Exact repair layers

The primary audit retains the THM-4178 pool-wall arrangement, integrates the
`q=50` comb into its exact labelled failure atoms, and subset-zeta sums every
deletion mask. Before anchor restriction, the three layers are

| depth `d` | candidates | repairs | equalities | minimum edge delta | maximum nonedge delta | edge-mask SHA-256 |
|---:|---:|---:|---:|---:|---:|:---|
|5|`142,506`|`3,017`|`0`|`186,258,371,040`|`-69,055,421,400`|`2c01c096640ff8690496a424fce357cdb38d0a29cafc81bacfe94eea64729f06`|
|6|`593,775`|`85,324`|`0`|`348,508,440`|`-2,639,430,360`|`503ffc5e72825f275a9b351c1ceab1b13c452c5082ec4009020bb858c44aafdd`|
|7|`2,035,800`|`821,737`|`0`|`348,508,440`|`-402,312,960`|`2a2b26ffd0c26000260ead13f507955b332d5c61927219418a49225b1a004b7d`|

All comparisons are therefore strict; there is no unresolved threshold
boundary at any load-bearing depth.

For a size-five anchor, `|V_A|=25`. The primary audit literally tests every
subset of `V_A` through size five against `E_5^A`, including every extension
of a smaller cover. Among the 297 rows:

```text
41 rows:  no E5 cover through size 5,
256 rows: at least one such cover,
all covers through size 5: 21,506,
cover-size histogram: 3:33, 4:1,229, 5:20,244.                 (12)
```

Filtering all `21,506` covers against `E_6^A` leaves exactly `15` labelled
presentations. They collapse to two bodies. Put

```text
C={170,252,286},
D={88,95,145,168,170,193,240,252,286},
B_16=D union {16},                 B_80=D union {80},

Gamma_16={{16,95},{16,145},{16,193},
          {88,95},{88,145},{88,193},
          {95,168},{145,168},{168,193}},
Gamma_80={{88,95},{88,145},{88,193},
          {95,168},{145,168},{168,193}}.                       (13)
```

The complete residual presentation list is

```text
(A,K)=(C union e, B_x\(C union e)),
x in {16,80},              e in Gamma_x.                       (14)
```

Thus `B_16` has nine anchor presentations and `B_80` has six. Their exact
depth-seven repair incidences are

| body | E7 repairs disjoint from the body | first lexicographic repair |
|:---|---:|:---|
|`B_16`|`736`|`{8,10,15,20,80,85,290}`|
|`B_80`|`514`|`{8,10,15,16,20,85,120}`|

Because `5+5<25`, every cover of `E_6^A` through size five occurs in the
complete ledger `(12)`. Because `6+5<25`, every cover of `E_7^A` through
size five occurs among the 15 presentations `(14)`. The displayed positive
repair incidences cross every one, proving the middle row of `(9)`. For these
15 hostile presentations, deletion depth seven is sharp within the
anchor-preserving filtration: their choices cover both `E_5^A` and
`E_6^A`, while an `E_7^A` repair misses them.

For a size-six anchor, `|V_A|=24`, and all covers of `E_5^A` through size
four are exhausted directly. Twenty-two rows have none. Only the following
four covers occur in the other two rows:

| anchor `A` | cover `K` | E6 repairs disjoint from `A union K` | first repair |
|:---|:---|---:|:---|
|`{42,63,88,132,286,290}`|`{80,85,168,193}`|`913`|`{8,10,15,16,95,145}`|
|same|`{85,168,193,240}`|`708`|`{8,10,15,16,95,145}`|
|same|`{95,168,193,240}`|`461`|`{8,10,16,80,85,145}`|
|`{42,63,132,176,286,290}`|`{95,168,193,240}`|`651`|`{8,10,15,85,88,143}`|

There are no smaller covers. Since `5+4<24`, any cover of `E_6^A` through
size four would be one of these four; every row is crossed by the displayed
positive incidence. This proves the last row of `(9)`. Depth six is sharp
for the four displayed hostile presentations.

The ordered row-ledger fingerprints, using standard bytewise FNV-1a-64 on
little-endian semantic `uint64` words, are

```text
size-five rows: 9dbf51ecd6f8c2d5,
size-six rows:  a4f776dc340f0260.                               (15)
```

Both implementations independently reproduce `(15)`.

## 4. Complete body partition and transfer

There are `binom(27,10)=8,436,285` zero-original ten-subsets. A direct gcd
and thirteen-divisibility-bit scan finds exactly `1,491,665` good bodies.
The anchor-presentation counts are

```text
|M_4| binom(23,6)=3,230,304,
|M_5| binom(22,5)=7,821,198,
|M_6| binom(21,4)=  143,640.                                  (16)
```

The primary path unions those presentations and obtains `(7)`. The
independent path instead scans each of the `8,436,285` bodies once, checks
goodness, and records its first contained anchor size. It obtains the
disjoint increments `(6)`. This independently proves both that the final
union is complete and that presentation overlap has not inflated the count.

Fix `B in G_10`. Choose a contained anchor `A in M_j` with `j in {4,5,6}`
and put `K=B\A`, so `|K|=10-j`. Equation `(9)` gives a repair
`R in E_d^A` disjoint from `K`, at depth at most seven. Since `R` also avoids
`A`,

```text
B union {50} subset (P union {50})\R,
mu(G_(B union {50}))
 >= mu(G_((P union {50})\R))
 >= 4/63.                                                       (17)
```

The body is primitive and divisor-complete by `(2)`. THM-4150 applied to
`(17)` proves `(10)`. **QED relative to THM-4150.**

The four already named, pairwise disjoint original-anchor-count strata at
fixed `q=50` now have exact core-body counts

```text
all three of A_0 (THM-4174):       888,030,
exactly two of A_0 (THM-4175):   6,660,225,
exactly one of A_0 (THM-4179):   1,071,961,
none of A_0 (this theorem):      1,491,665,
                                  ---------
total named transferred cores: 10,111,881.                    (18)
```

Equation `(18)` is a count of these four proved fixed-pool strata, disjoint
because the number of labels from `A_0` differs. It is not a count of all
primitive divisor-complete eleven-label cores in `P union {50}`.

## 5. Exact audit and reproduction

The primary path hash-pins the THM-4182 Python audit, whose own frozen
dependency is the THM-4178 pool-wall atom construction. It uses Python
big-integer incidence masks for cover tests but enumerates every candidate
combination through the target budget. Its output prints every size-five and
size-six row, all 19 hostile presentations, both residual body states, all
witness counts, and the full body hierarchy.

The independent path includes the hash-recorded THM-4179 C++20 joint-wall
geometry. It reconstructs

```text
denominator=91,205,797,082,400,
7,214 joint walls, 7,213 atoms,
labelled atom FNV-1a-64=0ccd305ae47c79ea.                       (19)
```

It uses signed-integer cell lengths, direct active-edge scans for every cover
candidate, and a body-first census. The C++ native delta normalization is
one twentieth of the primary normalization, with the same signs and zero
equalities. Its global E5/E6/E7 edge fingerprints are respectively

```text
6beebf7aa4701ef0, e9951e7b4a3fba86, 808ca531b29c5e83.          (20)
```

Reproduction from the repository root:

```bash
python3 04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy.py \
  > /tmp/lrc14_hierarchy.out
python3 -O 04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy.py \
  > /tmp/lrc14_hierarchy.opt.out
PYTHONHASHSEED=8731 \
  python3 04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy.py \
  > /tmp/lrc14_hierarchy.hash.out
cmp /tmp/lrc14_hierarchy.out /tmp/lrc14_hierarchy.opt.out
cmp /tmp/lrc14_hierarchy.out /tmp/lrc14_hierarchy.hash.out

clang++ -std=c++20 -O2 -Wall -Wextra -pedantic \
  04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy_independent.cpp \
  -o /tmp/lrc14_hierarchy_cpp
/tmp/lrc14_hierarchy_cpp > /tmp/lrc14_hierarchy_cpp.out

clang++ -std=c++20 -O0 -Wall -Wextra -pedantic \
  04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy_independent.cpp \
  -o /tmp/lrc14_hierarchy_cpp_o0
/tmp/lrc14_hierarchy_cpp_o0 > /tmp/lrc14_hierarchy_cpp_o0.out

clang++ -std=c++20 -O1 -g -Wall -Wextra -pedantic \
  -fsanitize=undefined -fno-sanitize-recover=all \
  04-computation/lrc14_q50_zero_original_minimal_anchor_hierarchy_independent.cpp \
  -o /tmp/lrc14_hierarchy_cpp_ubsan
/tmp/lrc14_hierarchy_cpp_ubsan > /tmp/lrc14_hierarchy_cpp_ubsan.out

cmp /tmp/lrc14_hierarchy_cpp.out /tmp/lrc14_hierarchy_cpp_o0.out
cmp /tmp/lrc14_hierarchy_cpp.out /tmp/lrc14_hierarchy_cpp_ubsan.out
```

All three Python streams have raw SHA-256
`f2ed373f70c073aa576c30aa97094eebed8687100784a925c8d1bbd2879196c9`.
All three C++ streams have raw SHA-256
`0d231f3cd589ded4b42f6d42da036f29fec92f68439f409cdf491dafd43e491c`.
The checked-in artifacts are the normal and O2 streams respectively.

## 6. Boundary and next question

The repair staircase is sharp at its displayed hostiles, but its endpoint is
still a fixed finite pool and the single newcomer `q=50`. Maximal-deletion
duality explains why no bounded-depth cover is a body obstruction, while the
new anchor hierarchy explains why every good zero-original ten-body enters:
the size-four family covers most bodies, and two successively larger anchor
layers cover exactly the holdouts.

The next structural question is whether this hierarchy persists under
newcomer variation or admits a pool-independent anchor-replacement theorem.
Neither exact audit addresses that parameter, arbitrary eleven-label cores,
mixed/even tails, or the full `LRC(14)` frontier.
