---
id: THM-4188
title: "Every-newcomer zero-original hierarchy and exact 23-resonance filtration"
status: >
  PROVED RELATIVE TO THM-4150/4170/4174/4175/4185 + VERIFIED-EXACT PRIMARY
  POOL-WALL AND INDEPENDENT JOINT-WALL AUDITS; LRC(14) OPEN. For every
  positive newcomer q outside the fixed thirty-label pool, all 1,491,665
  primitive divisor-complete zero-original ten-bodies have Haar mass at
  least 4/63 after adjoining q. The failure of q=50 edgewise inheritance is
  exactly the nested 15/19/23-label resonance filtration; every resonant row
  closes by a uniform shallower native repair hierarchy. A strict-limit
  component-discrepancy argument closes all q>=2587 independently of the
  finite census. A separate endpoint-erosion sidecar sharpens this sufficient
  cofinal range to q>=2479 but is not load-bearing for the theorem. Together
  with the disjoint THM-4174 and THM-4175 slices, this gives 9,039,920 named
  cores per newcomer.
source: codex-lrc-q-deformation-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4185-q50-complete-zero-original-minimal-anchor-hierarchy
related:
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
  - THM-4182-q50-zero-original-four-anchor-transfer-and-two-repair-obstruction-graph
script: 04-computation/lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp
output: 05-knowledge/results/lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.out
independent_audit_script: 04-computation/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp
independent_audit_output: 05-knowledge/results/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.out
sharpening_script: 04-computation/lrc14_q_deformation_endpoint_erosion_sidecar_thm4188.cpp
sharpening_output: 05-knowledge/results/lrc14_q_deformation_endpoint_erosion_sidecar_thm4188.out
script_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
output_sha256: 32d3d30cf5b7001e0d685d62a194265e4d0891185fc5f05b1cbfc495526161a5
independent_audit_script_sha256: 58817d5f5e1a8cc07384f3ea82a1feb221af37ab0907afde890ab4fbdd949137
independent_audit_output_sha256: 2e08ffe40ced44f9a33584af1d599e7acf76b1a48a2b0fbd6077d53377825528
sharpening_script_sha256: 323314b12e0baebad558adbf74c35ed3c4dea2d31345727fabf58e36aaacec70
sharpening_output_sha256: 2356a281a2977f4fa56ea7f64030b3e1b9adf23e0ff0d777cb2c911c4495c32a
sharpening_primary_dependency_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
hash_basis: raw LF bytes
primary_audit: >
  PASS. A warning-clean signed-integer C++20 pool-wall path rebuilds the
  complete q=50 E5/E6/E7 layers, proves their strict limiting surpluses and
  component-discrepancy cutoffs, scans every one of the 2,556 newcomers at
  q<=2586, and exhausts every labelled cover through the load-bearing
  budgets for all 23 resonance labels. Its standard bytewise FNV-1a-64
  little-endian semantic ledger is 188e752eabb43991.
independent_audit: >
  PASS. A separate warning-clean C++20 path explicitly constructs the joint
  P-union-{q} wall refinement in lcm(D,14q) coordinates, remeasures every
  joint atom, counts circular components by run transitions, and proves the
  native rows by recursive transversal search. Direct per-runner mass
  controls and exhaustive q=6 cover controls are positive. O2, O3, and
  UBSan streams byte-match; the independent bytewise ledger is
  6016229d0317ff06.
sharpening_audit: >
  VERIFIED-EXACT, NON-LOAD-BEARING. An endpoint-erosion sufficient
  certificate audits every q=50 E5/E6/E7 edge at Q=2478 and Q=2479. Exactly
  one E7 edge fails the erosion certificate at 2478, while every edge passes
  at 2479; monotonicity of the eroded lengths proves cofinal inclusion for
  q>=2479. The sidecar includes the pinned primary geometry, so it is not a
  second independent implementation and does not claim the true inclusion
  cutoff. O2 and UBSan outputs byte-match.
---

# THM-4188 -- every-newcomer zero-original hierarchy and resonance filtration

**PROVED RELATIVE TO THM-4150/4170/4174/4175/4185 + VERIFIED-EXACT; LRC(14)
REMAINS OPEN.**

## 1. Exact all-newcomer statement

Retain the fixed pool, the three original labels, and the zero-original
ground set

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143},                 U=P\A_0.              (1)
```

A nonempty `S subset U` is **good** when

```text
gcd(S)=1
and for every 2<=m<=14 some s in S satisfies m|s.        (2)
```

Put

```text
G_10={B in binom(U,10):B is good}.                       (3)
```

THM-4185 proves that `|G_10|=1,491,665` and that every `B in G_10`
contains an inclusion-minimal good anchor in one of the exact layers

```text
|M_4|=32,                 |M_5|=297,                 |M_6|=24. (4)
```

For every positive integer `q notin P` and `d in {5,6,7}`, define the
global repair layer

```text
E_d(q)={R in binom(P,d):
        mu(G_((P union {q})\R))>=4/63},                 (5)

G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                (6)
```

Every deletion in `(5)` is taken from `P` only: the newcomer `q` is never
deleted. For an anchor `A subset P`, write

```text
E_d^A(q)={R in E_d(q):R intersect A=empty}.              (7)
```

> **Every-newcomer zero-original theorem.** For every positive integer
> `q notin P` and every `B in G_10`,
>
> ```text
> mu(G_(B union {q}))>=4/63.                             (8)
> ```
>
> Consequently, for every positive integer `c` and every two distinct
> positive odd integers `a,b`, there is an `x in R/Z` such that
>
> ```text
> min_(v in 2c(B union {q}) union {a,b})||vx||>=1/14.    (9)
> ```

The body `B` is primitive and divisor-complete by `(2)`; goodness defines the
declared normalized slice and its count, not a hypothesis of THM-4150. That
theorem applies to any nonempty finite positive body once its complete safe
set has mass at least `4/63`. The arbitrary positive content in `(9)` follows
from Haar invariance under `y mapsto cy` (equivalently, THM-4150's inherited
content corollary), while distinct positivity and oddness are the tail
conditions.

Every newcomer supplies exactly `1,491,665` distinct labelled zero-original
cores. For `N>=290`, newcomer labels `1<=q<=N` supply

```text
(N-30)*1,491,665.                                        (10)
```

Combining the three pairwise-disjoint original-anchor-count slices proved by
THM-4174 (all three originals), THM-4175 (exactly two originals), and this
theorem (zero originals) gives exactly

```text
888,030+6,660,225+1,491,665=9,039,920                   (11)
```

named cores per newcomer. Only the exactly-one-original slice is not promoted
uniformly in `q`: at `q=50`, THM-4179 supplies its `1,071,961` cores, and
`9,039,920+1,071,961=10,111,881`, the four-slice total recorded in THM-4185.
Equation `(11)` neither duplicates nor replaces that `q=50` count.

## 2. Exact edgewise phase diagram

The `q=50` global layers have

```text
|E_5(50)|=3,017,       |E_6(50)|=85,324,
|E_7(50)|=821,737,     threshold equalities=0,0,0.       (12)
```

Define

```text
Q_5={6,24,25,48,70,96,105,110,128,140,192,206,210,256,366},

Q_6=Q_5 union {22,72,130,260},

Q_7=Q_6 union {100,186,220,294}.                         (13)
```

> **Exact resonance filtration.** For every positive `q notin P`,
>
> ```text
> E_d(50) subset E_d(q)    iff    q notin Q_d,
>                                                d=5,6,7. (14)
> ```

Thus `Q_5 subset Q_6 subset Q_7`, with cardinalities `15,19,23`. The set
`Q_7`, not merely a transversal-number failure set, is the complete global
exception set needed simultaneously on all three layers.

The stronger edgewise claim that `q=50` is worst is false. The first global
failure occurs at `q=6`:

```text
R={40,80,145,168,290} in E_5(50)\E_5(6).                (15)
```

Nor does anchor restriction repair the statement. At `q=48`, first exact
witnesses on the three THM-4185 rows are

```text
M_4: A={40,63,84,286},
     R={16,20,60,80,95,170,240} in E_7(50)\E_7(48),

M_5: A={8,42,60,63,286},
     R={10,15,16,80,120,170,240} in E_7(50)\E_7(48),

M_6: A={8,10,42,63,132,286},
     R={16,60,80,95,170,240} in E_6(50)\E_6(48).        (16)
```

Every displayed repair is disjoint from its anchor. The strongest survivor
is the dichotomy proved below: outside `Q_7`, global edgewise inclusion
transports the `q=50` hierarchy; at the 23 resonances, a strictly shallower
native hierarchy closes every anchor directly.

## 3. Analytic cofinal mechanism

Put

```text
D=lcm_(p in P)(14p)=18,241,159,416,480.                 (17)
```

For `R subset P`, let

```text
U_R=G_(P\R),       mu(U_R)=m_R/D,       c_R=#components(U_R). (18)
```

The component discrepancy lemma inherited from THM-4170 has a periodic
primitive `psi` with

```text
-3/49 <= psi <= 3/49.                                   (19)
```

On each interval component `[alpha,beta]` of `U_R`, integration against the
`q`-comb differs from `(6/7)(beta-alpha)` by
`(psi(q beta)-psi(q alpha))/q`. Summing the endpoints of the `c_R`
components gives

```text
|mu(U_R intersect G_{q})-(6/7)mu(U_R)| <= 6c_R/(49q).   (20)
```

Since

```text
(6/7)(m_R/D)-4/63
  =2(27m_R-2D)/(63D),                                   (21)
```

every edge with `27m_R-2D>0` belongs to `E_d(q)` whenever

```text
q >= ceil(27c_R D/[7(27m_R-2D)]).                       (22)
```

The exact geometry of every edge in `(12)` gives:

| layer | nonpositive surpluses | minimum `27m_R-2D` | attaining edge | maximum ceiling | attaining edge | `c_R` | `m_R` |
|:---|---:|---:|:---|---:|:---|---:|---:|
|`E_5(50)`|`0`|`5,601,528,091,116`|`{8,10,95,193,240}`|`2,297`|`{40,80,95,120,193}`|`184`|`1,559,969,234,178`|
|`E_6(50)`|`0`|`5,531,567,360,904`|`{15,80,84,95,176,240}`|`2,443`|same|`192`|`1,556,069,859,032`|
|`E_7(50)`|`0`|`5,524,132,150,368`|`{8,10,40,84,95,176,240}`|`2,587`|`{8,10,20,40,80,95,240}`|`204`|`1,556,695,841,942`|

Thus all three inclusions in `(14)` hold for every `q>=2587` by analysis,
not by extrapolating a finite census. Equality in an individual discrepancy
bound is harmless because repairs use the nonstrict threshold `>=4/63`; the
strict-limit surplus is the load-bearing positivity.

There is also an independently derived but computationally primary-dependent
sharpening. If a circular finite union is written as `U=disjoint_union_j I_j`,
endpoint erosion gives

```text
mu(U intersect G_q)
 >= (6/7) sum_j max(|I_j|-1/(7q),0).                    (23)
```

Thus an eroded length at least `2/27` at `Q` proves the threshold for every
`q>=Q`. The pinned erosion sidecar audits all `q=50` edges. At `Q=2478`
exactly one E7 edge,

```text
{8,10,15,20,40,95,145},
```

misses this sufficient certificate by `36,319,963,525,776` scaled units; at
`Q=2479` the same minimizing edge has positive margin
`3,018,582,994,392`, and every E5/E6/E7 edge passes. Hence `(14)` also follows
cofinally for every `q>=2479`. This endpoint result is sharp only for the
erosion certificate: it neither claims that inclusion fails at `q=2478`
(it does not) nor identifies a true inclusion cutoff. Because its executable
includes the pinned primary geometry rather than rebuilding a second path,
the theorem retains `(20)--(22)` and `q>=2587` as its load-bearing analytic
tail.

The endpoint sidecar also explains the last finite hostile. For

```text
D_*={8,15,145,193,290},
```

the fixed safe union has `170` components and hence `340` oriented endpoints.
At `q=366`, `gcd(366,D)=6`: the endpoint map collapses those endpoints onto
`244` residues while retaining oriented coefficient `ell^1` norm `340` and
cyclic boundary height `53`. The exact prefix-error term is
`-992,896,373,433,456`, giving negative repair delta
`-166,972,676,678,640`. At the adjacent coprime label `q=367`, all `340`
endpoint residues are distinct, cyclic height drops to `17`, the error is
`55,264,020,770,676`, and the repair delta is positive
`9,290,430,146,252,052`. This identifies the sixfold endpoint resonance at
the last hostile; injectivity or coprimality alone is not a uniform proof for
the other labels.

## 4. Finite exact resonance census

It remains to inspect `1<=q<=2586`, excluding the thirty labels in `P`: an
exact universe of `2,556` newcomers. The primary path fixes the `7,133`
positive pool-wall cells in denominator `D`, integrates each `q`-comb over
them by an exact integer safe-prefix formula, and subset-sums every deletion
through size seven.

The independent path instead uses

```text
D_q=lcm(D,14q)                                          (24)
```

and explicitly refines all pool walls and all newcomer walls

```text
(14k+1)/(14s), (14k+13)/(14s),       s in P union {q}.  (25)
```

On each joint atom it reclassifies `q` and records the exact integer length
`ell`; the repair comparison is literally

```text
63*sum(ell) >= 4D_q.                                    (26)
```

Both routes produce `(13)--(14)`. Across the finite universe there are zero
equality hits when the `q=50` E5/E6/E7 edges are tested at `(26)`. On the 23
resonances, the complete native ledger is:

|`q`|lost q50 E5/E6/E7 edges|affected M4/M5/M6 rows|native E5/E6 edges|
|---:|:---|:---|:---|
|6|`2/5/5`|`18/203/20`|`50,160/389,365`|
|22|`0/1/3`|`14/202/24`|`54,396/421,289`|
|24|`15/117/233`|`32/263/22`|`51,188/390,537`|
|25|`2/11/31`|`23/252/20`|`34,684/320,668`|
|48|`2/14/29`|`23/246/17`|`77,327/475,101`|
|70|`10/42/71`|`29/169/12`|`73,996/469,924`|
|72|`0/10/96`|`32/289/24`|`59,802/418,692`|
|96|`201/2,603/11,454`|`32/297/24`|`21,457/247,722`|
|100|`0/0/12`|`30/247/0`|`26,410/284,863`|
|105|`23/448/2,062`|`32/297/24`|`37,139/325,808`|
|110|`102/822/2,369`|`32/297/24`|`44,484/365,850`|
|128|`77/594/1,557`|`32/297/24`|`31,275/301,247`|
|130|`0/2/8`|`12/183/18`|`105,407/539,067`|
|140|`29/188/366`|`32/273/18`|`59,351/416,495`|
|186|`0/0/1`|`20/96/0`|`61,671/432,488`|
|192|`5/73/384`|`23/297/24`|`43,980/364,769`|
|206|`2/11/16`|`28/279/24`|`65,683/451,190`|
|210|`243/2,805/11,624`|`32/297/24`|`26,533/261,505`|
|220|`0/0/1`|`8/115/0`|`77,964/479,384`|
|256|`230/2,284/8,185`|`32/297/24`|`25,765/273,944`|
|260|`0/1/1`|`0/0/18`|`61,262/422,570`|
|294|`0/0/1`|`12/147/0`|`68,909/454,315`|
|366|`1/2/1`|`20/80/15`|`57,331/422,399`|

An affected row means that at least one lost `q=50` edge is disjoint from
the anchor; it is not a body-failure count. Native E5/E6 threshold equalities
are zero throughout the 23 labels.

Complete cover quantifiers, not edge density, give

```text
q in Q_7, A in M_4: tau(E_6^A(q))>6,
q in Q_7, A in M_5: tau(E_5^A(q))>5,
q in Q_7, A in M_6: tau(E_5^A(q))>4.                    (27)
```

The primary audit enumerates every candidate through the stated budgets.
The independent audit reaches the same zero-cover verdict by bounded
recursive transversal search, with positive and negative synthetic controls;
at `q=6` it additionally exhausts exact-budget candidates for one anchor of
each size. Any smaller cover extends outside the anchor to the displayed
budget, so this hostile control checks the full `<=` quantifier.

## 5. Why every body closes

Fix `B in G_10`. By THM-4185, choose `A subset B` in one of
`M_4,M_5,M_6`, and put `K=B\A`. Then `|K|` is respectively `6,5,4`.

If `q notin Q_7`, nesting in `(13)` and `(14)` gives every needed global
inclusion. Restriction to edges disjoint from `A` preserves inclusion, so
the THM-4185 hitting-number rows transport as

```text
A in M_4: tau(E_7^A(q)) >= tau(E_7^A(50))>6,
A in M_5: tau(E_7^A(q)) >= tau(E_7^A(50))>5,
A in M_6: tau(E_6^A(q)) >= tau(E_6^A(50))>4.             (28)
```

If `q in Q_7`, use the native shallower rows `(27)` instead. In either case,
`K` is not a transversal of the relevant restricted layer. Hence some
`R subset P` is disjoint from both `A` and `K`, so `R intersect B=empty`, and

```text
mu(G_(B union {q}))
 >= mu(G_((P union {q})\R))
 >= 4/63.                                               (29)
```

The first inequality is safe-set monotonicity: the ambient repaired pool
contains `B union {q}`. This proves `(8)`. Apply THM-4150 to the nonempty
finite positive body `B union {q}` and the two distinct positive odd tails;
Haar invariance under `y mapsto cy` supplies every positive common content
`c`, proving `(9)`. Primitivity and divisor-completeness remain the declared
slice conditions, not extra transfer hypotheses.

## 6. Audit contract and replay

The primary and independent paths agree on the fixed geometry, anchor counts,
all `q=50` layers, every analytic table entry, the complete `2,556`-newcomer
filtration, every native E5/E6 count, all restricted affected counts, zero
load-bearing threshold equalities, and all `353` native anchor rows per
resonance label. They do not share the newcomer integration or cover-search
algorithm.

Primary replay:

```bash
clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp \
  -o /tmp/lrc4188-primary
/tmp/lrc4188-primary \
  > 05-knowledge/results/lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.out
```

Independent normal and sanitizer replay:

```bash
clang++ -std=c++20 -O2 -pthread -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp \
  -o /tmp/lrc4188-independent
/tmp/lrc4188-independent \
  > 05-knowledge/results/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.out

clang++ -std=c++20 -O1 -g -pthread -fsanitize=undefined \
  -fno-sanitize-recover=all -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp \
  -o /tmp/lrc4188-independent-ubsan
/tmp/lrc4188-independent-ubsan > /tmp/lrc4188-independent-ubsan.out
```

The O2, O3, and UBSan independent outputs have the same SHA-256
`2e08ffe40ced44f9a33584af1d599e7acf76b1a48a2b0fbd6077d53377825528`.

Non-load-bearing endpoint-erosion sharpening:

```bash
c++ -std=c++20 -O2 -DNDEBUG -Wall -Wextra -pedantic \
  04-computation/lrc14_q_deformation_endpoint_erosion_sidecar_thm4188.cpp \
  -o /tmp/lrc4188-erosion
/tmp/lrc4188-erosion \
  > 05-knowledge/results/lrc14_q_deformation_endpoint_erosion_sidecar_thm4188.out
```

Its normal and UBSan output SHA-256 is
`2356a281a2977f4fa56ea7f64030b3e1b9adf23e0ff0d777cb2c911c4495c32a`.

## 7. Scope and firewalls

1. `q notin P` is essential to the newcomer universe and distinct-core
   count. No claim here is made for duplicate labels `q in P`.
2. The body universe is exactly the `1,491,665` good members of
   `binom(U,10)`. Arbitrary ten-bodies and other body sizes remain outside the
   statement.
3. Only the zero-original slice is deformed here. THM-4174 separately supplies
   the all-three-original slice and THM-4175 supplies the exactly-two-original
   slice for every newcomer. The exactly-one-original slice remains
   `q=50`-only through THM-4179.
4. The conclusion is the safe-set Haar bound and universal doubled-body
   transfer for positive common content and two distinct positive odd tails.
   It does not establish physical entry, mixed/even-tail branches, or an
   arbitrary thirteen-speed theorem.
5. This theorem does not prove LRC(14).
