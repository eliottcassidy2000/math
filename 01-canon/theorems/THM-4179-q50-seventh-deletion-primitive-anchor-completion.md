---
id: THM-4179
title: "q=50 seventh-deletion primitive-anchor completion and maximal-blocker duality"
status: >
  PROVED RELATIVE TO THM-4150/4178 + VERIFIED-EXACT PYTHON POOL-WALL
  BLOCKER-DESCENT AND INDEPENDENT C++20 JOINT-WALL AUDITS; LRC(14) OPEN.
  The two THM-4178 depth-six seven-cover rows have depth-seven transversal
  number exactly eight, so every one of the six primitive divisor-complete
  pool anchors transfers at q=50 by depth at most seven. The exactly-one-
  original-anchor union has 1,071,961 distinct core bodies, 591,261 beyond
  THM-4178's single genuine two-anchor row.
source: codex-lrc-post4178-exchange-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
script: 04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179.py
output: 05-knowledge/results/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179.out
independent_audit_script: 04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.out
script_sha256: 3310b12f80aac7c8b774c8af5ba356bdb28a8a95e64a154541ea9b18165f3d89
output_sha256: ea3fb40c8c613e25a1ff563cd82ae88bb9c11587092faf98c2acd9163d1b10eb
independent_audit_script_sha256: 47c0d8f1107d35b770d63348a21fcac04ffc08dc4a611b46b221f0be0ac5aa3b
independent_audit_output_sha256: 1f4ea2ca6e454075c9916f84a9cdbe3859c3243ef4070b426c6d3e397a94a3c0
hash_basis: raw LF bytes
primary_audit: >
  PASS. The frozen THM-4178 pool-wall atoms are subset-zeta summed at depths
  six and seven. A new Python-big-int incidence recursion exhausts every
  depth-six seven-cover (one and three respectively), and a depth-seven edge
  misses each. Explicit eight-covers, four direct maximal-body masses, and a
  literal labelled-body overlap census are checked. Normal, optimized, and
  hash-seeded streams byte-match.
independent_audit: >
  PASS. A warning-clean C++20 path independently inserts q=50 into the joint
  wall arrangement, tests all 888,030 seven-deletions in each hostile row,
  and directly rejects cover budgets zero through seven before verifying a
  cover at eight. O2, O0, and UBSan streams byte-match. The paths share the
  fixed pool and exact threshold, but use different wall geometries,
  arithmetic representations, storage, and transversal proofs.
---

# THM-4179 -- seventh-deletion completion and maximal-blocker duality

**PROVED RELATIVE TO THM-4150/4178 + VERIFIED-EXACT; LRC(14) REMAINS OPEN.**

## 1. Two monotone deletion lemmas

Let `V` be a finite labelled set, and let `g(R)` be a real-valued function on
`R subset V` which is monotone under deletion:

```text
R subset R'  implies  g(R)<=g(R').                    (1)
```

At a threshold `theta`, put

```text
E_d={R in binom(V,d):g(R)>=theta}.                     (2)
```

This abstraction has two exact consequences.

> **One-step blocker descent.** Suppose `d+k<|V|`. Every transversal of
> `E_(d+1)` of size at most `k` is a transversal of `E_d`.

Indeed, if `T` misses some `R in E_d`, choose
`x in V\(R union T)`. Such an `x` exists by the strict cardinality
hypothesis. Monotonicity makes `R union {x}` a member of `E_(d+1)`, and that
edge also misses `T`, a contradiction. Consequently, after all size-`k`
blockers of `E_d` have been enumerated, proving that every one misses an
edge of `E_(d+1)` rules out every transversal of `E_(d+1)` through size `k`.

There is also an all-arity endpoint. In the THM-4178 notation, fix an anchor
`A`, newcomer `q`, optional set `V=P\A`, and body choice `K subset V`. Write

```text
g_q(R)=mu(G_((P union {q})\R)).                        (3)
```

> **Maximal-deletion duality.** A threshold repair `R subset V\K` exists if
> and only if the maximal deletion `V\K` is itself a threshold repair. In
> that case
>
> ```text
> g_q(V\K)=mu(G_(A union {q} union K)).                (4)
> ```

The forward direction is `(1)`; the reverse direction takes `R=V\K`.
Equation `(4)` is the literal identity of the two remaining speed sets. Thus
a shallow transversal is only a failure of a bounded-deletion certificate.
It is not a body obstruction: after all deletion depths are allowed, the
repair criterion is exactly the body's own Haar threshold.

These lemmas are the mechanism of the theorem. The closest inherited result
is THM-4178's exact depth-six profile. Its two seven-covers are the canonical
hostiles. The corrected near miss is “a minimum shallow cover blocks the
body”; maximal-deletion duality refutes that implication. The least-used
sidecar is the complete labelled blocker, retained here across the next
deletion operation.

## 2. Exact depth-seven crossing

Retain the THM-4178 pool and newcomer

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},                       q=50.           (5)
```

For a primitive divisor-complete anchor `A subset P`, put `O_A=P\A` and

```text
E_d^A={R in binom(O_A,d):
       mu(G_((P union {50})\R))>=4/63}.                (6)
```

THM-4178 proves there are exactly six such anchors. Four already have
`tau(E_6^A)=8`. The two remaining rows have depth-six transversal number
seven. Their next layers are exactly:

| anchor `A` | `|E_6^A|` | all depth-six seven-blockers | `|E_7^A|` | `tau(E_7^A)` | one eight-cover |
|:---|---:|:---|---:|---:|:---|
|`{80,143,252}`|`33,636`|`{16,88,95,145,168,193,240}`|`298,279`|`8`|`{8,85,88,95,145,168,193,240}`|
|`{143,240,252}`|`31,234`|`{16,85,88,95,145,168,193}`; `{16,80,88,95,145,168,193}`; `{16,88,95,145,168,170,193}`|`286,291`|`8`|`{8,80,85,88,95,145,168,193}`|

There are zero threshold equalities in both depth-seven rows. The complete
blocker descent is witnessed by the following missed depth-seven repairs:

| anchor | depth-six blocker `K` | one `R in E_7^A` disjoint from `K` |
|:---|:---|:---|
|`{80,143,252}`|`{16,88,95,145,168,193,240}`|`{8,10,15,20,85,120,290}`|
|`{143,240,252}`|`{16,85,88,95,145,168,193}`|`{8,10,15,80,120,190,290}`|
|same|`{16,80,88,95,145,168,193}`|`{8,10,15,20,85,120,290}`|
|same|`{16,88,95,145,168,170,193}`|`{8,10,15,20,80,85,290}`|

Here `|O_A|=27`, `d=6`, and `k=7`, so `d+k<27`. Section 1 shows that any
seven-cover of `E_7^A` would be one of the displayed depth-six blockers; the
right column rules out every one. The displayed eight-covers prove the
matching upper bounds. Hence both transversal numbers are exactly eight.
Since THM-4178 gives transversal number seven at depth six, deletion depth
seven is sharp for both rows in this uniform anchor-preserving mechanism.

## 3. Transfer for all six anchors

Let `A` be any of the six primitive divisor-complete anchors classified in
THM-4178, let `K in binom(O_A,7)`, and put

```text
H_(A,K)=A union {50} union K.                          (7)
```

For four anchors, THM-4178 gives a depth-six repair disjoint from every `K`.
For the two rows above, `tau(E_7^A)=8>7` gives a depth-seven repair disjoint
from every `K`. In either case some `R` satisfies

```text
H_(A,K) subset (P union {50})\R,
mu(G_(H_(A,K)))>=mu(G_((P union {50})\R))>=4/63.       (8)
```

Therefore, for every positive integer `c` and every two distinct positive
odd integers `a,b`, THM-4150 supplies an `x in R/Z` with

```text
min_(v in 2cH_(A,K) union {a,b}) ||vx|| >= 1/14.       (9)
```

Every `H_(A,K)` is primitive and divisor-complete because `A` is. The labels
are distinct because `50 notin P`, and doubling separates the even body from
the odd tails. The content map preserves Haar measure exactly, as in
THM-4178. This proves `(9)` for all six anchors. **QED relative to THM-4150.**

The four depth-six blockers are visibly not body obstructions. Direct maximal
deletion gives their exact body masses, in the order displayed in Section 2:

```text
152056046161/894174481200,
3868553773463/22801449270600,
152056046161/894174481200,
7530508986601/45602898541200.                         (10)
```

All four are strictly above `4/63`; the smallest is about `0.16513`. This is
the hostile control promised by maximal-deletion duality: a depth-six
seven-cover can coexist with a body mass more than twice the transfer
threshold.

## 4. Distinct bodies, not anchor presentations

Let the original anchor be

```text
A_0={120,126,143}.                                    (11)
```

A core in the present primitive-anchor lane has exactly one member of `A_0`
precisely when it contains `143`, contains neither `120` nor `126`, and its
nine labels from `P\A_0` include

```text
252 and at least one of {40,80,240}.                  (12)
```

The three choices in `(12)` correspond to the anchors
`{40,143,252}`, `{80,143,252}`, and `{143,240,252}`. Each anchor has
`binom(25,7)=480,700` presentations. Their pairwise intersections each have
size `binom(24,6)=134,596`, and their triple intersection has size
`binom(23,5)=33,649`. Hence the exact number of distinct labelled core bodies
is

```text
3 binom(25,7)-3 binom(24,6)+binom(23,5)
 =binom(26,8)-binom(23,8)
 =1,071,961.                                           (13)
```

The literal presentation-multiplicity histogram is

```text
one anchor:   735,471 bodies,
two anchors:  302,841 bodies,
three anchors: 33,649 bodies.                         (14)
```

Thus the `1,442,100=3*480,700` anchor presentations in `(13)` must not be
summed as bodies. THM-4178 had already transferred the `{40,143,252}` row;
the completed union adds exactly

```text
1,071,961-480,700=591,261                             (15)
```

distinct scale-one cores beyond that row. Every core in `(13)` has exactly
one original anchor, whereas the fixed-`q=50` THM-4174/4175 anchor families
have three or two. Thus `(13)` is their disjoint exactly-one-original-anchor
slice. This is a labelled fixed-pool core count, not a claim of global
novelty against every canonical LRC family. Primitivity makes its positive
content towers internally distinct.

## 5. Exact audits

The primary Python path verifies the frozen THM-4178 source hash, rebuilds its
`7,134`-wall pool geometry and complete labelled failure atoms, and evaluates
every depth-six and depth-seven deletion by exact subset-zeta sums. It stores
edge incidence as Python big integers. At each state, the first uncovered
edge must contribute a chosen vertex; branching over that edge and
deduplicating the chosen mask is an exhaustive cover enumeration. It finds
exactly one and three depth-six seven-blockers, checks the four missed
depth-seven repairs above, and verifies both eight-covers. It also computes
the four maximal body masses and literally enumerates the histogram `(14)`.

The independent C++20 path instead inserts the newcomer walls into the joint
arrangement

```text
D=91,205,797,082,400,
7,214 walls and 7,213 open atoms.                     (16)
```

It directly sums signed-128-bit joint-atom lengths for all
`2*binom(27,7)=1,776,060` deletion masks. Its direct cover recursion searches
budgets zero through eight, without using the depth-six blocker list. It
rejects budget seven after `31,578` and `46,670` search nodes and accepts the
displayed eight-covers. The labelled joint-atom FNV-1a-64 control remains
`0ccd305ae47c79ea`; the depth-seven edge hashes are respectively
`ad0ce6d3119dda70` and `2fcc89f5b9d8a238`.

The two paths agree on both edge counts, zero equality counts, exact
transversal numbers, covers, and body-union arithmetic. They share the fixed
pool and threshold but not wall geometry, arithmetic representation, edge
storage, or the proof that no seven-cover survives.

## 6. Scope, connection contract, and replay

```text
source:       two q=50 depth-six seven-cover repair rows
target:       all six primitive divisor-complete anchor families at q=50
map:          blocker descent -> d=7 missed repair -> THM-4150
preserved:    q, anchor labels, exact 4/63 threshold, content, divisor pins
destroyed:    canonical safe phase and selected repair after existence
sidecar:      complete labelled d=6 blocker list through the d=7 operation
positive:     both d=7 hypergraphs have transversal number exactly eight
hostile:      every d=6 blocker body itself has Haar mass above 0.165
decisive test: 1,776,060 d=7 masks plus independent exact cover proofs. (17)
```

This theorem is `q=50`-only and fixed-pool-only. It treats primitive
divisor-complete anchor triples, uniform deletions, doubled bodies, and
distinct positive odd tails. It proves no all-`q` third-anchor theorem, no
arbitrary-body classification, no physical entry, no mixed/even-tail result,
and no necessity of the Haar threshold. In particular, it does not prove
LRC(14).

Primary replay:

```bash
python3 -B 04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179.py
python3 -B -O 04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179.py
PYTHONHASHSEED=4179 python3 -B \
  04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179.py
```

Independent replay:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp \
  -o /tmp/lrc4179-o2
/tmp/lrc4179-o2
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp \
  -o /tmp/lrc4179-o0
/tmp/lrc4179-o0
g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined \
  04-computation/lrc14_q50_seventh_deletion_primitive_anchor_completion_thm4179_independent_audit.cpp \
  -o /tmp/lrc4179-ubsan
/tmp/lrc4179-ubsan
```

All declared streams byte-match their frozen outputs. **QED.**
