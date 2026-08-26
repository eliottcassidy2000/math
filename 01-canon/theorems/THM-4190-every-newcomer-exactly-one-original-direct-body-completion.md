---
id: THM-4190
title: "Every-newcomer exactly-one-original direct-body completion"
status: >
  PROVED RELATIVE TO THM-4150/4174/4175/4178/4179/4188/4191 + VERIFIED-EXACT PRIMARY
  POOL-WALL FAILURE-SUBSET ZETA AND INDEPENDENT JOINT-WALL SAFE-SUPERSET
  ZETA AUDITS; LRC(14) OPEN. For every positive newcomer q outside the fixed
  thirty-label pool, all 1,071,961 bodies in THM-4179's exactly-one-original
  slice have complete safe-set Haar mass at least 4/63. THM-4188 transports
  the correctly matched q=50 E6/E7 anchor layers outside its 23 resonances;
  direct maximal-deletion coordinates close all 33,168,300 anchor-labelled
  resonance presentations, strictly and without equality. Together with
  THM-4174/4175/4191 this gives 17,056,501 pairwise-disjoint named cores per
  newcomer, not LRC(14).
source: lrc-one-original-deformation-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
  - THM-4191-complete-all-body-zero-original-newcomer-haar-transfer
related:
  - THM-4185-q50-complete-zero-original-minimal-anchor-hierarchy
script: 04-computation/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.cpp
output: 05-knowledge/results/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.out
independent_audit_script: 04-computation/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.cpp
independent_audit_output: 05-knowledge/results/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.out
script_sha256: 34ce4372e72ffac34b0e4f127b9b326ef85e7b3147243c6506cd134963dfa5a6
output_sha256: cf326163c535d645cab6ce9128165c176f40199325fd1be9f763bf0ecaa634f8
independent_audit_script_sha256: 550dd2ae3f5b9b0cbea294e7cc71b36393012d6fbfe0a444fd7adc50452d5e37
independent_audit_output_sha256: f95de26dd0c5ba8c49503fc92ebe3110bd853d67e23cbf8b0cb0c765cdee3ac1
hash_basis: raw LF bytes
primary_audit: >
  PASS. A warning-clean signed-integer C++20 path retains the 7,133 fixed
  pool cells, integrates each newcomer comb by an exact prefix primitive,
  projects failure masks to each 25-label optional ground, and uses a dense
  subset zeta to inspect all C(25,7) complements. It checks all 69
  anchor/resonance rows, 33,168,300 presentations, the unique-body formulas,
  zero failures, zero equalities, and a locked bytewise semantic ledger.
independent_audit: >
  PASS. A separate warning-clean C++20 path explicitly constructs each full
  P-union-{q} joint wall refinement, reclassifies every pool label on every
  atom, projects safe masks, and applies the dual superset zeta. It directly
  re-sweeps every one of the 69 reported minimizers runner by runner. Both
  paths reproduce every minimum and the common ledger 72db7e3090d5626e.
---

# THM-4190 -- every-newcomer exactly-one-original direct-body completion

**PROVED RELATIVE TO THM-4150/4174/4175/4178/4179/4188/4191 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Exact theorem and body universe

Retain the fixed pool and original anchor

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143}.                                      (1)
```

For a finite positive set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                (2)
```

The exactly-one-original slice completed by THM-4179 at `q=50` is

```text
H_1={B in binom(P,10):
     B intersect A_0={143},
     252 in B,
     B intersects {40,80,240} is nonempty}.             (3)
```

Equivalently, after fixing `{143,252}`, choose eight labels from the remaining
twenty-six and require at least one of `{40,80,240}`. Hence

```text
|H_1|=binom(26,8)-binom(23,8)=1,071,961.                (4)
```

> **Every-newcomer theorem.** For every positive integer `q notin P` and
> every `B in H_1`,
>
> ```text
> mu(G_(B union {q}))>=4/63.                            (5)
> ```
>
> Consequently, for every positive integer `c` and every two distinct
> positive odd integers `a,b`, some `x in R/Z` satisfies
>
> ```text
> min_(v in 2c(B union {q}) union {a,b})||vx||>=1/14.   (6)
> ```

The three primitive divisor-complete anchor triples used below are

```text
A_40 ={40,143,252},
A_80 ={80,143,252},
A_240={143,240,252}.                                    (7)
```

To avoid mistaking anchor presentations for bodies, choose an anchor by the
deterministic rule

```text
40 first; if 40 is absent, 80; otherwise 240.           (8)
```

The three disjoint priority classes have sizes

```text
binom(25,7)=480,700,
binom(24,7)=346,104,
binom(23,7)=245,157,                                    (9)
```

which sum to `(4)`. Auditing all three full anchor rows instead produces

```text
3 binom(25,7)=1,442,100                                 (10)
```

presentations. Equations `(4)` and `(10)` are deliberately different.

## 2. Inheritance and the nonresonant proof

For `d in {6,7}` define the global repair layer

```text
E_d(q)={R in binom(P,d):
        mu(G_((P union {q})\R))>=4/63},                 (11)

E_d^A(q)={R in E_d(q):R intersect A=empty}.             (12)
```

For `q notin P`, THM-4188 proves the exact resonance filtration

```text
E_6(50) subset E_6(q)  iff  q notin Q_6,
E_7(50) subset E_7(q)  iff  q notin Q_7,
Q_6 subset Q_7,                                        (13)

Q_7={6,22,24,25,48,70,72,96,100,105,110,128,130,140,
     186,192,206,210,220,256,260,294,366}.              (14)
```

The q=50 anchor rows must be matched to their proved depths:

```text
tau(E_6^A40(50))=8                         (THM-4178),
tau(E_7^A80(50))=tau(E_7^A240(50))=8       (THM-4179).  (15)
```

Their edge counts are respectively

```text
44,562,                 298,279,                 286,291. (16)
```

Fix `q notin P union Q_7`, `B in H_1`, choose `A` by `(8)`, and put `K=B\A`.
Then `|K|=7`. If `A=A_40`, use the first inclusion in `(13)` and
`Q_6 subset Q_7`; otherwise use the second inclusion. Restricting an edge
inclusion to repairs disjoint from `A` preserves it, so `(15)` implies that
`K` cannot meet every q-repair in the matched layer. Some `R in E_d^A(q)`
is disjoint from `K`, hence from all of `B`, and

```text
B union {q} subset (P union {q})\R,
G_((P union {q})\R) subset G_(B union {q}).             (17)
```

Safe-set monotonicity and `(11)` prove `(5)` for every nonresonant newcomer.
This imports THM-4188's finite-exact and analytic cofinal control without
claiming a new edgewise inclusion theorem.

## 3. Direct maximal-deletion coordinates at the 23 resonances

The closest inherited mechanism is THM-4188's edge-inclusion filtration. Its
canonical hostile is already `q=6`, where lawful q=50 repair edges disappear.
The corrected near miss is that loss of a bounded-depth repair edge makes a
body unsafe. THM-4179's maximal-deletion duality says the all-depth endpoint
is the body's own Haar mass. The least-used sidecar is therefore the complete
labelled body, retained here instead of another shallow repair count.

Put

```text
D=lcm_(p in P)(14p)=18,241,159,416,480.                 (18)
```

The primary audit partitions the circle at the `7,134` pool walls. On each
of the `7,133` open cells the complete failure mask `F_P` is constant. For a
wall tick `t`, write `qt=kD+r`, `0<=r<D`, and use the exact q-safe prefix

```text
J_q(t)=12kD+clamp(14r-D,0,12D).                        (19)
```

Thus each cell receives a nonnegative integer q-safe weight with total
`12qD`, on denominator `14qD`.

For one anchor `A` let

```text
W_A=P\(A union {120,126}),                    |W_A|=25. (20)
```

Discard cells on which an anchor label fails. Project each remaining failure
mask to `W_A`, and aggregate the exact cell weights as `w_(q,A)(F)`. For every
`K in binom(W_A,7)`, the body numerator is exactly

```text
N_(q,A)(K)=sum_(F subset W_A\K) w_(q,A)(F),             (21)

mu(G_(A union K union {q}))=N_(q,A)(K)/(14qD).          (22)
```

Equation `(21)` is one dense failure-mask subset zeta on `2^25` states. The
threshold comparison is the signed integer

```text
delta_(q,A)(K)=9N_(q,A)(K)-8qD,                        (23)

delta>=0  iff  mu(G_(A union K union {q}))>=4/63.       (24)
```

The exact audit evaluates every one of the

```text
23*3*binom(25,7)=33,168,300                            (25)
```

anchor-labelled presentations. All `69` rows have zero negative deltas and
zero equalities. The minimum over the complete resonance audit occurs at

```text
q=140,
A=A_80,
K={10,42,60,85,95,145,264},
B={10,42,60,80,85,95,143,145,252,264},                 (26)

N=5,459,802,690,013,840,
14qD=35,752,672,456,300,800,
delta=28,708,125,663,666,960,                          (27)

mu(G_(B union {140}))=0.152710337854804... .           (28)
```

This is a strict positive control, not a claimed equality or a global sharp
constant. Equations `(21)--(25)` prove `(5)` for every `q in Q_7`. Together
with Section 2 they prove `(5)` for every `q notin P`.

## 4. Odd-tail transfer and uniform four-slice synthesis

Every `B in H_1` contains one of the primitive divisor-complete anchors in
`(7)`. More importantly, THM-4150 applies to every nonempty finite positive
body once `(5)` is known; primitivity is not an extra transfer hypothesis.
Haar invariance under multiplication by `c` and THM-4150 prove `(6)`.
**QED relative to THM-4150.**

THM-4174, THM-4175, this theorem, and THM-4191 now give four pairwise-disjoint
fixed-pool slices for every newcomer, distinguished by the number of labels
from `A_0`:

```text
all three originals (THM-4174):          888,030,
exactly two originals (THM-4175):      6,660,225,
THM-4179 exactly-one slice (this):      1,071,961,
zero originals (THM-4191):             8,436,285,
                                       ----------
uniform four-slice total:              17,056,501.      (29)
```

THM-4188's `1,491,665` good zero-original bodies remain a distinguished
primitive divisor-complete subcount of THM-4191's `8,436,285`; they are not
added again in `(29)`. Before this theorem, the exactly-one row and hence the
four-slice total were `q=50`-only. Equation `(29)` is now uniform for every
positive `q notin P`.

## 5. Independent audit

The independent program does not use the fixed-cell prefix integral. For each
`q in Q_7` it constructs the full joint wall denominator

```text
D_q=lcm(D,14q)                                         (30)
```

and merges every pool and newcomer wall. It reclassifies all thirty pool
labels on every joint atom. After retaining q-safe and anchor-safe atoms, it
projects the **safe** label mask `S subset W_A` and aggregates its length as
`v_(q,A)(S)`. The dual superset zeta gives

```text
M_(q,A)(K)=sum_(S superset K) v_(q,A)(S),              (31)
```

the same body mass in joint-wall coordinates. This is independent of the
primary failure-mask formula `(21)`. The program then remeasures every one of
the `69` reported minimizers by a literal joint-atom sweep, testing `q` and
all ten body labels afresh. Every direct mass equals its zeta value.

The implementations agree on every row minimum, minimizing labelled body,
zero negative/equality counts, and the standard bytewise FNV-1a-64
little-endian semantic ledger

```text
72db7e3090d5626e.                                      (32)
```

Primary O3, O0, and UBSan streams byte-match; independent O2, O0, and UBSan
streams byte-match.

## 6. Connection contract and scope firewall

```text
source:       q=50 anchor repair layers plus the complete labelled bodies
target:       every q outside P in THM-4179's exactly-one slice
map:          nonresonant edge inclusion; at resonances, failure-mask
              disjointness zeta / maximal-deletion body mass
preserved:    q, anchor and body labels, exact 4/63 threshold, content,
              divisor pins, and odd-tail quantifiers
destroyed:    a selected repair and safe phase in the direct resonance route;
              the anchor presentation after taking the body union
sidecar:      deterministic 40/80/240 priority anchor and complete masks
positive:     all 33,168,300 resonance presentations have strict surplus
hostile:      q=6 loses q=50 edges; q=140 carries the smallest body mass
decisive test: two dual zeta transforms plus 69 literal minimizer sweeps. (33)
```

This theorem completes exactly the body universe `(3)`. It does not certify
every ten-body containing exactly one member of `A_0`, bodies outside this
fixed pool, physical entry of a hypothetical counterexample, mixed/even tail
branches, or necessity of the Haar threshold. The `17,056,501` cores in
`(29)` are named pairwise-disjoint fixed-pool families, not all possible
thirteen-speed sets. This theorem does not prove LRC(14).

## 7. Exact artifacts and replay

Primary:

```text
04-computation/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.cpp
sha256 34ce4372e72ffac34b0e4f127b9b326ef85e7b3147243c6506cd134963dfa5a6

05-knowledge/results/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.out
sha256 cf326163c535d645cab6ce9128165c176f40199325fd1be9f763bf0ecaa634f8
```

Independent:

```text
04-computation/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.cpp
sha256 550dd2ae3f5b9b0cbea294e7cc71b36393012d6fbfe0a444fd7adc50452d5e37

05-knowledge/results/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.out
sha256 f95de26dd0c5ba8c49503fc92ebe3110bd853d67e23cbf8b0cb0c765cdee3ac1
```

Replay:

```bash
clang++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.cpp \
  -o /tmp/lrc4190-primary
/tmp/lrc4190-primary | diff -u \
  05-knowledge/results/lrc14_all_newcomer_one_original_direct_body_completion_thm4190.out -

clang++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.cpp \
  -o /tmp/lrc4190-independent
/tmp/lrc4190-independent | diff -u \
  05-knowledge/results/lrc14_all_newcomer_one_original_joint_wall_audit_thm4190.out -
```

**QED.**
