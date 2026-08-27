---
id: THM-4240
title: "Fixed-fifty four-petal cofinal atlas and twenty-two-label chart"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191/4211/4229/4234 +
  FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14) OPEN. All 495 four-petal
  limiting charts are strict over 4,241,160 labelled profiles, giving a
  uniform cofinal tail for r>=608. The selected four-petal set
  {63,132,176,264} has all three newly exposed triple finite heads and its
  quadruple finite head exhausted, proving chi_50>=22. The independent
  atlas has 990 direct extremizer replays; the selected chart has 27,159,132
  literal comparisons and twelve direct replays. Arbitrary entry and LRC(14)
  remain OPEN. MISTAKE-524 withdraws a prior forward claim about reserved
  THM-4242; no result here depends on that empty stub.
source: codex-lrc14-breakthrough-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4229-fixed-fifty-nineteen-label-petal-haar-charts
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
related:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4238-forty-small-label-uniform-r590-haar-tail-closure
primary_script: 04-computation/lrc14_fixed_fifty_four_petal_haar_chart_thm4240.cpp
atlas_independent_script: 04-computation/lrc14_fixed_fifty_all_quadruple_limiting_independent_audit_thm4240.cpp
selected_independent_script: 04-computation/lrc14_fixed_fifty_selected_four_petal_independent_audit_thm4240.cpp
atlas_primary_output: 05-knowledge/results/lrc14_fixed_fifty_all_quadruple_limiting_primary_thm4240.out
atlas_independent_output: 05-knowledge/results/lrc14_fixed_fifty_all_quadruple_limiting_independent_audit_thm4240.out
selected_primary_output: 05-knowledge/results/lrc14_fixed_fifty_selected_four_petal_primary_thm4240.out
selected_independent_output: 05-knowledge/results/lrc14_fixed_fifty_selected_four_petal_independent_audit_thm4240.out
primary_script_sha256: e42ecf5a3c248cdab6c7f87b528674be2025293729bd8a972eadc548fefa3a52
atlas_independent_script_sha256: b9b07ba98b781f1e5a0af1ffc8d6fb301c9b5a18ec2741d5407a16fd2aacd9b8
selected_independent_script_sha256: 1ad9abf86efdb4ba05119beae95fdca5ef92f7934bda6be912ffeeb850410855
atlas_primary_output_sha256: 5a8aaa169c56e37849b735544d0047daf388e6e386fda194af1121832a956225
atlas_independent_output_sha256: b238ebaf3a957adb5ed2f0ec3b14161fc1d6e0ab3335f30366ad3ff7cc49d7fe
selected_primary_output_sha256: 448421efd60631545b5289d2baef3b492f5f2ade624130b6ca0b76720b8cc349
selected_independent_output_sha256: cf73b7bd3b46b2c70530973cc12caf7cee7558beca38f09111db6196c2f75263
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary midpoint-cell/core-subset-zeta engine and the
  independent grouped-event/petal-batched complement-scatter engine agree on
  all 495 consequence-bearing atlas records. Eleven displayed cutoff bodies
  differ because tied maximizers are traversed in opposite orders, and their
  body-dependent component counts differ too. Every consequence-bearing field
  agrees. A second grouped-event implementation agrees field-for-field on the
  four selected finite heads and directly replays every minimum, cutoff, and
  closest literal body.
---

# THM-4240 -- fixed-fifty four-petal cofinal atlas and twenty-two-label chart

**PROVED RELATIVE TO THM-4150/4156/4170/4191/4211/4229/4234 +
FINITE-EXACT + INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

For a finite positive set `S`, write

```text
G_S={x in R/Z:min_(s in S)||sx||>=1/14},
alpha=4/63.                                             (1)
```

Retain the THM-4234 partition

```text
C={8,16,40,42,80,84,85,88,95,
   120,126,143,145,168,193,240,252,286},

O={10,15,20,30,60,63,132,170,176,190,264,290},
P=C union O.                                            (2)
```

> **Four-petal cofinal atlas.** For every `T in binom(O,4)`, every integer
> `r>=608`, and every `K in binom(C union T,9)`,
>
> ```text
> mu(G_(K union {50,r}))>=alpha.                        (3)
> ```
>
> **Selected universal chart.** Put
>
> ```text
> U={63,132,176,264}.                                   (4)
> ```
>
> Equation `(3)` holds with `T=U` for every positive
> `r notin P union {50}`.

Consequently the fixed-fifty chart number of THM-4211 satisfies

```text
chi_50>=22.                                             (5)
```

The closest proved mechanism is THM-4234's complete triple cofinal atlas and
its one universal triple. The canonical hostile is THM-4207: individually
lawful newcomer marginals need not compose. The corrected near miss is
MISTAKE-520: inherited zero-newcomer faces come from full-pool heredity, not
from an unrelated cardinality certificate. The least-used decisive sidecar is
the circular component count of each literal base safe set.

The live concept board was

```text
triple hypergraph | quadruple faces | base surplus | component tariff
finite phase address | Pascal inheritance.             (6)
```

## 2. Exact tail mechanism

Fix `T in binom(O,4)` and `L in binom(C,5)`. Put

```text
V_(T,L)=G_(L union T union {50}),
mu(V_(T,L))=m_(T,L)/D,
D=91,205,797,082,400,
c_(T,L)=# positive-length circular components.         (7)
```

Endpoint-only components have measure zero and are discarded. THM-4170
gives, for every positive integer `r`,

```text
mu(V_(T,L) intersect G_r)
 >=(6/7)(m_(T,L)/D)-6c_(T,L)/(49r).                   (8)
```

Define

```text
delta_(T,L)=54m_(T,L)-4D,
kappa(T,L)=ceil(54c_(T,L)D/(7delta_(T,L)))             (9)
```

when `delta_(T,L)>0`. Then `(8)--(9)` imply the target for every
`r>=kappa(T,L)`. The ceiling is a sufficient analytic cutoff, not the first
literal safe outsider.

The new four-petal universe is exactly

```text
binom(12,4)binom(18,5)
=495*8,568
=4,241,160.                                            (10)
```

Every profile has positive `delta`. The complete atlas has

```text
strict quadruples                  495 / 495
nonstrict profiles                  0
largest sufficient cutoff         608
cutoff quadruple              {170,190,264,290}
cutoff core                 {168,193,240,252,286}.      (11)
```

At the extremal profile in `(11)`,

```text
m=18,169,518,380,550,
c=532,
delta=616,330,804,220,100.                            (12)
```

The maximin-surplus quadruple is

```text
{15,20,30,60},
min_L delta=723,857,556,028,932,
max_L kappa=320.                                       (13)
```

This is a maximin label, not an assertion that its finite head is uniquely
best. The already-audited zero-, one-, two-, and three-petal layers have
cutoffs at most `589`; hence `(11)` proves the full statement `(3)`, not only
its four-petal layer.

## 3. The selected twenty-two-label chart

For `U` in `(4)`, the triple `{132,176,264}` is universal by THM-4234. The
only proper faces not already covered by THM-4211, THM-4229, THM-4234's pair
theorem, and that triple are

```text
{63,132,176}, {63,132,264}, {63,176,264}, U.            (14)
```

Their exact ledgers are:

| required petals | core size | cutoff | literal rows | checks | closest `r` |
|:---|---:|---:|---:|---:|---:|
| `{63,132,176}` | 6 | 477 | 445 | 8,260,980 | 48 |
| `{63,132,264}` | 6 | 427 | 395 | 7,332,780 | 35 |
| `{63,176,264}` | 6 | 475 | 443 | 8,223,852 | 96 |
| `{63,132,176,264}` | 5 | 422 | 390 | 3,341,520 | 256 |

Across `(14)` there are

```text
limiting profiles       3binom(18,6)+binom(18,5)=64,260,
literal rows                                      1,673,
literal checks                               27,159,132,
failures / equalities                              0 / 0. (15)
```

The closest quadruple row is

```text
r=256,
L={8,80,143,168,193},
63mu(G_(L union U union {50,256}))-4
 =9,384,395,987,520,792 / 1,459,292,753,318,400>0.    (16)
```

Thus every nine-body in `C union U` is safe with `50` and every permitted
outsider. The internal layer identity is

```text
sum_(j=0)^4 binom(4,j)binom(18,9-j)
=48,620+175,032+190,944+74,256+8,568
=497,420=binom(22,9).                                  (17)
```

Together with the inherited zero- and one-newcomer layers, every eleven-face
of `(C union U) union {50,r}` is safe:

```text
binom(22,11)+2binom(22,10)+binom(22,9)
=705,432+1,293,292+497,420
=2,496,144=binom(24,11).                               (18)
```

Equations `(15)--(18)` prove `(5)`. THM-4150 supplies universal distinct
odd-tail completions after common positive scaling; it does not supply
physical entry of an arbitrary LRC(14) row.

## 4. Independent audit architecture

The primary implementation constructs the exact common wall grid, evaluates
midpoint safety on each open cell, indexes failure masks on the eighteen-label
core, and applies an eighteen-bit subset-zeta transform separately for each
petal requirement.

The atlas referee reverses both choices. It groups simultaneous wall events,
tracks exact safe/failing state through event toggles, batches petal lanes, and
uses complement scatter. It replays the minimum and cutoff profile of every
quadruple directly, for `990` body-local controls. It agrees with every
consequence-bearing field. Eleven cutoff-body displays choose the other member
of a tied maximizing fibre, so their body-dependent component counts also
differ. The cutoffs and all consequence-bearing fields agree.

The selected-chart referee is a third grouped-event implementation. It audits
the three triple lanes and one quadruple lane simultaneously, then directly
replays two base extremizers and one closest literal body per lane. It agrees
field-for-field and obtains zero failures and zero equalities.

## 5. Reproduction and scope firewall

From the repository root:

```bash
g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_four_petal_haar_chart_thm4240.cpp \
  -o /tmp/lrc14-thm4240-primary
/tmp/lrc14-thm4240-primary --all-petal-base 4
OMP_NUM_THREADS=4 /tmp/lrc14-thm4240-primary --selected-four-petal

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_all_quadruple_limiting_independent_audit_thm4240.cpp \
  -o /tmp/lrc14-thm4240-atlas-independent
OMP_NUM_THREADS=4 /tmp/lrc14-thm4240-atlas-independent

g++ -std=c++20 -O3 -DNDEBUG -fopenmp \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_fifty_selected_four_petal_independent_audit_thm4240.cpp \
  -o /tmp/lrc14-thm4240-selected-independent
OMP_NUM_THREADS=4 /tmp/lrc14-thm4240-selected-independent
```

Compare each stream with its named frozen output. No program uses floating
point or sampling.

This theorem does not prove that `608` is a minimal literal tail, that
`chi_50=22`, replacement of the fixed center `50`, physical entry, or
LRC(14). MISTAKE-524 records and withdraws the former sentence treating the
still-reserved THM-4242 as a strengthening. **QED.**
