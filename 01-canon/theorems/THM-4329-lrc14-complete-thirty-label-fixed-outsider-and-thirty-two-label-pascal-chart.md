---
id: THM-4329
title: "Complete thirty-label fixed-outsider and thirty-two-label Pascal chart"
status: >
  PROVED RELATIVE TO THM-4156/4191/4326 + INDEPENDENTLY AUDITED.
  For every two distinct positive outsiders q,r of the fixed thirty-label
  pool, all C(32,11)=129,024,480 eleven-subsets of P union {q,r} have
  1/14-safe Haar mass at least 4/63. Equivalently, every fixed outsider has
  the maximum possible thirty-label universal nine-body chart; in the
  existing notation chi_50=30. THM-4150 transfers every eleven-face to all
  distinct positive odd-tail pairs after arbitrary common scaling. This is a
  fixed-chart completion, not arbitrary-row entry or LRC(14).
source: root + entry_scout / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
related:
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
  - THM-4240-fixed-fifty-twenty-two-label-four-petal-haar-chart
audit: >
  PASS / ACCEPT. An independent quantifier audit checked the THM-4211 chart
  definition, outsider exclusions, every zero/one/two-outsider layer,
  collision boundaries, binomial counts, current MISTAKE entries, and all
  claimed nonconsequences. The proof is a formal corollary of the displayed
  dependencies and uses no new numerical experiment.
---

# THM-4329 -- complete fixed-outsider and Pascal charts

**PROVED RELATIVE TO THM-4156/4191/4326 + INDEPENDENTLY AUDITED.
LRC(14) REMAINS OPEN.**

## 1. Statement

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
alpha=4/63.                                            (1)
```

> **Complete thirty-two-label Pascal chart.** Let `q,r` be distinct positive
> integers outside `P`. For every
>
> ```text
> A in binom(P union {q,r},11),
> ```
>
> one has
>
> ```text
> mu(G_A)>=alpha.                                      (2)
> ```

For a fixed positive outsider `q`, define the pool-relative universal chart
number

```text
chi_q(P)=max{|C|:C subset P and
  mu(G_(K union {q,r}))>=alpha
  for every K in binom(C,9) and every positive
  r notin P union {q}}.                                (3)
```

Then

```text
chi_q(P)=30                    for every q notin P.    (4)
```

At `q=50`, equation `(3)` is exactly THM-4211's definition of `chi_50`.
Consequently

```text
chi_50=30.                                             (5)
```

This is the largest value allowed by the definition, because every chart is
a subset of the thirty-label pool.

Finally, for every `A` in `(2)`, every positive integer `c`, and every two
distinct positive odd integers `a,b`, there is an `x in R/Z` such that

```text
min_(v in 2cA union {a,b})||vx||>=1/14.                (6)
```

Equation `(6)` is a structured thirteen-speed family. It is not an entry
theorem for an arbitrary thirteen-speed row.

## 2. The three complete outsider layers

Fix `A in binom(P union {q,r},11)` and put

```text
j=|A intersect {q,r}|.                                (7)
```

There are exactly three cases.

### Layer zero

If `j=0`, then `A subset P`. THM-4156 proves the hereditary strict bound

```text
mu(G_A)>=mu(G_P)>alpha.                                (8)
```

No cardinality or anchor condition is needed for this hereditary
consequence.

### Layer one

If `j=1`, write

```text
A=B union {s},       B in binom(P,10),
s in {q,r}.                                            (9)
```

The label `s` is a positive outsider of `P`. THM-4191 gives

```text
mu(G_(B union {s}))>=alpha.                           (10)
```

### Layer two

If `j=2`, write

```text
A=B union {q,r},     B in binom(P,9).                 (11)
```

The outsiders are distinct and positive, so THM-4326 applies and gives

```text
mu(G_(B union {q,r}))>=alpha.                         (12)
```

Equations `(8)`, `(10)`, and `(12)` exhaust every eleven-face and prove
`(2)`.

The layer counts are

```text
j=0: binom(30,11)   = 54,627,300,
j=1: 2binom(30,10) = 60,090,030,
j=2: binom(30,9)   = 14,307,150.                     (13)
```

Thus

```text
54,627,300+60,090,030+14,307,150
 =129,024,480
 =binom(32,11),                                       (14)
```

the Vandermonde decomposition of the complete Pascal layer. No face is
discarded or counted twice.

## 3. Exact chart-number consequence

Fix `q notin P`. For every positive `r notin P union {q}` and every
`K in binom(P,9)`, THM-4326 gives

```text
mu(G_(K union {q,r}))>=alpha.                         (15)
```

Hence the choice `C=P` is admissible in `(3)`, so `chi_q(P)>=30`. The
reverse inequality `chi_q(P)<=30` is tautological from `C subset P`. This
proves `(4)` and `(5)`.

In particular, THM-4326 closes the finite-head body-safety obligations in
the earlier fixed-fifty petal hierarchy all at once. The lower bounds
`chi_50>=18,19,20,22` remain true historical stages, while `(5)` is the
current exact chart number. Their component, petal, cutoff, and repair
atlases remain independent structural results; only their chart-number
frontier is superseded.

## 4. Odd-tail transfer

For positive `c`, circle multiplication by `c` preserves normalized Haar
measure, so

```text
mu(G_(cA))=mu(G_A)>=alpha.                             (16)
```

Apply THM-4150 to the eleven-element body `cA`. It supplies `(6)` for every
distinct positive odd pair. All body speeds in `2cA` are even and pairwise
distinct, while the two tails are odd and distinct, so `(6)` has exactly
thirteen nonzero speeds.

## 5. Inheritance, mechanism, and scope

The closest proved mechanism is THM-4326's fixed-pool arbitrary-pair
closure. The canonical hostile remains its exact rank-two minimizer
`(q,r)=(50,70)` with body mask `031c7400`. The corrected near misses are
MISTAKE-524, which forbids forward use of a reservation, and MISTAKE-532,
which requires normalized rather than raw cross-grid ticks. Neither affects
the quantified safety statement used here. The least-used decisive sidecar
is the **complete layer boundary**: THM-4156, THM-4191, and THM-4326 cover
exactly zero, one, and two outsiders, respectively.

The connection contract is

```text
source:       three proved Haar-safe layers over P
target:       every eleven-face of P union {q,r}
map:          A |-> |A intersect {q,r}| in {0,1,2}
preserved:    labels, cardinality, outsider identity, and Haar threshold
destroyed:    safe-component addresses and theorem-specific repair witnesses
sidecar:      the outsider count and distinct-outside-P hypotheses
hostile:      (50,70), body 031c7400 at the rank-two layer
decisive test: every j-layer matches one theorem with the same quantifiers.
                                                               (17)
```

This theorem does not map an arbitrary speed row into the displayed pool,
does not normalize arbitrary parity, and proves no pool maximality among
other choices. It identifies no minimum full safe mass and supplies no
owner, arrival, or terminating descent. LRC(14) remains open. **QED.**
