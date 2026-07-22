---
id: THM-2076
title: "Guard-capacity terminal rank floor for dyadic safe-child towers"
status: >
  PROVED by Haar measure and compact/open separation. If an s-speed closed
  delta-safe set is contained in one strict 2delta-eligibility set, then
  s>(1-4delta)/(2delta). At delta=1/14 this forces s>=6. Therefore every
  nontrivial THM-2073 dyadic tower has terminal core size at least six and
  depth at most five, sharpening the elementary depth-eight cardinality
  bound. Together with THM-2075, its top safe set has measure strictly below
  2^(1-r)/7. This is a structural reduction, not LRC(14).
source: codex-2026-07-21-LRC-guard-capacity-rank-floor
depends_on:
  - THM-2073
  - THM-2075
related:
  - THM-2072
  - THM-775
  - THM-2061
  - THM-2080
---

# THM-2076 -- guard-capacity terminal rank floor

## 1. General guard-capacity lemma

Fix

```text
0<delta<1/4.                                             (1)
```

Let `Q` be a finite set of `s` nonzero integer speeds and let `h` be a
nonzero integer. Put

```text
G_Q(delta)={t in R/Z:||qt||>=delta for every q in Q},
E_h(2delta)={t in R/Z:||ht||<2delta}.                    (2)
```

If `G_Q(delta)` is nonempty and

```text
G_Q(delta) subset E_h(2delta),                           (3)
```

then

```text
s>(1-4delta)/(2delta).                                   (4)
```

### Proof

For each nonzero integer `q`, multiplication by `q` preserves Haar
probability measure on the circle. Hence its strict danger set has measure

```text
measure{t:||qt||<delta}=2delta.                          (5)
```

The union bound gives

```text
measure(G_Q(delta))
 =1-measure(union_(q in Q){||qt||<delta})
 >=1-2s delta.                                          (6)
```

Similarly,

```text
measure(E_h(2delta))=4delta.                             (7)
```

The inequality is strict after imposing (3). Indeed `G_Q(delta)` is compact
and `E_h(2delta)` is a proper open union of intervals. A nonempty compact
subset of it misses a positive-width collar of its boundary, so

```text
measure(G_Q(delta))<4delta.                              (8)
```

Combining (6) and (8) gives

```text
1-2s delta<4delta,
```

which is exactly (4). QED.

The strict/open convention is load-bearing. Replacing `E_h` by a closed set
would give only a weak inequality and would lose the integer jump in the
application below.

## 2. LRC(14) specialization

At `delta=1/14`, formula (4) becomes

```text
s>5,    hence s>=6.                                     (9)
```

Equivalently:

> No closed `1/14`-safe set of five or fewer speeds can be contained in the
> strict `1/7`-eligibility teeth of one nonzero guard.

The endpoint `s=5` illustrates why compactness matters. The union bound gives
`measure(G_Q)>=2/7`, exactly the total measure of the open guard teeth. A
closed nonempty safe set cannot fill that open set up to equality because it
must stay a positive distance from the tooth boundaries.

## 3. Sharpened depth of the THM-2073 tower

Let

```text
C=Q_0,
Q_i=2Q_(i+1) union {h_i},    0<=i<r,                    (10)
```

be a nontrivial THM-2073 tower. Its terminal guard property says

```text
G_(Q_r)(1/14) subset E_(h_(r-1))(1/7).                  (11)
```

The safe set is nonempty by settled lower-dimensional LRC; more directly,
THM-2073's safe-child construction starts from it. Apply (9) to `Q_r`:

```text
|Q_r|>=6.                                               (12)
```

But every descent removes one speed, so

```text
|Q_r|=11-r.                                             (13)
```

Therefore

```text
r<=5.                                                   (14)
```

This sharpens THM-2073's purely hereditary/cardinality bound `r<=8`. Thus the
only possible nontrivial terminal sizes and depths are

```text
terminal size: 10,9,8,7,6,
tower depth :  1,2,3,4,5.                               (15)
```

THM-2080 subsequently resolves the equality-pressure rank-six lane by an
exact unequal-comb overlap floor. It improves (12) to `|Q_r|>=7`, hence
`r<=4`. The present Haar lemma remains the general first capacity gate.

Combining (7)--(8) with THM-2075's exact measure scaling gives the additional
top-level tax

```text
measure(G_C)
 =2^(-r) measure(G_(Q_r))
 <2^(-r)*(2/7)=2^(1-r)/7.                               (16)
```

This is strongest on deep towers; it records an exact cost of each safe-child
level even though it does not by itself contradict the folded-diamond cap in
THM-2061.

## 4. Scope and next target

The theorem removes terminal sizes three, four, and five, and hence the three
deepest valuation patterns allowed by THM-2073. It does not rule out terminal
sizes six through ten. The union-bound proof cannot do so uniformly: for six
speeds it supplies only `measure(G_Q)>=1/7`, leaving room inside a guard set of
measure `2/7`.

A further rank floor must use information discarded by scalar measure, such
as pairwise danger-set overlap, terminal endpoint owners from THM-2075,
hereditary primitivity, or the divisor-completeness constraints through `14`.
This identifies a concrete next inequality:

```text
prove measure(G_Q)>=2/7 for every hereditarily primitive,
divisor-complete six- through ten-core Q,
or prove that any smaller safe set contains an antipodal pair.              (17)
```

Either conclusion would conflict with the strict guard containment.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that depth should be bounded by counting speeds
or 2-adic valuations. The sharper state variable is Haar capacity: a single
guard has only `2/7` of strict eligibility space, while `s` core dangers can
remove at most `s/7` in the first union bound. The quotient preserves exactly
the containment predicate (3), but discards component locations and owner
labels, which are needed beyond the rank-six floor.

There is no meaningful tournament here. Treating danger sets as vertices and
orienting by measure or containment gives a transitive scheduling order but
cannot encode a many-set union bound. The exact carrier is the hypergraph of
danger sets with Haar weights and the strict guard set as its target. QED.
