---
id: THM-2898
title: "Unique maximum-gate five-parity matching closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.  The
  unique THM-2896 root with least hitting-gate size 25 is closed through
  five singleton-complement parity branches, 277 literal H4-pair
  residuals, one exact disjoint-matching repair, and scalar cascades.
source: root-2026-07-29
depends_on:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
related:
  - THM-2900-flag-conditioned-rank-selective-partition-closure
verification:
  - 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898.py
  - 05-knowledge/results/lrc14_j6_k25_five_parity_matching_closure_thm2898.out
  - 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898_independent_audit.py
  - 05-knowledge/results/lrc14_j6_k25_five_parity_matching_closure_thm2898_independent_audit.out
---

# THM-2898 -- unique maximum-gate five-parity matching closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.**

## 1. Statement

Let

```text
E=(1,8,10,11,12,13,14)
```

and let `G_E` be its gap-`1/14` good set.  No six distinct external speeds
`w>=15` have danger combs whose union covers `G_E`.

This is the unique root on which the least THM-2896 hitting gate has size
`25`.  It is distinct from the four roots closed by THM-2895, so those two
theorems close five pairwise distinct seven-body/six-slot roots.  The result
is not a uniform seven-body closure and does not prove LRC(14).

## 2. The ordered gate schedule

THM-2896 gives the exact hitting gate

```text
(23,27,19,46,18,17,25,34,38,100,63,156,29,
 125,130,44,37,50,92,168,72,32,110,54,182).
```

Every hypothetical six-cover meets this set.  Assign it to the first gate
label in the following locked proof order:

```text
PARITY  23
PARITY  27
PARITY  19
PARITY  46
PARITY  18
SCALAR  168,182
MATCH   17
SCALAR  25,156,125,44
SCALAR  130,72,54
SCALAR  34,100,92,32,110
SCALAR  38,63,29,37,50.
```

At an apex `a` with prior prefix `P`, the other five cover labels lie
outside `P union {a}` by THM-2893's ordered-suffix refinement.  Every
displayed scalar activation is certified by the strict inequality

```text
q_1(P)+...+q_5(P)<|G_E minus D_a|.                         (1)
```

The exact scalar-closed-state search has minimum seed number `6`, uses
`68,407` closed states, `323,776` edges and `68,407` closure calls, and
returns the seed bank

```text
(23,27,19,46,18,17).
```

This minimum is a scalar-bootstrap workload statistic.  The matching
certificate below replaces the sixth parity proof; it is not a claim that
five is a globally minimum number of nonscalar arguments.

## 3. Five parity branches

For apices

```text
(23,27,19,46,18)
```

on their literal current prefixes, the strict THM-2895 hypothesis

```text
q_1<3h/7
```

holds.  Their globally sealed `H_4` core sizes and pair counts are

```text
apex        23   27   19   46   18
|H_4|       13   13   10   11    7
pairs       78   78   45   55   21.
```

All `277` literal pair residuals are nonempty and admit no three-cover.
The fresh child top-three singleton sum closes `276`; the union of that
test with the exact child `B_2+B_1` partition cap closes all `277`.
Exactly `2,549` finite child-pair evaluations are paid, and no recursive
`H_2` row survives.

The literal carrier, the actual excluded prefix, and the pair are retained
in every ledger row.  A later gate label is never used as an earlier
exclusion.

## 4. The apex-17 matching repair

After the fifth parity branch, scalar closure adds `168,182`.  The live
prefix for apex `17` is therefore

```text
P=(23,27,19,46,18,168,182).
```

The primary verifier proves the stronger certificate in the larger allowed
universe obtained by excluding only

```text
P_0 union {17}=(23,27,19,46,18,17).
```

Monotonicity under further deletion then transfers it to the actual
prefix.  On the carrier of mass

```text
h=189841/1021020
```

the exact reference-prefix statistics are

```text
q_1 = 1023641/23823800      at 25,
q_5 = 81891/2290750         at 125,
B_2 = 33906683/440740300    at {25,37},
T   = h-q_5                 = 26834677/178678500,
L   = T-B_2                 = 121070701/1652776125.
```

The THM-2897 finite matching gap is

```text
delta=L-q_1-h/7=1492091/400673000>0,
```

so every edge of weight at least `L` has both endpoints at most `1844`.
The exact intrinsic undirected `L`-heavy threshold graph has three edges:

```text
{25,37}: 33906683/440740300,
{25,45}: 3985831/53603550,
{37,50}: 65302143/881480600.
```

Ties are retained by the construction; none occur in this three-edge
ledger.  Its only vertex-disjoint edge pair is

```text
{25,45}, {37,50},
```

whose total weight is `47104891/317333016`, below `T` by

```text
69186919/39666627000>0.                                  (2)
```

The other two edge pairs exceed `T` but share a vertex and therefore are
not matchings.

This threshold graph is sufficient: if two disjoint edges had total weight
at least `T`, then each would have weight at least `L`, because the other
edge is bounded by `B_2`.  They would consequently occur in the enumerated
graph, contradicting `(2)`.  Hence

```text
q_5+M_(2,2)<h,
```

and THM-2897 excludes the apex-17 five-cover.  The number in `(2)` is the
gap of the unique disjoint pair in the threshold graph; no claim is made
that it is the exact global margin `h-(q_5+M_(2,2))`.

Scalar closure then activates the remaining seventeen gate labels in the
four displayed rounds.  Every possible earliest gate branch is excluded,
contradicting the THM-2896 hitting conclusion.  This proves the statement.

## 5. Exact verification

The primary verifier pins THM-2895, THM-2896, THM-2897 and all imported
exact artifacts.  Its aggregate is

```text
(parity branches, rank-pair rounds, rank-pair activations, matching,
 H4 pairs, child top3, child B2+B1 union, closed, hard, paid child pairs)
=
(5,0,0,1,277,276,277,277,0,2549).
```

The standalone independent audit imports only the separately implemented
THM-2895 rational-interval engine.  It reconstructs the root and gate,
all apex profiles, the scalar-closed BFS, the five parity branches, and the
matching graph directly at the actual prefix.  It independently finds
`40` singleton-pruned matching candidates, the same three heavy edges,
the same unique disjoint pair and the same strict gap `(2)`.

Reproduction:

```bash
python3 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898.py
python3 -O 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898.py
python3 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898_independent_audit.py
python3 -O 04-computation/lrc14_j6_k25_five_parity_matching_closure_thm2898_independent_audit.py
```

Both ordinary/optimized pairs are byte-identical.  Neither verifier uses a
Python `assert`.  SHA-256:

```text
primary source
5625b512ea8de6e14ee046b5c38909e1b4291f9b7a50ac3b17b2711b7554b104

primary output
41f5e443f6d1ee2553c332da7709bd0c89f400b9ca154ddb6047f8ca724c6a40

independent source
289a63e4392613bdb82ede770411be51d074ea38e9e4b0a8129d687d238e6b69

independent output
382e5ae119f67610b098ff5a393c1b52fa37832d95a0cbd12fa50ec78cdf8007

primary canonical ledger
79c37deea9d831456e5ed913adadef7c465ff454a690fed5ff4aad7f0ed911a2

independent canonical ledger
4eb7bb8f0d4fdb02cb6c2dc26eacce23b2e41749aefdbc05eac5506b89efc970
```

The matching object is an undirected weighted graph with a disjointness
invariant.  Orienting it as a tournament would add a gauge, destroy ties,
and lose the theorem-bearing matching constraint. ∎
