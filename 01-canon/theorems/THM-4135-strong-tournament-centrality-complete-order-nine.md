---
id: THM-4135
title: "Complete strong-tournament Johnson centrality at order nine"
status: >
  FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. All 178,133 strong
  tournament isomorphism classes of order nine have only central maximizers
  t=+-1 of the THM-4128 rational support floor and the THM-4123 exact-coset
  floor. The minimum central-versus-outer coset margin is 90. Actual response
  maxima are noncentral-only in 3,248 classes. THM-4133 refutes the all-order
  extension at order twelve; orders ten and eleven remain unclassified.
source: codex-frontier-synthesis-creative-20260825s
depends_on:
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4131-strong-tournament-centrality-through-order-eight
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
script: 04-computation/tournament_strong_centrality_order_nine_thm4135.py
output: 05-knowledge/results/tournament_strong_centrality_order_nine_thm4135.out
independent_audit_script: 04-computation/tournament_strong_centrality_order_nine_thm4135_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_strong_centrality_order_nine_thm4135_independent_audit.out
script_sha256: ab247ff7da23d2659efd53dfc2f0e7a8841134734004eafa04c13ad596b28986
output_sha256: 24412327b2d56362ddcdad2a042192ff04df56d0001e8456aba53698c037d46c
independent_audit_script_sha256: 5ed81bd2ab3054d1a05471c42fe724c14eb0df278dbc3ee3b7d6d38253fff530
independent_audit_output_sha256: b217807222b0b8cd2b2234f37a9bf40129469671ff64d6ea28dd0fcd88657b04
semantic_sha256: 0b3d2c65723e0ecc78cbf02d1735320794a6bd6e5f7c3a371ee94432e19c5b49
primary_profile_sha256: f000c999cab374ffcdc03cffaa08675f5e97152340a511c3d8de745230141746
independent_profile_sha256: c01572932582bf16071dfdedef810ff845197c4b2154518a31e251c04dc59655
hash_basis: raw LF bytes
primary_audit: >
  PASS. Homebrew nauty gentourng -q -c 9 supplies one representative of every
  strong order-nine class. The pinned and previously independently audited
  THM-4131 Start/End/exposed-gap evaluator recomputes every response, exposure
  packet, rational floor, layer lattice, exact-coset floor and actual maximum.
  Both raw universes and every ordered exact profile are fingerprinted.
independent_audit: >
  ACCEPT. A clean-room C++ contracted-good/reversed-block implementation
  imports no primary code, evaluates the full stream in one and eight OpenMP
  threads with identical output, hashes all 512 responses and every exact
  layer per class, and reproduces all counts, margins, histograms, worst
  packets and three named hostile controls by literal child DP.
---

# THM-4135 -- complete strong centrality at order nine

**FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4131 proves rational and exact-coset Johnson centrality through order
eight. The complete next census gives the last currently classified positive
order and supplies the structured extremum from which THM-4133's eventual
counterexample was found.

## 1. Statement

For a strong tournament `T` of order nine, retain

```text
F_T(S)=H(T+x_S),
J_m=E_m F_T+Var_m(F_T)/(E_m F_T-H),
L_m=ceil_(a_(T,m)+d_(T,m)Z)(J_m),
t=9-2m.                                                    (1)
```

The central layers are `t=+-1`.

> **Finite order-nine theorem.** For every strong tournament of order nine:
>
> 1. every layer maximizing `J_m` is central;
> 2. every layer maximizing `L_m` is central; and
> 3. the best central `L_m` is strictly larger than the best noncentral
>    `L_m`.

The exact universe has

```text
all tournament classes:       191536,
strong tournament classes:    178133.                     (2)
```

There are zero rational failures and zero exact-coset failures. The smallest
central-minus-outer coset margin is

```text
90.                                                            (3)
```

## 2. Exact order row

The rational optimizer histogram is

```text
t=(-1):       87736,
t=(-1,1):      2661,
t=(1):        87736.                                      (4)
```

Exact-coset rounding changes the optimizer tuple in `10,772` classes but
never leaves the central pair:

```text
t=(-1):       82350,
t=(-1,1):     13433,
t=(1):        82350.                                      (5)
```

The worst normalized tilt occurs in one reversal pair. Its packets are

```text
(H,W,D_4,C_hd)=(405,4665,6029505,+-1478700),
theta=+-32860/44663,
rho=16430/44663.                                          (6)
```

Thus the order-nine extremum uses less than thirty-seven percent of the
rational centrality budget. The negative-tilt member is the substitution of
a regular `R_5` block into the quotient later used by THM-4133; deleting one
vertex from larger cyclic blocks exposes the order-ten near-boundary and the
order-twelve failure.

## 3. Actual maxima remain a different coordinate

Actual layer maxima have histogram

```text
t=(-3):          1624,        t=(-3,-1):       83,
t=(-1):         85051,        t=(-1,1):      4617,
t=(1):          85051,        t=(1,3):         83,
t=(3):           1624.                                      (7)
```

Therefore `3,248` strong classes have no central actual maximizer. For the
negative member of `(6)`, the support floors optimize at `t=-1` while the
actual maximum occurs at `t=-3`. This preserves the THM-4128/4131 firewall:
a certified support floor does not locate the actual maximizing ear.

## 4. Universe and redundant audits

The raw nauty streams have SHA-256 fingerprints

```text
all:    4f7d6c43cfed87e1e5293dc751736efe2d7bc1554946cdc83f4026a575fbbbf8
strong: 3073d5ecf5f34345aa5f35c349c51b35f4c244e687db25c16627fcb602b019a1. (8)
```

The primary audit inherits THM-4131's Start/End/exposed-gap DP and hashes its
complete ordered exact profiles as

```text
f000c999cab374ffcdc03cffaa08675f5e97152340a511c3d8de745230141746. (9)
```

The independent C++ audit contracts each good arc block together with its
reversed one-defect block, reconstructs the full cut response independently,
and includes all `512` mask values in its profile serialization. Its digest is

```text
c01572932582bf16071dfdedef810ff845197c4b2154518a31e251c04dc59655. (10)
```

The distinct serialization explains why `(9)` and `(10)` differ; every named
mathematical invariant agrees. Nonstrong codes `2,140` and strong code `20`
reproduce the THM-4128 rational, coset-reordering and support/actual hostiles.

## 5. Scope and replay

This theorem is complete only at order nine. THM-4133 now refutes all-order
strong rational and exact-coset centrality at order twelve. Orders ten and
eleven remain unclassified, and no actual-maximizer or interval-completeness
claim is made.

Run the primary audit with

```text
python3 -B 04-computation/tournament_strong_centrality_order_nine_thm4135.py
```

and the independent audit with

```text
clang++ -O3 -std=c++17 -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/openssl/include \
  04-computation/tournament_strong_centrality_order_nine_thm4135_independent_audit.cpp \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  -L/opt/homebrew/opt/openssl/lib -lcrypto \
  -o /tmp/thm4135_independent
gentourng -q -c 9 | env OMP_NUM_THREADS=8 /tmp/thm4135_independent
gentourng -q -c 9 | env OMP_NUM_THREADS=1 /tmp/thm4135_independent
```

The one- and eight-thread independent outputs are byte-identical and match
the frozen output. **QED.**
