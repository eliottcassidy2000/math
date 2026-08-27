---
id: THM-4273
title: "Endpoint 542,732 minimum one-atom carrier augmentation"
status: >
  PROVED RELATIVE TO THM-4254/4261 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  At the historical hostile pair (542,732), THM-4261's 3,227-mask active
  subdeck misses exactly two nine-bodies, whose intersection has eight labels.
  Exactly 2,172 pair-active rank-eight repairs close both failures, and these
  are precisely all one-atom augmentations that close every labelled body.
  The least is {8,16,42,63,84,88,120,126}; zero atoms fail and one succeeds,
  so the augmentation number is exactly one relative to the frozen carrier.
  THM-4266 already owns this edge; no current proof-graph deletion, band,
  physical entry, or LRC(14) follows.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
  - THM-4261-semantic-endpoint-band-prefix-union-lift
related:
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4270-uniform-four-primitive-outsider-rays-common-deck-closure
primary_script: 04-computation/lrc14_endpoint_542_732_minimum_one_atom_primary_thm4273.cpp
primary_output: 05-knowledge/results/lrc14_endpoint_542_732_minimum_one_atom_primary_thm4273.out
primary_script_sha256: d976add2ce9dd6d07cf4163d9fa9f968aaba9e18a21b6f820f21de6fc9045966
primary_output_sha256: 3dae51b727fabcca0e2f7d05318113227c93a1b085db2b446da7b5127248fb4b
independent_script: 04-computation/lrc14_endpoint_542_732_minimum_one_atom_literal_wall_audit_thm4273.cpp
independent_output: 05-knowledge/results/lrc14_endpoint_542_732_minimum_one_atom_literal_wall_audit_thm4273.out
independent_script_sha256: 0ac6d0b3297b2627842a1ec1b0de43e896c53faead1e0e6ccb3d470d37b19b65
independent_output_sha256: fc3fd82468b06c1cbc10e1fc242a57fe247e03ae9d30d130879b595e8423b4ba
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary recomputes signed endpoint-cocycle atoms at
  primitive pair 271:366 with scale two. The clean-room route instead builds
  literal joint walls on a fresh common denominator. Both reconstruct the
  inherited carrier, the same two failures, all 125,970 common-complement
  candidates, the 2,172-mask active ledger, the exact canonical mass, and the
  zero-failure augmented scan. O2/O3 transcripts agree after newline
  normalization.
---

# THM-4273 -- endpoint `(542,732)` minimum one-atom carrier augmentation

**PROVED RELATIVE TO THM-4254/4261 + FINITE-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Frozen carrier and exact statement

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Let `C` be THM-4254's ordered first-occurrence union of 4,675 rank-eight
repair masks. At the fixed pair `(q,r)=(542,732)`, put

```text
D={R in C:mu(G_((P\R) union {542,732}))>=4/63}.        (1)
```

Then `|D|=3,227`. Exactly two bodies in `binom(P,9)` are not disjoint from a
member of `D`:

```text
A={80,85,95,132,145,168,193,252,286},   mask 0x151a5400,
B={80,95,132,145,168,170,193,252,286},  mask 0x153a4400. (2)
```

Their intersection has eight labels and their union has ten. Define

```text
Q={R in binom(P,8):
     R intersect (A union B)=empty and
     mu(G_((P\R) union {542,732}))>=4/63}.              (3)
```

The exact classification is

```text
|Q|=2,172,                 FNV(Q)=829ae906d6b54c9a,     (4)
```

where masks are ordered increasingly and serialized as little-endian `u64`
words. There are no activation equalities and `Q intersect C=empty`.

Moreover, `Q` is exactly the set of pair-active rank-eight repairs `R` for
which `D union {R}` covers every member of `binom(P,9)`. Consequently the
minimum number of pair-active rank-eight repairs needed to augment `D` to a
cover is exactly one.

## 2. Canonical least atom and exact mass

The least member of `Q` in unsigned-mask order, equivalently in zero-based
colex order, is

```text
R*={8,16,42,63,84,88,120,126},
mask=0x1aa89,              colex_rank_0=21,647.         (5)
```

Its exact mass is

```text
mu(G_((P\R*) union {542,732}))
 =2395416707526053/37693075789228860
 >4/63.                                                   (6)
```

Appending `R*` to `D` closes all `binom(30,9)=14,307,150` bodies. Thus, for
every `K in binom(P,9)`,

```text
mu(G_(K union {542,732}))>=4/63.                         (7)
```

This is an alternate fixed-edge proof and a minimum-augmentation
classification. THM-4266 already closes `(542,732)` using its broader
round-two learned carrier, so `(7)` removes no edge from the current proof
residual.

## 3. Why the classification is exact

The active inherited deck `D` already covers every body except `A` and `B`.
Therefore a single pair-active repair closes the deck if and only if it is
disjoint from both exceptional bodies. Since `|A union B|=10`, the complete
candidate universe for that one atom has size

```text
binom(30-10,8)=binom(20,8)=125,970.                     (8)
```

The primary enumerates every mask in `(8)` and evaluates its exact mass by
signed endpoint-cocycle atoms. It finds precisely `(4)`. This proves the
“only if” direction as well as existence; `R*` is one of 2,172 successful
atoms, not a unique witness.

For any body `K` disjoint from an active repair `R`,

```text
K subset P\R,
G_((P\R) union {542,732}) subset G_(K union {542,732}). (9)
```

Thus the repair mass lower-bounds the target in the required direction.
The inherited deck supplies `(9)` for all nonexceptional bodies and `R*`
supplies it for both bodies in `(2)`. Zero added atoms leave those two bodies
uncovered; one succeeds. This proves the minimum exactly within the declared
rank-eight active-repair language.

## 4. Two exact implementations

The primary reduces `(542,732)` to primitive pair `271:366` with scale `2`,
constructs the signed cocycle atoms on primitive grid `1,388,604`, and
rebuilds all 59 THM-4254 prefix transcripts. Its frozen work ledger is

```text
carrier/active/failures       4675 / 3227 / 2,
inherited body checks         474,614,827,
common candidates/active      125,970 / 2,172,
augmented failures/checks     0 / 474,614,829.          (10)
```

The independent checker never calls that primitive-prefix evaluator. It
constructs literal safe intervals for `542` and `732` on grid

```text
301,544,606,313,830,880,
```

forms `9,653` joint cells and `2,545` nonzero failed-mask atoms, then
recursively enumerates all rank-nine bodies. It reproduces `(2)`, `(4)`,
`(6)`, and `(10)`. Its raw threshold numerator for `R*` is
`1,111,595,337,807,192`; the primary numerator is exactly 168 times this,
`186,748,016,751,608,256`, as required by their denominator conversion.

## 5. Reproduction and scope

```bash
g++ -std=c++20 -O3 -DNDEBUG -pthread \
  04-computation/lrc14_endpoint_542_732_minimum_one_atom_primary_thm4273.cpp \
  -o /tmp/thm4273-primary
/tmp/thm4273-primary \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254

g++ -std=c++20 -O3 -DNDEBUG -pthread \
  04-computation/lrc14_endpoint_542_732_minimum_one_atom_literal_wall_audit_thm4273.cpp \
  -o /tmp/thm4273-literal
/tmp/thm4273-literal \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band
```

Cardinality one is minimal only relative to augmenting the frozen active deck
`D` by pair-active rank-eight repairs. This theorem gives no globally minimum
proof language, endpoint band, neighboring edge, learned-carrier descent,
physical entry, or LRC(14). **QED.**
