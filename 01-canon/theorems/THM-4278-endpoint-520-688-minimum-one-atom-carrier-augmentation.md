---
id: THM-4278
title: "Endpoint 520,688 minimum one-atom carrier augmentation"
status: >
  PROVED RELATIVE TO THM-4266 + FINITE-EXACT + INDEPENDENTLY AUDITED;
  STRUCTURAL ONLY. At the former top edge (520,688), THM-4266's frozen
  8,319-mask carrier has 5,934 active masks and misses exactly two
  nine-bodies. Exactly 72 pair-active rank-eight masks are disjoint from
  both failures, and these are precisely all one-atom completions of the
  frozen active deck. The least is {8,40,42,63,80,84,120,143}; zero atoms
  fail and one succeeds, so the augmentation number is exactly one in this
  fixed proof language. THM-4271 already owns this edge and the current
  descent, so the current proof-graph contribution is zero. No band,
  physical entry, or LRC(14) follows.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4266-three-round-learned-carrier-endpoint-descent
related:
  - THM-4271-fourth-round-learned-carrier-endpoint-descent
  - THM-4273-endpoint-542-732-minimum-one-atom-carrier-augmentation
  - THM-4276-six-atom-endpoint-671-augmentation-and-one-layer-descent
primary_script: 04-computation/lrc14_endpoint_520_688_minimum_one_atom_primary_thm4278.cpp
primary_output: 05-knowledge/results/lrc14_endpoint_520_688_minimum_one_atom_primary_thm4278.out
primary_script_sha256: d5b64aa6157a0db7020f40ffd958de558356d7dcf1bb4906cbe225ee96306e2a
primary_output_sha256: 4ae51e53de00505389f80d2bfb62b06f4d61ada598146647a1057d8c8b71a647
literal_script: 04-computation/lrc14_endpoint_520_688_minimum_one_atom_literal_wall_audit_thm4278.cpp
literal_output: 05-knowledge/results/lrc14_endpoint_520_688_minimum_one_atom_literal_wall_audit_thm4278.out
literal_script_sha256: 2e9ee3361938f1ad36d0acdae481a93c8e4f3ac7e8544d7558bd896b43320c5d
literal_output_sha256: b51281cdf13ffa8e60195a6432fe7f91959ce892b7f7f5cb89e92c4bee84b09d
event_script: 04-computation/lrc14_endpoint_520_688_minimum_one_atom_event_sweep_referee_thm4278.cpp
event_output: 05-knowledge/results/lrc14_endpoint_520_688_minimum_one_atom_event_sweep_referee_thm4278.out
event_script_sha256: 2c48cb50a352a37093412427e817df198ef93bdd4e6cfc8558b51d183389a9e1
event_output_sha256: c1efea110ddd901d1e8aeeb156932b801a6be0b8584314cec4cc929ada8e7b5e
family_script: 04-computation/lrc14_one_atom_common_active_family_survey_thm4278.cpp
family_outputs:
  - 05-knowledge/results/lrc14_one_atom_common_active_family_520_688_thm4278.out
  - 05-knowledge/results/lrc14_one_atom_common_active_family_416_704_hostile_thm4278.out
family_script_sha256: d520628745aa3e1c878a89d41ef15632b12c828094201abaaf65da2251ac90b1
family_output_sha256:
  - 31a7f43fd194872ce71f2d94c9e2fd6425b12b4266b0ed6fd3643b78d6825b7f
  - 8f1655965b6be8397d7aff31826c30cdad0c206058768b061eebfcf0f6591003
postprocess_script: 04-computation/lrc14_endpoint_520_688_one_atom_overlap_postprocess_thm4278.py
postprocess_output: 05-knowledge/results/lrc14_endpoint_520_688_one_atom_overlap_postprocess_thm4278.out
postprocess_script_sha256: 88203ba6443408a7d5f8a5afadf78ca96537fb5547cfc8086453e06ecf5f0934
postprocess_output_sha256: 342a7fedb270eebd9fdcf7683732c483b31880b3bd166617ee779935c435faa8
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary signed endpoint-cocycle path and two detached
  literal-wall/event-sweep paths reconstruct the same 8,319-mask carrier,
  5,934-mask active subdeck, two failed bodies, complete 75,582-candidate
  universe, 72 successful repairs, canonical exact mass, and zero-failure
  augmented scan. The event implementation's O0, O3, and
  _GLIBCXX_ASSERTIONS transcripts are byte-identical after newline
  normalization and have empty stderr. The family survey supplies a sharp
  disconnected hostile to any universal Johnson-connectivity claim. The
  portable postprocessor verifies zero overlap with THM-4271's full 5,398-mask
  discovery prefix in normal, optimized, and fixed-hash-seed runs.
---

# THM-4278 -- endpoint `(520,688)` minimum one-atom carrier augmentation

**PROVED RELATIVE TO THM-4266 + FINITE-EXACT + INDEPENDENTLY AUDITED;
STRUCTURAL ONLY; CURRENT PROOF-GRAPH CONTRIBUTION ZERO AFTER THM-4271;
LRC(14) REMAINS OPEN.**

## 1. Frozen carrier and exact classification

Retain the labelled pool and target

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14},   alpha=4/63.       (1)
```

Let `C` be THM-4266's frozen ordered carrier, with

```text
|C|=8,319,                     FNV(C)=e08b227730f6793c.     (2)
```

At `(q,r)=(520,688)`, let its active subdeck be

```text
D={R in C:mu(G_((P\R) union {520,688}))>=4/63}.            (3)
```

Then `|D|=5,934`. It covers every member of `binom(P,9)` except exactly

```text
H1={16,85,88,95,145,168,193,240,252},   mask 0x07187008,
H2={16,85,88,168,176,193,240,252,290},  mask 0x27503008.   (4)
```

Their intersection has seven labels and their union has eleven:

```text
H1 intersect H2={16,85,88,168,193,240,252},
H1 union H2={16,85,88,95,145,168,176,193,240,252,290}.    (5)
```

Define the common-disjoint active family

```text
Q={R in binom(P,8):
     R intersect (H1 union H2)=empty and
     mu(G_((P\R) union {520,688}))>=4/63}.                (6)
```

The complement of `H1 union H2` has nineteen labels, so `(6)` has the
complete candidate universe

```text
binom(19,8)=75,582.                                        (7)
```

Exact enumeration gives

```text
|Q|=72,                         FNV(Q)=ed1bfbaf6eaa06a3,
activation equalities=0,        Q intersect C=empty.       (8)
```

Masks in `(8)` are in increasing unsigned-mask order and the FNV is FNV-1a-64
over little-endian unsigned `u64` words.

The family `Q` is precisely the set of all single pair-active rank-eight
repairs that complete `D`. In particular, the minimum number of such repairs
needed to complete this frozen active deck is exactly one.

## 2. Failure-hypergraph lemma and canonical minimum repair

The exact mechanism is a general finite set-cover reduction. Let `U` be any
finite body universe, let `D_0` be a frozen active deck, and put

```text
F={B in U:no R in D_0 is disjoint from B}.              (L1)
```

For each eligible new active repair `R`, define its response

```text
C(R)={B in F:R intersect B=empty}.                      (L2)
```

Then an augmentation `A` closes `U` if and only if

```text
union_(R in A) C(R)=F.                                  (L3)
```

Indeed, `D_0` already closes `U\F`, while `(L3)` is exactly the statement
that the new repairs close every remaining body. Thus the minimum augmentation
number relative to the frozen deck is the set-cover number of this response
hypergraph. In the one-atom case, `C(R)=F` if and only if `R` is disjoint from
`union F`; for rank `k` repairs drawn from a pool `P`, the complete unfiltered
candidate universe therefore has size

```text
binom(|P\union F|,k).                                   (L4)
```

after which activity remains a separate exact filter. This lemma is
deck-relative: it does not assert global deck minimality.

For `(3)--(6)`, `F={H1,H2}`, so the response quotient has a full-response
class exactly when `(6)` is nonempty. Equation `(8)` says that class has 72
realizers and hence set-cover number one. At the next canonical resistance,
THM-4276 applies the same quotient to 27 row-labelled obligations: 5,852,925
repairs collapse to 330 nonempty response patterns, and a six-obligation
packing plus a six-pattern cover proves minimum six. This is a proved
mechanism shared across the two stages, not an assumption that activity or a
chosen repair transfers between pairs.

The least repair in `(8)`, equivalently the zero-based colex-minimum, is

```text
R*={8,40,42,63,80,84,120,143},
mask=0x00048ec1,                  colex_rank_0=51,083.       (9)
```

Its exact mass and surplus are

```text
mu(G_((P\R*) union {520,688}))
 =1559831620541/24511557965895
 =4/63 + 3542225881/24511557965895.                       (10)
```

Appending `R*` gives `5,935` active masks and covers all
`binom(30,9)=14,307,150` bodies. The ordered augmented scan performs
`486,545,454` disjointness checks and finds zero failures. Hence, for every
`B in binom(P,9)`,

```text
mu(G_(B union {520,688}))>=4/63.                          (11)
```

To prove the classification rather than mere existence, observe that `D`
already covers every body other than `H1,H2`. A single active repair completes
the deck only if it is disjoint from each missed body, and any active repair
with that property completes both. Thus the successful one-atom
augmentations are exactly `(6)`. With no added atom, the two failures in `(4)`
remain; `(9)` succeeds. This proves the augmentation number one in the stated
fixed rank-eight carrier language.

For every covered body, the target implication is the monotone inclusion

```text
B subset P\R,
G_((P\R) union {520,688}) subset G_(B union {520,688}).    (12)
```

No converse to `(12)` is used.

## 3. Structure of the 72-repair family

Let `J(Q)` join two repairs when one is obtained from the other by one label
swap, equivalently when their symmetric difference has size two. The exact
family survey gives

```text
union support size                 16,
common core size                    0,
allowed labels unused by Q          {20,30,60},
J(Q) components                     1 (all 72 vertices),
J(Q) edges                          318,
J(Q) degree range                   [2,19].               (13)
```

The transversal number of the hypergraph `Q` is two. Its thirteen minimum
transversals are

```text
{8,40}, {8,42}, {8,80}, {8,120}, {8,190}, {40,42},
{42,63}, {42,80}, {42,120}, {42,132}, {42,190},
{42,264}, {80,264}.                                      (14)
```

These are structural facts about this exact active family, not a general
description of active repairs. In particular, simple arithmetic summaries do
not classify the 72 masks:

```text
gcd profile                         {1:70,2:2},
number of odd labels                {0:2,1:40,2:29,3:1},
sum of labels modulo 14             all fourteen residues occur. (15)
```

The postprocessor also finds

```text
Q intersect (THM-4271 full 5,398-mask discovery prefix)=empty. (16)
```

Hence `Q` is disjoint in particular from THM-4271's 199 novel masks. The
targeted one-atom route and the learned-prefix route close the same seed by
genuinely different certificates.

## 4. Sharp boundary: connectivity is not universal

The connected graph and transversal number in `(13)--(14)` do not transfer to
other hostiles. At THM-4266's earlier seed `(416,704)`, let

```text
H0={85,88,168,190,193,240,252,264,290}.                  (17)
```

Form the family of all pair-active rank-eight masks disjoint from `H0`. It has
`5,872` members, but its one-swap graph has two components of sizes

```text
5,871 and 1.                                               (18)
```

Its transversal number is four. The isolated active repair is

```text
I={8,15,16,30,60,120,143,286},       mask 0x1004812d.   (19)
```

On the literal common grid its threshold numerator is `932,285,998,476>0`,
whereas its best allowed one-swap neighbor is

```text
{8,16,30,60,120,143,145,286},        mask 0x100c8129,
threshold numerator=-1,799,610,100,092.                  (20)
```

Thus the first failed implication in a universal connectivity guess is exact:
activity need not survive even the best adjacent label swap. The strongest
survivor is family-specific—Johnson components and transversals are useful
coordinates, but they must be recomputed on each arithmetic fibre.

## 5. Three exact paths and reproduction

The primary reduces `520:688` to primitive pair `65:86` with scale eight and
rebuilds THM-4266's carrier from its frozen stage data. The first detached
audit constructs literal joint walls on denominator
`784,369,854,908,640`. The referee implementation independently event-sweeps
those walls, reconstructs all five carrier stages, re-enumerates the complete
candidate family, and scans every body before and after augmentation. They
agree on every load-bearing ledger in `(2)--(11)` and on the full ordered list
of 72 masks.

```bash
g++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_endpoint_520_688_minimum_one_atom_primary_thm4278.cpp \
  -o /tmp/thm4278-primary
/tmp/thm4278-primary \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254 \
  05-knowledge/results/lrc14_three_round_learned_carrier_thm4266

g++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_endpoint_520_688_minimum_one_atom_literal_wall_audit_thm4278.cpp \
  -o /tmp/thm4278-literal
/tmp/thm4278-literal \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band \
  05-knowledge/results/lrc14_three_round_learned_carrier_thm4266

g++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_endpoint_520_688_minimum_one_atom_event_sweep_referee_thm4278.cpp \
  -o /tmp/thm4278-event
/tmp/thm4278-event \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band \
  05-knowledge/results/lrc14_three_round_learned_carrier_thm4266

g++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_one_atom_common_active_family_survey_thm4278.cpp \
  -o /tmp/thm4278-family
/tmp/thm4278-family 520 688 0x27587008
/tmp/thm4278-family 416 704 0x2f903000

python3 -B \
  04-computation/lrc14_endpoint_520_688_one_atom_overlap_postprocess_thm4278.py
```

## 6. Scope and current accounting

This theorem classifies a minimum augmentation only relative to THM-4266's
frozen active deck at one pair and the declared rank-eight active-repair
language. It does not prove a globally minimum deck, a neighboring edge, a
uniform endpoint band, or a physical entry theorem.

Most importantly, THM-4271 already closes `(520,688)` and owns the complete
2,419-edge descent from that seed. At the post-THM-4271 stage the residual was
`174,904`, with `(256,671),(384,671)` alone on its top layer. THM-4276 has
since closed those two rows and 161 endpoint-670 rows, leaving the current
residual `174,741`, topped by `(256,670),(384,670)`.

THM-4278 causes none of either theorem's deletions. Its current proof-graph
contribution is therefore **zero**, and it changes no residual ledger. It
supplies a compact alternate seed certificate, a complete one-atom
classification, the reusable response-hypergraph mechanism, and a sharp
structural hostile only. **LRC(14) remains open. QED.**
