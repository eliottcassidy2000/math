---
id: THM-4267
title: "Uniform four-to-five outsider-ray common-deck closure"
status: >
  PROVED RELATIVE TO THM-4231/4266 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  For every integer g>=73 and every nine-body in the fixed thirty-label pool,
  adjoining (4g,5g) leaves Haar-safe mass at least 4/63. A 9,806-repair deck
  is active at all 81 bridge scales 73 through 153 and closes every labelled
  body; THM-4231 supplies every g>=154. The exact post-THM-4266 residual
  contribution is 63 edges, leaving 177,522. No other ratio, below-boundary
  scale, physical entry, or LRC(14) follows.
source: root/cross-frontier-bridge/2026-08-27
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4266-three-round-learned-carrier-endpoint-descent
related:
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
  - THM-4262-uniform-three-four-outsider-ray-common-deck-closure
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
primary_script: 04-computation/lrc14_ratio45_common_deck_bridge_thm4267.cpp
primary_output: 05-knowledge/results/lrc14_ratio45_common_deck_bridge_thm4267.out
primary_script_sha256: 10e8f07589dfdcf42b4fa9f43769d3addd830d5e85d2b76a09e704e2aa09797e
primary_output_sha256: 1a49e73f68d029d72cb5face9dff01a434eb740eb944c5653b38dde2ea54e98c
independent_script: 04-computation/lrc14_ratio45_common_deck_literal_independent_audit_thm4267.cpp
independent_output: 05-knowledge/results/lrc14_ratio45_common_deck_literal_independent_audit_thm4267.out
independent_script_sha256: 8b6f5c19b0afc373afd15e0b55ac5c4697831e8f977562c49bfec7111d7e9563
independent_output_sha256: 411fbdefb173f0e75288f63cfec1430ca5a45e87bc593474b670ff3aa2c563de
proof_graph_script: 04-computation/lrc14_ratio45_common_deck_proof_graph_independent_audit_thm4267.py
proof_graph_output: 05-knowledge/results/lrc14_ratio45_common_deck_proof_graph_independent_audit_thm4267.out
proof_graph_script_sha256: f13b553470af4425ec802511c8800f82ff2af4176dc04ad70819be7cc864db8d
proof_graph_output_sha256: 211740c11a47d2a5f4da757a266b31523794a0570b362cb686810fcf08dd5001
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. A detached checker constructs literal joint walls at every
  scale and audits 2,654,370 repair-scale cells, all 14,307,150 bodies, both
  hostile 32K bodies and the targeted repair. O0/O3/ASan+UBSan transcripts
  agree. A separate semantic proof-graph reconstruction recovers the exact
  three-edge overlap with THM-4266, 63 novel edges and final residual.
---

# THM-4267 -- uniform four-to-five outsider-ray common-deck closure

**PROVED RELATIVE TO THM-4231/4266 + FINITE-EXACT + INDEPENDENTLY AUDITED;
LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Then

```text
for every integer g>=73 and every B in binom(P,9),
mu(G_(B union {4g,5g}))>=4/63.                         (1)
```

This is exactly the strict-above-pool range: `4g>max(P)=290` iff `g>=73`.
The exact finite bridge below covers `73<=g<=153`. THM-4231 proves every
distinct outsider pair with maximum at least `770`, so it supplies `(1)` for
`g>=154`, because `5g>=770`. No theorem failure or optimality is claimed
below `73`.

## 2. Primitive quotient with scale retained

Let `A=G_{4,5}`. On denominator

```text
N=14*4*5=280,
```

its positive components are exactly

```text
[5,52], [60,65], [75,108], [116,135],
[145,164], [172,205], [215,220], [228,275].            (2)
```

Their total length is `208`, so `mu(A)=26/35`. Under circle multiplication
`m_g(x)=gx`,

```text
G_{4g,5g}=m_g^(-1)(A).                                (3)
```

Let

```text
D=lcm(14s:s in P)=18,241,159,416,480.
```

For `R in binom(P,8)`, let `U_R=G_(P\R)`, with nonwrapping components
`[a_i/D,b_i/D]`, and define

```text
I(z)=280D integral_0^(z/D) 1_A(t)dt.                  (4)
```

Exact lifting gives

```text
mu(U_R intersect G_{4g,5g})
 =sum_i(I(gb_i)-I(ga_i))/(280Dg).                     (5)
```

Thus `R` is active at scale `g` exactly when

```text
63 sum_i(I(gb_i)-I(ga_i))-4*280*D*g>=0.               (6)
```

Quotienting `(4g,5g)` by its gcd destroys `g`, endpoint phase and repair
activation. Equation `(6)` is the required scale sidecar; primitive-ratio
truth is not copied blindly along the fibre.

## 3. Exact targeted common-deck bridge

Order all `binom(30,8)=5,852,925` repair masks by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).        (7)
```

Among the first `32,768` candidates, exactly `9,805` repairs satisfy `(6)`
at every one of the `81` scales `73..153`. Their ordered FNV is
`5ff9df10a934b071`. This common deck misses exactly two bodies:

```text
{16,88,95,145,168,176,193,240,252},
{16,88,95,145,168,193,240,252,290}.                   (8)
```

The later common-active repair

```text
R_*={63,80,85,126,143,170,264,286},
mask=0x18251600, candidate rank=35,908                  (9)
```

is disjoint from both. Appending it gives the theorem deck

```text
count=9,806,                  FNV=c2e18c7ecbcb9c72.    (10)
```

Exhaustive enumeration of all `binom(30,9)=14,307,150` labelled bodies has
zero failures, `540,608,482` ordered repair/body checks and maximum prefix
`9,806`.

For a disjoint pair `B,R`,

```text
B subset P\R,
G_((P\R) union {4g,5g}) subset G_(B union {4g,5g}).   (11)
```

The active repair mass therefore lower-bounds the target in the correct
direction. Equations `(6)`, `(10)` and `(11)` prove `(1)` on the bridge;
THM-4231 supplies the tail.

As a nested hostile control, the common counts at candidate budgets
`4,096, 8,192, 16,384, 32,768` are

```text
1,245, 2,417, 4,942, 9,805,
```

and the corresponding body-failure counts are `1,377, 219, 30, 2`. At
budget `65,536`, all `19,763` common repairs close every body independently;
the smaller targeted deck `(10)` is the promoted certificate.

## 4. Independent literal-wall audit

The referee does not evaluate the primitive prefix `(4)`. For each scale it
constructs the literal safe intervals of speeds `4g` and `5g`, sweeps their
joint walls against the fixed-pool geometry, and evaluates the first 32,768
candidates plus the targeted and primary-minimum controls. Its exact matrix
has

```text
cells=2,654,370, nonnegative=2,351,848, equalities=0,
FNV=94b4f38a25821fe3.                                  (12)
```

It reproduces the `9,805+1` deck, both bodies in `(8)`, and the zero-failure
full scan. The smallest checked literal margin occurs at `g=104` for

```text
R={10,42,85,95,120,170,240,264},
mu=68128192903/1073009377440,
mu-4/63=204461/357669792480>0.                         (13)
```

Multiplying the detached integer margin `1,313,866,386` by the exact
normalization `14,560` gives the primary margin `19,129,894,580,160`.
Ordinary, optimized and ASan+UBSan literal transcripts are byte-identical.
The sanitizer emits no diagnostics. A full 65,536-candidate primary replay
also byte-matches its frozen transcript.

## 5. Exact proof-graph contribution

Semantic selection from the post-THM-4256 residual finds `66` strict `4:5`
ray edges at

```text
g=73..131,133,134,135,137,138,139,141.                (14)
```

Their ledgers are

```text
FNV=7b63466b6d0274a1,
SHA256=b91b9eb54db8bd77649a1fcaf3a653e0f108008c01ba196b90878381e8eb6c31.
```

THM-4266 already removes exactly the three scales `138,139,141`, namely
`(552,690),(556,695),(564,705)`. The remaining `63` are the theorem's new
proof-graph contribution:

```text
FNV=2280174481469c0d,
SHA256=a63b3309dc30775c82ab7533472f4523ae4bda69525bc262ac06e749011dde23.
```

Removing only those edges from THM-4266's residual gives

```text
count=177,522,
FNV=33142f955cc93379,
SHA256=d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087,
maximum endpoint=688, unique top edge=(520,688).       (15)
```

## 6. Reproduction and scope

From the repository root:

```bash
clang++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_ratio45_common_deck_bridge_thm4267.cpp \
  -o /tmp/thm4267-primary
/tmp/thm4267-primary > /tmp/thm4267-primary.out
cmp /tmp/thm4267-primary.out \
  05-knowledge/results/lrc14_ratio45_common_deck_bridge_thm4267.out

clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_ratio45_common_deck_literal_independent_audit_thm4267.cpp \
  -o /tmp/thm4267-literal
/tmp/thm4267-literal > /tmp/thm4267-literal.out
cmp /tmp/thm4267-literal.out \
  05-knowledge/results/lrc14_ratio45_common_deck_literal_independent_audit_thm4267.out

python3 -B \
  04-computation/lrc14_ratio45_common_deck_proof_graph_independent_audit_thm4267.py \
  --repo . > /tmp/thm4267-graph.out
cmp /tmp/thm4267-graph.out \
  05-knowledge/results/lrc14_ratio45_common_deck_proof_graph_independent_audit_thm4267.out
```

This theorem proves only the fixed pool, primitive ratio `4:5`, and
strict-above-pool scales. It proves no different ratio, smaller scale,
physical entry of an arbitrary row into this pool, or LRC(14). **QED.**
