---
id: THM-4262
title: "Uniform three-to-four outsider-ray common-deck closure"
status: >
  PROVED RELATIVE TO THM-4231/4256 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  For every integer g>=97 and every nine-body in the fixed thirty-label pool,
  adjoining (3g,4g) leaves Haar-safe mass at least 4/63. A 3,751-repair deck
  is active at all 96 bridge scales 97 through 192 and closes every labelled
  body; THM-4231 supplies every g>=193. The exact post-THM-4256 residual
  contribution is 72 edges. No other ratio, below-boundary scale, physical
  entry, or LRC(14) follows.
source: root/cross-frontier-overnight/2026-08-27
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
primary_script: 04-computation/lrc14_ratio34_common_deck_bridge_thm4262.cpp
primary_output: 05-knowledge/results/lrc14_ratio34_common_deck_bridge_thm4262.out
primary_script_sha256: 71fcc341deb7dfc207b3f5069c2108047b354b789dc7ba3e8905cb9e45808406
primary_output_sha256: 50618b93f7cbb583c4fe1aff7b2f57d52444e6aec8bde308e41fded4a4b95d93
independent_script: 04-computation/lrc14_ratio34_literal_joint_wall_audit_thm4262.cpp
independent_output: 05-knowledge/results/lrc14_ratio34_literal_joint_wall_audit_thm4262.out
independent_script_sha256: 90e533cbaa0ba93777380d42c3baaaabd65808f52eb8abb897d4aef83c16fb29
independent_output_sha256: 6c0a88952b3d580d373da056303841d7ab7435641d0e8b295913ebd6ae4dce6d
postprocess_script: 04-computation/lrc14_ratio34_residual_postprocess_thm4262.py
postprocess_output: 05-knowledge/results/lrc14_ratio34_residual_postprocess_thm4262.out
postprocess_script_sha256: 141d20aed4822d35b09fa95b6f13cc1c018d98451b2b3e735456cbbe3212dbb2
postprocess_output_sha256: 4760d9ce8304b38b105dafe57b6203a7253d02c370423546b233e8ae2f9bb0da
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary checker uses the exact primitive 3:4 prefix
  integral. A clean-room implementation instead builds literal joint safe
  intervals at every scale, streams the candidate order through a priority
  queue, and audits all 786,432 repair-scale cells. They agree on the
  3,751-mask common deck, minimum margin, hostile truncation, and full body
  scan. O0/O3 and assertion-enabled outputs agree after newline normalization.
---

# THM-4262 -- uniform three-to-four outsider-ray common-deck closure

**PROVED RELATIVE TO THM-4231/4256 + FINITE-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Then

```text
for every integer g>=97 and every B in binom(P,9),
mu(G_(B union {3g,4g}))>=4/63.                         (1)
```

The lower endpoint is the first strict-above-pool scale:
`3g>max(P)=290` iff `g>=97`. No theorem failure or optimality is claimed for
smaller `g`.

The proof order is finite plus cofinite. The exact bridge below covers
`97<=g<=192`. THM-4231 proves every distinct outsider pair with maximum at
least `770`; hence it supplies `(1)` for `g>=193`, since `4g>=772`.

## 2. Primitive quotient with the scale sidecar retained

Let `A=G_{3,4}`. On denominator

```text
N=14*3*4=168,
```

its positive-length components are exactly

```text
[4,39], [45,52], [60,81], [87,108], [116,123], [129,164]. (2)
```

Their total length is `126`, so `mu(A)=3/4`. Under circle multiplication
`m_g(x)=gx`,

```text
G_{3g,4g}=m_g^(-1)(A).                                (3)
```

Let

```text
D=lcm(14s:s in P)=18,241,159,416,480.
```

For `R in binom(P,8)`, let `U_R=G_(P\R)`, whose positive components have
nonwrapping lifts `[a_i/D,b_i/D]`. Define

```text
I(z)=168D integral_0^(z/D) 1_A(t)dt.                   (4)
```

Exact lifting gives

```text
mu(U_R intersect G_{3g,4g})
 =sum_i(I(gb_i)-I(ga_i))/(168Dg).                     (5)
```

Thus `R` is active at scale `g` exactly when

```text
63 sum_i(I(gb_i)-I(ga_i))-4*168*D*g>=0.               (6)
```

Quotienting `(3g,4g)` by its gcd preserves the ratio but destroys `g`, the
endpoint phase, and repair activation. Equation `(6)` is the required scale
sidecar; no primitive-ratio truth is copied blindly along the fibre.

## 3. Exact common-deck bridge

Order all `binom(30,8)=5,852,925` repair masks by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).        (7)
```

The first 8,192 candidates have FNV `60148ca1fc61dbcb`. Retain a repair only
if `(6)` holds at every one of the 96 integer scales from 97 through 192.
The common deck has

```text
count=3,751,
FNV=d708c4fa029105e1.                                  (8)
```

There are no activation equalities. The least exact numerator margin is

```text
1,171,377,341,568 at g=122,
R={20,60,88,168,190,252,264,290}.                      (9)
```

Exhaustive enumeration proves that every one of the
`binom(30,9)=14,307,150` labelled bodies is disjoint from a repair in `(8)`.
The scan uses `467,450,372` ordered repair/body checks, has maximum prefix
`2,522`, and has zero failures.

For a disjoint pair `B,R`,

```text
B subset P\R,
G_((P\R) union {3g,4g}) subset G_(B union {3g,4g}).    (10)
```

The active repair mass therefore lower-bounds the target in the correct
direction. Equations `(6)`, `(8)`, and `(10)` prove `(1)` on the finite bridge;
THM-4231 supplies the cofinite tail.

### Hostile candidate budget

Among only the first 4,096 candidates, the 1,937 common-active repairs miss
exactly two bodies. The first is

```text
{8,80,88,143,168,190,193,252,264}.                    (11)
```

The full deck closes both. This is a sharp failure of the smaller common-deck
certificate, not a counterexample to `(1)`.

## 4. Independent literal-wall audit

The independent checker does not evaluate the primitive prefix `(4)`. At each
scale it constructs the literal safe intervals of speeds `3g` and `4g` on

```text
L_g=lcm(D,168g),                                       (12)
```

intersects their sorted interval lists, and sweeps the joint walls against
all 7,133 pool atoms. Candidate selection is independently streamed through
a priority queue before sorting only the retained 8,192 masks.

It evaluates all

```text
8,192*96=786,432 repair-scale cells                   (13)
```

and reproduces

```text
candidate FNV                 60148ca1fc61dbcb,
common count/FNV              3751 / d708c4fa029105e1,
full mass-matrix FNV          1e32f8e31b79a257,
nonnegative cells             716,014,
body failures/checks          0 / 467,450,372.         (14)
```

Its reduced-denominator minimum is `3,486,242,088` at the same scale and
repair; multiplying by the exact conversion factor `336` gives `(9)`. It also
finds exactly the two hostile bodies from the 4,096-candidate truncation.

## 5. Exact proof-graph contribution

Semantic selection from the post-THM-4256 residual finds exactly 72 new ray
edges, at

```text
g=97..166,172,180.                                     (15)
```

Their ledgers are

```text
FNV=512cbba28e2235fd,
SHA256=b1c89073d9b82351b663e97c18c807f0
       3f3fd2d40ddcfafe038d8cad0535cb2c.                (16)
```

Removing only these proof-graph-new edges gives

```text
count=180,919,
FNV=9fae8a515ea17db3,
SHA256=a44da1b8dd6daa484ae2133f6fdd4c4c8
       d753c25292280231b4d1aefd5af9a0e.                 (17)
```

The maximum endpoint remains 754. THM-4261's 297-edge band is disjoint from
this 72-edge ray contribution; applying both leaves 180,622 residual edges.

## 6. Reproduction and scope

```bash
g++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_ratio34_common_deck_bridge_thm4262.cpp \
  -o /tmp/thm4262-primary
/tmp/thm4262-primary

g++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_ratio34_literal_joint_wall_audit_thm4262.cpp \
  -o /tmp/thm4262-direct
/tmp/thm4262-direct

python3 -B 04-computation/lrc14_ratio34_residual_postprocess_thm4262.py
python3 -B -O 04-computation/lrc14_ratio34_residual_postprocess_thm4262.py
```

This theorem proves only the fixed pool, the primitive ratio `3:4`, and
strict-above-pool scales. It proves no different ratio, smaller bridge deck,
physical entry of an arbitrary row into this pool, or LRC(14). **QED.**
