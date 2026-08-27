---
id: THM-4269
title: "Uniform five-to-six outsider-ray common-deck closure"
status: >
  PROVED RELATIVE TO THM-4231/4267 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  For every integer g>=59 and every nine-body in the fixed thirty-label pool,
  adjoining (5g,6g) leaves Haar-safe mass at least 4/63. A targeted
  3,163-repair deck is active at all 70 bridge scales 59 through 128 and
  closes every labelled body; THM-4231 supplies every g>=129. The exact
  post-THM-4267 contribution is 53 edges, leaving 177,469. No other ratio,
  below-boundary scale, physical entry, or LRC(14) follows.
source: root/lrc-ratio56-cleanroom/2026-08-27
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4267-uniform-four-five-outsider-ray-common-deck-closure
related:
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
  - THM-4262-uniform-three-four-outsider-ray-common-deck-closure
  - THM-4266-three-round-learned-carrier-endpoint-descent
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
common_header_sha256: 03420c21edbd95647503396bcb81d4e34645acbc9fb829d364495cf213e6637d
primary_script: 04-computation/lrc14_ratio56_common_deck_primary_thm4269.cpp
primary_output: 05-knowledge/results/lrc14_ratio56_common_deck_primary_thm4269.out
primary_script_sha256: 9fcb2d076f77043dad8bedf6109e2c381ff9d05c4912dd39d8e4c644528c96c6
primary_output_sha256: 40116a4894f9f12c83dbbde4a49d4c526811567547a41c655743b354114d488e
independent_script: 04-computation/lrc14_ratio56_literal_joint_wall_independent_audit_thm4269.cpp
independent_output: 05-knowledge/results/lrc14_ratio56_literal_joint_wall_independent_audit_thm4269.out
independent_script_sha256: 32408f12eb4a22edaf6640586d693adc3c48ad60121ee31eef4468b08e7fb2d4
independent_output_sha256: 56a77f225e090d0a0a4f6b197e2c4210ad7c7428389ba6ac5052dd99bec81eb2
independent_sanitizer_stderr: 05-knowledge/results/lrc14_ratio56_literal_joint_wall_independent_audit_thm4269_san.err
independent_sanitizer_stderr_sha256: e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
proof_graph_script: 04-computation/lrc14_ratio56_proof_graph_postprocess_thm4269.py
proof_graph_output: 05-knowledge/results/lrc14_ratio56_proof_graph_postprocess_thm4269.out
proof_graph_script_sha256: 7ab8fcb457875917506e19844ded912329667087161258aa2bbe5ee2ee968107
proof_graph_output_sha256: cd04010317dd2a639cfc40fb5a7a5da0c3c84dce41fdbdb9bb99e2af16356728
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. A clean-room checker streams the candidate prefix, constructs
  literal joint walls at every bridge scale, audits all 1,146,880
  repair-scale cells and all 14,307,150 bodies, and reproduces the primary
  primitive, common-deck, hostile-body and targeted-witness ledgers. O0/O3/
  ASan+UBSan literal transcripts agree and sanitizer stderr is empty. A
  portable semantic graph replay reads the canonical THM-4266 residual and
  canonical THM-4267 audit, then recovers the 53 novel edges and final
  residual independently.
---

# THM-4269 -- uniform `5:6` outsider-ray common-deck closure

**PROVED RELATIVE TO THM-4231/4267 + FINITE-EXACT + INDEPENDENTLY AUDITED;
LRC(14) REMAINS OPEN.**

## 1. Statement and boundary

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

Then

```text
for every integer g>=59 and every B in binom(P,9),
mu(G_(B union {5g,6g}))>=4/63.                         (1)
```

This is the strict-above-pool range: `5g>max(P)=290` iff `g>=59`. The exact
finite bridge below covers `59<=g<=128`. THM-4231 proves every distinct
outsider pair with maximum at least `770`, so it supplies `(1)` for `g>=129`,
because `6g>=774`. At `g=58`, the smaller label is the pool element `290`;
no failure or optimality is claimed there.

## 2. Primitive quotient with exact scale sidecar

Let `A=G_{5,6}`. On denominator

```text
N=14*5*6=420,
```

its positive components are exactly

```text
[6,65], [75,78], [90,135], [145,162], [174,205],
[215,246], [258,275], [285,330], [342,345], [355,414]. (2)
```

Their total length is `310`, so `mu(A)=31/42`. The centered primitive prefix

```text
C(t)=420*length(A intersect [0,t])-310*t
```

has exact extrema `-4630` and `4630`, hence range `9260`. Under circle
multiplication `m_g(x)=gx`,

```text
G_{5g,6g}=m_g^(-1)(A).                                (3)
```

Let

```text
D=lcm(14s:s in P)=18,241,159,416,480.
```

For `R in binom(P,8)`, let `U_R=G_(P\R)`, with nonwrapping components
`[a_i/D,b_i/D]`, and define

```text
I(z)=420D integral_0^(z/D) 1_A(t)dt.                  (4)
```

Exact lifting gives

```text
mu(U_R intersect G_{5g,6g})
 =sum_i(I(gb_i)-I(ga_i))/(420Dg).                     (5)
```

Thus `R` is active at scale `g` exactly when

```text
63 sum_i(I(gb_i)-I(ga_i))-4*420*D*g>=0.               (6)
```

The primitive quotient forgets `g`, endpoint phase and repair activation.
Equation `(6)` is the exact scale sidecar; primitive-ratio truth is not copied
blindly along the fibre.

## 3. Exact targeted common deck

Order all `binom(30,8)=5,852,925` repair masks by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).        (7)
```

The first `16,384` masks have FNV `adf20f0ef1cadc1f`. The common-active
prefix ledgers across all seventy bridge scales are

```text
candidate budget   common count   common FNV             body failures
4,096              1,630          e0ffc2e129f24678        116
8,192              3,158          2c493b7ae8a8a847          6
16,384             6,416          9749384b49c2bba0          0. (8)
```

The six bodies missed by the `8,192` budget are

```text
{80,85,88,95,132,143,145,168,252},
{80,85,88,95,132,145,168,193,252},
{16,85,88,95,143,145,168,240,252},
{16,85,88,95,132,168,193,252,290},
{16,85,88,143,168,190,193,252,290},
{16,85,88,95,168,193,252,264,290}.                   (9)
```

Appending the first later common repair for each still-uncovered body adds
exactly five repairs, at candidate ranks `8300,8514,10169,15634,11267` in
that discovery order. The theorem deck has

```text
count=3,163,                    FNV=4e8a91621047de5c.  (10)
```

Exhaustive enumeration of all `binom(30,9)=14,307,150` labelled bodies has
zero failures, `486,631,847` ordered repair/body checks and maximum prefix
`3,163`. Every repair in `(10)` is strictly active on every bridge scale. The
smallest primary numerator margin is

```text
36,651,492,783,360 at g=64,
R={8,30,63,85,88,126,143,286}.                        (11)
```

For a disjoint pair `B,R`,

```text
B subset P\R,
G_((P\R) union {5g,6g}) subset G_(B union {5g,6g}).   (12)
```

The active repair mass therefore lower-bounds the target in the correct
direction. Equations `(6)`, `(10)` and `(12)` prove `(1)` on the bridge;
THM-4231 supplies the tail.

## 4. Independent literal-wall audit

The referee does not use the primitive antiderivative `(4)`. At each scale it
constructs the literal safe intervals of speeds `5g` and `6g` on

```text
L_g=lcm(D,420g),
```

intersects their sorted lists, and sweeps the joint walls against all `7,133`
pool atoms. It streams the first `16,384` candidates through a priority queue
instead of sorting the full repair universe. The full exact matrix has

```text
cells=1,146,880, nonnegative=987,824,
FNV=a039beab75bfbc2c.                                  (13)
```

The literal numerator is

```text
63*(safe ticks on L_g)-4*L_g,
```

and the primary numerator is larger by the exact integer factor `420Dg/L_g`.
The independent audit reproduces every candidate/common/targeted FNV, hostile
body, witness, work count and minimum. Its least theorem-deck literal margin
is

```text
10,908,182,376 at the same g=64 and repair as (11).    (14)
```

Ordinary, optimized and ASan+UBSan literal transcripts are byte-identical;
the sanitizer emits no diagnostics. Primary `-O2/-O3` transcripts also agree.

## 5. Exact proof-graph contribution

The canonical THM-4267 residual has

```text
count=177,522,
FNV=33142f955cc93379,
SHA256=d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087.
```

Semantic selection of the strict `5:6` ray finds exactly the 53 edges

```text
g=59..111.                                             (15)
```

Their ledgers are

```text
FNV=09a794c92e5ab5c9,
SHA256=94a21039c9e379433890e370a3e6985cf5d85e2f61054bc6897d778bffb34823.
```

The `4:5` and `5:6` contributions are disjoint. Removing only `(15)` gives

```text
count=177,469,
FNV=4d1feae0c1e653d5,
SHA256=289cede32347b364123827e7dea02d728b71e8c87d079a9892d3e0492b4a08ae,
maximum endpoint=688, unique top edge=(520,688).       (16)
```

The portable graph audit reads the canonical THM-4266 edge list and verifies
the canonical THM-4267 proof-graph transcript before deriving `(15)`--`(16)`.

## 6. Reproduction and scope

From the repository root after placing the listed artifacts at their canonical
paths:

```bash
clang++ -std=c++20 -O3 -pthread \
  -I 04-computation/lrc14_two_three_outsider_ray_thm4256 \
  04-computation/lrc14_ratio56_common_deck_primary_thm4269.cpp \
  -o /tmp/thm4269-primary
/tmp/thm4269-primary | cmp \
  05-knowledge/results/lrc14_ratio56_common_deck_primary_thm4269.out -

clang++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_ratio56_literal_joint_wall_independent_audit_thm4269.cpp \
  -o /tmp/thm4269-literal
/tmp/thm4269-literal | cmp \
  05-knowledge/results/lrc14_ratio56_literal_joint_wall_independent_audit_thm4269.out -

python3 -B \
  04-computation/lrc14_ratio56_proof_graph_postprocess_thm4269.py \
  --repo . | cmp \
  05-knowledge/results/lrc14_ratio56_proof_graph_postprocess_thm4269.out -
```

This theorem proves only the fixed pool, primitive ratio `5:6`, and
strict-above-pool scales. It proves no different ratio, smaller scale,
physical entry of an arbitrary row into this pool, or LRC(14). **QED.**
