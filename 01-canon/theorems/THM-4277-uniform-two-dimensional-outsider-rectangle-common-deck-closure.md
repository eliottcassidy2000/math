---
id: THM-4277
title: "Uniform two-dimensional outsider-rectangle common-deck closure"
status: >
  PROVED FOR THE DISPLAYED FINITE RECTANGLE + FINITE-EXACT + INDEPENDENT
  LITERAL-WALL AUDIT PASS. For every 450<=q<=499, 600<=r<=650, and every
  labelled nine-body in the fixed thirty-label pool, adjoining q,r leaves
  Haar-safe mass at least 4/63. One 5,257-mask deck is active for all 2,550
  pairs and disjoint-covers all 14,307,150 bodies. The exact post-THM-4276
  contribution is 2,419 residual edges, leaving 172,322 with top layer
  (256,670),(384,670). No neighbouring pair, physical entry, or LRC(14)
  follows.
source: root/lrc-parametric-channel-envelope/2026-08-27
depends_on:
  - THM-4276-six-atom-endpoint-671-augmentation-and-one-layer-descent
related:
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
  - THM-4266-three-round-learned-carrier-endpoint-descent
  - THM-4267-uniform-four-five-outsider-ray-common-deck-closure
artifact_root: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277
artifact_manifest: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/SHA256SUMS
artifact_manifest_sha256: fe4ef5c042274d0d37abb70c48a138b5533842702b174174af2851212b9fc860
common_header: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/lrc_rectangle_common.hpp
common_header_sha256: b0764f8e66fb7909755eb8ff3f56806d4cb776ec501ff5b168865892dcd4de4b
primary_script: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/primary_endpoint.cpp
primary_output: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/primary_endpoint_O3.out
primary_script_sha256: 371da2c67b9079c254ece3f6e37508c5c250db4ca7649ec4a0d231a1d761a163
primary_output_sha256: 33992e5f8e0a92a9a04c6974e432c0bf3e47cdde47d95c122824fc10942574d0
independent_script: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/independent_literal.cpp
independent_output: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/independent_literal_O3.out
independent_script_sha256: edab6ee879681fc2cbbf249067ddb8ed169efbd9d8734fab9e4b384dd2e7f42c
independent_output_sha256: eb6e3c5f8f96c6eb47f0940c9e9ccdafc5ccc3349d1d634faa49350b2ed27f30
transcript_audit_script: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/verify_transcripts.py
transcript_audit_output: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/verify_transcripts.out
transcript_audit_script_sha256: d8a6e31c4f653ada95e43a068736e0d7d773ce57a839f03dcd5ffe41e7ccd70d
transcript_audit_output_sha256: 9b23fbddb2225f9afb0579357ecfb1d612e816a8cd1d78652cfb8c036a8d86ba
comparator_audit_script: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/u256_comparator_audit.py
comparator_audit_output: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/u256_comparator_audit.out
comparator_audit_script_sha256: 6343d641803c0f4b9c05bf518997755f4404d63d490ff041551184549ad5d98c
comparator_audit_output_sha256: d7985b78d3b8f7ec9769c6f3e8ec7760c0025d776b7b490889fc99a8f704f48b
postprocess_script: 04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277/postprocess_current.py
postprocess_output: 05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277/postprocess_current.out
postprocess_script_sha256: 3a7c8aa6c79d26d391125725b470c1b455353ef1a994b023c80d29959561fe95
postprocess_output_sha256: aed1a6329115fdb8da5ad511e221653c814b2517c202aa5e0dbb2083819b09da
hash_basis: raw LF bytes; C++ semantic outputs omit only the SECONDS line
audit: >
  PASS / ACCEPT. Corrected primary endpoint and detached literal-wall
  engines agree on every reported pair/candidate/matrix/deck/body ledger and
  on the exact weakest pair, repair, mass, and reduced positive gap. O3,
  O3+NDEBUG, and O1+UBSan semantic outputs are byte-identical within each
  engine; sanitizer stderr is empty. A detached Python-bigint comparator audit
  passes 81 boundary products and 100,000 deterministic comparisons. Normal
  and optimized Python replays agree for the transcript verifier, comparator
  oracle, and exact post-THM-4276 proof-graph reconstruction.
---

# THM-4277 -- uniform two-dimensional outsider-rectangle common-deck closure

**PROVED FOR THE DISPLAYED FINITE RECTANGLE + FINITE-EXACT + INDEPENDENT
LITERAL-WALL AUDIT PASS; LRC(14) REMAINS OPEN.**

## Statement

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290},
G_A={x in R/Z:min_(a in A)||ax||>=1/14}.
```

For every pair of integers

```text
450 <= q <= 499,   600 <= r <= 650,
```

and every `B in binom(P,9)`, one has

```text
mu(G_(B union {q,r})) >= 4/63.                         (1)
```

This is a genuinely two-dimensional finite family of `50*51=2,550` pairs.
It contains `2,500` distinct primitive ratios: `2,474` occur at only one
scale; 16 ratios occur twice, seven occur three times, two occur five times,
and one occurs 13 times. Thus `(1)` is not a finite list on one or a few rays.

## Exact pair sidecars and activation predicate

Let

```text
D=lcm(14s:s in P)=18,241,159,416,480.
```

For a pair in the rectangle, put

```text
g=gcd(q,r),   u=q/g,   v=r/g,   N=14uv,   A_(u,v)=G_{u,v}.
```

The family has `77` distinct gcd values, with minimum `1` and maximum `162`.
Neither the primitive ratio nor the scale may be discarded. Let the positive
components of `U_R=G_(P\R)`, for `R in binom(P,8)`, be the nonwrapping lifts
`[a_i/D,b_i/D]`, and define

```text
I_(u,v)(z)=ND integral_0^(z/D) 1_(A_(u,v))(t)dt.       (2)
```

Splitting at primitive walls gives the exact identity

```text
mu(U_R intersect G_q intersect G_r)
 =sum_i(I_(u,v)(g b_i)-I_(u,v)(g a_i))/(NDg).         (3)
```

Hence `R` is active for `(q,r)` exactly when

```text
63 sum_i(I_(u,v)(g b_i)-I_(u,v)(g a_i))-4NDg >= 0.   (4)
```

The primary implementation builds the primitive wall arcs for every one of
the 2,550 typed triples `(u,v,g)` and evaluates `(4)` with signed 128-bit
integer arithmetic. Consequence-side comparisons between two positive
normalized gaps use an exact unsigned 256-bit product, not signed-128
cross-multiplication. It does not infer one pair from a neighbour, quotient by
scale, or use floating point.

## One common-active deck

Order all `binom(30,8)=5,852,925` repairs by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).        (5)
```

The first `16,384` candidates have FNV
`adf20f0ef1cadc1f`. Intersect their active subdecks over every pair in the
rectangle. The resulting common deck is

```text
count=5,257, FNV=60f329212844f8ac.                    (6)
```

The primary checks all

```text
2,550*16,384=41,779,200
```

activation cells. Exactly `36,657,425` are nonnegative, none is an equality,
and the pair-major, candidate-minor activation-sign byte stream has FNV
`9d3e995e23a7695a`.

The smallest invariant normalized gap among common-deck cells occurs at

```text
(q,r)=(462,626),
R={8,16,30,63,88,143,240,252},
mass=16477591782853/259521949879920,
mass-4/63=7663493/259521949879920 > 0.                 (7)
```

## Repair hypergraph consequence

Let `E_rect` be the deck in `(6)`. Exact enumeration of every
`binom(30,9)=14,307,150` labelled body proves

```text
for every B there exists R in E_rect with B intersect R=empty. (8)
```

For such a repair,

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}).        (9)
```

Every `R in E_rect` is active for every rectangle pair, so `(3)`, `(8)`, and
`(9)` prove `(1)`. The exhaustive body scan has zero failures,
`508,103,822` ordered checks, and maximum checked prefix `4,129`. The latter
occurs at

```text
B={16,85,88,95,145,168,193,240,252},
R={63,80,120,143,190,264,286,290}.
```

There is also a universal, activity-independent combinatorial lower bound on
any deck satisfying `(8)`. Delete duplicate repair masks, then complement each
distinct eight-repair to a 22-subset of the thirty-label pool. Condition `(8)`
says that these 22-blocks cover every 9-subset, so the number of distinct masks
is at least the covering number `C(30,22,9)`. In a `(v,k,t)` cover, the blocks
through each point cover all `(t-1)`-subsets of the remaining `v-1` points.
Summing these point degrees counts each block `k` times and gives the
Schoenheim recursion

```text
C(v,k,t) >= ceiling((v/k) C(v-1,k-1,t-1))
```

Iterating from the innermost bound outward yields

```text
C(22,14,1)>=2, then 4,6,9,13,19,27,38,52.
```

Consequently every such deck has at least `52` distinct eight-subset repair
masks.
No equality case, attainability of `52`, or minimality of the present
5,257-repair deck is claimed.

## Deterministic half-budget hostile

The first `8,192` candidates have FNV `60148ca1fc61dbcb`. Their common-active
subdeck has

```text
count=2,572, FNV=44a01e1ab114723e,
```

but misses exactly seven bodies. The first is

```text
{16,85,88,95,143,145,168,193,252}.                   (10)
```

Thus a large common-active intersection does not by itself imply the needed
transversal bound. Doubling the deterministic budget supplies the missing
repair addresses and closes the hypergraph. No minimal common deck or sharp
candidate budget is claimed.

## Independent literal-wall route

For each literal pair, the independent source forms

```text
L_(q,r)=lcm(D,14q,14r).                                (11)
```

It constructs the safe intervals for the two actual speeds `q,r` directly on
this grid, intersects the sorted interval lists, and integrates each fixed-pool
repair component after the exact scaling `L_(q,r)/D`. It never builds the
primitive quotient, calls `(2)`, or uses the primary full-sort candidate
selection: the first 16,384 candidates are selected by a bounded heap. The two
engines intentionally share the fixed-pool cell/repair-geometry builder and
basic mask/FNV utilities. Pair integration, candidate selection, and complete
body enumeration are detached.

The frozen transcript cross-check obtains byte agreement on the reported pair,
candidate, matrix, common-deck, and body-scan ledger lines, including the FNV
of the 41,779,200 activation-sign bytes. It also obtains exact agreement on the
pair, repair address, mass, and reduced gap in `(7)`. The raw signed margins
differ, as expected, because the two grids use different exact integer
normalizations.

## Arithmetic audit of the normalized-gap selector

Activation, deck intersection, and both body scans use only the sign of the
signed-128 quantity in `(4)`; the weakest-cell selector is consequence-side
diagnostic information. An earlier exploratory selector compared two such
positive rational gaps by multiplying signed-128 factors. Those products can
exceed `2^127-1`, so that exploratory weakest-cell line is not retained.

Both frozen engines instead multiply two nonnegative signed-128 factors into
eight little-endian 32-bit limbs. In each schoolbook inner step,

```text
(2^32-1)^2 + (2^32-1) + (2^32-1) = 2^64-1,
```

so the unsigned-64 accumulator is exact; lexicographic comparison from the
highest limb gives the exact order of the two 256-bit products. A detached
Python-bigint oracle checks 81 boundary products and 100,000 deterministic
random comparisons. The final replays agree on the corrected pair, repair
address, mass, and reduced positive gap. Their raw margins need not agree
because the primitive and literal grids use different exact integer
normalizations.

## Current proof-graph consequence after THM-4276

The canonical post-THM-4276 residual has

```text
count=174,741,
FNV=f13745b05320f83c,
SHA256=51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0.
```

Its intersection with the rectangle is exactly

```text
count=2,419,
FNV=67b373dc22ac870d,
SHA256=e9ee49675d5345b06157e64a060506be3f6dd6a94835cc92f6eb5138f346cffb. (12)
```

This set is disjoint from the already proved `4:5`, `5:6`, `3:5`, `7:8`,
`8:9`, and `11:12` outsider-ray components, from THM-4271's `r>=671`
descent, and from THM-4276's 163-row `r>=670` deletion. Removing only `(12)`
from the post-THM-4276 residual gives

```text
count=172,322,
FNV=30b2a7e597ac548c,
SHA256=7228658eae4067db4bbcceb6c9b1ccebf1bd3e6f128e202ea184854acc53f309,
maximum endpoint=670,
top layer={(256,670),(384,670)}.                       (13)
```

The inequality theorem is absolute for the displayed rectangle; only the
novelty accounting in `(12)--(13)` depends on the current proof graph.

## Connection ledger and scope

```text
source:  the pair-indexed active hypergraphs E(q,r)
target:  E_rect = intersection_(q,r in rectangle) E(q,r)
map:     retain the same labelled repair mask in every fibre
preserved predicate: exact activation and body/repair disjointness
destroyed data: pair-specific repairs, slack profile, and failure address
sidecars retained: (q,r,u,v,g), endpoint phase, exact signed margin
cheapest decisive hostile: the seven failed bodies at budget 8,192
```

The closest proved mechanism is THM-4267's single-ray common deck with a
scale sidecar. The corrected near miss is THM-4254/4266's pair-labelled
carrier: a union of masks is not a common-active deck. The present theorem
retains both primitive ratio and scale and intersects activity before testing
the hypergraph.

This proves only the fixed thirty-label pool and the displayed finite integer
rectangle. It proves no neighbouring pair, monotonicity in either coordinate,
cofinite two-dimensional wedge, physical entry of an arbitrary LRC row into
this pool, or LRC(14).
