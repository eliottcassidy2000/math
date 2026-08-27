---
id: THM-4261
title: "Semantic endpoint-band prefix-union lift"
status: >
  PROVED RELATIVE TO THM-4231/4254/4256 + FINITE-EXACT + INDEPENDENTLY
  AUDITED. For every current post-THM-4256 residual pair with second endpoint
  from 733 through 754, and every nine-body in the fixed thirty-label pool,
  adjoining the pair leaves Haar-safe mass at least 4/63. There are exactly
  297 such pairs. The 4,675-mask set-theoretic union inherited from THM-4254
  supplies 1,344,725 pair-active incidences and closes all 4,249,223,550
  labelled body cases. Removing the band leaves 180,694 residual edges with
  maximum endpoint 732. The same carrier misses exactly two bodies at the
  boundary pair (542,732), so no endpoint-732 uniformity is claimed.
source: root/cross-frontier-overnight/2026-08-27
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
  - THM-4256-uniform-two-three-outsider-ray-endpoint-cocycle-closure
related:
  - THM-4252-exact-signed-endpoint-cocycle-residual-edge-closure
primary_script: 04-computation/lrc14_prefix_union_semantic_band_exact_verifier_thm4261.cpp
primary_output: 05-knowledge/results/lrc14_prefix_union_semantic_band_exact_thm4261.out
primary_script_sha256: 8f59d03785fed46bdf1f6d131b2a90aec06bc80449b0d634dfb603b4ff41470b
primary_output_sha256: 32c2a8c4f6695bf4bc7e3603865356e4ae8e939fa70cbfa5cd54cd2c3bd2ae33
independent_script: 04-computation/lrc14_prefix_union_semantic_band_direct_wall_audit_thm4261.cpp
independent_output: 05-knowledge/results/lrc14_prefix_union_semantic_band_direct_wall_audit_thm4261.out
independent_script_sha256: b55c4fb98e514402cd1900666656e863f4bd3e5c3fd35cac5e2de4452ca7eccb
independent_output_sha256: 22a2840f5c62eec3dbaecdcfcbabd3f7a29f502dff4a87dd0642fb4a6d245031
postprocess_script: 04-computation/lrc14_prefix_union_semantic_band_residual_postprocess_thm4261.py
postprocess_output: 05-knowledge/results/lrc14_prefix_union_semantic_band_residual_postprocess_thm4261.out
postprocess_script_sha256: 8891e1b6ba38fc70cdd46001249b8b053a3f091d83c4c09095181c5303bd0d7e
postprocess_output_sha256: c1de6e52865f729e05725cda9a6e1d4b041ecfd1948fb85b9303877758190e05
hash_basis: raw bytes
audit: >
  PASS / ACCEPT. The primary checker rebuilds every candidate mass from
  signed primitive endpoint-cocycle atoms. The clean-room checker instead
  constructs fresh literal joint walls, aggregates direct safe-cell widths,
  and independently regenerates every nine-body. They agree on every active
  incidence and on 125,275,100,068 ordered repair/body checks. The endpoint
  732 hostile is reproduced by both implementations. No sampling or floating
  point participates.
---

# THM-4261 -- semantic endpoint-band prefix-union lift

**PROVED RELATIVE TO THM-4231/4254/4256 + FINITE-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and semantic universe

For a finite positive set `A`, put

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
alpha=4/63,
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}.
```

Let `C` be the set of pairs in the exact post-THM-4256 residual whose second
endpoint lies in `[733,754]`. Then

```text
for every (q,r) in C and every B in binom(P,9),
mu(G_(B union {q,r}))>=4/63.                            (1)
```

This is a semantic selection made before any certificate outcome. Exact
reconstruction of the inherited residual gives

```text
|C|=297,
FNV-1a-64=e923d1494185b820,
SHA256(lines q,r\n)=745ef7c8809335e6d6e9623314beff9
                     17edc71cfaaaa88e7210ede9dcd97d11b. (2)
```

The endpoint-layer profile from `733` upward is

```text
23,23,17,19,13,18,17,16,15,16,12,16,8,9,10,14,12,12,9,6,8,4. (3)
```

The postprocessor freezes the complete pair list and both cryptographic and
noncryptographic ledgers; `(2)` is not a sample or a hand-selected band.

## 2. The inherited carrier and typed lift

THM-4254 has 59 pair-labelled short repair prefixes in the same fixed pool.
Taking distinct repair masks in first-occurrence order gives the inherited
carrier

```text
U subset binom(P,8),        |U|=4,675,
FNV(U)=ce4e76ec11df057c.                                  (4)
```

Crucially, `U` is only a set of **candidates**. It is not asserted that one
repair is active for every pair. For each target `(q,r)`, define

```text
U(q,r)={R in U:mu(G_((P\R) union {q,r}))>=4/63}.        (5)
```

The primary verifier rebuilds the primitive pair `(q/g,r/g)`, the exact
signed endpoint-cocycle atoms, and every mass in `(5)`. It finds no equality
at the threshold. Across the 297 targets,

```text
sum_(q,r)|U(q,r)|=1,344,725,
min |U(q,r)|=3,210 at (704,744),
max |U(q,r)|=4,675 at (714,733).                        (6)
```

Thus the transfer from THM-4254 is lawful: it preserves the repair masks and
their order, but recomputes the pair-dependent active predicate instead of
transporting truth across a changed pair.

## 3. Repair hypergraph consequence

For every `(q,r) in C`, exact enumeration proves

```text
for every B in binom(P,9), there exists R in U(q,r) with B intersect R=empty.
                                                                  (7)
```

For such `B,R`,

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}).         (8)
```

The inclusion in `(8)` points from the more constrained certified set to the
target. Equations `(5)`, `(7)`, and `(8)` prove `(1)`. The exact workload is

```text
candidate pair-mask tests          1,388,475,
labelled body cases              4,249,223,550,
ordered repair/body checks       125,275,100,068.       (9)
```

The pair-active incidence ledger is `247b8bb8de3112c0`; the ledger of
`(q,r,|U(q,r)|)` is `c2657c090f345f51`.

## 4. Independent literal-wall audit

The independent implementation does not use a primitive cocycle, endpoint
atoms, the discovery superset transform, or the primary body generator. For
each pair it:

1. constructs the joint wall arrangement of `P union {q,r}` on an exact
   integer grid;
2. evaluates every open cell at its midpoint and aggregates literal cell
   widths by the pool-failure mask;
3. thresholds all 4,675 carrier masks by direct width summation; and
4. recursively generates all `binom(30,9)=14,307,150` bodies.

Across the band this rebuilds

```text
joint cells                         2,944,191,
literal pair-safe mask atoms          770,215.          (10)
```

It agrees exactly with `(6)` and `(9)`, including both ledgers, and has zero
body failures in every row. This checks the cocycle formula and the body
quantifier through materially different inner loops.

## 5. Sharp carrier boundary

Apply the same inherited carrier at the next residual layer to `(542,732)`.
Both checkers obtain

```text
|U(542,732)|=3,227,
uncovered nine-bodies=2,
ordered repair/body checks=474,614,827.                (11)
```

The primary checker's first uncovered body mask is `0x151a5400`. Hence the
fixed carrier does **not** uniformly close endpoint 732. This is a sharp
boundary only for this carrier and semantic-band proof: it is not a Lonely
Runner counterexample, not proof that `(542,732)` is intrinsically hard, and
not a claim that the full repair universe fails.

## 6. Proof-graph consequence and scope

Removing `C` from the post-THM-4256 residual gives

```text
count=180,694,
FNV=50c911cf48e3f50a,
SHA256=19906a8f773517f0b29767cb16b3b6450
       2fa38c7d03dfbbbfaeff87ba71c702c,
maximum second endpoint=732.                           (12)
```

The top layer contains exactly 24 edges, frozen by the postprocessor. The
cumulative removal from the post-THM-4242 residual is `432` edges: three from
THM-4252, 59 from THM-4254, 73 from THM-4256, and 297 here.

This theorem proves neither all pairs in a numerical rectangle nor any
endpoint-732 row, arbitrary physical entry, a globally minimal repair deck,
or LRC(14).

## 7. Reproduction

From the repository root, compile the two checkers and feed each the semantic
pair stream:

```bash
g++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_prefix_union_semantic_band_exact_verifier_thm4261.cpp \
  -o /tmp/thm4261-primary
python3 -B \
  04-computation/lrc14_prefix_union_semantic_band_residual_postprocess_thm4261.py \
  --pairs | /tmp/thm4261-primary \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254

g++ -std=c++20 -O3 -pthread \
  04-computation/lrc14_prefix_union_semantic_band_direct_wall_audit_thm4261.cpp \
  -o /tmp/thm4261-direct
python3 -B \
  04-computation/lrc14_prefix_union_semantic_band_residual_postprocess_thm4261.py \
  --pairs | /tmp/thm4261-direct \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band
```

The primary and independent transcripts each have 302 lines and end in the
stated PASS verdict. **QED.**
