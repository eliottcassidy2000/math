---
id: THM-3585
title: "LRC common A4 channel plane and centering complement"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Over the pinned prime
  755373809845391722745761, the seventeen-row nested source/current tensor
  and the 169-row two-current tensor have the same rank-four raw image in
  Fun(V4 x F13,F_p).  Their state/relation-centred images are also the same
  rank-four plane.  The common raw and common centred planes intersect only
  in zero, so their stack has rank eight.  This is a static finite-field
  identity only: no chronology, current, physical entry, row exclusion,
  characteristic-zero lift, or LRC(14) conclusion is proved.
source: kps-s188 / delegated LRC channel comparison, 2026-08-21
audit: >
  An independent reconstruction checked both typed tensors entrywise, used a
  separately written modular eliminator to recover every raw, centred, stack,
  and four-axis rank and canonical basis, verified all 186 margin rows, the
  prime/congruence hypotheses, both replays, and every pinned digest.  The
  finite-field/static-only and ambient-complement scope restrictions survived.
two_current_parent: 04-computation/lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_probe_20260816.py
two_current_parent_sha256_lf: 3dab580e479e4ba7ac8801c1e5d8523018e0b3dc1c2176c072e7c609033eb6c8
source_current_parent: 04-computation/lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_probe_20260816.py
source_current_parent_sha256_lf: 1188df8aa2a7a84c1e8ada5fc3cc8d3b839ece70298b94f1d94c9d440caa88f3
exact_script: 04-computation/lrc_common_a4_channel_plane_centering_complement_thm3585.py
exact_output: 05-knowledge/results/lrc_common_a4_channel_plane_centering_complement_thm3585.out
exact_script_sha256_lf: dad78dd317a25958f42b3df44f395d4da344330b07ea4f5041a7e838e46dd0e8
exact_output_sha256_lf: c0dbb1792545af29c3cfbd16b22efe24e5c32fecb25a3d78321a9dc93a22929f
hash_basis: LF-normalized bytes
---

# THM-3585 -- LRC common A4 channel plane and centering complement

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Set

```text
p = 755373809845391722745761,
H = Fun(V4(state) x F13(relation),F_p),
dim H = 52.
```

Reconstruct the following two tensors from their hash-pinned parent programs.

```text
T_source(channel,state,relation),       17 x 4 x 13,
T_two(state,r0,r1,relation),             4 x 13 x 13 x 13.
```

The words `source` and `current` in these names are inherited artifact labels;
they do not assert a temporal order.  Flatten the first tensor by `channel`
and the second by `(r0,r1)`, with common target-coordinate order
`(state,relation)`.  This gives row families

```text
R_source in H^17,       R_two in H^169.
```

Their row spaces are equal and four-dimensional:

```text
A_raw := rowspan(R_source) = rowspan(R_two),
dim A_raw = 4.                                             (1)
```

Let

```text
C4  = I4  - J4/4,
C13 = I13 - J13/13,
C   = C4 tensor C13.                                      (2)
```

Since `p` divides neither `4` nor `13`, this is a lawful finite-field
projection.  Apply `C` on the common target coordinates of every row.  The
two centred images again coincide and have dimension four:

```text
A_ctr := rowspan(R_source C) = rowspan(R_two C),
dim A_ctr = 4.                                             (3)
```

Finally,

```text
dim(A_raw + A_ctr) = 8,
A_raw intersect A_ctr = {0}.                              (4)
```

Thus `A_raw` and `A_ctr` are complementary **to each other inside their
eight-dimensional sum**.  Neither is claimed to be a complement in all of
the 52-dimensional target.

## 2. Pinned parent reconstructions

The verifier first checks the LF-normalized hashes of the two imported parent
programs:

```text
two-current parent =
  3dab580e479e4ba7ac8801c1e5d8523018e0b3dc1c2176c072e7c609033eb6c8,

source/current parent =
  1188df8aa2a7a84c1e8ada5fc3cc8d3b839ece70298b94f1d94c9d440caa88f3.
                                                               (5)
```

Those parents in turn pin the lower input programs.  In typed order

```text
(two-current owner,
 source/current source-sheet,
 source/current inverse-owner),
```

their hashes are

```text
(ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348,
 592aa0bce31f2da5d5e2ddff7f3ffe6f1398f3a07b5ce927e0d97c9fe309ae3b,
 ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348).
                                                               (6)
```

All thirteen character slices are recomputed by each parent.  No stored
rank or previously printed row-space digest is accepted in place of the
tensors.  The resulting reconstruction hashes, in order

```text
(T_two, two-current tau core,
 source/current gamma bank, T_source),
```

are

```text
(53eb6e618d0669bdb27841a1800c46e32456bb6c6d3698c590ae0d5e68822033,
 fd1b837d9de3e4f9e586d29b69ed6726364ce97535d2e48f441c6ccd694250de,
 9d0aa6823ed9bb83b338350d9da3d86c52e552a844b4334e3571bf2393f04cd1,
 53e1d8e6de69dcd1cf5740acce7b245f7d9a388cfa2b27da4aefb00f3075df14).
                                                               (7)
```

Equation (7) is the repaired parent gate.  In particular, the source/current
tensor is checked against the explicit fourth constant in (7); the broken
reference to a nonexistent parent `EXPECTED_DIGESTS` attribute is absent.

## 3. Canonical row-space proof

For each row family the verifier performs Gaussian elimination over `F_p`,
normalizes every pivot to one, clears each pivot column above and below, and
retains the nonzero rows.  This produces the unique reduced-row-echelon basis,
so equality below is equality of actual subspaces in the same 52 coordinates,
not an inference from equal dimensions.

For the raw families the exact record is

```text
(rank R_source, rank R_two, rank [R_source;R_two]) = (4,4,4),
canonical RREF SHA-256 =
  1d9293d05fa3551b785a1537e78bc8be585fcc43dbb5172036c9b32546ca8560.
                                                               (8)
```

Both canonical bases in (8) agree entrywise, proving (1).  After applying
the double-centering projection (2), the record is

```text
(rank R_source C, rank R_two C,
 rank [R_source C;R_two C]) = (4,4,4),
canonical RREF SHA-256 =
  0cfa2e3330f92ab59fd183e5664715c490a702df2ad74491a8180793cae4a21e.
                                                               (9)
```

Every centred row is separately checked to have zero state margins and zero
relation margins.  The inherited two-current four-way ANOVA tensor has the
same canonical basis as (9), with individual and joined ranks `(4,4)`.

Finally, stacking raw and centred rows gives

```text
(source stack rank, two-current stack rank, all-row stack rank) = (8,8,8),
canonical RREF SHA-256 =
  9f4ec33d31337b0100a55871f6284443c6e5cfc0f5133a1493f807d563670821.
                                                               (10)
```

The dimension formula applied to (8)--(10) gives

```text
dim(A_raw intersect A_ctr) = 4 + 4 - 8 = 0,
```

which proves (4).  The complete semantic surface comprising the prime,
parent hashes, reconstruction hashes, rank records, and three canonical RREF
hashes has digest

```text
861d4d95fd834e62cec842a5fc548e779554b17c79fabcf79cf7338e9db29848.
                                                               (11)
```

## 4. Scope and corrected interpretation

This theorem does not use Fourier-support fullness.  MISTAKE-417 shows that
a separable delta-cell lift can have full Fourier support while retaining
matrix rank one.  THM-3585 instead compares canonical bases of two explicitly
reconstructed row spaces in one typed target.

The conclusion also does not identify the common four-plane with the pointed
six-dimensional relation carrier `P6`.  The two numbers are ranks of different
tensor flattenings:

```text
A4: rows indexed by character/address, target V4 x F13;
P6: pointed relation carrier, target F13 after its typed marginal.
```

Equality of the two `A4` presentations therefore supplies no chronology,
source-to-arrival map, current, physical-entry theorem, closure edge, scalar
row exclusion, or LRC(14) result.  It also proves nothing in characteristic
zero: the statement is exactly over the pinned finite field.

## 5. Reproduction

Run

```text
python -B 04-computation/lrc_common_a4_channel_plane_centering_complement_thm3585.py
python -B -O 04-computation/lrc_common_a4_channel_plane_centering_complement_thm3585.py
```

Fresh normal and optimized replays both pass and have byte-identical captured
stdout with SHA-256

```text
56e6853e2edc6bab4b34694915ed3282982f10bcd15f1203d24af099a739ab88.
```

The maintained result file is LF-normalized and has SHA-256

```text
c0dbb1792545af29c3cfbd16b22efe24e5c32fecb25a3d78321a9dc93a22929f.
```

**End of proof.**
