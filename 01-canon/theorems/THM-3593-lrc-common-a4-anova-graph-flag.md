---
id: THM-3593
title: "LRC common A4 ANOVA graph flag"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Over the pinned field of THM-3585, the common raw A4 plane is the
  graph of an isomorphism between a four-dimensional additive ANOVA plane
  and its four-dimensional interaction image.  The additive plane in turn
  projects isomorphically onto a pure-relation four-plane and has only a
  rank-two state/constant correction, with a two-dimensional pure-relation
  kernel.  This is a static finite-field flag only: no chronology, current,
  physical entry, characteristic-zero lift, or LRC(14) conclusion is proved.
source: kps-s188 / LRC ANOVA graph continuation, 2026-08-21
audit: >
  Independently accepted by agent Godel on 2026-08-21 after reconstructing
  all ANOVA identities on the 52 coordinate vectors, all component and
  intersection bases, both graph directions, the rank-two correction and
  kernel, the field primality certificate, normal/optimized replay, and every
  pinned digest.  The audit requested the explicit kernel/surjectivity steps
  now inserted in Sections 2--3; no mathematical or scope failure remained.
depends_on:
  - THM-3585-lrc-common-a4-channel-plane-and-centering-complement
script: 04-computation/lrc_common_a4_anova_graph_flag_thm3593.py
output: 05-knowledge/results/lrc_common_a4_anova_graph_flag_thm3593.out
script_sha256: 3ee8486833d539599a4c8add304172b72929c69f3859e9cc1ce22e3018199516
output_sha256: 71327c65fc83628c117c2ed0bd1cbc5880f2d52277c73b6e135ffdea68ef550f
hash_basis: LF-normalized bytes
---

# THM-3593 -- LRC common A4 ANOVA graph flag

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This resolves the internal ANOVA geometry of THM-3585's common
four-plane.  It does not promote that static plane to a dynamical carrier.

## 1. Four ANOVA idempotents

Work over

```text
F=F_p,             p=755373809845391722745761,
H=Fun(V4(state) x F13(relation),F),        dim H=52.   (1)
```

Let `J_q/q` denote averaging on a `q`-point coordinate.  Define the four
pairwise complementary ANOVA projections

```text
P_0 =(J4/4) tensor (J13/13),
P_s =(I4-J4/4) tensor (J13/13),
P_r =(J4/4) tensor (I13-J13/13),
P_i =(I4-J4/4) tensor (I13-J13/13).                    (2)
```

They are idempotent, have pairwise disjoint images, and sum to the identity.
Their ambient image dimensions are `1,3,12,36`.

Let `A` be the common raw four-plane proved in THM-3585.  Put

```text
C0=P_0 A,             S=P_s A,
R=P_r A,              I=P_i A,
B=(P_0+P_s+P_r)A,     E=(P_0+P_s)A.                   (3)
```

Exact row-space elimination gives

```text
dim(A,C0,S,R,I,B,E)=(4,1,2,4,4,4,2).                 (4)
```

Thus the relation projection already has the full rank of `A`, while the
state/constant correction has rank only two.

## 2. First graph: raw plane versus additive and interaction planes

Every `a in A` decomposes uniquely as

```text
a=b+i,              b=(P_0+P_s+P_r)a in B,
                    i=P_i a in I.                     (5)
```

The exact stack ranks are

```text
dim(A+B)=8,          dim(A+I)=8,
dim(B+I)=8.                                             (6)
```

Therefore

```text
A intersect B=A intersect I=B intersect I=0.           (7)
```

Put `P_a=P_0+P_s+P_r`.  The two restricted projections have kernels

```text
ker(P_a|A)=A intersect I=0,       ker(P_i|A)=A intersect B=0.
```

They both have rank four by `(4)`, hence are isomorphisms.  Equivalently,

```text
T=(P_i|A) composed with (P_a|A)^(-1):B->I
```

is the unique isomorphism such that

```text
A={b+T(b):b in B},                    A+I=B direct_sum I. (8)
```

This explains THM-3585's raw/centred stack rank eight: it is the graph stack
of one four-dimensional map, not two unrelated four-planes.  The canonical
basis of `A+I` equals that of `B+I` entrywise and has digest

```text
9f4ec33d31337b0100a55871f6284443c6e5cfc0f5133a1493f807d563670821. (9)
```

## 3. Second graph: the additive plane is relation-led

The projection is onto because

```text
P_r(B)=P_r(P_a A)=P_r A=R.                            (10)
```

Since `dim R=dim B=4`,

```text
P_r|B:B->R
```

is an isomorphism.  Hence `B` is the graph of a unique linear map

```text
U=(P_0+P_s) composed with (P_r|B)^(-1):R->E,
im U=E,                                  rank U=2.       (11)
```

The rank in `(11)` is exactly `dim E=2`.  Its kernel is the part of the
additive plane that is already pure relation:

```text
K_2:=B intersect im(P_r)=ker U,             dim K_2=2. (12)
```

Thus the exact flag is

```text
K_2 subset R_4 subset im(P_r),          dimensions 2<4<12,
B=graph(U:R_4->E_2),
A=graph(T:B_4->I_4).                                    (13)
```

The state correction is genuinely tied to the constant coordinate.  The
separate images have dimensions `dim C0=1`, `dim S=2`, while their paired
image `E=(P_0+P_s)A` still has dimension two.  Projection `E->S` is an
isomorphism because `P_s(E)=P_s(A)=S` and both spaces have dimension two.
Thus `E` is itself the graph of one linear functional `S->C0`.

## 4. Hostile intersections and interpretation

Let `H_s^0`, `H_s`, `H_r^0`, `H_r` denote the standard pure-state,
state-only, pure-relation, and relation-only spaces, of dimensions
`3,4,12,13`.  The exact intersections are

```text
dim(B intersect H_s^0)=0,
dim(B intersect H_s)=0,
dim(B intersect H_r^0)=2,
dim(B intersect H_r)=2.                                (14)
```

So `B` is neither a pure-state chart nor a pure-relation subspace.  It is a
relation-led graph with exactly two corrected directions and two already
pure-relation directions.  This is the missing coordinate behind the coarse
phrase "rank-four amplitude quotient."

The relation four-plane `R` is still not identified with the pointed
six-dimensional relation carrier `P6`.  The next lawful comparison must map
the flag `(13)` into that typed carrier and preserve its two-dimensional
kernel; equality of dimensions or Fourier support cannot supply the map.

## 5. Exact verification and scope

The companion hash-pins and reconstructs the THM-3585 source tensor, applies
all four projections in `(2)`, and verifies idempotence component by
component.  It pins the eight component ranks and canonical basis digests,
the seven stack ranks, the six recorded standard-subspace intersections plus
`B intersect I` through the additive--interaction stack rank, and graph
equality `(9)`.  Its semantic digest is

```text
44d2d062a447eb68fc58b33eee8d23fb8092e1829534fb7c4f9e9f082e037f76. (15)
```

Reproduce with

```text
python -B 04-computation/lrc_common_a4_anova_graph_flag_thm3593.py
python -B -O 04-computation/lrc_common_a4_anova_graph_flag_thm3593.py
```

This theorem is exactly over the pinned finite field.  It proves no temporal
order, transfer law, current, arrival ancestry, row exclusion, physical
entry, characteristic-zero lift, or LRC(14) conclusion.

**End of proof.**
