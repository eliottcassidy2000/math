---
id: THM-2887
title: "Quaternionic Arf lift of the semantic V4 and global carry no-go"
status: >
  PROVED + VERIFIED-EXACT.  The THM-2884 semantic V4 carries the
  anisotropic Arf-one quadratic form q(h)=1_(h nonzero), with determinant
  polarization.  Its unique central-extension class having determinant
  commutator and nontrivial square on all three nonzero directions is
  Q8.  The marked source QB acts by inner conjugation with
  sign det(QB,h), recovering exactly the THM-2886 q3/q11/q7
  selector-pair XOR parity (0,1,1); it does not supply the origin
  occupancy or signed orientation.  Thus the source mark, not the
  unpointed projective line, supplies this parity coordinate.  The local
  polarization also equals the
  base-thirteen carry on the 8+9 triangle.  This cannot globalize to a
  base-independent C13 increment label: carry mod two has the unique
  normalized primitive h mod2, forcing every even increment to the zero
  semantic direction and contradicting the required labels at 8 and 4.
  Aut(Q8)=S4, Inn(Q8)=V4 and Out(Q8)=S3; every semidirect action
  C13 or C169 -> Aut(Q8) is trivial, so that particular joint candidate
  is the direct product C169 x Q8 of size 1352.  This does not classify
  arbitrary extensions or groups of order 1352.  No physical Q8 sheet
  lift, canonical Q8-valued current, row exclusion, or LRC(14) proof
  follows.
source: root/lrc-quaternionic-arf-lift-2026-07-29
depends_on:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2884-macro-semantic-diagonal-horn-carrier-and-origin-even-boundary
  - THM-2886-stepped-origin-provenance-transport-on-the-v4-horn
related:
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2814-projective-allocation-square-holonomy-and-idempotent-provenance-no-go
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2878-endpoint-factor-exit-carry-transducer-and-flat-vertex-boundary
  - THM-2882-event-twisted-all-q-coefficient-carry-lift
script: 04-computation/lrc14_quaternionic_arf_semantic_v4_carry_no_go_thm2887.py
output: 05-knowledge/results/lrc14_quaternionic_arf_semantic_v4_carry_no_go_thm2887.out
script_sha256: 4fa39842cf0a7f407bd7d436315ac85e132161933c10d2c078644d91dafd91b8
output_sha256: b2edbfdb4964573810a7147e8676dd152ba60aafe7953dde7afd3a2a65da3d5a
hash_basis: LF-normalized bytes
---

# THM-2887 -- quaternionic Arf lift of the semantic V4 and global carry no-go

**PROVED + VERIFIED-EXACT.**

The diagonal carrier of THM-2884 is not merely a four-element support
set.  Its nonzero-support parity is the anisotropic quadratic form on a
symplectic plane over `F2`, so the associated central sign has a forced
quaternionic model.  Keeping THM-2886's distinguished source then
resolves a seeming symmetry obstruction: the horn's selector-pair XOR
parity is the inner-conjugation character of that source.

This section proves both the positive local structure and its sharp
global boundary.  The quaternionic lift explains the transported
`(0,1,1)` parity, but a single base-thirteen increment labelling cannot
realize it on all carry triangles.

## 1. The semantic plane and its quadratic form

Write the THM-2884 carrier as

```text
H = {Q0, QA, QB, QAB} ~= F2^2,
QA=(1,0), QB=(0,1), QAB=(1,1).                         (1)
```

The corresponding physical triples are

```text
Q0 =(0,0,0),  QB =(1,0,1),
QA =(1,1,0),  QAB=(0,1,1).                            (2)
```

Define

```text
q_H(u,v) = u + v + uv = 1_((u,v) != (0,0)).            (3)
```

For all `x,y in H`,

```text
q_H(x+y)+q_H(x)+q_H(y) = det(x,y).                     (4)
```

Thus `q_H` is the anisotropic, Arf-one quadratic refinement of the
determinant pairing.  On `(Q0,QA,QB,QAB)` its values are `(0,1,1,1)`.
By contrast, the THM-2779 refinement `q_D(u,v)=uv` has values
`(0,0,0,1)` and Arf invariant zero.  The two forms have the same
polarization but different square laws; this is the exact distinction
between the quaternionic and dihedral boundaries.

## 2. The forced quaternionic extension

For `x=(x1,x2)` and `y=(y1,y2)`, put

```text
alpha(x,y) = x1 y1 + x2 y2 + x1 y2                    (5)
```

and define a central extension on `H x F2` by

```text
(x,e)(y,f) = (x+y, e+f+alpha(x,y)).                    (6)
```

The cocycle identity is exact, and

```text
(x,e)^2 = (0,q_H(x)).                                  (7)
```

Consequently the group has order census

```text
1^1, 2^1, 4^6,                                         (8)
```

so it is `Q8`.  In the section

```text
i=(QA,0), j=(QB,0), k=(QAB,0),
```

one has

```text
i j = -k,                 j i = k.                     (9)
```

Replacing the final term `x1 y2` in (5) by `x2 y1`
changes the section by the coboundary of `q_H` and reverses which ordered
product pays the central sign.  The commutator and square law are
unchanged.

This lift is unique with the stated labelled properties.  Exhaustion of
all normalized functions `H x H -> F2` gives:

```text
normalized 2-cocycles                                      16
with determinant commutator                                 8
quadratic refinements                         4, twice each
normalized coboundaries                                     2. (10)
```

The eight determinant-commutator cocycle tables split into four
cohomology classes, one for each quadratic refinement.  There are eight
normalized one-cochains but only two distinct normalized coboundaries:
the four linear characters form the kernel.  Thus each refinement has
two normalized cocycle-table representatives related by the unique
nonzero coboundary.  Three classes are Arf-zero labelled `D8` forms; the
only class whose three nonzero directions all square to the nonidentity
central element is the Arf-one `Q8` class above.  The uniqueness here is
among normalized central `F2`-extensions with determinant commutator and
this all-nonzero square law; it is not a claim that the LRC sheet has
already been lifted.

## 3. The carry triangle and the marked-source character

The two successive horn edges and their composite have semantic
directions

```text
q3 -> q11 : QA,             q11 -> q7 : QB,
q3 -> q7  : QAB.                                      (11)
```

Their base-thirteen increments are `8`, `9`, and `4 mod 13`.  Therefore

```text
det(QA,QB)
  = q_H(QA)+q_H(QB)+q_H(QAB)
  = 1
  = floor((8+9)/13).                                  (12)
```

The equality is gauge invariant: the local carry event is exactly the
polarization of the two semantic edge directions.

More importantly, retain THM-2886's marked source `QB`.  Conjugation by
its lift `j` acts on every section lift `s(h)` by

```text
j s(h) j^(-1) = (-1)^det(QB,h) s(h).                  (13)
```

It fixes `Q0` and `QB`, and negates `QA` and `QAB`.  On the transported
target path

```text
q3 : Q0,                  q11 : QA,                  q7 : QAB
```

the sign character is exactly

```text
(det(QB,Q0), det(QB,QA), det(QB,QAB)) = (0,1,1),      (14)
```

which is the selector-pair XOR parity independently proved in THM-2886.
This identifies the algebraic source of that one coordinate: it is a
**pointed** symplectic character, not an arbitrary orientation on the
three nonzero elements.  The conjugation sign is independent of the
choice of central lift of `QB` or `h`.  It does not reconstruct
THM-2886's empty/full origin occupancy `(0,1)` or its signed orientation
`-1`; those remain separate transported data.

The full `S3=GL(2,F2)` action preserves
`det(source,target)` when both entries move.  After `QB` is marked, its
stabilizer has order two and preserves (13).  If the source is forgotten,
however, exhaustive equivariance gives no map

```text
P1(F2) -> F2
```

intertwining the `S3` action with its sign character; for the trivial
action only the two constant maps remain.  Thus the unpointed projective
line cannot choose this seam parity coordinate, while the already
physical fixed source determines that coordinate canonically.

## 4. The global base-thirteen no-go

Let

```text
epsilon(h,k)=floor((h+k)/13),       0 <= h,k <= 12,     (15)
```

with `h+k` reduced modulo `13` in coboundary expressions.  If
`r(0)=0` and

```text
r(h)+r(k)+r(h+k mod 13)=epsilon(h,k) mod 2             (16)
```

for every pair, then exhaustive solution, equivalently induction around
the cyclic order, gives the unique primitive

```text
r(h)=h mod 2.                                          (17)
```

Suppose a base-independent semantic increment map

```text
ell : {0,...,12} -> H,       ell(0)=Q0
```

globalized the local quadratic/carry match, so that

```text
delta(q_H o ell)=epsilon mod 2.                         (18)
```

Equations (16)--(17) force `q_H(ell(h))=h mod 2`.
Because `q_H` is anisotropic, every even `h` must satisfy `ell(h)=Q0`.
But the proved horn requires

```text
ell(8)=QA,                 ell(4)=QAB,                  (19)
```

both nonzero.  This is impossible.  The one distinguished triangle in
(12) is genuine, but it cannot be promoted to a translation-invariant
labelling of all `C13` carry pairs.

## 5. Automorphisms and the joint-state invoice

Exact enumeration of the `Q8` multiplication table gives

```text
|Aut(Q8)|=24,          element-order census 1^1 2^9 3^8 4^6,
Inn(Q8)=V4,            Out(Q8)=S3.                     (20)
```

Equivalently, automorphisms are the orientation-preserving signed
permutations of the three imaginary axes; their faithful action on the
four body diagonals identifies

```text
Aut(Q8) ~= S4,          S4/V4 ~= S3.                   (21)
```

Every nontrivial subgroup of `Q8` contains its center.  Hence a faithful
permutation action needs a free orbit and has minimum degree eight.

Since `gcd(13,24)=gcd(169,24)=1`, every homomorphism

```text
C13 -> Aut(Q8),             C169 -> Aut(Q8)            (22)
```

is trivial.  Therefore the semidirect-product construction in which the
THM-2874 `C169` address acts on the present quaternionic sign cannot be
nontrivial.  That specific group-level candidate is

```text
C169 x Q8,                   size 169*8 = 1352.         (23)
```

This is an invoice for that structured semidirect candidate, not a
classification of arbitrary extensions or all groups of order `1352`,
and not a construction of the object.

## 6. Exact verification and boundary

The companion pins the exact THM-2779, THM-2884, and THM-2886 inputs.  It
checks both extension laws on all `8^3=512` triples, all `16` normalized
central cocycles, and separately all `8` normalized one-cochains yielding
`2` distinct coboundary tables.  It also checks all `24` automorphisms,
all `96` pointed determinant covariance cells, all `13^2` carry pairs,
and all `2^12` normalized carry primitives.  Normal, `python -O`, and
stored outputs agree byte for byte.

What is now proved is:

1. the THM-2884 parity has a forced Arf-one/quaternionic interpretation;
2. the retained `QB` source recovers exactly the THM-2886 selector-pair
   XOR parity by inner conjugation, but not its origin occupancy or
   signed orientation; and
3. a normalized base-independent `ell:C13->V4` satisfying
   `delta(q_H o ell)=epsilon` and the horn values is impossible, as is a
   nontrivial semidirect action `C169->Aut(Q8)`.

What remains open is the decisive physical step: construct, or exclude,
a sheetwise `Q8` lift and a canonical two-channel or `Q8`-valued
observable compatible with the Prony transport.  Nothing here converts
the marked raw current into the unmarked scalar current, excludes any of
the `165` residual rows, or proves LRC(14).
