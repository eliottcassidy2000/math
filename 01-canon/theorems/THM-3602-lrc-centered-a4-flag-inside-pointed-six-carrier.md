---
id: THM-3602
title: "LRC centered A4 flag inside the pointed six carrier"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
  AUDIT.  Over the pinned field, state averaging transports the THM-3593
  relation flag K2<R4 into Fun_0(F13), where it lies inside the centered
  pointed six-carrier Q6 with dimensions 2<4<6<12.  The raw carrier P6 is the
  graph of one nonzero augmentation functional over Q6.  Its intersections
  with R4 and K2 have dimensions three and two.  Among all 156 affine
  reindexings of F13, only the identity retains either intersection.  This is
  a static finite-field target-space theorem only: no difference-fibre lift,
  chronology, current, characteristic-zero statement, or LRC(14) conclusion
  is claimed.
source: kps-s188 + agent Godel / THM-3593 pointed-carrier continuation, 2026-08-21
audit: >
  Author/implementer exact audit.  The standard-library companion hash-pins
  and independently reconstructs the THM-3593 relation flag and the clean-room
  pointed carrier, verifies every rank, intersection, augmentation kernel,
  parent digest, and all 156 affine target reindexings.  Normal and optimized
  runs are byte-identical.  Independent hostile audit remains pending.
depends_on:
  - THM-3593-lrc-common-a4-anova-graph-flag
related:
  - THM-3585-lrc-common-a4-channel-plane-and-centering-complement
script: 04-computation/lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.py
output: 05-knowledge/results/lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.out
script_sha256: c6954206f6a3187a98772fdf2313f0fd91d7df12f38ec7438d39be7bad7235e8
output_sha256: bdb1a7cb2b85c1dce1ad1dfce54c9e02b4c935ff810101768fd08030d6531d56
semantic_sha256: fdce64910b2eb0c200ef9761cd6518f65c6b3295083c161d00c76a5c8c83d095
affine_atlas_sha256: f074cff98c370dbebfefc998f96ef761641ba28657a21c1b5c20f95b790678a6
hash_basis: LF-normalized bytes for files; canonical JSON for semantic ledgers
---

# THM-3602 -- LRC centered A4 flag inside the pointed six carrier

**PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
AUDIT.**  This supplies the typed comparison requested by THM-3593.  It also
isolates exactly what that comparison forgets: the difference coordinate and
all temporal data.

## 1. Shared field and typed relation map

Work over

```text
F=F_p,       p=755373809845391722745761,
H=Fun(V4(state) x F13(relation),F),
Z=Fun_0(F13,F).                                           (1)
```

The pinned primitive thirteenth root in both constructions is

```text
zeta_13=266737884585332483769981.                         (2)
```

THM-3593 places its pure-relation spaces in `im(P_r) subset H`.  Every
`h in im(P_r)` is independent of the state coordinate and has zero relation
sum.  Consequently

```text
T:im(P_r)->Z,
T(h)(t)=1/4 sum_(s in V4) h(s,t)                          (3)
```

is a typed isomorphism, with inverse given by repeating the relation profile
on all four state slices.  Transport the THM-3593 flag through `(3)` and retain
the names

```text
K2 subset R4 subset Z,                    dim(K2,R4,Z)=(2,4,12).   (4)
```

This state average is not a heuristic marginal: on `im(P_r)` all four slices
are exactly equal, and the companion verifies the equality before transport.

## 2. Raw and centered pointed carriers

Independently reconstruct the clean-room pointed tensor at the six ordered
states

```text
(0,0), (1,0), (1,6), (3,6), (3,12), (2,12).             (5)
```

For each pointed state, sum over the difference coordinate and retain its
thirteen-entry relation profile.  Let `P6 subset Fun(F13,F)` be the span of
these six raw profiles.  Let

```text
C13(f)=f-(1/13 sum_t f(t)) 1,
Q6=C13(P6) subset Z.                                    (6)
```

Exact row reduction gives

```text
dim P6=6,               dim Q6=6,               ker(C13|P6)=0.    (7)
```

Thus centering is an isomorphism `P6->Q6`; it is not an inclusion of the raw
carrier in `Z`.

## 3. The centered flag and the raw defect

The complete centered containment is

```text
K2 subset R4 subset Q6 subset Z,
dimensions       2 < 4 < 6 < 12.                        (8)
```

The canonical row-space digests are

```text
K2  193c82c6337f88c9b6c1bb2464198808336f42cc7036f5c3ae30af449e147357
R4  28bca0c67b0d94431fe3c43181a46d7b55e6cfce98fdc232f837acfde0ccec8c
P6  6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4
Q6  ae86ad2a3fd03bea95c823b2454b78f244581aa048cb2da63a03a75f484cc596. (9)
```

The raw carrier is not interchangeable with its centered image.  Since
`C13|P6` is an isomorphism, there is a unique linear functional

```text
lambda:Q6->F
```

such that

```text
P6={q+lambda(q)1:q in Q6}.                               (10)
```

The functional has rank one, so

```text
dim(P6 intersect Q6)=dim ker(lambda)=5.                  (11)
```

Its restrictions to the nested flag have ranks

```text
rank(lambda|R4)=1,             rank(lambda|K2)=0.        (12)
```

Therefore the sharper raw/centered flag is

```text
K2 subset (R4 intersect P6) subset R4 subset Q6 subset Z,
dimensions       2       <       3       < 4  < 6 < 12. (13)
```

Equivalently, centering adds one relation direction to the part of `R4`
already visible in the raw pointed carrier, while the kernel `K2` is completely
insensitive to the augmentation sidecar.

The remaining stack and projection ranks are

```text
dim(R4+P6)=7,                    dim(R4+Q6)=6,
dim(K2+P6)=dim(K2+Q6)=6,
rank(R4 -> Fun(F13)/P6)=1,
rank(R4 -> Z/Q6)=rank(K2 -> Z/Q6)=0.                    (14)
```

These identities distinguish equality after centering from equality in the
raw target.

## 4. Affine-coordinate hostile test

For every affine permutation

```text
g_(a,b):t |-> a t+b,                a in F13^x, b in F13,          (15)
```

reindex the target coordinates of `Q6`.  There are `12*13=156` choices.  The
exact intersection histogram is

```text
(dim(R4 intersect gQ6), dim(K2 intersect gQ6))
       (4,2): 1 occurrence,
       (0,0): 155 occurrences.                             (16)
```

The sole nonzero case is `(a,b)=(1,0)`.  Thus the containment `(8)` depends on
the exact shared relation order and shared root `(2)`; it is not explained by
dimension, Fourier support, or an unrecorded affine gauge choice.

This is a hostile control on coordinate alignment, not a uniqueness theorem
among arbitrary permutations or linear automorphisms of `Z`.

## 5. Verification and boundary

The companion hash-pins two disjoint parents: the THM-3593 exact source and the
clean-room pointed-six reconstruction.  It rebuilds both spaces rather than
importing stored row bases, verifies primality and shared root data through the
parents, and pins every canonical intersection basis.  Reproduce with

```bash
python -B 04-computation/lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.py
python -B -O 04-computation/lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.py
```

This theorem lives entirely in the static relation target over the one pinned
finite field.  The pointed construction has already summed out the difference
coordinate.  Nothing here supplies a lift back to difference fibres, a legal
time order, current conservation, source ancestry, an integer or
characteristic-zero transfer, a forbidden LRC row, or an LRC(14) conclusion.
The next lawful bridge must restore at least one of those missing coordinates;
the flag alone cannot do so.

**QED.**
