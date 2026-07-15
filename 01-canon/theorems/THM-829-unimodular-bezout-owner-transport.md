---
id: THM-829
title: Unimodular continued-fraction matrices transport inverse owners by the contragredient Bezout-row action
status: PROVED (all-size GL2(Z) transport, normalization, reflection conjugacy, and centralizer) + FINITE-EXACT (q=13 and THM-808 tables)
source: codex-2026-07-15-S13-postjoin
depends_on: [THM-778, THM-808, THM-812, THM-813, THM-819]
related: [THM-825, HYP-6880]
verification:
  - 04-computation/unimodular_bezout_owner_transport_codex_S13_postjoin.py
  - 05-knowledge/results/unimodular_bezout_owner_transport_codex_S13_postjoin.out
  - 05-knowledge/results/unimodular_bezout_owner_transport_codex_S13_postjoin.json
---

# THM-829 — the inverse-owner stalk is a contragredient Bezout row

Let

```text
v=(a,q)^T,       0<a<q,       gcd(a,q)=1,               (1)
```

be a primitive witness centre.  Write

```text
b=(j,k),         j a+k q=1,                              (2)
```

where the canonical choice has `1<=j<q`; then `j=a^(-1) mod q` is THM-819's
right inverse owner and `q-j` is its reflected left owner.  For every
unimodular matrix `A in GL_2(Z)`, put

```text
v'=A v,          b'=b A^(-1).                            (3)
```

Then

```text
b'v'=1.                                                   (4)
```

After the elementary sign and numerator normalization which writes
`v'=(a',q')^T` with `q'>1` and `0<a'<q'`, the first coordinate of the
normalized row `b'` is exactly `a'^(-1) mod q'`.  Thus inverse-owner transport
is a linear contragredient action on the full Bezout row.  It is not a
function of the old inverse residue alone when the denominator changes.

Reflection does not generally commute with this action.  If

```text
R=[-1 1; 0 1],       R(a,q)^T=(q-a,q)^T,                 (5)
```

then the reflection-compatible action on the reflected branch is

```text
A^R=R A R.                                               (6)
```

Indeed `A^R R=R A`.  A single `A` commutes with reflection only for

```text
A in {I,-I,R,-R}.                                       (7)
```

Consequently the reflected high/low substitutions retained by THM-812/813
are structurally necessary: no nontrivial positive continued-fraction matrix
can act identically on both branches.

The result was first targeted at THM-824, then moved past live THM-824--828
reservations.  A live-main check reserved THM-829 before the proof artifacts
were added.

## 1. Proof of the transport law

Since `det A=+/-1`, the inverse `A^(-1)` is integral and `A` preserves
primitive columns.  Equations (2)--(3) give immediately

```text
b'v'=b A^(-1) A v=bv=1,
```

proving (4).  If the second coordinate of `v'` is negative, multiply both
`v'` and `b'` by `-1`.  If

```text
v'_raw=(a'+tq',q')^T,
```

then numerator reduction uses the shear

```text
v'_raw=[1 t;0 1] (a',q')^T.
```

Therefore the corresponding row becomes

```text
(j',k') -> (j',k'+t j').                                (8)
```

Finally reducing `j'` modulo `q'` and making the compensating change in `k'`
gives the canonical Bezout row.  Primitivity ensures that a normalized
numerator cannot be zero when `q'>1`.  Denominators `q'=0,+/-1` are the cusp
or unit-denominator boundary and do not carry a THM-819 interval witness.

This proves the normalized action for both determinant signs.

## 2. Reflection swaps the two owners and conjugates the action

Because `R^2=I`, the reflected source row is

```text
b R=(-j,j+k).                                            (9)
```

Its first coordinate modulo `q` is `q-j`, exactly swapping THM-819's inverse
owner pair.  For an arbitrary action,

```text
(R A R)(R v)=R A v,                                     (10)
```

so applying `RAR` on the reflected source gives the reflection of the target,
including the owner swap.

For completeness, write `A=[alpha beta;gamma delta]`.  The equality `AR=RA`
forces

```text
gamma=0,       delta=alpha+2 beta.                       (11)
```

Unimodularity then says `alpha(alpha+2 beta)=+/-1`, so both factors are
`+/-1`.  The four solutions are precisely (7).  This proves the centralizer
claim, rather than merely checking the displayed matrices.

## 3. The centered ratios from THM-812 and THM-813

Use the standard convergent-matrix convention

```text
C(d)=[d 1;1 0],       A_[d0;...;dr]=C(d0)...C(dr).       (12)
```

The reduced ratios behind the first two centered schedules give

```text
3/2=[1;2],       A_32=[3 1;2 1],
4/3=[1;3],       A_43=[4 1;3 1],
A_43 A_32=[14 5;11 4].                                  (13)
```

One must distinguish the Euclidean digit list from the centered increment
word: THM-813's increment word is `(1,2,1)`, while the Euclidean digits of
`4/3` are `(1,3)`.  Feeding the increment letters into (12) would be a
category error.

The reflection-conjugate matrices are

```text
R A_32 R       =[ 1 -1;-2  3],
R A_43 R       =[ 1 -1;-3  4],
R A_43 A_32 R  =[ 3 -4;-11 15].                          (14)
```

None belongs to (7), so none acts by the same matrix on the two reflection
branches.  This arithmetic fact is the primitive-witness counterpart of
THM-813's geometric rule that reflection-orbit lines, rather than projected
node-pair cells, are the safe action carrier.

THM-812/813 supply the centered ratios and their tiling coordinate-copy
actions; they do not yet identify a tournament tiling with a primitive
witness column.  Equations (12)--(14) are the standard arithmetic actions to
place in that future fibre product, not a claim that the existing coordinate
copies already induce them on LRC witnesses.

## 4. Exact `q=13` witness table

Apply (13) to all twelve primitive witnesses `(a,13)`.  Each table entry is
`q',j'_+`; the other owner is `j'_-=q'-j'_+`.

```text
a  old j+     A_32       A_43       A_43 A_32
1      1      15, 1      16, 1       63, 4
2      7      17, 9      19,10       74,39
3      9      19,13      22,15       85,58
4     10      21,16      25,19       96,73
5      8      23,14      28,17      107,65
6     11      25,21      31,26      118,99
7      2      27, 4      34, 5      129,19
8      5      29,11      37,14      140,53
9      3      31, 7      40, 9      151,34
10     4      33,10      43,13      162,49
11     6      35,16      46,21      173,79
12    12      37,34      49,45      184,169.             (15)
```

For `A_32` and `A_43`, numerator normalization leaves `a'=a` and changes the
denominator to `2a+13` and `3a+13`, respectively.  For the composite it gives
`a'=3a+13`, `q'=11a+52`.  The inverse owner changes nontrivially even in the
first two columns; retaining only the old `j` without `(a,q)` cannot transport
it across moduli.

The verifier checks 831 exact action records from all primitive columns
`2<=q<=30` under the three matrices.  There are zero Bezout, normalization,
reflection, or owner-swap failures.  On this bank `(old inverse,matrix)` has
5,577 target-conflict pairs and even `(old Bezout row,matrix)` has 2,160:
the row is a stalk *over the primitive column*, not a replacement for it.
`(a,q,matrix)` and the full column-row stalk have zero conflicts.

## 5. Compatibility boundary with THM-808

THM-808 works in one fixed prime sheet `F_p`; its owner step is
`w_a^(-1) mod p`.  The fixed-denominator shear

```text
T_m=[1 m;0 1],       (a,p)->(a+mp,p)                    (16)
```

has

```text
(j,k)T_m^(-1)=(j,k-mj),                                 (17)
```

so it leaves the inverse owner `j` unchanged.  This proves that THM-808's
root translation

```text
d' = d-sum_a c_a j_a mod p                              (18)
```

is invariant under changing integer representatives of the owner speeds.
Using the exact owner-count blocks from its ten-wall movie, the verifier
reproduces the five translations

```text
4,3,3,4,3 mod 7                                         (19)
```

from the Bezout first coordinates, with zero failures under shears
`-3<=m<=3`.

A general matrix in (3) changes the denominator and therefore leaves the
common `F_p` sheet on which the sum (18) lives.  It does not define a new root
translation until all owners have again been placed over one specified
modulus.  This is the precise boundary between THM-808's affine root cocycle
and the present cross-modulus witness action.

## 6. Metric witness intervals and preservation boundary

At denominator `q`, THM-819's primitive interval has the labelled arithmetic
stalk

```text
(a,q; j_+,j_-; delta_q/j_+, delta_q/j_-),
delta_q=1/(q(q+1)).                                      (20)
```

Equations (3), (8), and (15) transport `(a,q;j_+,j_-)` exactly.  At the new
denominator, (20) then canonically recomputes the two new extents using
`delta_(q')`.  This is a rebase to the THM-819 witness chart at `q'`; it is
not the geometric image of the old interval under the Möbius map, and it does
not preserve the total measure or a fixed LRC speed set.  In particular the
nontrivial matrices in (13) leave the `q=13` chart relevant to LRC(14).

Tournament Analysis uses `inverse only`, `inverse+matrix`, `Bezout+matrix`,
`column+matrix`, and the full stalk as arithmetic-carrier vertices.  Its
pairwise observable is the number of target-distinct finite-bank records
separated.  Retention and retention per partition cell are the switches and
flip nine of ten pair orientations.  The exact carriers are the primitive
column plus matrix, or its full Bezout-row lift; runners and bare tournament
nodes are not the relevant vertices.

The theorem preserves primitivity, orientation after normalization, both
inverse owners, the full Bezout identity, the CF matrix, and reflection after
conjugating the branch action.  It destroys the original denominator chart,
metric interval image, owner-block chronology, prime-sheet token assignment,
metagraph line/node identity, LRC loneliness predicate, and continuation.
Coupling to THM-813 must therefore be a fibre product with its literal
reflection-orbit carrier `Q`, and coupling to THM-808 must retain the common
sheet modulus and owner-count block. ∎
