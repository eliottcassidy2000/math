---
id: THM-2868
title: "Two-origin endpoint projector and projective Kummer carry descent"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the THM-2847 common 42-cell horn, a two-origin endpoint difference
  and a variable-offset two-sample Prony projector give a physical
  frequency-dual C13 atlas from 26 actual 91-unit multiplier samples.
  Its split-left branch is literally character 3, its split-right branch
  is trivial, and the resulting multiplier-by-target table has unique
  diagonal-invariant channel (3,10).  The projective branch ratio is a
  primitive degree-thirteen Kummer coordinate over Q(zeta_91).  This is a
  signed coefficient atlas, not one positive current or a physical
  ancestry co-shift; common E3 support and same-triangle transport remain
  open.
source: root/lrc-holotopy-allocation-2026-07-28
depends_on:
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
related:
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2863-endpoint-prony-splitter-and-carry-character-three-intertwiner
  - THM-2870-prime-power-convolution-versus-physical-diagonal-intertwiner-obstruction
script: 04-computation/lrc14_two_origin_endpoint_projective_kummer_thm2868.py
output: 05-knowledge/results/lrc14_two_origin_endpoint_projective_kummer_thm2868.out
script_sha256: 3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5
output_sha256: ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9
hash_basis: LF-normalized bytes
---

# THM-2868 -- Two-origin endpoint projector and projective Kummer carry descent

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

## 1. Statement

Let `P` be the literal THM-2806 left endpoint coefficient on the retained
unit sheet at the distinguished source origin.  On every cell of the
THM-2847 `q3/q11` common horn, all six remaining factor predicates are
one on the same source and target intervals.  This horn has exactly 42
cells, and `P` is nonzero in both certified endpoint fields.

At multiplier `m`, put

```text
Y(m)=12+26m
```

and evaluate the right endpoint coefficient at physical frequency
`-Y(m)`.  Let `R_a^m(q)` denote this coefficient at target origin `a`.
For every multiplier used below, exact replay gives

```text
R_00^m(3)=R_00^m(11)=R_12^m(11)=C_m != 0,
R_12^m(3)=0,
R_00^m(7)=R_12^m(7)=0.                              (1)
```

Thus the signed full-current coefficient

```text
S_m=P [R_00^m(3)-R_12^m(3)]=P C_m                  (2)
```

is nonzero at `q3`; the analogous `q11` origin difference cancels and the
`q7` difference is zero.

Define formal frequency sections

```text
m_r=1+42r,                 0<=r<13,
e_r=1 if r=4,  2 if r=8,  0 otherwise,
n_r=m_r+e_r.                                           (3)
```

Every one of the 26 raw multipliers `n_r,n_r+1`, and every corresponding
right frequency `Y(n_r),Y(n_r+1)`, is coprime to `91`.  From the two
actual measurements `S_(n_r),S_(n_r+1)`, the endpoint projector followed
by the known node transport reconstructs nonzero coefficients `U_r,V_r`
at the desired formal section `m_r`.  They satisfy cyclically

```text
U_(r+1)=omega^3 U_r,       V_(r+1)=V_r.               (4)
```

Consequently the frequency-dual table

```text
H(r,q)=U_r delta_3(q)                                  (5)
```

has normalized Fourier support

```text
supp Hhat={3} x F_13,
Hhat(3,b)=13^(-1) U_0 omega^(-3b).                     (6)
```

The unique channel invariant under the formal diagonal translation of
this finite multiplier-by-target table is `(3,10)`.

The projective ratio `t_r=U_r/V_r` cancels `P`.  It is a primitive
degree-thirteen Kummer coordinate over `F=Q(zeta_91)` and obeys

```text
t_(r+1)=omega^3 t_r.                                   (7)
```

Equations (4)--(7) are a physical **frequency-dual** `C13/chi_3` fibre:
the index `r` labels actual multiplier-frequency measurements plus their
explicit local Prony transport.  It is not a physical ancestry action,
one raw affine pair, or one positive current.

## 2. Two-node endpoint projector

Let `[L,R)` be the unique live endpoint interval.  In the relevant
cyclotomic factor, with primitive root `xi`,

```text
C_m=A_L lambda_L^m-A_R lambda_R^m,
A_L=xi^(12 R_dil L),        lambda_L=xi^(26 R_dil L),
A_R=xi^(12 R_dil R),        lambda_R=xi^(26 R_dil R). (8)
```

The nodes are distinct.  Therefore `S_m=P C_m` obeys the same order-two
recurrence as `C_m`.  The four-window `m=1,2,3,4` proves and checks the
recurrence, and two adjacent values split the oriented branches:

```text
U_m=(S_(m+1)-lambda_R S_m)/(lambda_L-lambda_R)
    =P A_L lambda_L^m,

V_m=(lambda_L S_m-S_(m+1))/(lambda_L-lambda_R)
    =-P A_R lambda_R^m.                                (9)
```

This split keeps the full common source coefficient; it does not divide
by `P` or discard the rest of the current.

The naive adjacent windows at the desired `m_r` are not all lawful:

```text
m_4=169=0 mod 13,
m_8+1=338=0 mod 13.                                   (10)
```

These are the only two index-zero sections.  The offsets in (3) repair
them.  Apply (9) at the actual unit pair `n_r,n_r+1`, obtaining
`U_(n_r),V_(n_r)`, and set

```text
U_r=lambda_L^(-e_r) U_(n_r),
V_r=lambda_R^(-e_r) V_(n_r).                          (11)
```

Then (9) gives identically

```text
U_r=P A_L lambda_L^(m_r),
V_r=-P A_R lambda_R^(m_r).                            (12)
```

This is why a variable local offset loses no coefficient information.
The mod-7 residues of the raw pairs are `(1,2)` except `(2,3)` at `r=4`
and `(3,4)` at `r=8`, so the deepest coefficient clearing is valid as
well.  The atlas consists of 26 unit samples; the exact companion also
replays `m=3,4` for the independent four-term recurrence, hence 28 total
sample entries.

## 3. Literal character-three frequency cycle

The exact endpoint nodes satisfy

```text
lambda_L^42=xi^546=omega^3,
lambda_R^42=xi^7098=1.                                (13)
```

Since `m_(r+1)-m_r=42`, equations (12) and (13) prove (4), including the
cyclic wrap from `r=12` to `r=0`.  In particular all 13 values are
nonzero.  Applying the normalized two-axis transform

```text
Hhat(a,b)=13^(-2) sum_(r,q) H(r,q)omega^(-ar-bq)       (14)
```

to (5) gives (6).  Under the formal diagonal action
`(r,q)->(r+t,q+t)`, a Fourier channel gains
`omega^(-(a+b)t)`, so precisely `(a,b)=(3,10)` among the support in (6)
is invariant.

This is an invariant of the finite **multiplier-section by target-label
table**.  Calling it a physical ancestry co-shift would be a type error.
At the representation level its one-dimensional `chi_3` summand matches
the character in THM-2851 and the abstract intertwiner of THM-2863; the
normalized map between these one-dimensional representations is unique.
No common-support ancestry action or `E3` contraction is supplied here.
This does not evade THM-2870's same-mask diagonal/convolution obstruction:
the construction is an explicit frequency-indexed Prony line with local
scalar transport, not an invertible intertwiner between a Boolean physical
mask and its convolution kernel.  The missing ancestry comparison is
precisely the coefficient/character reference isolated there.

There is a suggestive but deliberately non-consequential seam alignment.
The canonical section and the two repaired sections are

```text
{0,4,8} = {q3-q3, q7-q3, q11-q3},
q7-q11=9 mod 13.                                      (15)
```

Numerically, (15) comes from `m_4=169` and `m_8+1=338=2*169`.  Thus the
frequency charts form the same affine triangle of labels as the
`q3/q7/q11` sidecar in THM-2851.  This is only a typed Cech-seam
coincidence: (1) still says that the signed `q11` row cancels and the
`q7` row is zero.

## 4. Projective Kummer coordinate and the non-rational scale

The common scale `P` is not `F=Q(zeta_91)`-rational.  The cyclotomic
automorphism with exponent `u=547` fixes `zeta_91`, since `547=1 mod 91`,
but exact evaluation gives

```text
field 352341050142921841:
  P=254455016269350867 -> 248769464494275593,

field 956354278959359281:
  P=318932490657369324 -> 646245011020437966.          (16)
```

Thus a pure-character claim for the full branch under the full-field
Galois action would be false.  The frequency reindexing in Section 3 is
different: it keeps the same physical source coefficient `P` fixed while
varying only the endpoint multiplier.

Projectivization removes this scale:

```text
t_m=U_m/V_m=-L_m/R_m=xi^(1111-156m).                  (17)
```

For `m=1,2,3,4`, the exponents are

```text
955, 799, 643, 487,
```

and their thirteenth-power exponents are

```text
585, 923, 1261, 1599.
```

For the 13 formal sections `m_r`, equation (17) becomes

```text
t_r=xi^(955+546r),
exponents mod 2366:
955,1501,2047,227,773,1319,1865,45,591,1137,1683,2229,409.  (18)
```

They are distinct and each is coprime to `2366`.  Hence every `t_r`
generates `K=Q(xi)`.  Since `[K:F]=13`, `t_r^13` belongs to `F`, and
`t_r` does not, its minimal polynomial over `F` is

```text
X^13-t_r^13.                                          (19)
```

Equations (18) and (13) give the rotation (7).

## 5. Representative gauge and exact evidence

Under the marked-sheet representative change

```text
(ell,h_source,h_target) -> (ell+W,h_source+1,h_target-1),
```

the source and target factors in (2) are unchanged for every entry of
the 28-sample replay bank in both certified fields.  The corresponding
phase increments

```text
-12 R_dil (T/13),       (12+26m) R_dil (T/13)
```

vanish modulo the endpoint conductor.

Reproduce the exact check from the repository root:

```text
python 04-computation/lrc14_two_origin_endpoint_projective_kummer_thm2868.py
python -O 04-computation/lrc14_two_origin_endpoint_projective_kummer_thm2868.py
```

Both runs equal the stored 16-line transcript.  The companion checks:

- the exact 42 common horn cells and the full 28-entry endpoint bank;
- all 26 raw multiplier and right-frequency unit conditions;
- the support/value pattern (1), representative gauge, and four-term
  Prony recurrence in both certified fields;
- all 13 transported split pairs, their cyclic `chi_3`/trivial laws,
  the 169 multiplier-target Fourier coefficients, and the unique
  `(3,10)` invariant channel;
- all projective exponents, Kummer conditions, and the `u=547`
  non-rationality witnesses.

LF-normalized SHA-256:

```text
script 3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5
output ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9
```

## 6. Exact boundary and connection contract

The first remaining failure is semantic/common-support transport, not
multiplier support, representative gauge, or omitted endpoint
factorization.  The signed selector works because the `q3` interval is
present at origin `00` and `E3`-absent at origin `12`; `q11` is present
at both origins and cancels; `q7` is `E3`-absent at both.  Therefore this
theorem does not transport its nonzero split branch to the physical
THM-2851 `q11-to-q7` ancestry carry.

```text
source:
  THM-2847 two-origin 42-cell horn and the literal endpoint coefficients;

map:
  take the q3 origin difference, split each 91-unit adjacent pair by the
  two-node projector, and transport the two repaired local windows back
  to the formal sections m_r=1+42r;

preserved:
  all 42 horn cells, the full common source coefficient P, endpoint
  orientation, marked-sheet gauge, q3 target label, all 13 frequency
  sections, the literal chi_3/trivial branch cycle, the unique finite-table
  invariant channel (3,10), and the projective Kummer coordinate;

destroyed / missing:
  positivity, realization as one raw affine pair or one physical current,
  one common marked triangle across the variable local offsets, q11 mass,
  q7 E3 support, a common ancestry action, E3/complement contraction,
  and row closure;

cheapest next test:
  construct a lawful common-support E3 mapping cone that intertwines the
  frequency chi_3 atlas with the physical q3/q7/q11 ancestry triangle,
  or prove that any such intertwiner requires an additional positive
  co-support/reference channel.
```
