---
id: THM-3076
title: "Finite affine-plane line-quotient tomography and the p=2 three-view law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over any field in
  which p is invertible, the p+1 line-
  quotient averages on F_p^2 split into one constant channel and p+1
  pairwise-annihilating centered direction channels.  Every nonempty labelled
  s-view bank has rank 1+s(p-1), exact common-mean image, canonical observable
  projector, and sharp omitted-line hostiles; all p+1 views reconstruct.  In
  characteristic p the cleared line norms are square-zero and the full deck
  has rank p(p+1)/2 instead.  At p=2 this is exactly the fibrewise 2/3/4 rank
  law behind THM-3072, and reduction of PSL2(Z)=C2*C3 modulo three identifies
  its generators with translation and the ternary A4 chart cycle after one
  relabelling.  No physical quartic, Farey, Keller, or LRC intertwiner follows.
source: root-finite-plane-tomography-2026-08-01
audit: >
  An immutable hostile audit rederived the projector algebra, every good- and
  bad-characteristic rank, the sharp missing-line witnesses, the joint-label
  versus separate-table distinction, and the generatorwise PSL2(F3)=A4
  relabelling.  It independently reproduced the normal transcript and LF
  hashes; the author normal and optimized replays both match stored output.
depends_on: []
related:
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
  - THM-850-chi7-face-charge-cyclotomic-law
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss
script: 04-computation/finite_affine_plane_line_quotient_tomography_thm3076.py
output: 05-knowledge/results/finite_affine_plane_line_quotient_tomography_thm3076.out
script_sha256: 170e44f7aac923d04d3e68e35bca0d5e1c742aeca35a704b2bd63b220792c7b2
output_sha256: 4a744564d04eac25472b92a1d6f82203911902ba0c43ced85fcaded28e6b1084
hash_basis: LF-normalized bytes
---

# THM-3076 -- finite affine-plane line-quotient tomography and the p=2 three-view law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The affine-plane projector algebra

Let `p` be prime, let

```text
V=F_p^2,                                                     (1)
```

and let `k` be a field with `char(k) != p`.  Translation by `v in V` acts on
`k^V` by `T_v f(x)=f(x+v)`.  For every line `L <= V`, put

```text
P_L=(1/p) sum_(ell in L) T_ell,
P_V=(1/p^2) sum_(v in V) T_v.                               (2)
```

There are `p+1` lines, indexed by `P^1(F_p)`.  The following identities hold:

```text
P_L^2=P_L,
P_L P_M=P_V                         for L != M,
sum_(L in P^1(F_p)) P_L=I+pP_V.                            (3)
```

Consequently, with

```text
E_L=P_L-P_V,                                                (4)
```

one has

```text
E_L^2=E_L,                 E_L E_M=0 for L != M,
sum_L E_L=I-P_V,
k^V=im(P_V) direct_sum direct_sum_L im(E_L),               (5)
rank(P_V)=1,               rank(E_L)=p-1.                  (6)
```

**Proof.**  The first identity in `(3)` is the Reynolds identity for the
order-`p` translation subgroup `L`.  Distinct lines satisfy `V=L direct_sum M`,
so every `v in V` occurs exactly once as `ell+m`; this proves the second.
Every nonzero vector lies on a unique line, while zero lies on all `p+1`
lines.  Counting translation coefficients gives the third identity.

Also `P_LP_V=P_V`.  Substitution into `(4)` proves the idempotence and
annihilation laws in `(5)`, and summing `(4)` proves the resolution of the
identity.  Finally, `P_L` projects onto the `p`-dimensional space of functions
constant on `L`-cosets, whereas `P_V` projects onto constants.  This proves
`(6)` and the direct decomposition.  QED.

The `p+2` operators

```text
{P_L : L in P^1(F_p)} union {P_V}                          (7)
```

are linearly independent.  Indeed, on the translation coefficient of a
nonzero vector in `L`, only `P_L` and `P_V` occur; comparing different lines
forces all line coefficients, and the zero coefficient then forces the last
one.  This is the precise operator span in which the reconstruction
coefficients below are unique.

## 2. Every labelled subbank and its exact image

For a labelled set `S` of `s>=1` lines, consider the separate-table map

```text
O_S:k^V -> product_(L in S) im(P_L),
O_S(f)=(P_L f)_(L in S).                                   (8)
```

Its canonical observable projector is

```text
Q_S=P_V+sum_(L in S)E_L
   =sum_(L in S)P_L-(s-1)P_V.                              (9)
```

The exact ledger is

```text
rank(O_S)=rank(Q_S)=1+s(p-1),
ker(O_S)=direct_sum_(L notin S) im(E_L),
dim ker(O_S)=(p+1-s)(p-1).                                (10)
```

The image in `(8)` has an equally concrete description.  A tuple
`(y_L)_(L in S)`, with `y_L in im(P_L)`, occurs if and only if all its global
means agree:

```text
P_V y_L=c                         independently of L.      (11)
```

In that case

```text
f_S=sum_(L in S)y_L-(s-1)c                               (12)
```

realizes the tuple and equals `Q_S f` for every preimage `f`.

**Proof.**  In the splitting `(5)`, the table `P_L f` consists of the common
constant component plus the single centered component `E_Lf`.  Thus `(8)`
sees exactly `P_V` and the `s` selected centered summands, proving `(9)--(11)`.
Equation `(12)` is then immediate.  QED.

Within the scalar operator span `(7)`, `(9)` is the unique decoder: testing a
putative linear combination on each selected `E_L` sector forces every view
coefficient to be one, and testing constants forces the `P_V` coefficient to
be `1-s`.  Arbitrary linear decoders are not unique, because the common mean
can be extracted from any displayed table.

For the full bank, `(9)` is the identity and gives the reconstruction law

```text
f=sum_L P_Lf-pP_Vf.                                       (13)
```

The empty bank is separate: it has rank zero.  The expression
`1+s(p-1)` must not be extrapolated to `s=0`, because no displayed table then
supplies the common mean.

## 3. Sharp missing-direction and aggregation boundaries

If `M notin S` and `a+M,b+M` are distinct cosets, then

```text
h_(M,a,b)=1_(a+M)-1_(b+M)                                (14)
```

has mean zero, is fixed by `P_M`, and is killed by every `P_L` with `L != M`.
For fixed `M`, these differences span the full `(p-1)`-dimensional missing
channel `im(E_M)`.  Thus every kernel dimension in `(10)` is attained by
literal integral hostiles.

There is a second, logically distinct boundary.  For two different lines
`L,M`, the point-label map

```text
V -> V/L x V/M                                            (15)
```

is bijective, because `L intersect M={0}`.  Nevertheless the two separately
aggregated quotient tables have rank only

```text
2p-1                                                     (16)
```

and blind dimension `(p-1)^2`.  Joint labels remember which two cosets meet;
separate marginals forget that coupling.  This is the finite-plane form of
the joint-label/separate-table boundary in
[THM-3072](THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md).

## 4. The characteristic boundary is exactly `char(k)=p`

There is no extra exclusion when `char(k)` divides `p+1`.  As long as `p` is
invertible, the proof of `(3)--(13)` is unchanged.  In a field where
`p+1=0`, the raw sum `sum_LP_L` kills constants, while the correction
`-pP_V` in `(13)` restores them.  This is a sharp witness that the mean
sidecar is load-bearing.

In characteristic `p`, the normalized operators `(2)` do not exist.  Clear
denominators and write

```text
N_L=sum_(ell in L)T_ell,
N_V=sum_(v in V)T_v.                                      (17)
```

Then in `F_p`-coefficients

```text
N_L^2=0,
N_LN_M=N_V                         for L != M,
sum_LN_L=N_V.                                             (18)
```

The stacked cleared observation deck has exact rank

```text
p(p+1)/2.                                                 (19)
```

**Proof of `(19)`.**  In the group algebra,

```text
F_p[V]=F_p[X,Y]/(X^p,Y^p),                                (20)
```

with augmentation ideal `I=(X,Y)`.  If `g_L` generates `L`, then

```text
N_L=sum_(t=0)^(p-1)g_L^t=(g_L-1)^(p-1).                  (21)
```

For the finite slopes and the vertical line, the degree-`p-1` initial forms
are

```text
(X+aY)^(p-1)  (a in F_p),                 Y^(p-1).        (22)
```

All binomial coefficients in `(22)` are nonzero modulo `p`; the resulting
Vandermonde matrix shows that these `p+1` forms span the full homogeneous
degree-`p-1` piece.  Hence

```text
sum_L F_p[V]N_L=I^(p-1).                                  (23)
```

The monomials `X^iY^j` with `0<=i,j<p` and `i+j>=p-1` number
`p+(p-1)+...+1=p(p+1)/2`.  Since every line-sum matrix is symmetric, this
ideal dimension equals the row rank of the stacked observation deck.  QED.

Thus the bad-characteristic object is an augmentation filtration with
nilpotent line norms, not projector tomography.  Equations `(18)` do not
rescue the lost semisimple direction splitting.

## 5. Why `p=2` and `p=3` meet at `A4`

For `p=2`, `V` has three lines and the nonzero linear map `rho` of order
three cycles them.  The ranks of one, two, and three views are

```text
2,3,4.                                                     (24)
```

On each fixed-direction four-point fibre in THM-3072, its three involutions
are precisely translations by these three lines.  Taking the direct sum over
the three direction fibres multiplies `(24)` by three and gives THM-3072's
`6,9,12` table ranks and reconstruction

```text
I=P_0+P_1+P_2-2P_V.                                      (25)
```

The order-two/order-three modular interface is literal.  Put

```text
S=[0 -1;1 0],                   T=[1 1;0 1],
R=ST.                                                      (26)
```

Reduction modulo three acts on `P^1(F_3)`.  In
`PSL_2(F_3)`, the permutations induced by `S,R,SR` have orders `2,3,3`, so

```text
PSL_2(F_3)=Delta(2,3,3)=A4.                               (27)
```

After the explicit relabelling `(1,0,2,3)`, `S` becomes translation by
`(1,0)` on `F_2^2`, while `R` becomes `rho`.  Hence

```text
A4=F_2^2 semidirect C3                                   (28)
```

generator by generator, not merely by comparing group orders.

This explains one precise co-occurrence of the repository's binary and
ternary grammars: the three binary quotient directions form the `V4`
translation torsor, and the modular order-three generator cycles them.  The
extra relation `(SR)^3=1` is also the stopping rule.  It is imposed when

```text
PSL_2(Z)=C2*C3 -> A4,                                     (29)
```

so free modular normal forms, Bass--Serre origins, and binary/ternary tree
ancestry are not retained.  Moreover `p=3` affine-plane tomography has four
line directions; the word *ternary* in `(26)--(29)` refers to the subgroup
order, not to a three-view theorem.

THM-2768 and THM-2632 already own the abstract modular quotient and its CRT/
Farey boundaries.  Equation `(28)` is recorded here only to identify the
finite-plane operators with THM-3072's three `C2` charts.

## 6. Relation to the existing finite Radon theorem

[THM-850](THM-850-chi7-face-charge-cyclotomic-law.md), Section 7.3, already
proves over characteristic zero that `s` of the eight line directions in
`F_7^2` have rank `1+6s`, and that all eight recover the plane signal by
Fourier slices.  The present theorem generalizes that antecedent to arbitrary
prime `p` and arbitrary coefficient field with `char(k) != p`; its additional
content is the intrinsic projector algebra, the exact common-mean image and
decoder, sharp omitted-channel bases, the good/bad characteristic boundary,
and the generatorwise `p=2` modular bridge.

This theorem still supplies no map from affine-plane points to:

```text
Farey flanks or THM-2056 Gram defects;
quartic roots, V4 resolvent channels, or Keller fibres;
LRC speeds, owners, endpoint currents, or relation addresses.             (30)
```

In particular, `(15)` does not turn two aggregated physical observables into
a common atom, and `(28)` does not identify a quartic `V4` torsor with an LRC
parity torsor.  Those remain separate realization and coupling debts.

## 7. Exact evidence and truth surface

Run

```bash
python 04-computation/finite_affine_plane_line_quotient_tomography_thm3076.py
python -O 04-computation/finite_affine_plane_line_quotient_tomography_thm3076.py
```

The companion uses rational and finite-field arithmetic, explicit exceptions,
and no truth-bearing Python assertions.  It checks every labelled nonempty
subbank for `p=2,3,5,7`; every projector, centered channel, pair product,
joint-label bijection, and omitted-line hostile; good-characteristic controls
`(p,ell)=(2,3),(3,2),(5,2),(7,2)` where `ell` divides `p+1`; bad-characteristic
ranks `3,6,15,28`; and the generatorwise `A4` conjugacy.  Normal and optimized
transcripts LF-normalize byte-for-byte to the stored output.

```text
PROVED IN THE CANDIDATE:    equations (3)--(29), including every subbank rank,
                            exact reconstruction, sharp hostiles, and the
                            characteristic-p augmentation boundary.

NOT PROVED:                 a binary/ternary tree intertwiner; a physical
                            quartic, Farey, Keller, or LRC realization; any
                            row exclusion or ledger decrement.             (31)
```

QED.
