---
id: THM-3524
title: "Folded C7/K4 common-drift transporter classification and rank obstruction"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.  For the fixed raw
  THM-2594 7-by-13 response and the fixed THM-3514 U_full 4-by-13 Walsh
  response, the Paley-C7 negation quotient is a weighted looped symmetric
  four-block digraph with off-diagonal K4 support.  The source has exact rank
  three and all 91 split Fourier coordinates nonzero; the target has rank
  four and all 52 Walsh/drift Fourier coordinates nonzero.  Complete
  48/192/2401 labelled allocation universes admit no common scalar
  convolution even modulo independent amplitude and cyclic-origin gauges.
  More generally, fixed channel mixing followed by any one common right
  operator cannot raise rank three to four.  The independent projective
  systems have annihilator-only nullity.  Four channel-dependent invertible
  circulants survive only as a marked split-field equivalence.  No ancestry,
  Boolean-current, physical-current, absolute-H1, bispectrum, grouped-row,
  scalar-exclusion, or LRC(14) conclusion follows.
source: codex/folded-c7-k4-canon-audit/2026-08-16
depends_on:
  - THM-2594-realized-theta-slaved-contraction-at-the-r5-window
  - THM-3514-ufull-owner-boundary-k4xf13-endpoint-factorization-and-walsh-spectrum
related:
  - THM-3518-ufull-all-role-phase-normalized-dual-contraction-and-cut-cycle-certificate
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
script: 04-computation/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.py
output: 05-knowledge/results/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.out
script_sha256: 5fe3f696f122869462ff73b1b4ebdb957fa8ca7ee3692c25ef94f0f7efae81cf
output_sha256: b97db9c8b0582d868c8386be5f5d8dff72968062a9b5ae79956678cd86c0b98a
semantic_sha256: bae3749e01361ebbaf6c9cc0e7160d2e4d22d6df420d4512ba94f778718af3a1
independent_audit_script: 04-computation/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.py
independent_audit_output: 05-knowledge/results/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.out
independent_audit_script_sha256: 0bf56db1acfd83d6f161dd95b5d428df90f8cb5e6ec24e729081f24a86f287f2
independent_audit_output_sha256: c20c59264590cd42b108363c91736df8c9e921ac94b09c1f5fbb5d977f300859
independent_audit_semantic_sha256: a6d6ad7a9891af8bddcba6fe7cf9f2554214731043844098f8556d1003885222
hash_basis: LF-normalized UTF-8 bytes; exact rational and split-field ledgers
---

# THM-3524 -- the folded `C7` source cannot linearly carry the four `K4` Walsh rows

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**

This theorem compares two already proved, fixed finite objects.  It classifies
one linear transporter class; it does not construct a new LRC object.

## 1. The fixed source and target

Let

```text
A(ell,theta),       ell in F_7, theta in F_13,          (1)
```

be the raw common-base response table materialized in THM-2594.  Its only
possibly nonzero columns are `theta=0,1,2`.  Over a field containing primitive
seventh and thirteenth roots, define the seven septimal source profiles and
their drift spectra by

```text
X_b(theta)=sum_ell A(ell,theta) zeta_7^(-b ell),
Xhat_b(k)=sum_theta X_b(theta) zeta_13^(-k theta).       (2)
```

Let `D_ij(d)` be the fixed THM-3514 `H-q5` endpoint bridge bucket, with
`(i,j)` in `F_2^2` the two owner-boundary chamber bits and `d in F_13` the
relative endpoint drift.  Its four Walsh rows are

```text
Y_0=D_00+D_01+D_10+D_11,
Y_L=D_00+D_01-D_10-D_11,
Y_R=D_00-D_01+D_10-D_11,
Y_X=D_00-D_01-D_10+D_11.                               (3)
```

All comparisons below use the common certified split field

```text
P=572252886246508880869,
zeta_7 =353818603057943120846,
zeta_13=505438565698892403012.                          (4)
```

Nonzero reduction in `(4)` certifies the corresponding characteristic-zero
cyclotomic quantity because all pinned denominators are units modulo `P`.

## 2. The exact Paley fold is weighted and looped

Orient the Paley tournament `T_7` by the quadratic residues `{1,2,4}` and
fold its vertices under negation:

```text
O_0={0}, O_1={1,6}, O_2={2,5}, O_3={3,4}.              (5)
```

The directed arc-count matrix between the four blocks is exactly

```text
0 1 1 1
1 1 2 2
1 2 1 2
1 2 2 1.                                              (6)
```

Since `-1` is a quadratic nonresidue modulo seven, simultaneous negation
reverses every Paley arc.  This forces the symmetry in `(6)`.  The diagonal
ones record the internal arc in each nontrivial two-point orbit.  Therefore
the quotient is a **weighted, looped, symmetric four-block digraph** whose
off-diagonal support is a bidirected or undirected `K4`; it is not a
four-vertex tournament.

Multiplication by two cycles the three nonzero blocks,

```text
O_1 -> O_2 -> O_3 -> O_1,                              (7)
```

while the three nontrivial characters of `F_2^2` cut the four Walsh states
into the three perfect matchings

```text
01|23, 02|13, 03|12.                                   (8)
```

This is a common three-matching grammar, not a canonical identification.
The rational convolution algebras retain different descent data:

```text
Q[C_7]^inv = Q x Q(zeta_7+zeta_7^-1),
Q[F_2^2]   = Q^4.                                      (9)
```

The cubic factor in `(9)`'s first line has minimal polynomial
`x^3+x^2-2x-1`, discriminant `49`, and no rational root, hence is a cyclic
cubic field.  The two algebras have respectively two and four rational simple
factors, so they are not isomorphic over `Q`.  They both split as four copies
of the field only after a marked splitting extension; that extension forgets
the rational descent datum.

## 3. Full scalar spectra coexist with a rank defect

The exact two-dimensional transform in `(2)` satisfies

```text
Xhat_b(k) != 0       for every (b,k) in F_7 x F_13.     (10)
```

Thus the source has split Fourier support `91/91`.  Independently, reductions
at `547`, `1093`, and `2003` reproduce both support `91/91` and rank three.
The raw rational table has a nonzero `3 x 3` minor on rows `(1,3,4)`:

```text
48132400318030844942585467137428417511789621695304247254245705613312000.
                                                               (11)
```

Because only three columns can be nonzero, `(11)` proves

```text
rank(X)=3                                                (12)
```

already in characteristic zero.

THM-3514 gives, and both present implementations reconstruct,

```text
Yhat_chi(k) != 0    for all 4*13=52 pairs (chi,k),
rank(Y)=4.                                               (13)
```

The first four target drift columns have determinant
`154930576847978146912 mod P`, so `(13)` is witnessed by an explicit
nonzero minor.  Equations `(10)`--`(13)` separate two notions that must not be
conflated: every individual source row is a cyclic `C_13` vector, but the
seven-row source bank spans only three channel directions.

## 4. Complete allocation census

The negation classes of source modes are

```text
{0}, {1,6}, {2,5}, {3,4}.                              (14)
```

For an allocation `f` from the four target channels to labelled source modes,
define its four Fourier multiplier rows

```text
mu_chi(k)=Yhat_chi(k)/Xhat_(f(chi))(k).                 (15)
```

All entries in `(15)` are units by `(10)` and `(13)`.  Two multiplier rows
are equivalent under independent nonzero amplitude and cyclic-origin gauges
exactly when

```text
mu_1(k)=c zeta_13^(ks) mu_2(k)                         (16)
```

for some `c != 0` and `s in F_13`.

The primary and independent programs exhaust the following universes by
different constructions:

| universe | exact predicate | size |
|---|---|---:|
| trivial-preserving folded lifts | channel `0` receives mode `0`, and the other channels receive one labelled member of each nonzero pair | `3! 2^3=48` |
| all folded labelled bijections | the four channels receive the four classes in `(14)` in arbitrary order, with one labelled member from each nonzero pair | `4! 2^3=192` |
| unrestricted labelled selections | every function from four channels to seven modes, repetitions allowed | `7^4=2401` |

In every one of these allocations, the four rows `(15)` occupy **four**
distinct equivalence classes under `(16)`.  Consequently there is

```text
no exact common multiplier;
no common multiplier after independent amplitude gauges;
no common multiplier after independent cyclic-origin gauges.  (17)
```

The `48` and `192` universes are marked split-field lifts, not rational
folded-algebra maps.  The `2401` universe is complete for labelled row
selection, not for arbitrary linear combinations; the next section handles
arbitrary fixed channel mixing.

## 5. The general fixed-linear obstruction

Regard the full source bank as a `7 x 13` matrix `X` and the target bank as a
`4 x 13` matrix `Y` over the common split field.  A fixed channel mixer and
one common drift operator would have the form

```text
Y=M X C,                                                (18)
```

where `M` is `4 x 7` and the same `13 x 13` operator `C` acts on every
channel.  It need not be circulant or invertible for the following inequality:

```text
rank(M X C) <= rank(X)=3 < 4=rank(Y).                  (19)
```

Thus `(18)` is impossible.  This is the mechanism behind the finite censuses:
three active deep windows cannot linearly carry four independent Walsh rows.

The independent audit also keeps every projective multiplier rather than
eliminating it.  For each of the eight choices of one member from each
nonzero pair in `(14)`, it solves

```text
M x_k=lambda_k y_k,       k in F_13,                   (20)
```

with sixteen entries of a `4 x 4` matrix `M` and thirteen scalars `lambda_k`.
Every resulting `52 x 29` system has

```text
rank 25, nullity 4, multiplier-projection rank 0.       (21)
```

Every null vector has all `lambda_k=0`, and its matrix component annihilates
the rank-three source.  Equivalently, the multiplier-eliminated projective
system has rank twelve and its four-dimensional nullspace is
annihilator-only.  Hence no nonzero projective transport is hidden in the
rank count.

## 6. Sharp formal survivor

For one chosen source row and one target row, `(10)` and `(13)` make the
multiplier `(15)` a unit in the split `C_13` group algebra.  Its inverse DFT
is therefore an invertible circulant carrying that source row to that target
row.  For the calibration

```text
b=0,1,2,3  ->  Y_0,Y_L,Y_R,Y_X,                        (22)
```

the four channel-dependent kernels each have support `13/13`; the four kernel
rows have rank four.  Their full and augmentation digests are

```text
9f03166d1b8d730673aeccf7bff9196d2de42e0ff9fa7c08680d38a772dc041f,
351976b6d1ac51915d686b32e1cc0f3f3dfa53d1f33f7c8831f7a74a495e3025. (23)
```

Full transport minus centred transport is a constant-kernel margin in each
channel.  The four kernels are pairwise inequivalent under `(16)`.  Their
existence is the sharp positive survivor, but it is automatic for pairs of
cyclic vectors and supplies no common ancestry map.

## 7. The THM-3518 sidecars do not evade this class

THM-3518's factor `zeta_13^(q_t a_L)` lives on the U_full endpoint guard
sheets, after the source ancestry coordinates have been forgotten.  No proved
map identifies `a_L` with a THM-2594 root or deep coordinate.  A single common
phase/right action on the source would preserve rank three; channel-dependent
phases leave the common-operator class and recover only the formal situation
of Section 6.

Likewise, THM-3518's thirteen-edge arrays are gradients.  Their incidence
rank is seven, cycle rank is six, and all `56,592` certified cycle pairings
vanish.  A zero absolute-cycle packet adds no fourth source response.  The
nonzero weighted tree polynomial is a function of a chosen coboundary
representative, not a hidden absolute `H^1` coordinate.  These facts explain
why THM-3518 is a relevant hostile sidecar but not a dependency for the rank
proof `(19)`.

## 8. Independent verification

The primary program imports the pinned THM-3514 all-role parent and parses the
pinned THM-2594 aggregate.  The independent program imports neither the
submitted transporter nor that U_full parent.  It instead rebuilds the target
from the independently audited THM-3514 atom tables by direct `tau` slices,
recomputes the source transform directly, generates the `48/192` universes by
predicate-filtering all `7^4` choices, and uses the augmented systems
`(20)`--`(21)`.

Run from the repository root:

```text
python -B 04-computation/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.py
python -B -O 04-computation/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.py
python -B 04-computation/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.py
```

All four executions reproduce their stored outputs exactly.  The two semantic
digests are the values pinned in the frontmatter.

## 9. Exact scope and failure boundary

The proved connection contract is

```text
source:    fixed THM-2594 raw 7 x 13 common-base response;
target:    fixed THM-3514 U_full 4 x 13 endpoint Walsh response;
map:       marked split-field row allocation or fixed linear channel mixing
           followed by one common right operator;
preserved: character labels, drift coordinate, scalar Fourier support;
lost:      rational descent under a marked split, positivity, atoms, owner,
           word, horizon, chronology, and common ancestry;
obstruction: source row rank 3 versus target row rank 4;
needed linear repair: at least one independent lawful source response;
broader live repairs: nonlinear, channel-dependent, frequency-dependent, or
           unmarginalized common-stalk relations.                         (24)
```

The needed extra source response is necessary only for this fixed-linear
class and is not sufficient for a physical transplant.  This theorem does
not exclude the broader repairs in `(24)` and proves no U_full ancestry
relation, Boolean coupling, lawful or physical current, nonzero absolute
`H^1`, LRC `7 x 13` bispectrum bridge, grouped coefficient, all-unit
projector, scalar-row exclusion, or LRC(14). **QED.**
