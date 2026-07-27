---
id: THM-2592
title: "Fallback-rail digit-diagonal pullback and primitive Bockstein"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED core;
  PROVISIONAL EXACT ADDENDUM UNDER INDEPENDENT HOSTILE AUDIT.
  On the canonical typed row, the literal route-two pullback of the 84
  THM-2586 theta-zero rail cells against all thirteen normalized THM-2585
  target sections is positive in exactly 39 of 1092 cases: precisely the
  three fallback cells (s4,ell4,v,t)=(7,4/5/6,6,12), for every q.  After one
  global primitive-content reduction, 38 of the 39 positive slices have a
  nonzero first Bockstein in all six owner colours and 37 have unit septimal
  slice polynomial.  More sharply, all 81 primary rails are orthogonal to the
  entire unconditioned delayed word, so no future-digit selector or affine
  digit offset can attach them; the fallback digit support is exactly
  h=1,...,11.  Boolean q-section aggregates and one-section-per-chart
  transversals have the finite no-cancellation spectra stated below.  The
  construction is an actual common-x fibre product, but does not identify the
  two clock labels, a named old head, or a canonical semantic endpoint.
source: codex-2026-07-28-joint-rail-projector-pullback
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-2586-depth-five-arrival-to-future-root-diagonal
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2590-boolean-bockstein-and-theta-selector-incidence-spectrum
script: 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
output: 05-knowledge/results/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.out
script_sha256: 87aa2407d7c3417e093e4ab3d336c704dd8bc4ccf1dcc8e357c6671a60f8cce1
output_sha256: 002ed091c79a96c084dc579ab149c815705a092e7a7eadadc381a0c15f8e5c33
hash_basis: LF-normalized bytes
---

# THM-2592 -- fallback-rail digit-diagonal pullback and primitive Bockstein

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED core.**

**PROVISIONAL EXACT ADDENDUM UNDER INDEPENDENT HOSTILE AUDIT:** the
union-first, affine-offset, Boolean-section, and cross-chart transversal
statements in Sections 3 and 4 below.

## 1. The common physical carrier

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),                    (1)
```

with

```text
T=297836897838480,   d=13^5,   R=13^6.                 (2)
```

THM-2584 writes its `sigma={b}`, depth-five packet in the route-two
physical coordinate `x` as

```text
13 U(x)V(x-s4/13),                                     (3)
```

where

```text
U=P_d(1_Qb P_169 1_E),   V=P_d 1_E.                    (4)
```

Before integration retain

```text
ell4 : frac(169x) lies in the ell4-th owner window,
v    = floor(13x),
t    = floor(26x) mod 13.                              (5)
```

For `s4!=0`, THM-2586's deterministic theta-zero choice is

```text
(v,t)=(0,0),
```

except at

```text
(s4,ell4)=(7,4),(7,5),(7,6),                           (6)
```

where `(v,t)=(6,12)`.  In every case `7t=v mod 13`.

THM-2585 has a different word and a different target action.  Let

```text
F_(ell5,s5)(x)
```

be its present `sigma={a}` source/target packet, let

```text
Delta_r(x)=d(c3 x-r/13),                               (7)
```

and let `Q^+_(ell5,h)` denote its delayed word restricted to root digit
`h`.  Thus

```text
Q^+_(ell5,h)(Rx)
 =Q^+_ell5(Rx) 1_(floor(13 frac(Rx))=h).                (8)
```

For `q in F_13` put `s5=-q`.  The joint tensor is

```text
J_(s4,ell4,q;ell5,r)
 =13 int U(x)V(x-s4/13)
       1_(ell4,v,t)(x)
       F_(ell5,-q)(x) Delta_r(x)
       Q^+_(ell5,v)(Rx) dx,                            (9)
```

where `(v,t)` is the selected cell above and `1_(ell4,v,t)` denotes all
three restrictions in (5).

Although `U,V` are Perron averages, (9) is still a finite Boolean
fibre-product bank.  With THM-2586's preimages

```text
X_a=(x+a)/d,
Z_(a,beta)=(X_a+beta)/169,
Y_e=({x-s4/13}+e)/d,
```

its exact expansion is

```text
J=13/(169d^2) sum_(a,beta,e)
    mu(S_(a,beta,e) intersection
       F_(ell5,-q) intersection Delta_r intersection
       {Q^+_(ell5,v)(Rx)=1}).                          (9a)
```

Here `S_(a,beta,e)` is THM-2586's literal Boolean ancestry sheet with the
three root/owner indicators already inserted.  Thus positivity of any entry
of (9) really supplies a positive Boolean pullback sheet; no probabilistic
coupling is being chosen after the fact.

Equation (9) is the required common carrier.  It is formed **before** the
THM-2584 rail marginal and before the THM-2585 target/deep transforms.
Moreover the two occurrences of `v` give the exact time identity

```text
floor(13x)=v=floor(13 frac(13^6 x)).                    (10)
```

Thus (9) is an actual fibre product over the same physical `x`, not a
composition inferred from equal marginal labels.  This is precisely the
distinction required by THM-2315.

## 2. One global integer normalization

Write the integer numerator profiles of (4) as

```text
U=N_U/(169d),   V=N_V/d.                                (11)
```

The exact prefix formula for the last factor in (9) gives integers
`A_(s4,ell4,q;ell5,r)` such that

```text
J_(s4,ell4,q;ell5,r)
 =13 A_(s4,ell4,q;ell5,r)/(169 d^2 R T).                (12)
```

No cell is normalized separately.  Across the complete
`84*13*7*13` bank, with the identically zero `r=0` cells retained, the
nonzero `A` have global content

```text
gcd(A)=25465440
      =2^5*3*5*7*11*13*53.                             (13)
```

Including the route-two factor `13`, the physical raw content is

```text
g=331050720
 =2^5*3*5*7*11*13^2*53.                                (14)
```

Define the one globally primitive tensor

```text
a=A/25465440.                                           (15)
```

Then every joint mass has the common scalar form

```text
J=(2/202345849311155871624914247) a.                    (16)
```

The serialization digest of the full primitive bank, ordered by
`(s4,ell4,q,ell5,r)` in the companion's deterministic rail order, is

```text
1286b1f2baa1f05299d93a1074db19afc3f99eb3d44500fd1de39c97c62c300c. (17)
```

This single normalization is load-bearing.  Re-primitivizing the 1092
pair slices separately would define different mod-13 objects and is not
used here.

## 3. Exact support: only the fallback rail attaches

The delayed word has the exact digit gate

```text
Q^+_(ell5,0)=empty                                      (18)
```

for every `ell5=0,...,6`, while every `Q^+_(ell5,6)` has positive
measure.  This is checked directly as a half-open interval identity on the
grid `T`; it is stronger than an integral cancellation.

Every one of the 81 primary THM-2586 cells has `v=0`.  Equations (8) and
(18) therefore make its complete pullback (9) identically zero for all
thirteen `q`.  The three fallback cells have `v=6`, and exact integration
gives

```text
sum_(ell5,r) J_(s4,ell4,q;ell5,r)>0

iff

s4=7, ell4 in {4,5,6}, (v,t)=(6,12), q arbitrary.       (19)
```

Consequently exactly

```text
39/1092                                                  (20)
```

rail-by-target-section pullbacks are positive.  Their unscaled `A` totals
range from

```text
163706661702634491201720960
```

to

```text
315647624785807202076027840.                             (21)
```

At the fine `(ell5,r)` level there are exactly

```text
1404/91728                                               (22)
```

positive entries.  Among the 39 positive pair slices, the number of
positive fine entries has histogram

```text
24:7, 36:25, 48:7.                                      (23)
```

Thus the attachment is genuine but maximally localized: the same digit
gate that proves positivity on the fallback cell annihilates the 81-cell
primary rail.

### 3.1 The primary obstruction is the whole word, not digit zero

The digit-zero explanation is true for the diagonal (9), but is not sharp.
Replace `Q^+_(ell5,v)` in (9) by the whole unconditioned delayed word
`Q^+_ell5`, and call the resulting nonnegative entries `J^all`.  Exact
common-`x` integration gives

```text
sum_(ell5,r) J^all_(s4,ell4,q;ell5,r)=0                (23a)
```

for every one of the `81*13=1053` primary-rail/target-section pairs.
Every summand and every underlying Boolean sheet in (9a) is nonnegative.
Therefore (23a) implies that each fine primary entry is zero almost
everywhere on the support of the **entire** delayed word.

The thirteen digit pieces in (8) are a disjoint nonnegative partition of
that word.  It follows at once that every digit selector vanishes on every
primary rail.  This includes selectors depending arbitrarily on
`(s4,ell4,q,ell5,r)`, arbitrary subsets of future digits, and arbitrary
linear mixtures supported on those pieces.  No re-selection of the future
digit inside the same delayed word can attach a primary rail.

### 3.2 Complete affine-offset spectrum

For `c in F_13`, define the affine graph section by replacing the last factor
of (9) with

```text
Q^+_(ell5,v+c)(Rx),                                    (23b)
```

where the sum is modulo thirteen.  On all 81 primary rails it vanishes by
(23a), for every `c`.  On the three fallback rails `v=6`; the exact positive
pair census by literal future digit `h` is

```text
h=0: 0,   h=1,...,11: 39 each,   h=12: 0.              (23c)
```

Consequently the number of positive rail-by-`q` pairs for (23b) is

```text
39,  c in {0,1,2,3,4,5,8,9,10,11,12};
 0,  c in {6,7}.                                       (23d)
```

Whenever it is `39`, the support is again precisely the three cells (6),
for all thirteen `q`.  Thus no fixed affine offset attaches even one primary
rail.  In particular `c=6`, the most tempting translation from the old
digit-zero hole to the observed fallback digit, instead sends the fallback
digit `6` to the empty future digit `12` and kills the complete bank.

Condition (23b) is a lawful literal affine correspondence between the two
chosen `F_13` root-digit charts on the common physical coordinate.  It is
not, without a further transported-head or root-gauge sidecar, equality of
semantic old and future heads, a named deck transition, a relation residue,
or a THM-2305 endpoint identification.

## 4. The globally primitive first Bockstein

For a positive pair slice define

```text
Y_ell5=sum_(r=1)^12 a_(s4,ell4,q;ell5,r) r^(-1) mod 13,

Y(z)=sum_(ell5=0)^6 Y_ell5 z^ell5
       in F_13[z]/(Phi_7).                              (24)
```

As in THM-2585, the deep augmentation identity follows from the literal
`r=0` zero.  If `kappa in F_7^*`, the first Bockstein of the corresponding
primitive deep cycle is

```text
beta(D^(kappa))=Omega Y(zeta_7^kappa),

Omega=(1,2,...,12) mod 13.                              (25)
```

The computation of (24) is made anew from the primitive tensor (15); no
THM-2585 marginal Bockstein is inherited through conditioning.

Exactly one positive slice has zero Bockstein:

```text
(s4,ell4,v,t,q)=(7,6,6,12,10),
(Y_0,...,Y_6)=(0,0,0,0,0,0,0).                         (26)
```

Every other positive slice has nonzero (25) for every one of the six
nonzero owner colours.  Hence

```text
38/39 positive slices,
228/234 positive-slice owner profiles,                  (27)
```

survive the first Bockstein.

Multiplication by `Y` in `F_13[z]/(Phi_7)` is invertible on 37 of the 39
positive slices.  Besides the zero slice (26), the unique nonzero nonunit
is

```text
(s4,ell4,v,t,q)=(7,4,6,12,11),
(Y_0,...,Y_6)=(9,11,0,10,0,0,0),                       (28)
```

whose multiplication determinant is zero mod 13.  Its six individual
owner-colour Bocksteins are nevertheless all nonzero.

### 4.1 Boolean sections and the three-chart transversal atlas

The one global primitive normalization (15) permits coefficient sums across
the thirteen literal `q` sections without reprimitivizing a slice.  For one
fallback chart `ell4` and a Boolean subset `S subset F_13`, put

```text
Y_(ell4,S)=sum_(q in S) Y_(ell4,q) in R_7,
R_7=F_13[z]/(Phi_7).                                   (28a)
```

The complete `2^13` census on each actual common-`x` fallback carrier is

| `ell4` | distinct images | zero subsets | max fibre | units | nonunits |
|---:|---:|:---|---:|---:|---:|
| 4 | 864 | `empty` | 40 | 8118 | 74 |
| 5 | 516 | `empty` | 88 | 8029 | 163 |
| 6 | 328 | `empty,{10}` | 110 | 8114 | 78 |

The unit/nonunit columns count all 8192 subsets, so the zero class is included
among the nonunits.  Thus every nonempty Boolean aggregate on charts 4 and 5
has nonzero first Bockstein in all six owner colours.  On chart 6 the only
nonempty zero is the already known singleton `S={10}` from (26); adjoining
any other target section repairs it.

For completeness, this owner-colour conclusion is not an inference from a
unit determinant.  For every `kappa in F_7^*`, substitution
`z -> z^kappa` is an automorphism of `R_7`.  Moreover, with
`u=zeta_13-1`, the factor in (25) is

```text
Omega=-u^11 in F_13[u]/(u^12),                         (28b)
```

so multiplication by `Omega` embeds `R_7` into the corresponding socle
summand.  Hence every nonzero class in the table remains nonzero under every
labelled owner substitution and after multiplication by `Omega`, including
the displayed nonunit zero divisors.

There is also no cancellation among chart transversals.  If exactly one
section is selected on each of two distinct charts, all `13^2=169` sums are
nonzero.  Their `(distinct images,max fibre)` data are

```text
(ell4,ell4')=(4,5): (47,15),
              (4,6): (53,30),
              (5,6): (50,20).                         (28c)
```

All `13^3=2197` three-chart sums are nonzero; they have

```text
295 distinct images, max fibre 103,
2173 units and 24 nonunits.                             (28d)
```

Thus every nonempty one-section-per-chart transversal is charged except the
single chart-6 section `q=10`.  This is a finite common-carrier
coefficient/cosheaf no-cancellation statement.  It does not assert that the
chosen sections share one positive Perron ancestry sheet, nor can it cross
the primary whole-word orthogonality (23a).

## 5. Independent exact controls

The companion uses only integer and `Fraction` arithmetic and contains the
following separate checks.

1. THM-2584's tensor is rebuilt by its root-chart and physical-`x` routes;
   all entries agree before this pullback is formed.
2. The selected physical rail pieces re-sum exactly to the stored
   `K_(ell4)(s4,v,t)` numerators.
3. Three deep-comb restrictions are computed both by direct nonwrapping
   window intersection and by the inherited complement construction.
4. Eight positive delayed overlaps are computed both by one weighted
   prefix sweep and by grouping equal integer profile levels and invoking
   THM-2585's original unweighted sweep term by term.
5. Three unconditioned cells recover the original THM-2585 delayed-overlap
   formula before the rail factor is inserted.
6. All 1053 primary pairs are tested against the complete delayed word before
   digit splitting.  The thirteen exact digit pieces re-sum to that word;
   all thirteen affine offsets are then checked on the fallback cells.
7. All 24,576 within-chart Boolean subsets, all 507 two-chart singleton
   transversals, and all 2197 three-chart singleton transversals are enumerated
   in `R_7`, including unit, fibre, zero, and owner-colour checks.
8. Normal and optimized execution reproduce the stored transcript exactly.

An independent hostile audit rederived the common-`x` and digit-time typing,
the Boolean-sheet expansion behind the two Perron profiles, the single global
content and its 13-adic valuations, the complete support census, and both
exceptional Bockstein slices.  It also replayed normal and optimized modes
against the stored transcript and recomputed both declared LF hashes.  No
theorem defect remained.

The positive support census, both contents, primitive scalar, primitive
digest, support histogram, Bockstein census, the two exceptional slices
(26)--(28), whole-word zero, affine spectrum, and Boolean/transversal atlases
are executable assertions rather than printed approximations.

Reproduce with

```text
python 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
python -O 04-computation/lrc14_fallback_rail_digit_diagonal_pullback_thm2592.py
```

## 6. Exact scope

The theorem supplies a literal common-`x` span joining three THM-2586
fallback rail cells to every normalized THM-2585 target section, and a
globally primitive nonzero Bockstein on 38 of those 39 sections.  It is the
first actual attaching carrier between these two tensors; THM-2315's
pair-marginal hostile does not apply after (9) has been formed.

It does **not** prove any of the following.

- The other 81 theta-zero rail cells attach.  They are orthogonal to the
  entire delayed word by (23a), not merely to digit zero.  Thus no
  future-digit selector, affine offset, uniform 84-cell, or Cech filler
  follows.
- The owner clocks `ell4` and `ell5` are the same coordinate.  Both are
  retained, not identified.
- The future digit in (10) is the preselected named `k_a` boundary, the old
  semantic head, or a THM-2305 endpoint phase.
- The globally primitive Bockstein selects the same hidden Perron Boolean
  sheet as any separately chosen positive mass witness.
- A relation residue, all-coordinate unit address, all-165 conclusion, row
  decrement, or LRC(14) closure follows.

The sharp next problem is therefore not transport across one missing digit.
It is to change the physical rail/word support or supply a genuine
semantic/clock/root intertwiner: every resection of the existing delayed word
already lies inside the primary orthogonality obstruction (23a).

QED.
