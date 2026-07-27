---
id: THM-2436
title: "Punctured ninety-one-stalk mixed mode and repeated-step closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Normalize
  the common quotient stalk by the guard
  G={0,...,25}. In target row h its two guard sites are {h,h-1}.
  For one or two exact-depth top blockers, every hole of the
  guard-plus-five-ordinary multiplicity lies on a blocker source
  column. Therefore an h-independent zero-sum defect is necessarily
  zero: every source occurs in a guard pair. Every essential parent
  consequently has a genuine mixed F_13 x F_7 character, hence a
  unit character modulo 91. THM-2435's essential-parent invoice
  yields one fixed mixed quotient mode on parent mass at least 1/252
  uniformly, with typed k=2 floor 1/168. THM-2506 later proves all
  72 primitive modes are nonzero pointwise and upgrades those floors
  to 2/7 and 3/7. Exact enumeration gives
  2,629 one-source cover assignments and 38,750 distinct-two-source
  assignments; every cover repeats an unsigned ordinary step. The
  one-source repeated-pair cell bank has maximum 75/91, excluding
  the positive-depth shape against its 6/7 parent floor. Fixed-step
  refinement then caps the exceptional M=0 one-source spectra by
  57/91<5/7 and the exceptional two-source spectra by 67/91<6/7.
  Thus all three deep-c_3 residual shapes are empty and every
  hypothetical counterexample has nu_7(c_3)<=M. The theorem does not
  remove a scalar row, transport the quotient mode to a canonical
  physical endpoint current, close the c_3<=M noncirculant-graft
  branch, or prove LRC(14).
source: codex-2026-07-26-punctured-common-root-atlas
depends_on:
  - THM-2430-guard-top-common-ninety-one-root-tiling-spectrum
  - THM-2431-repeated-step-rounding-exclusion-of-guard-top-zero-blocker-types
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
related:
  - THM-2421-all-clock-septimal-ancestry-endpoint-event-detector
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
  - THM-2439-cyclic-marker-replica-degree-and-homometric-gram-boundary
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
script:
  - 04-computation/lrc14_punctured_91_stalk_mixed_mode_thm2436.cpp
  - 04-computation/lrc14_punctured_91_fixed_spectrum_referee_thm2436.py
output:
  - 05-knowledge/results/lrc14_punctured_91_stalk_mixed_mode_thm2436.out
  - 05-knowledge/results/lrc14_punctured_91_fixed_spectrum_referee_thm2436.out
script_sha256:
  - 4992f1ab58f07a3699f6c04f877f387cd48e8e4206750645d536eb11edbf1fd7
  - b4f9550249dac0ab53d8d1992440731b7ca4dfac5687d90bf65ad877f8c1e8fe
output_sha256:
  - c2659364e9d312a670f2606cde8eeb8f9ef5d82c062cfb44837970047f717f8b
  - bba34d1c67d68b992bf274a22c0e92f29b1c80c58ce15d9a7ca631b0cceff584
hash_basis: working-tree bytes (LF)
---

# THM-2436 -- a punctured stalk cannot stay vertical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2435 leaves a very specific stopping boundary. Its signed
`13 x 7` row defect either carries a unit character modulo `91`, or
is independent of the target row. The second alternative is in fact
incompatible with the normalized guard:

```text
every positive defect site is a blocker source
  + every blocker source occurs inside some guard row
  -> no nonzero vertical defect
  -> every essential parent has a mixed unit character.            (1)
```

The same finite atlas retains the repeated-step coordinate discarded
by the defect argument. With the sharpened parent image in THM-2435,
a fixed-step refinement of its formerly noncontradictory cell banks
removes all three residual shapes.

## 1. The punctured common-root model

Retain the three deep-`c_3` shapes entering THM-2435:

```text
M=0:   (k,t,b,W)=(1,5,1,8), (2,5,2,9),

M>0:   (k,t,b,W)=(2,5,1,8).                                    (2)
```

Write

```text
P=T_(7^(M+1))(A)
```

for its quotient-parent image. THM-2435 proves

```text
mu(P)>=(4+k)/7.                                                  (3)
```

For every generic `Y in P`, the top guard, five top ordinary words,
and the `b` exact-depth top blockers cover the full common quotient
stalk

```text
{(Y+s)/91:s in Z/91Z}.                                          (4)
```

Normalize the guard by an affine automorphism of `Z/91Z` to

```text
G={0,1,...,25}.                                                  (5)
```

Use CRT coordinates

```text
h=s mod 13,                    r=s mod 7.                        (6)
```

Each target row contains the two guard integers `h` and `h+13`.
Since `13==-1 mod 7`,

```text
G_h={h,h-1}                 in F_7.                              (7)
```

Every unit-step ordinary progression of length thirteen meets each
target row once. Every exact-depth top blocker has normalized speed
`13w`, with `7 not|w`, and therefore occupies one vertical source
column

```text
B_j={r_j}                    in every F_13 row.                  (8)
```

Two blocker labels may have the same source column. The finite cover
only sees the source multiset through its support, but the physical
packet must retain the labels separately.

The affine guard gauge has the form

```text
s -> epsilon u_0 s+kappa(Y),             epsilon in {+-1}.      (8a)
```

Its multiplier is a unit modulo `91`. It therefore sends unit AP13
supports to unit AP13 supports and every vertical modulo-seven
blocker coset to another vertical blocker coset. The translation
`kappa(Y)`, including the common integer carry from the quotient-root
coordinate, translates all intrinsic centres equally and cancels
from every centre difference. Reversal is retained by using directed
symmetric centre banks. Thus normalization loses neither a blocker
source possibility nor a repeated-pair centre difference.

## 2. The exact hole-to-duplicate flow

Let `m_Y(h,r)` be the multiplicity of the guard and five ordinary
unit masks, before adding the top blockers. Row incidence is exactly

```text
sum_(r in F_7)m_Y(h,r)=2+5=7.                                  (9)
```

Define

```text
d_Y(h,r)=1-m_Y(h,r).                                           (10)
```

Then `d_Y(h,-)` is an integral zero-sum vector. Because the complete
top cover can repair a hole only with an exact-depth blocker,

```text
d_Y(h,r)>0        implies        r in S_Y,                       (11)
```

where `S_Y` is the set of one or two blocker source columns.

For one source `r_0`, equation (9) gives the exact dichotomy

```text
d_Y(h,-)=0,

or

d_Y(h,-)=e_(r_0)-e_(r_1(h)),       r_1(h)!=r_0.                 (12)
```

Thus every essential row has an intrinsic arrow from its actual
blocker source to its unique duplicate.

For two distinct sources, a row has zero, one, or two holes. The
exact atlas proves the stronger finite law

```text
m_Y(h,r) in {0,1,2}.                                           (13)
```

Hence:

- zero holes means a flat row;
- one source hole is paired with one duplicate; and
- two holes are the two distinct sources and are balanced by two
  distinct duplicate sites.

No triple collision survives. At most two rows in one stalk have two
holes. That upper bound is sharp only for source pairs `{2,5}` and
`{2,6}` in the normalized gauge.

If the two blocker labels have coincident sources, the base defect is
again (12). One blocker quantum repairs its possible hole and the
second is a retained co-owner/surplus label; it must not be collapsed
into a fictitious second arrow.

## 3. The vertical alternative is empty

Suppose

```text
d_Y(h,r)=b(r)                         for every h.                (14)
```

If `b` were nonzero, its zero sum would force `b(r_0)>0` for some
`r_0`. Equation (11) makes `r_0` a blocker source. But row

```text
h=r_0 in {0,...,6}
```

has `r_0 in G_h` by (7). Therefore

```text
m_Y(h,r_0)>=1,

b(r_0)=1-m_Y(h,r_0)<=0,                                      (15)
```

a contradiction. Equivalently, off the source columns positivity is
forbidden, while each source is nonpositive on one guard row; hence
`b<=0`, and its zero sum gives `b=0`.

Thus

```text
d_Y is h-independent        iff        d_Y=0.                   (16)
```

This proof is uniform for one source, two distinct sources, and two
coincident labelled sources. It does not use the finite enumeration.

## 4. A genuine mixed quotient character

Fix primitive roots `zeta_13,zeta_7` and define the discrete stalk
transform

```text
dtilde_Y(alpha,beta)
 =sum_(h in F_13)sum_(r in F_7)
   d_Y(h,r) zeta_13^(-alpha h) zeta_7^(-beta r).                (17)
```

The row-sum law makes every coefficient with `beta=0` vanish. If all
coefficients with

```text
alpha!=0,                    beta!=0                            (18)
```

also vanished, Fourier inversion in the `h` coordinate would make
each nontrivial `F_7` coefficient constant in `h`. Inversion in `r`
would then make `d_Y` independent of `h`. Equations (16)--(18)
therefore give

```text
d_Y!=0
  implies
dtilde_Y(alpha,beta)!=0
for some alpha!=0,beta!=0.                                    (19)
```

Under CRT duality, the `72=(13-1)(7-1)` pairs in (18) are exactly the
unit characters modulo `91`.

THM-2435 proves

```text
mu(P minus mathcal E)>=(1+k)/7,                                (20)
```

where `mathcal E` is the flat six-unit tiling locus. Every parent in
the difference has nonzero defect. The union of the `72` measurable
nonvanishing loci in (19) contains that difference. Hence one fixed
mixed quotient character is nonzero on parent mass at least

```text
(1+k)/(7*72).                                                   (21)
```

The typed floors are

| shape entering the theorem | fixed mixed-mode parent mass |
|:---|---:|
| `M=0,(1,5,1,8)` (`k=1`) | `>=1/252` |
| `M=0,(2,5,2,9)` (`k=2`) | `>=1/168` |
| `M>0,(2,5,1,8)` (`k=2`) | `>=1/168` |

Thus the uniform floor is `1/252`.

**Current upgrade (THM-2506).** Rational Galois transitivity makes all `72`
mixed modes nonzero at every essential parent.  Thus each preassigned mode has
nonvanishing locus exactly `P minus mathcal E`, with floors `2/7` for `k=1`
and `3/7` for `k=2`.  Equations (21) and the table remain valid but are
strictly superseded.

This is a **discrete signed quotient-stalk mode**. It is linear in
the guard/ordinary incidence defect and retains both prime
coordinates. It is not yet the Fourier coefficient of a fixed
physical Boolean packet, a transported Abel endpoint, or a canonical
relation current. At `M>0`, physical pullback still has the
`7^M`-ancestor kernel identified by THM-2435.

## 5. Complete one-source atlas

Write

```text
A(s,d)={s+jd mod 91:0<=j<=12},             gcd(d,91)=1.          (22)
```

There are exactly `3,276` unoriented supports. For each fixed source
column `r_0`, enumerate unordered five-support families for which

```text
G union A_1 union ... union A_5 union B_(r_0)=Z/91Z.             (23)
```

The exact results are:

| `r_0` | covers | nonflat | defect words | arrow-count histogram | affine arrow words | step spectra |
|---:|---:|---:|---:|:---|---:|---:|
| 0 | 264 | 202 | 81 | `0:62,1:103,2:72,3:24,4:3` | 232 | 21 |
| 1 | 257 | 195 | 74 | `0:62,1:99,2:69,3:24,4:3` | 232 | 18 |
| 2 | 229 | 167 | 64 | `0:62,1:92,2:57,3:16,4:2` | 218 | 19 |
| 3 | 257 | 195 | 74 | `0:62,1:99,2:69,3:24,4:3` | 232 | 18 |
| 4 | 264 | 202 | 81 | `0:62,1:103,2:72,3:24,4:3` | 232 | 21 |
| 5 | 679 | 617 | 232 | `0:62,1:200,2:230,3:137,4:44,5:6` | 481 | 22 |
| 6 | 679 | 617 | 232 | `0:62,1:200,2:230,3:137,4:44,5:6` | 481 | 22 |

A defect word records, in target-row order, either a flat symbol or
the target `r_1(h)` in (12). `Affine arrow words` means that this
partial function agrees with `ah+b mod 7` on its domain; it is a
classification statistic, not an affine symmetry assertion.

The strongest one-source refinement is

```text
h -> r_1(h) is injective on the essential rows.                 (24)
```

No duplicate target repeats. Thus a `b=1` defect is a partial
injection into the six-spoke star centred at `r_0`, not an arbitrary
row tournament.

Let `chi_7` be the quadratic character. The three positive and three
negative values split `F_7^*` into two triples. Consequently, if
there are at least four arrows, injection forces the differences

```text
r_1(h)-r_0
```

to contain both `chi_7` signs. This is an honest Fano/`chi_7`
edge-address refinement. It does not assert Fano-plane symmetry or
orient rows with ties.

The guard stabilizer is generated by

```text
s -> 25-s,                    r -> 4-r,             h -> 12-h.  (25)
```

Across all seven sources there are

```text
2,629 cover assignments,
15 reflection-fixed assignments,
1,322 reflection orbits.                                      (26)
```

Every source contains the same `62` flat THM-2430 tilings.

## 6. Repeated steps and the positive-depth contradiction

Every one of the `2,629` one-source covers repeats an unsigned
ordinary step. Since the five normalized speed residues are fixed
before `Y` varies, their complete unsigned step multiset is fixed.
The existence of any parent cover therefore shows that this fixed
multiset has a repeated value. Choose two fixed physical labels in
one repeated class once and for all. If the class has more than two
labels, the computational banks below include every pair in that
class and hence overestimate the bank of the chosen pair.

For a repeated unsigned step `d`, let `S_d` be the union, over all
seven source atlases, of directed intrinsic-centre differences.
Normalize by `d` as in THM-2431:

```text
Q_d=d^(-1)S_d mod 91.                                          (27)
```

The sharp rounding-error locus is the union of numerator cells

```text
C_d=Q_d union (Q_d-1).                                         (28)
```

The exact banks are:

| `d` | `|S_d|` | `|C_d|` | coarse three-error cells |
|---:|---:|---:|---:|
| 1 | 72 | 75 | 76 |
| 2 | 52 | 57 | 62 |
| 3 | 14 | 20 | 26 |
| 4 | 12 | 15 | 18 |
| 5 | 18 | 22 | 26 |
| 44 | 12 | 14 | 16 |
| 45 | 16 | 20 | 24 |

The source column can vary with `Y`. This is why `S_d` is the union
over all seven normalized columns; equation (8a) proves that these
are all, and only, the possible blocker gauges. The union is taken
before the measure bound, so every parent is covered without choosing
a source or orientation parentwise.

The derivation of THM-2431 equations (20)--(23) uses only:

- the fixed signed speed relation for the repeated pair;
- its intrinsic centre difference belonging to the finite bank; and
- nearest-integer rounding.

It does not use disjointness of the other masks after the centre bank
has been supplied. Therefore it applies to the punctured atlas and
gives

```text
mu(P)<=|C_d|/91<=75/91.                                        (29)
```

For completeness, write the fixed signed normalized steps as

```text
a_j=sigma a_i,                       sigma in {+-1}.
```

The corresponding normalized physical speeds satisfy

```text
u_j=sigma u_i+91n.                                             (29a)
```

The atlas contains no repeated AP support. Hence `n!=0`: at `n=0`
the two norm-symmetric physical masks coincide on every common
stalk. With

```text
R(x)=floor(x+1/2),             rho_i(Y)=u_iY-R(u_iY),
```

the physical, pre-normalization rounding equation is

```text
R(u_jY)-sigma R(u_iY)
 =91nY+sigma rho_i(Y)-rho_j(Y),        |rho_i|,|rho_j|<1/2.     (29b)
```

After the common `kappa(Y)` cancels, membership of the directed centre
difference in `S_d` turns (29b) into

```text
dist_(R/91Z)(91{nY},Q_d)<1.                                   (29c)
```

The numerator cell is therefore in `Q_d` or `Q_d-1`. Since
`n!=0`, multiplication by `n` preserves Haar measure, proving (29)
directly in the physical parent coordinate.

For the positive-depth shape, `k=2`. Equations (3) and (29) say

```text
78/91=6/7<=mu(P)<=75/91,                                       (30)
```

impossible. Hence

```text
M>0,(k,t,b,W)=(2,5,1,8) is empty.                              (31)
```

The exact slack is

```text
6/7-75/91=3/91.                                                (31a)
```

This is near-sharp, not a coarse rounding accident. As positive
finite controls, every source atlas contains all `62` flat tilings
and hundreds of nonflat covers; the step-one union itself occupies
`75` numerator cells. None of these finite controls is asserted to
be a global scalar-cover packet.

For the one-source `M=0` shape, equation (3) gives

```text
mu(P)>=5/7=65/91.                                               (32)
```

If its fixed step spectrum has any repeated `d!=1`, the table above
already gives

```text
mu(P)<=57/91,
```

contradicting (32). It remains to condition on a fixed spectrum whose
only repeated value is one. This conditioning is lawful because the
entire five-step multiset is fixed speed-residue data. The exact
atlas leaves only:

| fixed unsigned step spectrum | cover assignments | step-one centre differences | sharp cells |
|:---|---:|---:|---:|
| `(1,1,1,1,1)` | 144 | 40 | 46 |
| `(1,1,1,1,30)` | 106 | 36 | 41 |
| `(1,1,1,1,45)` | 262 | 48 | 57 |

For each row, the centre bank includes every pair among all
step-one supports. It therefore overestimates the bank for any fixed
labelled pair when four or five labels have step one. Nevertheless
the largest sharp bank has only `57` cells. Thus (32) again
contradicts

```text
mu(P)<=57/91.                                                   (33)
```

Consequently

```text
M=0,(k,t,b,W)=(1,5,1,8) is empty.                              (34)
```

The exceptional one-source lane has exact slack

```text
65/91-57/91=8/91.                                             (34a)
```

## 7. Complete two-source atlas

For two distinct blocker source columns, the exact table is:

| sources | covers | exact-both-source covers | defect words | covers with a double-hole row | max total holes | step spectra |
|:---|---:|---:|---:|---:|---:|---:|
| 0,1 | 1019 | 560 | 448 | 10 | 5 | 21 |
| 0,2 | 602 | 171 | 258 | 12 | 6 | 23 |
| 0,3 | 717 | 258 | 315 | 10 | 5 | 24 |
| 0,4 | 912 | 446 | 460 | 0 | 6 | 25 |
| 0,5 | 2600 | 1719 | 1500 | 10 | 9 | 29 |
| 0,6 | 2168 | 1287 | 1021 | 41 | 6 | 25 |
| 1,2 | 944 | 520 | 404 | 0 | 5 | 19 |
| 1,3 | 583 | 131 | 234 | 18 | 6 | 21 |
| 1,4 | 717 | 258 | 315 | 10 | 5 | 24 |
| 1,5 | 1790 | 916 | 882 | 61 | 7 | 27 |
| 1,6 | 2096 | 1222 | 1042 | 56 | 7 | 27 |
| 2,3 | 944 | 520 | 404 | 0 | 5 | 19 |
| 2,4 | 602 | 171 | 258 | 12 | 6 | 23 |
| 2,5 | 1677 | 831 | 798 | 61 | 6 | 26 |
| 2,6 | 1677 | 831 | 798 | 61 | 6 | 26 |
| 3,4 | 1019 | 560 | 448 | 10 | 5 | 21 |
| 3,5 | 2096 | 1222 | 1042 | 56 | 7 | 27 |
| 3,6 | 1790 | 916 | 882 | 61 | 7 | 27 |
| 4,5 | 2168 | 1287 | 1021 | 41 | 6 | 25 |
| 4,6 | 2600 | 1719 | 1500 | 10 | 9 | 29 |
| 5,6 | 10029 | 8733 | 5098 | 82 | 8 | 32 |

`Exact-both-source` means that the union of hole residues over all
thirteen rows is exactly the displayed pair. A defect word records
the unordered row flow

```text
[hole set > duplicate multiset].                               (35)
```

The maximum of two double-hole rows is attained. A smallest canonical
witness for sources `{2,5}` has ordinary supports recorded by
`centre:unsigned-step`

```text
33:1, 47:1, 58:1, 72:1, 86:1
```

and double rows

```text
h=0: [25>03],                  h=1: [25>14].                    (36)
```

The `{2,6}` witness printed by the companion is its reflection.

There are

```text
38,750 distinct-two-source cover assignments,
100 reflection-fixed assignments,
19,425 reflection orbits.                                      (37)
```

Including the seven coincident-source configurations gives

```text
41,379 source-multiset assignments,
26,535 distinct ordinary support configurations,
115 reflection-fixed assignments,
20,747 reflection orbits.                                      (38)
```

Every two-source cover also repeats an unsigned step. The union
centre banks are:

| `d` | `|S_d|` | `|C_d|` |
|---:|---:|---:|
| 1 | 80 | 81 |
| 2 | 66 | 71 |
| 3 | 56 | 62 |
| 4 | 22 | 25 |
| 5 | 26 | 30 |
| 6 | 6 | 7 |
| 30 | 8 | 10 |
| 31 | 4 | 5 |
| 36 | 4 | 6 |
| 44 | 20 | 22 |
| 45 | 38 | 44 |

For the two-source `M=0` shape, `k=2` and

```text
mu(P)>=6/7=78/91.                                               (39)
```

Every repeated bank except `d=1` has at most `71` cells, so a
spectrum with such a repeated value contradicts (39). Condition on a
fixed spectrum whose only repeated value is one. This time coincident
as well as distinct source columns must be included; the two physical
blocker labels can occupy the same quotient column. The complete
exceptional atlas is:

| fixed unsigned step spectrum | cover assignments | step-one centre differences | sharp cells |
|:---|---:|---:|---:|
| `(1,1,1,1,1)` | 2681 | 60 | 66 |
| `(1,1,1,1,15)` | 85 | 30 | 35 |
| `(1,1,1,1,18)` | 13 | 26 | 38 |
| `(1,1,1,1,30)` | 1611 | 54 | 59 |
| `(1,1,1,1,45)` | 4314 | 62 | 67 |
| `(1,1,1,18,45)` | 16 | 10 | 16 |
| `(1,1,23,30,45)` | 144 | 6 | 8 |

Again every possible pair of step-one supports is included, so these
are safe overestimates for a fixed labelled pair. The maximum is

```text
mu(P)<=67/91<78/91,                                             (40)
```

contradicting (39). Hence

```text
M=0,(k,t,b,W)=(2,5,2,9) is empty.                              (41)
```

The exceptional two-source lane has exact slack

```text
78/91-67/91=11/91.                                            (41a)
```

## 8. Consequence and scope

Equations (31), (34), and (41) empty the complete residual entering
this theorem:

```text
nu_7(c_3)>M        is impossible.                               (42)
```

THM-2426 had already proved that every deep-`c_3` survivor lies in
the top-guard side classified by THM-2427, and THM-2431--2432 reduce
that side exactly to the three shapes (2). Therefore (42) is the
complete deep-`c_3` exclusion, not merely a subcase.

The mixed-mode theorem remains structurally informative: before the
rounding contradiction, every essential parent has a genuine unit
quotient character, and one fixed such character is nonzero on
parent mass at least `1/252`. That mode is not needed for the
repeated-step closure above. To reuse it outside the now-empty deep
branch would require a lawful transplant retaining the blocker label
and physical endpoint/current.

What remains is the complementary `nu_7(c_3)<=M`
noncirculant-graft branch. Its missing theorem is owner-conditioned
graft survival, or an alternative scalar-cover decrement retaining
the same owner data.

THM-2439's degree-ten marker boundary is orthogonal to this closure.
The repeated pair is chosen among fixed physical ordinary-speed
labels after conditioning on their fixed unsigned step spectrum; the
proof never reconstructs THM-2435's lexicographic blocker marker from
self-Gram data. Thus the punctured atlas bypasses, rather than solves,
the marker-replica obstruction.

The theorem does not make a tournament out of tied rows, identify a
unique owner in the two-blocker word, remove one of the `165` scalar
rows, close the complementary `nu_7(c_3)<=M` noncirculant-graft
branch, or prove LRC(14).

## 9. Exact companion and hostile controls

Compile and run:

```bash
c++ -O0 -std=c++20 \
  04-computation/lrc14_punctured_91_stalk_mixed_mode_thm2436.cpp \
  -o /tmp/lrc14_thm2436_O0

c++ -O3 -std=c++20 \
  04-computation/lrc14_punctured_91_stalk_mixed_mode_thm2436.cpp \
  -o /tmp/lrc14_thm2436_O3
```

Both executables must reproduce

```text
05-knowledge/results/lrc14_punctured_91_stalk_mixed_mode_thm2436.out
```

byte-for-byte. The companion:

- reconstructs all `3,276` unoriented unit AP13 supports;
- enumerates all `28` blocker-source multisets by an order-free
  missing-cell recursion;
- independently reconstructs the `62` flat tilings by exact
  Algorithm X;
- verifies exact equality under reflection and redundant
  source-superset reconstruction;
- reconstructs every row multiplicity, defect word, repeated step,
  centre bank, and double-hole witness;
- checks that no nonzero `h`-independent defect survives;
- checks one-source target injection and the four-arrow
  `chi_7` sign law; and
- raises on every count, orbit, bank, or structural drift.

Optimized and unoptimized transcripts are required to agree.

The independent hostile audit rebuilt the companion separately at
`-O0` and `-O3`; both transcripts byte-match the stored output and have
SHA-256

```text
c2659364e9d312a670f2606cde8eeb8f9ef5d82c062cfb44837970047f717f8b.
```

It also checked the two places where a finite atlas can silently lose
the physical quantifiers.

First, the cover recursion is complete without assuming disjointness:
at every depth it branches through every AP containing a missing pivot;
the row and cardinality tests are necessary conditions only; and at
depth four it includes every AP containing the whole remaining set
(or every AP when nothing remains). Sorting and set insertion remove
order but not a support multiset. The later zero repeated-support check
then proves that no physically relevant completion used a duplicate
support. The separately implemented flat Algorithm-X search, exact
source-superset reconstruction, and reflection transport are hostile
controls on this recursion rather than restatements of it.

Second, the repeated pair is selected before the parent phase varies.
The five unsigned normalized steps are fixed physical speed-residue
data. Once their fixed spectrum has a repeated class, choose two fixed
physical labels in that class. The centre bank is then formed by
unioning over every source gauge, every admissible cover, both
directions, and every same-step support pair. It is therefore an
overestimate for that fixed labelled pair even when the support
assignment varies with `Y`. The common translation `kappa(Y)` cancels,
while reversal is already in the directed bank. The resulting signed
speed relation has fixed nonzero integer lift `n`; multiplication by
`n` preserves Haar measure. Conditioning on the entire fixed step
spectrum in the two exceptional step-one arguments is consequently
lawful, and their `57/91` and `67/91` caps apply to the whole parent,
not to a parentwise-selected subcell.

The dependency-free Python referee independently reconstructs the
load-bearing conditioned banks without importing the C++ atlas or a
stored solution. It represents a step-one AP13 as a cyclic interval of
radius six, indexes every other AP by its centre, and uses a separate
canonical dynamic-pivot owner recursion. It reproduces the full
one-source exceptional table

```text
(1,1,1,1,1):   144 assignments, 40 differences, 46 cells,
(1,1,1,1,30):  106 assignments, 36 differences, 41 cells,
(1,1,1,1,45):  262 assignments, 48 differences, 57 cells.
```

On all `28` coincident/distinct two-label source multisets it also
reconstructs the five critical spectra with at least four step-one
labels. In particular the extremal `(1,1,1,1,45)` row has exactly
`4,314` assignments, `62` directed centre differences, and `67` sharp
cells. Normal and `-O` runs both byte-match

```text
05-knowledge/results/lrc14_punctured_91_fixed_spectrum_referee_thm2436.out
```

with output SHA-256

```text
bba34d1c67d68b992bf274a22c0e92f29b1c80c58ce15d9a7ca631b0cceff584.
```

Finally, the audit traced the residual list through hostile-audited
THM-2426, THM-2427, THM-2431, and THM-2432 and the parent floors through
hostile-audited THM-2435. Hence the three shapes in (2) exhaust the
deep-`c_3` branch. QED.

## Second independent one-source referee (kind-pasteur-2026-07-26-S131)

A second code path, written before the C++ companion landed and
committed after it, independently confirms the one-source layer:

```text
04-computation/lrc14_punctured_91_stalk_repeated_step_thm2436.py
05-knowledge/results/lrc14_punctured_91_stalk_repeated_step_thm2436.out
sha256 (LF): 6968b676edabe504d1587f82154827fa36036b0501b6e66e4ddc8cf1c7a602eb
             bd7ebd5c0bee18e75ffff7c7e8dc0cd6e13a9a0063ed9197cee6180fed547f13
```

It is a genuinely different algorithm (memoized existence set-cover
with dynamic uncovered-point pivots and waste budget, per
guard-stabilizer orbit representative `r in {0,1,2,5}`), not a port.
It confirms exactly: the `3276` unit-AP universe; the reflection laws;
the punctured targets `56/56/56/55`; that no one-source cover uses
five distinct unsigned steps (`distinct_cover=0`, the quantifier
behind the fixed repeated pair); that centre difference `0` never
occurs (`n != 0`); and the complete one-source repeated-pair bank --
all seven step classes with raw/sharp/coarse sizes identical to the
C++ table, giving the same `75/91` cap and
`parent_mass_78_over_91_surviving_steps=NONE`. Because the Python
bank admits any three completing masks (a formal relaxation), the
exact agreement is a tightness result, not a tautology.

Two caveats and two free consequences for future auditors:

1. Column order differs: the Python prints `step:raw/coarse/sharp`,
   the C++ prints `step raw sharp coarse`. A naive diff misreads
   every row.
2. The Python covers the one-source layer only; the two-source
   atlas, the fixed-spectrum exceptional tables, and the defect/row
   layer rest on the C++ and its promotion audit alone.
3. The exact d=1 structures are `Q_1 = {9,...,82} \ {22,69}` and
   `C_1 = {8,...,82}`, an interval of `75`; the `M>0` contradiction
   margin is exactly the 16-cell complement of that interval minus
   the 13-cell floor complement, i.e. `3` cells.
4. The depth-limited recursion in both engines proves, though
   neither states it: no cover of any punctured one-source stalk by
   four or fewer ordinary words exists -- five distinct supports are
   necessary as well as sufficient.
