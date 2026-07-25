---
id: THM-2287
title: "Anchored scalar rank-six Plucker flag and finite label atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every one of the 120
  interior first-depth-one scalar profiles has six relations of heights
  at most 9841, 4921, 7381, 11072, 16608, and 24912 whose reductions
  modulo thirteen have rank six. A nested six-by-six pivot contains c_1,
  avoids the primitive-pair partner, contains at least three guard/unit
  labels, and belongs to an exact 21-pattern label atlas. Its determinant
  is a thirteen-unit of absolute value at most
  15114636407967321147310080. The construction supplies 13^(6n)
  lossless relation addresses at every thirteen-adic depth, a
  three-column residue quotient, and dark-column/all-unit address
  alternatives. A separate unanchored six-row packet of total height
  1453 holds on all 165 profiles. No scalar profile is excluded, no
  Fourier amplitude or owner service is selected, and LRC(14) remains
  open.
source: codex-2026-07-25-anchored-scalar-rank-six-flag
depends_on:
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
  - THM-2279-shallow-blocker-anchored-relation-minor
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-2295-weighted-basis-safe-floor-and-scalar-rank-five-harvest
  - THM-2298-weighted-rank-five-facet-deficit-and-scalar-rank-six-harvest
related:
  - THM-2069-k-deletion-code-cogirth-crt-wheel
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
script: 04-computation/lrc14_anchored_scalar_rank_six_flag_thm2287.py
output: 05-knowledge/results/lrc14_anchored_scalar_rank_six_flag_thm2287.out
script_sha256: ba71f348879de95f1a749860b38a18dd812f8c0c54be385c3ec8104d4f47a3b0
output_sha256: cb43d4b51611b884fcc4e035171550338494ee02d17e0621129bc3e20f7e8dbe
hash_basis: working-tree bytes (LF)
---

# THM-2287 -- anchored scalar rank-six Plucker flag

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2266 supplies one primitive relation whose reduction modulo thirteen
is pinned to the shallow blocker `c_1`. THM-2295 supplies five-dimensional
rational relation rank at height `196`, and THM-2298 supplies a sixth
direction at height `526`. Neither rank theorem places its residues.

Applying THM-2284's centered pivot-extension lemma five times joins the
three inputs without replacing the anchor:

```text
primitive c_1 row
  + four rationally fresh height-196 seeds
  + one rationally fresh height-526 seed
  + centered residue cancellation and division by thirteen
  -> a nested c_1-anchored rank-six flag over F_13.             (1)
```

Every new pivot column can avoid the primitive-pair partner. This leaves an
exact `21`-pattern terminal atlas, a three-column quotient, and a finite
dark-column versus all-unit address alternative.

## 1. Scalar input and the primitive anchor

Use the scalar row

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3)                     (2)
```

in one of the `120` interior first-depth-one profiles

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

3<=b<=c-2,                    5<=c<=19.              (3)
```

Put

```text
Lambda_*={x in Z^9:x.w_*=0}.                         (4)
```

This kernel is saturated in `Z^9`.

THM-2266 gives a primitive relation

```text
p in Lambda_*                                        (5)
```

supported on

```text
{c_1,z},       z in {H,q_1,...,q_5},                 (6)
```

such that

```text
||p||_infinity<=9841,
0<|p_(c_1)|<=757,
13 does not divide p_(c_1),
13 divides p_z.                                      (7)
```

The scalar value at `z` is a thirteen-unit. Thus `p mod 13` is supported
exactly at `c_1` and is nonzero there.

The two bounded-rank inputs are

```text
dim_Q W_196^*(w_*)>=5,               [THM-2295]

dim_Q W_526^*(w_*)>=6,               [THM-2298]      (8)
```

where

```text
W_K^*(w_*)
 =span_Q{x in Lambda_*:||x||_infinity<=K}.           (9)
```

## 2. Five centered pivot extensions

THM-2284's pivot-extension lemma says that if

```text
rho_1,...,rho_d in Lambda_*                          (10)
```

have independent reductions modulo thirteen, heights at most
`H_1,...,H_d`, and `s_0` is a rationally new seed of height at most `S_0`,
then the new row can be chosen to raise the mod-thirteen rank and have
height at most

```text
max(S_0,ceil((H_1+...+H_d)/2)).                      (11)
```

Start with

```text
rho_1=p,                    H_1=9841.                (12)
```

For `d=1,2,3,4`, equation (8) supplies an actual height-`196` relation
outside the current rational `d`-plane: otherwise all bounded generators,
and hence the five-dimensional space `W_196^*(w_*)`, would lie in that
plane. Four applications of (11) give

```text
H_2=4921,
H_3=7381,
H_4=11072,
H_5=16608.                                          (13)
```

These are the rank-five prefix heights from the earlier candidate.

The first five rows span a rational five-plane. Equation (8) now supplies
a height-`526` seed outside it. One more application of (11) gives

```text
H_6
 =max(526,ceil((9841+4921+7381+11072+16608)/2))
 =24912.                                             (14)
```

The five prime-thirteen recurrence slacks are

```text
12H_2-6H_1=6,
12H_3-6(H_1+H_2)=0,
12H_4-6(H_1+H_2+H_3)=6,
12H_5-6(H_1+H_2+H_3+H_4)=6,
12H_6-6(H_1+...+H_5)=6.                             (15)
```

We obtain

```text
rho_1,...,rho_6 in Lambda_*                          (16)
```

with the stated heights and

```text
rank_(F_13)(rho_1,...,rho_6)=6.                     (17)
```

The output rows need not remain in `W_196^*` or `W_526^*`; centered
division preserves the explicit height boxes instead.

## 3. A nested pivot avoiding the pair partner

We construct nested pivot sets

```text
C_d={c_1,k_2,...,k_d},             1<=d<=6,          (18)
```

such that `z notin C_d` and the `d` by `d` minor of the first `d` rows on
`C_d` is a unit modulo thirteen.

For `d=1`, use `C_1={c_1}` and (7). Suppose `C_d` has been chosen for
some `d<6`. Reduce `rho_(d+1)` modulo thirteen against the first `d` rows
so that the residual row `t` vanishes on `C_d`. Since the row rank rises,

```text
t!=0.                                                (19)
```

The residual remains a relation modulo thirteen:

```text
t.w_*=0 mod 13.                                     (20)
```

It cannot be supported only at `z`, because the scalar value at `z` is a
thirteen-unit and (20) would force `t_z=0`. Hence a fresh label

```text
k_(d+1) notin C_d union {z}
```

has

```text
t_(k_(d+1))!=0.                                     (21)
```

Adjoining that column proves the induction. In particular,

```text
C_6={c_1,k_2,k_3,k_4,k_5,k_6}                       (22)
```

is a six-column unit pivot avoiding `z`.

## 4. The finite label and complement atlases

Remove `c_1` and the guard/unit partner `z` from the nine scalar labels.
The remaining bank is

```text
B=({H,q_1,...,q_5} minus {z}) union {c_2,c_3}.       (23)
```

It contains five guard/unit labels and two blocker labels.

The rank-five prefix chooses four labels from `B`. Its exact `35`-pattern
atlas is

```text
10  with c_2,c_3 and two unit labels,
20  with one blocker and three unit labels,
 5  with four unit labels.                           (24)
```

The terminal rank-six pivot chooses five labels from `B`. Its exact
`21`-pattern atlas is

```text
10=binomial(5,3)    with c_2,c_3 and three units,
10=2binomial(5,4)  with one blocker and four units,
 1=binomial(5,5)   with five unit labels.            (25)
```

Thus every terminal pivot contains at least three guard/unit labels outside
`c_1`.

The three nonpivot labels form the complementary atlas:

```text
10  contexts: {z and two unit labels},
10  contexts: {z, one unit label, one blocker},
 1  context:  {z,c_2,c_3}.                           (26)
```

These are exact candidate sets, not claims that every profile realizes
every pattern.

## 5. Determinant bounds

Let `Delta_d` be the nested minor on `C_d`. The actual row `p` is supported
on `{c_1,z}` and every `C_d` avoids `z`, so expansion along the first row
gives

```text
0<|Delta_d|
 <=(d-1)! |p_(c_1)| product_(i=2)^d H_i.             (27)
```

The exact bounds are

```text
0<|Delta_2|<=3725197,

0<|Delta_3|<=54991358114,

0<|Delta_4|<=1826592951114624,

0<|Delta_5|<=121344222928446701568,

0<|Delta_6|
 <=5!*757*4921*7381*11072*16608*24912
  =15114636407967321147310080.                       (28)
```

Every `Delta_d` is a unit modulo thirteen.

## 6. Smith, address, and three-column quotient sidecars

Let

```text
U=Zrho_1+...+Zrho_6 subset Z^9.                      (29)
```

The gcd of the six-by-six minors is the index of `U` in its rational
saturation. Since `Delta_6` is a thirteen-unit,

```text
U has no thirteen-primary Smith obstruction.         (30)
```

For every `n>=1`, the pivot remains invertible modulo `q=13^n`. Hence

```text
(Z/qZ)^6 -> (Z/qZ)^9,

(a_1,...,a_6) |-> sum_(i=1)^6 a_i rho_i             (31)
```

is injective and gives exactly

```text
13^(6n)
```

distinct relation addresses.

The scalar row-height sum is

```text
H_1+...+H_6=74735.                                   (32)
```

Centered representatives give every address an exact integer relation of
height at most

```text
B_n=74735(13^n-1)/2,                                (33)

B_1=448410.                                          (34)
```

There is also an exact residue quotient. For any coefficient vector

```text
v in (Z/qZ)^9,
```

there is a unique `a in (Z/qZ)^6` for which

```text
(v+aR)|_(C_6)=0,                                    (35)
```

where `R` is the six-by-nine matrix with rows `rho_i`. Thus every coset
modulo the relation-row module has a unique residue representative supported
on the three nonpivot columns in (26).

Choose integer lifts of `v` and `a`. The unreduced shift

```text
v -> v+aR
```

preserves the exact scalar frequency because every row of `R` annihilates
`w_*`. Its pivot coordinates are generally multiples of `q`, not literally
zero over the integers. Subtracting those multiples of `q e_i` to make the
integer support literal would change scalar frequency. The quotient is
therefore a three-column **residue** normal form.

Neither relation shifting nor residue normalization preserves factorwise
Fourier coefficients, convolution weights, cancellation, an owner label,
or service by THM-2296's common atom.

## 7. Rank-six Kakeya address alternatives

Let

```text
V_i=((rho_1)_i,...,(rho_6)_i) in F_13^6             (36)
```

be the nine coefficient columns. Call label `i` **dark** if `V_i=0`.
A pivot column cannot be dark, so every dark label lies among the three
nonpivots in (26).

Suppose no label is dark. For a row-combination

```text
t in F_13^6,
```

the coefficient at label `i` is `t.V_i`. Its zero set is one linear
hyperplane of `13^5` points. The raw union bound gives at least

```text
13^6-9*13^5=4*13^5                                  (37)
```

base combinations for which all nine coefficients are thirteen-units.
This is a valid floor, not a sharpness claim.

Every good base combination has exactly `13^(6n-6)` lifts modulo `13^n`,
and (31) makes all lifted addresses distinct. Thus the no-dark branch has
at least

```text
4*13^(6n-1)                                         (38)
```

all-nine-unit relation addresses at depth `n`, each with height bounded by
(33).

An anchored-normalized version is also useful. Impose

```text
(tR)_(c_1)=1 mod 13^n.                              (39)
```

The base affine space has `13^5` points. Each of the other eight coordinate
zero conditions is empty or removes at most `13^4` points, so at least

```text
13^5-8*13^4=5*13^4                                  (40)
```

base addresses are all-unit. At depth `n` this gives at least

```text
5*13^(5n-1)                                         (41)
```

all-unit addresses with `c_1` coefficient congruent to `1 modulo 13^n`.
The centered integer coefficient at `c_1` need not be literally `1`.

THM-2294 proves a sharper essential-arrangement count in rank three.
Equations (37)--(41) deliberately retain only elementary union-bound floors;
the sharper visible-rank-six pivot-torus counts are proved separately in
THM-2301 (`essential-affine-arrangement-and-visible-rank-six-address-bank`).

## 8. A repaired prescribed two-digit alternative

Assume the no-dark branch, so `V_(c_3)!=0`. The first coordinate of
`V_(c_1)` in the row ordering above is the unit `p_(c_1)`, while the first
coordinate of `V_(c_3)` is zero. Hence

```text
V_(c_1),V_(c_3) are linearly independent over F_13. (42)
```

By basis exchange, the rank-six pivot can be chosen to contain both
`c_1,c_3` and still avoid `z`.

Fix any unit

```text
m in Z/13^nZ
```

and impose

```text
(tR)_(c_1)=1,
(tR)_(c_3)=m                         mod 13^n.        (43)
```

The base fibre has `13^4` points. There is one genuine additional
obstruction. A nonzero labelled column may satisfy

```text
V_i=alpha V_(c_1)+beta V_(c_3),

alpha+beta m=0                       mod 13.          (44)
```

Then its coefficient vanishes on the entire prescribed base fibre even
though the column is not dark. Call (44) a **prescribed-digit collapsed
column**.

If no column satisfies (44), each of the other seven zero conditions is
empty or a proper affine hyperplane of at most `13^3` points. The union
bound leaves at least

```text
13^4-7*13^3=6*13^3                                  (45)
```

all-unit base addresses. At depth `n`, at least

```text
6*13^(4n-1)                                         (46)
```

all-unit relation addresses realize the two congruences in (43), with the
same height bound (33).

The collapsed-column branch is distinct from a dark column. Without it the
two-digit assertion is false: the allowed nonzero column

```text
V_i=V_(c_3)-mV_(c_1)
```

vanishes identically on (43).

## 9. The THM-2293 spectral augmentation

THM-2293 supplies on all `120` interior profiles a terminal spectral
coefficient row

```text
a_0=m e_(c_3),

0<|m|<=2533,                 gcd(m,13)=1.            (47)
```

This is not an integer relation:

```text
a_0.w_*=m c_3!=0.                                   (48)
```

If `c_3 in C_6`, the relation pivot is already anchored at both `c_1` and
the spectrally named deepest blocker. This is an annotation, not an
additional rank direction.

If `c_3 notin C_6`, append `a_0` as a seventh row and `c_3` as a seventh
pivot column. The augmented determinant is block triangular:

```text
Delta_7^augmented=+/-m Delta_6,

Delta_7^augmented!=0 mod 13,

0<|Delta_7^augmented|
 <=2533*15114636407967321147310080
  =38285374021381224466136432640.                    (49)
```

This is a rank-seven **relation/spectral** flag, not rank-seven relation
rank. Modulo `13^n`, the augmented rowspace has a two-nonpivot-column normal
form. Its first six row shifts preserve exact scalar frequency; shifting by
`a_0` does not. The augmentation neither preserves Fourier amplitude nor
proves that THM-2293's terminal atom is incident with THM-2296's common
ancestry atom.

Using the particular `m` from (47) in Section 8 gives a lawful
prescribed-digit alignment only on the noncollapsed branch. It still does
not supply atom incidence or amplitude preservation.

## 10. Exact base-phase hostile boundary

THM-2299 proves that even the owner, expiration clock, target blocker,
one-sheet root word, and an anchored unit minor do not determine the signed
terminal base phase. Its exact local witness has scalar speeds

```text
H=1,

(q_1,...,q_5)=(4,2,3,6,10),

(c_1,c_2,c_3)=(13,2197,742586).                     (50)
```

This row already has a `c_1`-anchored rank-six packet of height `13`.
Take

```text
r_i=e_(q_i)-q_i e_H,                1<=i<=5,

p_0=13e_(q_1)-4e_(c_1).                            (51)
```

All six rows are exact relations. On the columns

```text
(q_1,q_2,q_3,q_4,q_5,c_1)
```

their determinant is

```text
-4!=0 mod 13.                                       (52)
```

Nevertheless THM-2299's two-component source set and prescribed
`c_1 -> c_2` handoff satisfy the exact cancellations

```text
(1_F)_hat(4)=0,

(1_E)_hat(52)=0,

(W_4)_hat(0)=0.                                     (53)
```

The witness is not a global scalar cover and does not contradict the present
theorem. It is the sharp scope boundary: even an anchored rank-six unit
pivot of far smaller height than (14) cannot by itself align the terminal
base phase. A phase/component or amplitude-incidence sidecar remains
load-bearing.

## 11. Fixed-section original-row lift

THM-2203 lifts scalar relation coefficients on the fixed nine-coordinate
section by

```text
(x_H,x_rest) |->(2x_H,x_rest).                       (54)
```

The map is integral and injective. It preserves mod-thirteen rank because
`2` is a unit. The lifted row heights are at most

```text
(19682,9842,14762,22144,33216,49824).                (55)
```

At most one pivot column is `H`, so

```text
0<|Delta_6^original|
 <=30229272815934642294620160.                       (56)
```

The fixed-section address height is at most

```text
B_n^original=74735(13^n-1),

B_1^original=896820.                                 (57)
```

## 12. Optional anchored height-20 branch

THM-2275 supplies a scalar relation `r` of height at most `20`. Suppose,
as an additional branch hypothesis, that

```text
r in Qp.                                             (58)
```

THM-2279 then proves that the same primitive `c_1` anchor satisfies

```text
||p||_infinity<=20,
0<|p_(c_1)|<=20.                                    (59)
```

Four height-`196` extensions and the final height-`526` seed give

```text
(H_1,...,H_6)=(20,196,196,206,309,526),

H_1+...+H_6=1453.                                    (60)
```

The same partner-avoiding pivot and `21`-pattern atlas apply. The terminal
determinant satisfies

```text
0<|Delta_6|
 <=5!*20*196*196*206*309*526
  =3086987197593600.                                 (61)
```

There are `13^(6n)` addresses of scalar height at most

```text
1453(13^n-1)/2,

B_1=8718.                                            (62)
```

The fixed-section determinant is at most

```text
6173974395187200,
```

and the depth-one fixed address height is at most

```text
17436.                                               (63)
```

This sharpening is conditional on (58). It is not used in the uniform
anchored theorem.

## 13. All-165 unanchored contrast

THM-2301 (`essential-affine-arrangement-and-visible-rank-six-address-bank`)
proves a shorter unanchored packet on every one of the `165` live scalar
profiles. We record its exact packet and blocker atlas here only to expose
the contrast with the anchored construction above. THM-2275 supplies a
nonzero relation of height at most `20`.
Divide it by its integer content. The primitive row still has height at
most `20` and has nonzero reduction modulo thirteen.

THM-2295 then supplies four rationally fresh height-`196` seeds, and
THM-2298 supplies a sixth height-`526` seed. The same pivot-extension lemma
gives

```text
(20,196,196,206,309,526),            sum=1453,       (64)
```

with mod-thirteen row rank six on all `165` profiles. The packet has
`13^(6n)` addresses of height at most

```text
1453(13^n-1)/2,

B_1=8718,                    B_1^original=17436.     (65)
```

Without the anchored cofactor expansion, a terminal minor has the
Leibniz bound

```text
6!*20*196*196*206*309*526
 =18521923185561600.                                 (66)
```

This packet still retains one blocker label. Let `R` be its six-by-nine
matrix and reduce `Rw_*=0` modulo thirteen. The three blocker speeds
vanish, while all six guard/unit speeds are units, so the six unit columns
obey a dependence with all six scalar weights nonzero. Therefore the
all-unit six-column set cannot be a pivot. The exact possible pivot atlas is

```text
18=3binomial(6,5) pivots with one blocker and five units,
45=binomial(3,2)binomial(6,4) with two blockers and four units,
20=binomial(6,3) with three blockers and three units,

18+45+20=83.                                         (67)
```

The corresponding complement triples contain respectively

```text
two blockers and one unit,
one blocker and two units,
three units.                                         (68)
```

This all-`165` corollary loses the prescribed `c_1` anchor, the partner
`z`, and the `21`-pattern atlas. It complements rather than replaces the
anchored `120`-profile theorem.

## 14. Frontier effect and stopping boundary

The proved chain is

```text
120 interior scalar profiles
  -> primitive c_1 anchor with a unit mod-thirteen coefficient
  -> rational rank five at height 196
  -> rational rank six at height 526
  -> five centered pivot extensions
  -> nested anchored rank-six flag
  -> exact 21-pattern pivot / three-column complement atlas
  -> full mod-13^n quotient and address banks
  -> dark column or many all-unit relation words.    (69)
```

The connection and loss ledger is

```text
source:
  THM-2266's primitive shallow pair, THM-2295/2298's bounded
  rational rank, and THM-2284's pivot-extension lemma;

target:
  a bounded c_1-anchored rank-six Plucker flag, finite label atlas,
  residue quotient, Kakeya bank, and spectral augmentation;

map:
  repeatedly choose a rationally fresh bounded seed, cancel its centered
  old residue, divide the complete relation by thirteen, and choose a
  fresh pivot column away from the primitive partner;

preserved:
  integrality, scalar relation equation, rational independence,
  mod-thirteen rank, c_1 anchor, explicit heights, determinant bounds,
  exact scalar frequency under relation shifts, and all thirteen-adic
  relation addresses;

destroyed:
  original short seeds, descent lengths, exact outside labels within
  the finite atlas, signs, factorwise Fourier amplitudes, root history,
  expiration ancestry, and owner transition;

stopping objects:
  one of three labelled dark nonpivot columns, a prescribed-digit
  collapsed column, or a nonzero but amplitude-blind quotient word;

needed sidecar:
  prove that a quotient/address word has the same nonzero Fourier
  amplitude and owner service as THM-2296's prescribed-expiration atom. (70)
```

No scalar profile is excluded. The anchored theorem applies only to the
`120` interiors covered by THM-2266. Rank six, residue normal forms, and
the address banks do not prove a successor cover, an ancestry return, or
LRC(14).

Reproduction:

```bash
python3 04-computation/lrc14_anchored_scalar_rank_six_flag_thm2287.py
python3 -O 04-computation/lrc14_anchored_scalar_rank_six_flag_thm2287.py
```
