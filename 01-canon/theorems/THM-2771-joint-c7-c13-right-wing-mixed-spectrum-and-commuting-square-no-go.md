---
id: THM-2771
title: "Joint C7-by-C13 right-wing spectrum, target-uniformizer decoder, and commuting-square boundary"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  The physical 7-by-13 cofiber table, its primitive
  cyclotomic sectors, the intrinsic mod-13 coefficient Bockstein decoder,
  and the chart-by-target commuting-square boundary have exact companions
  replayed byte-identically under normal and optimized Python.  Until an
  independent audit promotes this file, no result may depend on it.
source: root/joint-clock-target-wing-holotopy-2026-07-28
depends_on:
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2751-root-zero-clutch-mayer-vietoris-wing-shear
related:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - THM-2754-diagonal-clock-81-label-root-zero-clutch-addendum
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
script: 04-computation/lrc14_joint_c7_c13_right_wing_mixed_spectrum_thm2771.py
output: 05-knowledge/results/lrc14_joint_c7_c13_right_wing_mixed_spectrum_thm2771.out
script_sha256: 39bced4e4b8633b003ae0aa09f79c0482ed49b684b1fa1685592b29182a88a27
output_sha256: 02b42b1f7d53c09792268b7c738ba11e4511aeb7711b36c2894c1d5e0ff2736e
secondary_script: 04-computation/lrc14_joint_c7_c13_commuting_square_no_go_thm2771.py
secondary_output: 05-knowledge/results/lrc14_joint_c7_c13_commuting_square_no_go_thm2771.out
secondary_script_sha256: fc927430f31022621ab9cd68ff32a74270e7956b84a00b79da419b0c9775c2b1
secondary_output_sha256: 4dd63ef4ae2f24da0d5d51e66756ddb96e3835a2ef0ad975ff393a0a5d917355
hash_basis: LF-normalized bytes
---

# THM-2771 -- joint clock-target wing spectrum and target-uniformizer boundary

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

The corrected THM-2751 right cofiber has substantially more structure than
either of its one-variable shadows suggests.  Retaining the physical
present-clock label `e in F_7` and target label `t in F_13` gives a primitive
`C_7 x C_13` coefficient with all seventy-two fully mixed characters alive
and invertible modulo `91`.  That mixed unit cannot itself pay the THM-2591
chart-augmentation invoice: every fully mixed character has nontrivial
`C_7` character and therefore zero chart augmentation.

The useful signal lies in the apparently worse pure-target sector.  Its
mod-`13` target collapse is a uniformizer of the local group algebra
`F_13[C_13]`.  An explicit target convolution decodes it to `u-1`; applied
rowwise to the intrinsic coefficient Bockstein of the **raw** physical table,
it gives the sparse seven-chart cochain

```text
(0,2,0,10,0,0,0),
```

whose uniform `C_7` convolution is pointwise `(-1,...,-1)`.  Algebraically
this cancels a constant `(+a)` chart obstruction after scaling by `a`.
What remains open is exactly the physical typing: no lawful common-ancestry
transition realizes this target convolution and uniform chart average.

Independently, chart pullback and target relabelling commute as literal
interval operations on all `169` labelled squares.  Their alternating
four-corner boundary is nonzero, but it is a mixed coboundary rather than a
commutator or holotopy `2`-class.

## 1. The physical `7 x 13` right-cofiber table

Use THM-2749's rail-eight physical carrier.  For each present-clock sheet
`e in F_7` and target label `t in F_13`, reconstruct the corrected THM-2751
objects

```text
A_(e,t),                         B_(e,t)=T_tau^(-1) A_(e,t),
M_(e,t)=A_(e,t) intersect B_(e,t),
R_(e,t)=B_(e,t) minus A_(e,t).                         (1)
```

Both relative-present complements, the source seam, and the actual delayed
`D^6 Q_(3,{1,2})` terminal prefix are retained before taking the constant-tail
coefficient.  The seven physical present-clock sheets are checked literally
disjoint.  Let `R_raw(e,t)` be the resulting integer table.

Its exact common content is

```text
g0=5905329039529920
  =2^6*3^4*5*7^2*11^2*13*43*53*1297.                 (2)
```

Thus `v_13(g0)=1`, and the raw physical table is zero modulo `91`.  Division
by the full `g0` below exposes a primitive coefficient shape; it is not a
physical rescaling.  Put

```text
H=R_raw/g0,                 k1=483303,
k2=483287=7^3*1409.                                  (3)
```

Only rows `e=1,2,3` survive, and their target polynomials are

```text
H_1(u)=k1[2(u^3+...+u^11)+121u^12],

H_2(u)=k2[2(u^2+...+u^9)+265u^12],

H_3(u)=k2[2(u^2+...+u^7+u^10+u^11)+254u^12].        (4)
```

The other four rows vanish.  The support has exactly `28` cells, and the
matrix ranks are

```text
rank_Q(H)=3,             rank_F7(H)=1,
rank_F13(H)=3.                                      (5)
```

In particular the joint table is not a rank-one clock-by-target tensor over
`Q`, even though its reduction modulo `7` has rank one.  At `t=3` the exact
corrected full-clock hostile is recovered:

```text
left =(0,0,0,12,0,0,0),
right=(0,9,2,2,0,0,0)                    in F_13^7. (6)
```

## 2. Cyclotomic sectors and the mixed-unit boundary

Identify `(e,t)` with its CRT residue `n mod91`.  The primitive CRT support is

```text
{2,3,8,9,10,16,17,22,24,29,30,31,36,38,43,44,45,
 50,51,57,58,59,64,71,72,80,85,86}.                 (7)
```

The exact polynomial gcd degrees against the corresponding cyclotomic
polynomial, in characteristics `(2,7,13)`, are

```text
order 7:  (0,0,0),
order 13: (12,0,1),
order 91: (0,0,0).                                  (8)
```

Consequently all six order-`7`, all twelve order-`13`, and all seventy-two
order-`91` complex characters survive.  The order-`7` and fully mixed
order-`91` sectors are units modulo `91`; the pure-target sector is not.  The
small exact resultants satisfy

```text
Res(Phi_7,H_clock)
 =2164833805035694204114534905543204737802313970779
 =50 mod91,

Res(Phi_13,H_target)
 =708762831620208775530213399682753305758900865474119127375721370505898041395766813306566259706602082304
 =78 mod91,                         v_13=1.           (9)
```

The primitive augmentation is

```text
sum_(e,t) H_(e,t)=333470254=5 mod7=0 mod13.          (10)
```

Hence `H` is neither a unit of the full `(Z/91)[C_91]` group ring nor of its
full augmentation quotient.  It is invertible only after projection to the
fully mixed sector.

This limitation is conceptual, not just numerical.  Every fully mixed
character has a nontrivial `C_7` component, so chart augmentation kills it.
Arbitrary combinations of the seventy-two mixed characters remain in the
`C_7` augmentation ideal.  They cannot directly supply THM-2591's nonzero
chart invoice.  Any such contribution must pass through the `C_7`-trivial
pure-target sector.

## 3. The pure-target sector is an augmentation uniformizer

Let

```text
S(u)=sum_e H_e(u).
```

Its integer coefficient vector is

```text
(0,0,1933148,
 2899754,2899754,2899754,2899754,2899754,
 1933180,1933180,1933180,1933180,309305616),         (11)
```

and modulo `13` it is

```text
S=(0,0,9,0,0,0,0,0,2,2,2,2,9).                     (12)
```

Put `epsilon=u-1`.  Since

```text
F_13[C_13]=F_13[epsilon]/(epsilon^13),
```

direct binomial expansion gives

```text
S=7 epsilon+8 epsilon^2+9 epsilon^3+12 epsilon^4
  +2 epsilon^5+4 epsilon^6+4 epsilon^7+7 epsilon^8
  +6 epsilon^9+7 epsilon^10+6 epsilon^11+9 epsilon^12.
                                                               (13)
```

Thus

```text
S=epsilon V,                         V(0)=7 !=0.       (14)
```

The local factor `V` is a unit.  Therefore

```text
(S)=(epsilon)=I,                                      (15)
```

where `I` is the full target augmentation ideal.  In particular

```text
[S]=7[u-1] in I/I^2.                                  (16)
```

Under the target displacement `u -> u^a`, this becomes
`7a[u-1]`; orientation reversal gives the algebraic `-7a` tangent invoice.

An explicit inverse factor is

```text
K=(7,3,3,3,12,2,9,10,9,8,10,4,0) in the u basis,    (17)
```

and exact cyclic convolution proves

```text
S*K=u-1                    in F_13[C_13].             (18)
```

This is an augmentation-ideal decoder, not an inverse of `S` in the whole
local algebra and not a full `(Z/91)[C_91]` inverse.

## 4. Intrinsic Bockstein and the sparse chart cochain

The full-gcd normalization in `(3)` is unnecessary for the transgression.
Equation `(2)` gives

```text
v_13(g0)=1,                       (g0/13) mod13=11.   (19)
```

Therefore the coefficient Bockstein

```text
beta(R_raw)=R_raw/13 mod13=11H                              (20)
```

is intrinsic.  Since `11^(-1)=6`, use

```text
K_beta=6K=(3,5,5,5,7,12,2,8,2,9,8,11,0).            (21)
```

Then the target collapse of `beta(R_raw)` still satisfies

```text
S_beta*K_beta=u-1.                                    (22)
```

Rowwise convolution gives the exact decoded table

```text
e=1: (2,6,7,7,9,2,8,12,6,5,7,11,6),
e=2: (0,11,1,7,6,9,9,4,12,8,8,5,8),
e=3: (10,10,5,12,11,2,9,10,8,0,11,10,12),           (23)
```

with the other four rows zero.  Its target-`0` chart column is

```text
C=(0,2,0,10,0,0,0),                 sum_e C_e=-1.    (24)
```

The target column sums are exactly `(-1,1,0,...,0)`, as required by `(22)`.
Let

```text
N_7=1+x+...+x^6.
```

Cyclic chart convolution now gives the stronger pointwise identity

```text
C*N_7=(sum_e C_e)N_7=-N_7.                           (25)
```

Thus `a(C*N_7)` is the constant seven-chart correction `(-a,...,-a)`.  If
the physical chart edge obstruction of THM-2542 is identified with the
present clock coordinate and equals `(+a,...,+a)`, `(25)` cancels it
pointwise.  This makes the admissible-theta selector irrelevant **inside the
coefficient algebra**.

The qualification is load-bearing.  Operations `(20)--(25)` divide a raw
coefficient by `13`, convolve its target labels, and uniformly average its
chart coordinate.  No proved physical map carries those operations through
one common-ancestry endpoint packet or transition square.  Uniform averaging
also destroys local chart provenance.

There are two different coordinate identifications here.  The physical
present clock `e` is the `c1` comb sheet of THM-2749, whereas THM-2542's chart
clock is the delayed Latin scheduler of THM-2517/2535.  This mismatch no
longer matters to the **final algebraic vector**: `-N_7` is fixed by every
permutation of seven labels.  It still matters physically, because canon has
no common-ancestry carrier on which the two clock actions coexist.

More importantly, `t` is THM-2742's second endpoint-dipole coordinate
`lambda=e_c3-e_q2`; it is not the physical root-deck character of THM-2542.
A common numeral in `F_13` is not an identification.  Using `(25)` as a deck
correction therefore also requires an affine target-role-to-root-deck
intertwiner `pi` that preserves the same ancestry.  Hence `(25)` is an
explicit coefficient candidate and isolates the two remaining maps; it is
not yet a THM-2591 payment or an LRC row exclusion.

## 5. The chart-by-target square commutes

Freeze physical clock `e=1`.  Let `sigma_delta` be lawful target relabelling
and let the chart operation be pullback by `tau` followed by the common rail
cut.  Literal interval reconstruction gives, for every `t,delta in F_13`,

```text
T_tau^(-1) sigma_delta A_t
  =sigma_delta T_tau^(-1) A_t
  =B_(t+delta).                                        (26)
```

All `169` object squares commute; there is no exceptional square.  The
primitive fixed-clock right cofiber is

```text
r=(0,0,0,2,2,2,2,2,2,2,2,2,121),                    (27)
```

and its target-difference matrix has rational rank `12`.  The exact
alternating four-corner boundary census is

```text
delta 0:       (0 nonzero cells, 0 mass),
delta 1,12:    (21,89595253440),
delta 2,11:    (35,90705938400),
delta 3,10:    (49,91816623360),
delta 4,...,9: (56,92186851680).                      (28)
```

Thus every nonzero target displacement has a nonzero mixed boundary.  But
`(26)` shows that this boundary is

```text
Delta_target Delta_chart A,
```

not the operation commutator

```text
[Delta_target,Delta_chart]A=0.                        (29)
```

It is a mixed coboundary, not a transition `2`-cocycle or holotopy class.

## 6. Nonoverlap and exact boundary

THM-2754 is proved and supplies a rank-one independent-clock coefficient bank
on `81` full-target labels.  It is nonoverlapping with the present theorem:
it does not construct the physical present-clock-by-target right cofiber,
the decoder `(21)`, or a common transition square.  THM-2772's independent
carrier-allocation pullback is also deliberately excluded from this file.

This candidate proves, subject to independent audit:

```text
PROOF-COMPLETE HERE:
  exact physical 7 x 13 raw and primitive cofiber tables;
  all order-7/order-13/order-91 character survival;
  mixed order-91 unit and pure-target nonunit boundary;
  pure-target uniformizer and explicit augmentation-ideal decoder;
  intrinsic raw/13 coefficient Bockstein;
  sparse chart cochain and pointwise uniform correction;
  literal chart-target commutation and nonzero mixed-coboundary census.

NOT PROVED HERE:
  a lawful physical realization of target convolution K_beta;
  a target-dipole-to-root-deck intertwiner on the same ancestry;
  preservation of common ancestry, endpoint phase, or semantic owner;
  a noncommuting transition or THM-2591 holotopy class;
  a semantic arm, carrier-allocation K4, row exclusion, or LRC(14).        (30)
```

## 7. Reproduction

From the repository root run

```bash
python 04-computation/lrc14_joint_c7_c13_right_wing_mixed_spectrum_thm2771.py
python -O 04-computation/lrc14_joint_c7_c13_right_wing_mixed_spectrum_thm2771.py

python 04-computation/lrc14_joint_c7_c13_commuting_square_no_go_thm2771.py
python -O 04-computation/lrc14_joint_c7_c13_commuting_square_no_go_thm2771.py
```

Each normal/optimized pair byte-matches its stored transcript.  The scripts
use explicit `require` gates, including exact raw content, primitive table,
support, ranks, cyclotomic gcd degrees, resultants, decoder coefficients,
rowwise Bockstein output, and the full four-corner census.  The LF-normalized
addresses are pinned in the front matter.

**Awaiting independent audit; not QED.**
