---
id: THM-2754
title: "Two-clock rank-one separation on the 81-label root-zero clutch"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  rail-8 two-sided root-zero common section, keep the physical c1-present clock
  e and delayed coefficient clock j independent.  For every one of the 81
  full-target labels common to both strict endpoint cylinders, all 49 exact
  source and target coefficients form the same nonzero rank-one matrix
  B_(e,j)(s,t)=a_e(s,t) 1_(j!=0).  The carrier itself is independent of j:
  E3 makes every inherited Present_(j,7) complement automatic, while the
  terminal-Q delayed prefix is empty at j=0 and identical at j=1,...,6.  The
  diagonal therefore equals the physical amplitude vector and gives a private
  unit at roots 12 and 1 for all 81 labels.  This strengthens, but does not
  replace, THM-2749's all-rail coindexed and fixed-clock families.  It remains
  a common-section coefficient theorem, not an endpoint current, row
  exclusion, or LRC(14).
source: root/two-clock-root-zero-clutch-separation-2026-07-28
audit: >
  semantic-present-wall-2026-07-28 (independent 567-carrier reconstruction,
  all-label unit/Fourier audit, full 3969-cell object and coefficient replay,
  one-sided shear provenance, exact margin/Bezout checks, and normal/-O/hash
  replay: ACCEPT)
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2742-full-two-target-present-sheet-deepest-source-semantic-current
  - THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice
  - THM-2750-arm-blind-clutch-no-go-and-minimal-marked-leakage
  - THM-2751-root-zero-clutch-mayer-vietoris-wing-shear
script: 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
output: 05-knowledge/results/lrc14_root_zero_full_target_semantic_clutch_20260728.out
script_sha256: 208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8
output_sha256: 4f795d055b2b06e46b4250c67d874444d437217c45f4435063791923e377bcd6
hash_basis: LF-normalized bytes
---

# THM-2754 -- the root-zero common section separates its two clocks

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2749 contains two differently typed specializations of the fully marked
root-zero clutch: an all-rail family in which the physical and delayed clocks
are coindexed, and a rail-8 family in which the physical section is frozen at
clock one while the delayed clock varies.  This theorem constructs the object
containing both.  On rail `8` and the whole `9 x 9` common target-label bank,
the two clocks are independent coordinates of an exact response matrix.  That
matrix separates as an outer product for an object-level reason.

## 1. The two-clock common carrier

Use

```text
p=13,             R=13^6=4826809,
T=297836897838480, tau=7/R.                            (1)
```

The strict adjacent source and target points are

```text
q_-=47850889647341/100360982066072,
q_+=47851035194197/100360982066072=q_-+tau.            (2)
```

Their full-target labels common to both inherited open cylinders are

```text
s in S=(0,1,2,3,8,9,10,11,12),
t in U=(3,4,5,6,7,8,9,10,11).                         (3)
```

The minimum physical distance from either point to any boundary of any
accepted `q1/c2` or `q2/c3` factor is

```text
1541619/100360982066072=1541619*Q_RADIUS.              (4)
```

Thus all `81` labels in `S x U` are stable on both whole open cylinders, not
only at the two displayed centres.

Keep two clock indices independent:

```text
e in F_7: physical c1-present clock in F_(e,s,t),
j in F_7: relative-present and delayed coefficient clock.    (5)
```

Let `K_j` be THM-2744's rail-8 equal-weight overlap carrier at delayed clock
`j`; it contains both fixed-label-7 relative-present complements and the open
right-root-12/left-root-1 seam.  Put

```text
S_(e,s,t)=E3 intersect F_(e,s,t),

A_(e,j,s,t)
 =K_j intersect S_(e,s,t) intersect T_tau^(-1)S_(e,s,t). (6)
```

The target carrier is exactly `T_tau A_(e,j,s,t)`, with the translated rail
weight retained separately.  Both endpoint masks in `(6)` are essential;
using only the mask owned by each endpoint gives two different carriers.

## 2. Object-level delayed-clock collapse

The apparent `j`-dependence in `(6)` disappears before integration.  Every
inherited packet `Present_(j,7)` contains the unshifted `c3`-safe factor,
whereas the exclusive source `E3` contains `c3`-danger.  Exact interval
intersection gives

```text
E3 intersect Present_(j,7)=empty              for j=0,...,6. (7)
```

Consequently the two `Present_(j,7)^c` factors in `K_j` are automatic after
the two copies of `S_(e,s,t)` are imposed.  For fixed `(e,s,t)`, all seven
weighted interval tuples `A_(e,j,s,t)` are literally equal.  This is equality
of endpoints and weights, not merely equality of mass or coefficients.

Now intersect every THM-2640 delayed prefix with the actual terminal fork
`Q_(3,{1,2})` on the invariant coordinate `y={Rx}` and rebuild it.  The exact
prefix objects obey

```text
Pi_0=empty,
Pi_1=Pi_2=...=Pi_6,                                    (8)
```

each nonempty kappa prefix having mass `206986279500`.  Define

```text
w=(0,1,1,1,1,1,1).                                    (9)
```

As a `Phi_7` coefficient vector, `w` is
`z+...+z^6=Phi_7(z)-1`, hence the constant unit `-1`; here its more important
role is the exact support mask of `(8)`.

## 3. Rank-one separation and translation

Let `B^-_(e,j)(s,t)` be the source-carry-12 delayed numerator on `(6)`, and
let `B^+_(e,j)(s,t)` be the independently recomputed target-carry-6 numerator
on its translated carrier.  Define the physical-clock amplitudes from the
common nonempty prefix by

```text
a^-_e(s,t)=B^-_(e,1)(s,t),
a^+_e(s,t)=B^+_(e,1)(s,t).                             (10)
```

Equations `(7)--(9)` and extensionality of the integral give

```text
B^-_(e,j)(s,t)=a^-_e(s,t) w_j,
B^+_(e,j)(s,t)=a^+_e(s,t) w_j.                        (11)
```

Translation by `tau` fixes `{Rx}`, sends predecessor carry `12` to `6`, and
maps the weighted source carrier exactly to the target carrier.  Therefore

```text
a^-_e(s,t)=a^+_e(s,t)=a_e(s,t),                       (12)
```

and both endpoint matrices equal

```text
B_(e,j)(s,t)=a_e(s,t) w_j.                            (13)
```

The exact companion constructs and integrates every one of the

```text
81 labels * 7 physical clocks * 7 delayed clocks = 3969 cells. (14)
```

All `3969/3969` carrier-object comparisons, translations, raw source/target
comparisons, and outer-product entries pass.  For every label the amplitude
vector is nonzero, so `(13)` has exact rank one over `Q`; all `2 x 2` minors
vanish.

## 4. The diagonal loses no information

The zero delayed column in `(9)` could make a diagonal forget `a_0`.  On this
bank a separate exact edge calculation gives

```text
a_0(s,t)=0                         for all 81 labels.   (15)
```

Since `w_e=1` for `e=1,...,6`, it follows that

```text
(B_(0,0),B_(1,1),...,B_(6,6))
   =(a_0,a_1,...,a_6).                                 (16)
```

Thus the earlier diagonal-clock vectors are not a lossy shadow here: they are
exactly the physical-clock amplitude vectors.  Their nonzero supports have
the complete census

```text
support (1):       4,
support (1,2):    21,
support (1,3):    21,
support (1,2,3):  35.                                 (17)
```

In particular physical clocks `0,4,5,6` vanish throughout this bank.  No
amplitude vector is zero, and there are exactly `15` distinct vectors.

For `(s,t)=(0,3)`, put

```text
C=339633525654239542165440,
D=750593782703678965571520,
E=719200126392878704654080.                            (18)
```

Then

```text
a=(0,C,D,E,0,0,0),
B=a^T w.                                               (19)
```

The diagonal in `(19)` is THM-2749's displayed rail-8 coindexed vector; the
same bridge is checked directly at THM-2749's uniform label `(0,4)`.  Freezing
`e=1,s=0` instead gives `C w` for every `t` in the nine-label common bank,
recovering THM-2749's frozen-section row.

This theorem does **not** repeat THM-2749's fourteen-rail `(0,4)` table.  It
adds the full off-diagonal response on rail `8` and all `81` common labels.
Conversely, it does not evaluate the unrestricted clock-coindexed `t=2`
section: `t=2` is outside the common two-cylinder universe `(3)`, not asserted
to be physically empty.

## 5. Content, private units, and marked target characters

Every diagonal/amplitude vector in `(16)` is divisible by the inherited
content `26` and is a private unit at source root `12` and target root `1`.
Every vector gcd has `13`-valuation one, with exact range

```text
5905329039529920
 <= gcd <=
302530703523944466130560.                              (20)
```

For `(0,3)`, the joint gcd is `41337303276709440`.  After content division
and root normalization, its two `Phi_7` profiles are

```text
root 12: (0,4,10,8,0,0),
root  1: (0,9, 3,5,0,0),                              (21)
```

which are negatives and both have determinant one.  The sign comes solely
from `12^(-1)=-1` versus `1^(-1)=1` in `F_13`.

Zero-extend only the **marked common-cylinder profile** outside
`t=3,...,11`.  For every common `s`, all twelve primitive endpoint-dipole
target characters are nonzero, already witnessed by physical clock `e=1`;
this gives `108/108` per-`s` characters and `12/12` after summing over `s`.
The aggregated clock-one coefficient is

```text
C_target * (u^3+...+u^11),
C_target=2554386600508776388555200.                    (22)
```

After factoring out `C_target`, the support polynomial

```text
W(u)=u^3+...+u^11                                     (23)
```

is an integral cyclotomic unit, since for `V=u^2+u^6+u^10`,

```text
W V-1=(u^9+u^5+u-1)Phi_13(u),
Norm(W)=1.                                             (24)
```

The scalar multiple `C_target W` is not itself asserted to be an integral
unit.  Equation `(24)` is a coefficient decoder modulo the uniform target
sector, not a physical packet transport.

## 6. Hostile controls and type boundary

If the pulled opposite-endpoint copy of `S_(e,s,t)` is omitted, the chosen
diagonal source and target vectors are

```text
(0,339633525654239542165440,
   750593782703678965571520,
   722054095148406001101120,0,0,0),

(0,345341652135823400016960,
   756301720214733558465600,
   724908063903933297548160,0,0,0),                    (25)
```

and are unequal.  The older one-sided shear also omitted the physical
`c1`-clock; restoring that clock gives exactly `(25)`.  Its terminal-Q prefixes
were already the same, so the discrepancy is not a Q effect.

MISTAKE-313 records the first typed repair.  If `e=1` is fixed on both the
one-sided and common carriers, the source coefficient already equals the
common coefficient `C`; its left wing is coefficient-null.  The target
coefficient is `345341652135823400016960`, and its right wing is
`5708126481583857851520`.  The source, target, and right-wing profiles are
`9`, `8`, and `4`; the folded target/source ratio is `11`, not `7`.

A later `FINITE-EXACT + VERIFIED` source-clocked one-sided sidecar extends this
comparison across the canonical endpoint bank: all `81` lawful labels have
normalized natural target/source gain `11`, while the target has one extra
global label `t=12`.  This is a fixed-`e=1` physical-chart coefficient bank,
not the two-sided common carrier of this theorem, and it awaits independent
hostile audit.  By contrast, the earlier `81`-sheet variable-gain/rank-three
reference census never inserts a physical `c1` clock comb.  Its exact data are
a clock-blind one-sided quotient hostile, not a coindexing `e=j` and not a
counterexample to `(16)`.

There is a second, forgotten-`e` repair.  The legacy clock-blind carrier is
the disjoint union of all seven physical sections.  If that universe is
retained intentionally, its matching common object is the full same-clock
union, not the frozen `e=1` slice.  At `(0,3)` the sections are pairwise
disjoint; after the two-sided cut every cross-clock piece with `e!=e'` is
empty.  The union of the same-clock pieces is exactly the literal unclocked
intersection `M_full=A intersect B` as a weighted interval object.  Its three
amplitudes are

```text
A=1812281403506324508838080,
M_full=1809427434750797212391040=C+D+E,
B=1826551436254490256030720.                           (26)
```

Thus the forgotten-`e` wings are

```text
L=A-M_full= 2853968755527296447040,
R=B-M_full=17124001503693043639680.                    (27)
```

The literal weighted differences `A minus M_full` and `B minus M_full` integrate to
the same values.  Before augmentation, their physical-present rows are

```text
L_e=(0,0,0,2853968755527296447040,0,0,0),

R_e=(0,
     5708126481583857851520,
     5707937511054592894080,
     5707937511054592894080,
     0,0,0).                                            (28)
```

After division by `26` and root normalization, the physical-present residue
rows are

```text
s=(0,0,0,12,0,0,0),
p=(0,9,2, 2,0,0,0).                                  (29)
```

As `Phi_7` profiles, `s=12z^3` has determinant `1`, while

```text
p=9z+2z^2+2z^3=2z(z-1)(z+2)                          (30)
```

has determinant `11`.  Both physical-clock rows are units.  More precisely,

```text
g=2z+2z^2+2z^3+2z^4+6z^5,
det(g)=11,                 g s=p.                      (31)
```

This is an exact coefficient-algebra ratio, but no physical `L -> R` carrier
map or cyclic clock operation has been constructed.

Augmentation has a different truth value: `epsilon(s)=12` whereas
`epsilon(p)=9+2+2=0`.  Equivalently, `v_13(L)=1` while `v_13(R)=2`.  Folding
first and then repeating over the delayed-clock mask gives reduced profiles
`(1,0,0,0,0,0)` and zero, with determinants `1` and `0`.  Hence no
**augmented scalar** wing gain exists.
Together, the two repairs refute the gain-`2` frame of the former THM-2751 body:
it mixed the clock-blind union with the frozen `e=1` common slice.  In the
fixed-`e` universe the left wing is null; in the forgotten-`e` universe the
matching intersection is `C+D+E` and the strongest survivor is `(31)`.
The current THM-2751 slug is a separately rebuilt, reserved fixed-clock
candidate; it and both one-sided sidecars supply no dependency here.

The terminal fork is genuinely present in `(8)`, although removing it leaves
the selected coefficient unchanged on this restricted carrier.  A separate
carry-cell construction followed by the older weighted-prefix numerator
reproduces all fourteen chosen source and target coefficients.

The two axes `(e,j)` are a physical-present clock and a delayed coefficient
clock.  They are not THM-2615's independent present and conjugated-bare
endpoint-translation axes, so `(16)` is not a Radon diagonal.  The marked
`t`-character is THM-2742's endpoint dipole
`lambda=e_c3-e_q2`; it is not a physical deck character, THM-2334 paired
left-relation current, or endpoint amplitude.  Rank-one separation therefore
does not attach the canonical endpoint current or an external semantic arm.

No global clutch functor, all-rail off-diagonal tensor, physical endpoint
current, scalar-ledger decrement, row exclusion, or LRC(14) conclusion follows.

## 7. Exact reproduction

Run

```bash
python3 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
python3 -O 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_root_zero_full_target_semantic_clutch_20260728.out.
```

The companion pins the THM-2744 overlap and THM-2742 full-target scripts by
LF-normalized hashes; reconstructs all `3969` cells; checks the object-level
carrier and prefix laws before coefficient separation; integrates every cell
at both endpoints; checks all rank-one minors, edge cases, content, units, and
marked characters; and retains the one-sided and independent carry-cell
hostiles.  It contains no truth-bearing Python `assert`.

QED.
