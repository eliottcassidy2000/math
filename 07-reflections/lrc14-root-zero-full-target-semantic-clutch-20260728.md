# The root-zero clutch separates its physical and delayed clocks

> **STATUS: FINITE-EXACT + VERIFIED on the canonical rail-8 overlap.**
> The next refinement stated after THM-2744 has a positive answer.  After
> retaining both the source and pulled-target copies of a lawful full-target
> section, and inserting the genuine terminal fork in the delayed prefix, the
> source-root-12 and target-root-1 raw coefficient vectors are identical,
> divisible by the canonical content `26`, and primitive units.  This holds
> for all `81` labels common to both displayed open endpoint cylinders.  It is
> stronger: all `81*49=3969` physical-present/delayed-clock cells form exact
> rank-one two-clock matrices.  It is not a global clutch action, endpoint
> current, row exclusion, or LRC(14) conclusion.

## 1. Inheritance and the missing factor

The closest proved mechanism is [THM-2744](../01-canon/theorems/THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch.md): on rail `8`, translation by

```text
tau=7/R,                    R=13^6,
```

maps the right-root-`12` overlap carrier to the left-root-`1` carrier and
gives equal raw vectors before root normalization.  The corrected near miss
is MISTAKE-310: right root zero is a forbidden label, not empty physical
support, because the right-zero and left-one half-tooth charts overlap on
`(1,13)/182`.  The least-used sidecar is THM-2742's full two-target sheet
`F_(ell,s,t)` and its genuine `E3 -> D^6 -> Q_(3,{1,2})` semantic current.

The earlier overlap companion verified a common label only at the strict
endpoint pair.  Its integral still contained neither the full-target factor
nor the terminal semantic fork.  There is a subtle but load-bearing typing
point: inserting only `F(x)` on the source and only `F(x')` on the target
does not define a common translated carrier.  One must retain both endpoint
conditions on both sides.

This calculation is complementary to both typed families already promoted in
THM-2749.  Its frozen family fixes the **physical present clock** at `e=1`
and `s=0`, then varies the independent delayed clock `j`; its supported rows
are `(0,C,...,C)`.  Its all-rail family instead coindexes the clocks `e=j`,
uses the uniform label `(0,4)`, and separately records the raw `s=0` target
support `t=2,...,11`.  The present calculation neither repeats the all-rail
table nor conflates those target banks.  It stays on rail `8`, retains the
`9*9` whole-cylinder common bank, and first makes `e` and `j` independent.
The old diagonal vectors are then recovered as a corollary of the full
two-clock tensor.

## 2. The symmetric full-target carrier

Let `C_ell` be THM-2744's rail-`8` overlap carrier in the source coordinate.
It already contains

```text
rho_8(x) and rho_8(x+tau),
Present_(ell,7)^c(x) and Present_(ell,7)^c(x+tau),
H^R_12(x) and H^L_1(x+tau).                            (1)
```

The `Present_(ell,7)` factor in `(1)` is the fixed one-target label `7`; it
is not the full two-target sheet.  For a lawful label `(s,t)`, put

```text
S_(ell,s,t)=E3 intersect F_(ell,s,t).                  (2)
```

The correct source-coordinate refinement is

```text
C^*_(ell,s,t)
 =C_ell intersect S_(ell,s,t)
        intersect T_tau^(-1)S_(ell,s,t).               (3)
```

The target-coordinate carrier is exactly `T_tau C^*`.  Thus `(3)` retains
both `F(x)` and `F(x+tau)`, rather than comparing two unrelated one-sided
restrictions.

For the terminal word, let `y={Rx}`.  Instead of materializing an
`R`-fold pullback, intersect each underlying THM-2640 delayed prefix with

```text
Q=Q_(3,{1,2})                                           (4)
```

on the `y` circle before rebuilding the prefix.  This is exactly the factor
`D^(-6)Q` in the physical `x` coordinate.

## 3. Translation equivariance is the equality mechanism

Write `S=R/13`, `z={Sx}`, and

```text
c=floor(13z),              y={13z}={Rx}.               (5)
```

Under `x -> x+tau`,

```text
z -> z+7/13 mod1,
c -> c+7 mod13,
y -> y.                                                        (6)
```

In particular source carry `12` goes to target carry `6`.  Every delayed
clock prefix, including its intersection with `(4)`, depends only on the
invariant coordinate `y`.  On rail `8`, the source and pulled-target step
weights agree.  Therefore the exact delayed-carry functional obeys

```text
L_(6,ell)(T_tau A)=L_(12,ell)(A)                       (7)
```

for every symmetrically cut weighted carrier `A` of the form `(3)`.  Raw
source/target equality is consequently structural; computation is needed for
positivity, content, and the unit norm.

### The missing two-clock object

Keep the two clocks independent.  Write `e in F_7` for the physical
`c1`-present clock in `F_(e,s,t)`, and `j in F_7` for the delayed-prefix and
relative-present clock.  Let `B_(e,j)(s,t)` be the raw coefficient obtained
from the two-sided carrier with these two choices.  Two exact object laws
collapse this apparently two-dimensional response.

First, the inherited one-target packet `Present_(j,7)` contains unshifted
`c3`-safe, while `E3` contains `c3`-danger.  Hence

```text
E3 intersect Present_(j,7)=empty             for every j.       (7a)
```

After intersecting with `E3 intersect F_(e,s,t)` at both endpoints, both
relative-present complements in the THM-2744 carrier are automatic.  The
resulting weighted carrier is exactly independent of `j`, not merely equal
in mass or coefficient.

Second, after inserting the terminal fork, the rebuilt delayed-prefix pairs
obey the exact tuple identity

```text
Pi_0=empty,
Pi_1=Pi_2=...=Pi_6.                                      (7b)
```

Put `w=(0,1,1,1,1,1,1)`.  Equations `(7a)--(7b)` and extensionality of the
delayed integral give the rank-one law

```text
B_(e,j)(s,t)=a_e(s,t) w_j                               (7c)
```

at both endpoints.  The script reconstructs and integrates all
`81*7*7=3969` cells, checks every carrier tuple identity and translation, and
checks every entry of `(7c)`.  Every matrix has exact rational rank one: all
its `2 x 2` minors vanish and its amplitude column is nonzero.

There are two edge cases.  The delayed column `j=0` vanishes because `Pi_0`
is empty.  Separately, `a_0(s,t)=0` on all `81` labels.  Therefore the
diagonal loses no amplitude:

```text
(B_(0,0),...,B_(6,6))=(a_0,...,a_6).                   (7d)
```

This explains the relation to THM-2749 exactly.  THM-2749 is the fixed row
`e=1,s=0`; this addendum's seven-vector is the diagonal `e=j` of the same
rank-one tensor.

## 4. The chosen label `(s,t)=(0,3)`

The two marked endpoints have common label banks

```text
s=(0,1,2,3,8,9,10,11,12),
t=(3,4,5,6,7,8,9,10,11).                              (8)
```

An independent factor-boundary audit proves this is a whole-cylinder, not
merely pointwise, statement.  Across both endpoints and every accepted
`q1/c2` and `q2/c3` factor in `(8)`, the least physical distance to a shifted
width-`1/14` boundary is

```text
1541619/100360982066072
 =1541619*Q_RADIUS > Q_RADIUS.                         (8a)
```

For `(s,t)=(0,3)`, the exact source/target piece counts over
`ell=0,...,6` are

```text
(0,0), (239,239), (526,526), (504,504),
(0,0), (0,0), (0,0).                                  (9)
```

After applying the `Q`-refined delayed-carry functional, both raw vectors are

```text
A^-=A^+
 =(0,
   339633525654239542165440,
   750593782703678965571520,
   719200126392878704654080,
   0,0,0).                                             (10)
```

Their joint gcd is

```text
41337303276709440,
v_13(gcd)=1.                                           (11)
```

Hence every coordinate in `(10)` is divisible by the inherited canonical
content `26=2*13`, and dividing by `26` does not erase the row modulo `13`.

Root normalization gives the reduced `Phi_7` classes

```text
source root 12: (0,4,10,8,0,0),
target root  1: (0,9, 3,5,0,0).                       (12)
```

These are negatives modulo `13`.  This is the formal mechanism: root
`12` normalizes by `12^(-1)=-1`, while root `1` normalizes by `+1`.
The multiplication determinant on the six-dimensional `Phi_7` algebra is
homogeneous of even degree six, so negation preserves it.  Both exact
determinants in `(12)` are

```text
1 mod13.                                               (13)
```

Thus one content calculation and one norm calculation prove both endpoint
unit tests; they are not unrelated numerical successes.

## 5. The whole common-label bank survives

The exact universe is the Cartesian product of the two banks in `(8)`, hence
`9*9=81` labels.  For every one of those labels, the companion verifies

```text
translated weighted carrier identity: 81/81,
raw source/target vector equality:      81/81,
content-26 divisibility:                81/81,
private unit at roots 12 and 1:         81/81.          (14)
```

No vector is zero.  There are exactly `15` distinct vectors.  Their nonzero
clock supports have the census

```text
support (1):       4,
support (1,2):    21,
support (1,3):    21,
support (1,2,3):  35.                                  (15)
```

Every vector gcd has `v_13=1`; the gcd range is

```text
5905329039529920
 <= gcd <=
302530703523944466130560.                              (16)
```

Thus the positive answer is not a lucky property of `(0,3)`.  It is uniform
over the entire full-target label bank common to both strict cylinders.

## 6. The fully marked target profile

For each common `s`, zero-extend the vector-valued coefficient profile by

```text
A_(s,t)=the seven-clock vector above,   t=3,...,11,
A_(s,t)=0,                              otherwise.      (17)
```

The zeros in `(17)` are imposed by the requirement that `t` be lawful on
**both displayed strict endpoint cylinders**.  They do not assert that an
unrestricted physical `F_(ell,s,t)` section is empty outside the common
bank.

For each of the nine common `s` labels and every `k=1,...,12`, exact
power-basis reduction in `Q(zeta_13)` gives

```text
sum_(t in F_13) A_(s,t) zeta_13^(kt) != 0.             (18)
```

All `9*12=108` vector-valued certificates are already witnessed by clock
coordinate `ell=1`.  After summing over the nine `s` labels, that coordinate
has the exact profile

```text
0,0,0,C,C,C,C,C,C,C,C,C,0,

C=2554386600508776388555200.                           (19)
```

Consequently all twelve primitive target characters also survive after
aggregation.  The algebraic reason is sharp: a rational coefficient
polynomial vanishing at one primitive thirteenth root must be a multiple of
`Phi_13`, so all thirteen entries would be equal; `(19)` has both zero and
positive entries.

More strongly, for

```text
W(u)=u^3+...+u^11,             V(u)=u^2+u^6+u^10,
```

exact multiplication gives

```text
W V-1=(u^9+u^5+u-1)Phi_13(u).                         (19a)
```

After factoring out the scalar amplitude `C`, the support window `W` is an
integral cyclotomic unit with positive inverse `V`.  The aggregate coefficient
is `C W`, not itself an integral unit.  The norm of `W` is one because
`W(zeta)=zeta^3(1-zeta^9)/(1-zeta)` and multiplication by `9` permutes the
nontrivial thirteenth roots.  In `Z[C_13]`, the same identity reads
`W*V=(3,2,...,2)`.  This is a coefficient decoder modulo the uniform target
sector, not a physical packet transport.

This is the endpoint-dipole `lambda=e_c3-e_q2` target character of THM-2742.
It is not a physical deck character, the paired left-relation coordinate of
THM-2334, or a pointwise endpoint amplitude.

## 7. Hostile controls and exact boundary

Three controls identify what is load-bearing.

First, omit the pulled-target copy from `(3)` and integrate only the lawful
one-sided sections.  The source and target vectors become

```text
source=(0,339633525654239542165440,
          750593782703678965571520,
          722054095148406001101120,0,0,0),

target=(0,345341652135823400016960,
          756301720214733558465600,
          724908063903933297548160,0,0,0),              (20)
```

which are unequal.  Both `F(x)` and `F(x+tau)` are therefore essential for
the clutch identity.

Second, the terminal fork is genuinely inserted by `(4)`, but on this
particular overlap carrier removing `(4)` leaves `(10)` unchanged.  This is a
restricted redundancy of `Q` on the surviving carry/rail support, not a claim
that `D^(-6)Q` is globally redundant.

Third, an independent replay explicitly cuts the predecessor carry cells
`n=12 mod13` and `n=6 mod13` on the canonical grid, then applies the older
`delayed_weighted_numerator` routine.  It reproduces all fourteen chosen
source/target coefficients, independently of `delayed_carry_pair`'s carry
range-difference implementation.

Fourth, the physical clock sections give a hostile audit of the proposed
Mayer--Vietoris wings.  There are two typed repairs.  If `e=1` is frozen on
both the one-sided and common sheets, the source already equals the common
coefficient `C`, so the left wing is coefficient-null.  The target differs by
`5708126481583857851520`; the source, target, and right-wing profiles are
`9`, `8`, and `4`, and the folded target/source ratio is `11`.

If the legacy clock-blind universe is retained instead, the physical sections
are pairwise disjoint, their union is the unclocked sheet, and after the
two-sided cut every `e!=e'` cross piece is empty.  Hence the literal full
intersection is the union of the three nonempty same-clock pieces, with
amplitude

```text
M=C+D+E=1809427434750797212391040.
```

For the unclocked one-sided amplitudes

```text
A=1812281403506324508838080,
B=1826551436254490256030720,
```

the forgotten-`e` wings are `A-M=2853968755527296447040` and
`B-M=17124001503693043639680`.  After content `26` and root normalization,
their physical-present residue rows are `(0,0,0,12,0,0,0)` and
`(0,9,2,2,0,0,0)`.  Both are `Phi_7` units.  The first has determinant `1`;
the second factors as `2z(z-1)(z+2)` and has determinant `11`.  More exactly,

```text
g=2z+2z^2+2z^3+2z^4+6z^5,    det(g)=11,    g s=p.
```

This is a coefficient-algebra ratio; no physical wing carrier map or cyclic
clock operation is known.  Augmentation destroys it: the two augmentations
are `12` and `0`, the target cancellation `9+2+2=0` raises its scalar
`13`-valuation from one to two, and the corresponding reduced delayed-clock
determinants are `1` and `0`.  Thus there is no augmented scalar wing gain.
Together with the fixed-`e` null left wing, this rejects the THM-2751
candidate's gain-`2` frame, which mixed the clock-blind union with the frozen
`e=1` common slice.

The proved chain is now

```text
THM-2744 open chart overlap
 + both endpoint copies of E3 intersect F_(ell,s,t)
 + Q-refined invariant delayed prefix
 -> translated common weighted carrier
 -> equal recomputed raw vectors
 -> content 26 with v_13=1
 -> private units at roots 12 and 1
 -> every marked primitive target character.          (21)
```

This closes the explicit semantic/full-target-in-the-integral refinement left
after THM-2744 and strengthens THM-2749 with the full `81`-label two-clock
separability law.  THM-2749 already supplies its own fixed-clock and
clock-coindexed target profiles; neither theorem may substitute one typing for
the other.
It does not construct a global clutch functor, identify the target
`t`-character with a physical deck character, or turn an aggregate unit into
a pointwise endpoint amplitude.  The next live obligation is a **physical
endpoint-current attachment**: retain the physical deck or paired
left-relation coordinate through `(21)` and prove a nonzero oriented current,
rather than another target-only unit census.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
python3 -O 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_root_zero_full_target_semantic_clutch_20260728.out.
```

LF-normalized SHA-256:

```text
script  208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8
output  4f795d055b2b06e46b4250c67d874444d437217c45f4435063791923e377bcd6
```

The companion pins the THM-2744 overlap script and THM-2742 full-target
script by LF-normalized evidence bytes, freezes the exact `81`-label universe,
constructs all `3969` physical/delayed-clock cells at both endpoints,
intersects the delayed prefixes with `Q`, checks exact carrier-object
independence and translation before the rank-one law, verifies both typed wing
repairs and the physical-unit/augmentation split, checks content and norms, reduces
every marked target character in the exact `Phi_13` power basis, and includes
the direct carry-cell replay and one-sided hostile.  It contains no
truth-bearing Python `assert`.
