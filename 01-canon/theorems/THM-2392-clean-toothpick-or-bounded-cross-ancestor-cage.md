---
id: THM-2392
title: "Clean toothpick or bounded cross-ancestor cage"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Assuming
  THM-2388's last-lane hypotheses and exact 36/343 hole/excess ledger, let
  delta be the mass of unit holes which are outside the quotient-blocker
  union. Every such parent has exactly one double-covered inverse root,
  guard root count four, and five ordinary root counts two. A deterministic
  distinguished top label q_* gives either a singleton or adjacent
  two-root exclusive word whose normalized nonzero target spectrum is
  everywhere nonzero, with exact square/fourth-power sums respectively
  (12/169,12/28561) or (22/169,62/28561). An owner/status/translate cell
  has mass at least delta/52. THM-2391 simultaneously gives a unique
  excess among the seven septimal siblings; joint owner resolution costs
  only 338 cells and yields a same-parent C_7 x C_13 tensor nonzero in
  every septimal colour and every nonzero target colour. Conversely, the
  blocker cage has mass at least 36/343-delta. Different 7-adic valuations
  force exact danger overlap 1/49; the last-lane hypotheses apply this to
  both high cross-ancestor pieces. Thus if delta<6/4459, one of the two
  remaining low cross pairs has reduced product at most
  1/[2(6/4459-delta)]. Every repeated-first profile forces
  delta>=1/26754. Among the 150 strict profiles, all 135 rows with b>=3
  force delta>=6042/9796423 and a fixed charged cell of mass at least
  3021/254706998. The fifteen b=2 rows either have positive clean mass or
  their compatible four-piece low-blocker union lies in an exact bank of
  ten oriented ratios, namely six up to interchange. THM-2391 reduces the
  no-clean bank further to (4,3) at septimal depth M=2 and makes it empty
  at M>=3. This is a coefficient-level root word, not a canonical
  expiration target, branch exclusion, row decrement, or proof of LRC(14).
source: codex-2026-07-26-clean-toothpick-cross-ancestor-cage
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-sidecar
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
script: 04-computation/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.py
output: 05-knowledge/results/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.out
script_sha256: d88d5231f6efbefbcec515e2e50f78d7853471fea7ffc0c17424963e2e84b1be
output_sha256: 26701b0968864af671ddb7b0faa6a9a470f65799513f3a4a8371d3bb666649d2
hash_basis: working-tree bytes (LF)
---

# THM-2392 -- a clean toothpick or a bounded blocker cage

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2388 finds exact signed mass in the final septimal
`k=2,(t,b)=(1,0)` lane, but it does not put that mass outside the
quotient blockers.  This theorem keeps that missing bit rather than silently
discarding it.  The result is the exact alternative

```text
positive clean-hole mass
  -> one double root + a top-labelled singleton/adjacent charged word;

small clean-hole mass
  -> a large quotient/original blocker cage
  -> a bounded cross-ancestor reduced ratio.                         (1)
```

The repeated-first profiles and all `135` strict profiles with middle depth
at least three fall quantitatively on the first side.  The last `15`
middle-depth-two profiles retain a finite exact cross-ancestor bank.

## 1. The clean-hole/cage split

Use THM-2388's notation.  Thus

```text
T(x)=13x,

K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i)),

c_j=13C_j,

B=D_(C_1) union D_(C_2) union D_(C_3),                (2)
```

and the scalar cover implies

```text
{K=0} subset T^(-1)B,

{K>=2} subset B.                                      (3)
```

In the last septimal lane put

```text
X=D_(q_*)^c intersection D_(c_3)^c,

Z={K=0}.                                               (4)
```

THM-2388 proves

```text
mu(Z intersection X)>=36/343,

Z intersection X
 subset T^(-1)(
   (D_(C_1) union D_(C_2)) minus D_(C_3)
 ).                                                    (5)
```

Define the clean-hole mass

```text
S=(Z intersection X) minus B,

delta=mu(S),                                           (6)
```

and the geometric cage

```text
Gamma
 =B intersection
   T^(-1)((D_(C_1) union D_(C_2)) minus D_(C_3))
   intersection X.                                    (7)
```

Splitting `Z intersection X` according to membership in `B` gives

```text
mu(Gamma)>=36/343-delta.                               (8)
```

This is the first correction to the tempting finish.  The exact
`36/343` in THM-2388 is not `delta`; all of it could a priori lie in
`Gamma`.

## 2. Every clean hole is a labelled toothpick

Let `y in S` be generic and write its inverse roots as

```text
x_r=(y+r)/13,                         r in F_13.       (9)
```

Because `y notin B`, all three original blockers are safe at every root.
The scalar cover gives

```text
K(x_r)>=1                            for every r.      (10)
```

THM-2388's root reflection is

```text
sum_r K(x_r)=14-K(y).                                 (11)
```

Since `K(y)=0`, equations (10)--(11) force exactly one root of
multiplicity two and twelve roots of multiplicity one.  In particular,
the double root has a unique unordered pair among the six labelled unit
masks.

The individual root counts are stronger.  Every unit mask vanishes at
`y`, so the ordinary and guard reflection formulas give

```text
#{r:x_r in E_H}=4,

#{r:x_r in D_(q_i)}=2                 for i=1,...,5.  (12)
```

Their total is `4+5*2=14`, exactly the incidence count in (11).  Thus the
complete labelled root partition has one of two multiplicity profiles:

```text
double pair {H,q_i}:
  singleton counts (H,q_i,other q's)=(3,1,2,2,2,2);

double pair {q_i,q_j}:
  singleton counts (H,q_i,q_j,other q's)=(4,1,1,2,2,2).
                                                               (13)
```

There are at most

```text
13[
  5*12!/(3!(2!)^4)
  +10*12!/(4!(2!)^3)
 ]=648648000                                         (14)
```

complete labelled root words.

There is also a blocker label at the parent.  Equation (5) and
`y in X` say that at least one of `D_(c_1),D_(c_2)` contains `y`, while
`D_(c_3)` does not.  Choosing the least active low blocker supplies a
deterministic labelled owner in two possibilities.  This owner need not be
exclusive; retaining the complete nonempty low-blocker word instead gives
the three possibilities `{1},{2},{1,2}`.

## 3. The distinguished top-labelled charged word

Use the already distinguished top septimal label `q_*`.  It has two
adjacent active roots in (12).  If `q_*` is outside the double pair, both
are exclusive singleton roots.  If `q_*` belongs to the double pair,
exactly one is double-covered and the other remains exclusive.  Let `A`
be the indicator of the exclusive `q_*` roots and put `W=1-A`.  Thus

```text
q_* outside the double pair:
  A=1_({r_0,r_0+eta}),                 eta in F_13^*;

q_* inside the double pair:
  A=1_({r_0}).                                         (15)
```

This removes the unnecessary choice of a label outside the double pair.
Partition `S` only by

```text
chosen low-blocker owner:       2 possibilities,

q_* double-pair status:         2 possibilities,

singleton/edge translate:      13 possibilities.      (16)
```

One literal top-labelled charged cell `Y` therefore has

```text
rho:=mu(Y)>=delta/(2*2*13)=delta/52.                  (17)
```

If the complete low-blocker word rather than one deterministic owner is
required, the corresponding bound is `delta/78`.

Use the normalized target DFT

```text
a_k=(1/13)sum_(r in F_13)A(r)zeta^(-kr),

w_k=(1/13)sum_(r in F_13)W(r)zeta^(-kr),

zeta=exp(2 pi i/13).                                   (18)
```

For every `k!=0`,

```text
w_k=-a_k.                                             (19)
```

In the singleton status,

```text
|a_k|^2=1/169,

sum_(k!=0)|a_k|^2=12/169,

sum_(k!=0)|a_k|^4=12/28561.                           (20a)
```

In the adjacent status,

```text

|a_k|^2
 =(2+2cos(2 pi k eta/13))/169
 >0,

sum_(k!=0)|a_k|^2=22/169,

sum_(k!=0)|a_k|^4=62/28561.                           (20b)
```

The strict inequality uses oddness of thirteen: `-1` is not a thirteenth
root of unity.  Thus every nonzero target colour survives in either
status.

Consequently, on the fixed cell `Y`,

```text
integral_Y a_k conjugate(w_k)
 =-rho |a_k|^2<0,

singleton:
  |integral_Y a_k|=rho/13,
  integral_Y |a_k|^2=rho/169;

adjacent:
  |integral_Y a_k|>=(2rho/13)sin(pi/26),
  integral_Y |a_k|^2
   >=(4rho/169)sin^2(pi/26).                          (21)
```

The sums in (20a)--(20b) acquire the factor `rho`.

Normalization is load-bearing.  With the unnormalized DFT
`A_k^sharp=sum_r A(r)zeta^(-kr)`, the singleton sums are `12,12` and
the adjacent sums are `22,62`, not their normalized values above.

### 3a. The same parent carries a nonzero `C_7 x C_13` tensor

There is a second, transverse word over every generic `y in S`.  For
`s in F_7`, put

```text
L_7(s)
 =1_(E_H)(y+s/7)
  +sum_(q_i!=q_*)1_(D_(q_i))(y+s/7)
  +1_(D_(c_1))(y+s/7)+1_(D_(c_2))(y+s/7).            (22a)
```

Both high masks `q_*` and `c_3` are divisible by seven, so their safe
status at `y` persists at all seven siblings.  The scalar cover gives
`L_7(s)>=1`.  THM-2391 puts all seven displayed lower labels in the
primitive septimal layer, and their total seven-root incidence is

```text
2+4+1+1=8.
```

Therefore there is a unique `d(y) in F_7` such that

```text
J_y(s):=L_7(s)-1=1_({d(y)})(s).                       (22b)
```

Its normalized septimal DFT is

```text
j_l=(1/7)exp(-2 pi i l d(y)/7),        l in F_7,      (22c)
```

and is nonzero for every septimal colour, including zero.

The digit `d` also resolves the low-blocker owner.  Since `s=0` is the
original parent and `K(y)=0`,

```text
L_7(0)=1_(D_(c_1))(y)+1_(D_(c_2))(y).
```

Thus

```text
d=0  iff both low blockers are active;

d!=0 iff exactly one is active.                       (22d)
```

There is one simultaneous-owner category and `6*2` exclusive
`(d,owner)` categories, hence `13` in all.  Combining these with the
`26` singleton/adjacent support categories from Section 3 partitions
`S` into only

```text
13*26=338
```

cells.  On one fixed owner-resolved cell `Y_7x13`,

```text
mu(Y_7x13)>=delta/338,                                (22e)
```

and the rank-one tensor coefficient

```text
j_l a_k
```

is nonzero for every `l in F_7` and every `k!=0 in F_13`.  Its magnitude
is exactly `1/91` in the singleton status and at least

```text
2 sin(pi/26)/91
```

in the adjacent status.                               (22f)

This repairs the pre-THM-2391 stopping statement that the septimal word
need not meet `S`: it meets every generic clean parent.  The two words
live on transverse sibling/predecessor fibres over the same parent, not
on one common physical root.  No target/endpoint intertwiner is inferred.

## 4. The general cross-ancestor alternative

The cage in (7) is contained in

```text
union_(i=1)^2 union_(j=1)^3
  (D_(c_i) intersection D_(C_j)).                     (23)
```

The two same-line intersections have an exact universal mass:

```text
mu(D_(c_i) intersection D_(C_i))
 =mu(D_(13C_i) intersection D_(C_i))
 =1/91.                                               (24)
```

Indeed, multiplication by `C_i` preserves Haar measure.  In the standard
coordinate, the central tooth of `D_13` has length `1/91` inside `D_1`;
the two neighboring teeth meet the endpoints `+/-1/14` only, and every
other tooth is disjoint.

There is a second exact overlap identity which is easy to miss.
If two positive speeds have different `7`-adic valuations, then

```text
rho(a,b)=1/49.                                        (24a)
```

To prove it, divide by the gcd.  Exactly one reduced coefficient is
divisible by seven.  That coefficient is congruent to its negative modulo
fourteen, so the two folded endpoint arguments in THM-1166 have the same
value.  The correction to `1/49` is therefore zero.

THM-2388's last-lane condition (26) says

```text
nu_7(c_i)<M<nu_7(c_3)=nu_7(C_3),       i=1,2,
```

because division by thirteen does not change a `7`-adic valuation.
Consequently

```text
rho(c_1,C_3)=rho(c_2,C_3)=1/49.                       (24b)
```

Subtract (24) and (24b) from (8).  Only the two low cross pairs
`(c_1,C_2)` and `(c_2,C_1)` remain.  The union bound shows that one of
them satisfies

```text
rho(c_i,C_j)
 >=(36/343-delta-2/91-2/49)/2

 =1/49+3/4459-delta/2.                                (25)
```

Whenever

```text
delta<6/4459,                                         (26)
```

the right side is strictly above `1/49`.  Divide the selected pair by its
gcd and write the coprime reduced speeds as `a,b`.  THM-1166 gives

```text
rho(a,b)
 =1/49+Delta/(196ab),                  Delta<=49.      (27)
```

Equations (25)--(27) imply the finite product bound

```text
ab<=1/[2(6/4459-delta)].                              (28)
```

Equivalently, for every parameter

```text
0<theta<6/4459,                                       (29)
```

one has the explicit alternative

```text
delta>=theta:
  a fixed top-labelled charged cell has mass at least theta/52;

delta<theta:
  one of the two low cross pairs has
    rho>1/49+3/4459-theta/2
  and
    ab<[2(6/4459-theta)]^(-1).                        (30)
```

For example `theta=3/4459` gives

```text
charged-cell mass>=3/231868,

or

low-cross reduced product ab<=743.                    (31)
```

At the extreme `delta=0`, equation (25) gives

```text
rho(c_i,C_j)>=94/4459=1/49+3/4459,

ab<=371.                                              (32)
```

The linear trade (25) is the exact output of the cage mass, the two
same-line identities, the two `7`-valuation identities, and the remaining
two-pair union bound.  No independence of the six intersections is assumed.

## 5. Repeated-first profiles force the charged branch

Now impose one of the fifteen repeated-first profiles

```text
(lambda_1,lambda_2,lambda_3)=(1,1,c),

5<=c<=19.                                             (33)
```

Then

```text
nu_13(c_1)=nu_13(c_2)=1,

nu_13(C_1)=nu_13(C_2)=0,

nu_13(C_3)=c-1.                                       (34)
```

THM-2263 gives the following sharp pair caps.

1. The same-line pairs in (24) equal `1/91`.
2. The two shallow cross pairs `(c_1,C_2),(c_2,C_1)` have gap one and
   mass at most `23/1092`.
3. The two high cross pairs `(c_i,C_3)` equal `1/49` by (24b).

Therefore the whole cage obeys

```text
mu(Gamma)
 <=2/91+2*(23/1092)+2/49
 =401/3822.                                           (36)
```

Combining (5)--(8) and (36) gives the unconditional clean-hole floor
within the THM-2388 candidate:

```text
delta
 >=36/343-401/3822
 =1/26754
 >0.                                                   (37)
```

The labelled-owner charged cell from (17) consequently satisfies

```text
rho>=1/1391208.                                       (38)
```

The owner-resolved transverse tensor from (22e) has the simultaneous
floor

```text
mu(Y_7x13)>=1/9042852.                                (38a)
```

On that one fixed cell, every nonzero target colour has the coefficient
and energy floors

```text
|integral_Y a_k|
 >=sin(pi/26)/9042852,

integral_Y |a_k|^2
 >=sin^2(pi/26)/58778538,                             (39)
```

while the exact summed floors are

```text
sum_(k!=0)integral_Y |a_k|^2
 >=1/19592846,

sum_(k!=0)integral_Y |a_k|^4
 >=1/3311190974.                                      (40)
```

Thus every repeated-first last-lane packet has a positive literal
top-labelled singleton/adjacent charged word.  This is not merely the
formal alternative (30).

## 6. The strict-profile cage collapses except at middle depth two

Now let

```text
(lambda_1,lambda_2,lambda_3)=(1,b,c),

2<=b<c,                         5<=c<=19.             (40a)
```

There are `150` such profiles.  For a positive thirteen-adic gap define
the sharp THM-2263 upper cap

```text
u(d)=
  1/49+6/(49*13^d),             d even,
  1/49+5/(588*13^d),            d odd.                (40b)
```

If `b>=3`, the two low cross pairs have positive gaps `b` and `b-2`.
Equations (24), (24b), and (40b) give

```text
mu(Gamma)
 <=U_b:=2/91+u(b)+u(b-2)+2/49.                        (40c)
```

Within each parity, the correction in `U_b` decreases strictly with `b`.
The odd maximum is at `b=3`, the even maximum at `b=4`, and

```text
6/49(13^(-4)+13^(-2))
 -5/588(13^(-3)+13^(-1))
 =85/1199562>0.                                       (40d)
```

Thus the unique worst middle depth is `b=4`—all admissible `c` at that
middle depth tie—and

```text
U_b<=U_4=146022/1399489.
```

Consequently all `135` strict profiles with `b>=3` satisfy

```text
delta
 >=36/343-U_4
 =6042/9796423,                                       (40e)

rho(labelled charged cell)
 >=(6042/9796423)/52
 =3021/254706998.                                     (40f)

mu(owner-resolved Y_7x13)
 >=(6042/9796423)/338
 =3021/1655595487.                                    (40f')
```

The `15` remaining strict profiles have `b=2`.  Their unresolved low
cross pair is `(c_1,C_2)`, whose thirteen-adic gap is zero.  If

```text
nu_7(c_1)!=nu_7(C_2)=nu_7(c_2),                       (40g)
```

then (24a) makes that overlap exactly `1/49`.  The other low cross pair
has gap two and cap `25/1183`, so

```text
mu(Gamma)<=864/8281,

delta>=36/57967,

rho(labelled charged cell)>=9/753571.                 (40h)

mu(owner-resolved Y_7x13)>=18/9796423.                (40h')
```

It remains to record the exact boundary when the two valuations in (40g)
are equal.  Divide the unresolved pair by its gcd and orient it as

```text
(C_2,c_1)=13h(a,b),

gcd(a,b)=1, gcd(ab,91)=1.                             (40i)
```

The four low speeds then have the correlated form

```text
(C_1,C_2,c_1,c_2)=h(b,13a,13b,169a).
```

Thus the union of all four low-low cage pieces is exactly

```text
U_(a,b)
 =(D_b union D_(13a))
  intersection
  (D_(13b) union D_(169a)).                           (40j)
```

The two pieces involving `C_3` contribute at most `2/49`, so the
extreme no-clean-hole case `delta=0` forces

```text
mu(U_(a,b))>=36/343-2/49=22/343.                      (40k)
```

Before using this compatible union, a coarse pair invoice makes the
search finite.  Since `169=1 mod 14`, the two low cross-pair folded
defects are the same, while their denominators differ by `169`.  The
sum of all six pair caps is

```text
2/91+4/49+(170/169) Delta(a,b)/(196ab).
```

Absorption of `36/343` implies

```text
Delta(a,b)/(ab)>=156/595.                             (40l)
```

Since `Delta<=49`, this first gives `ab<=186`.  Exact enumeration of

```text
a<=b, gcd(a,b)=1, gcd(ab,91)=1, ab<=186,

[F((a+b) mod 14)-F((b-a) mod 14)]/(ab)>=156/595      (40m)
```

leaves `52` unordered preliminary pairs, with `ab<=177`,
`min(a,b)<=10`, and `max(a,b)<=85`.  Their small-coordinate census is

```text
1:24, 2:9, 3:8, 4:4, 5:4, 8:1, 9:1, 10:1.           (40n)
```

Now evaluate (40j) for both orientations.  Exactly ten oriented ratios
survive (40k):

| `(a,b)` | `mu(U_(a,b))` |
|---|---:|
| `(1,1)` | `193/1183` |
| `(1,2)` | `114/1183` |
| `(2,1)` | `239/2366` |
| `(1,3)` | `263/3549` |
| `(3,1)` | `95/1183` |
| `(4,1)` | `331/4732` |
| `(2,3)` | `43/546` |
| `(3,2)` | `95/1183` |
| `(3,4)` | `491/7098` |
| `(4,3)` | `331/4732` |

Up to interchange, this is the six-ratio bank

```text
1:1, 1:2, 1:3, 1:4, 2:3, 3:4.                       (40o)
```

THM-2391 sharpens the bank by septimal depth.  When `M>=2`, its two
blocker progressions have steps of one common absolute value, either
one or two.  Hence

```text
C_2/C_1=13a/b=+/-1 mod 7^M,

13a=+/-b mod 7^M.                                    (40p)
```

At modulus `49`, only `(a,b)=(4,3)` from the table survives.  At modulus
`343`, none survives.  Consequently the middle-depth-two no-clean
boundary is

```text
M=1:  one of the ten oriented ratios above;

M=2:  only (a,b)=(4,3);

M>=3: impossible, so delta>0.                         (40q)
```

The ten orientations also have a canonical first septimal digit.  Choose
the balanced nonzero residue

```text
sigma in {-3,-2,-1,1,2,3},

13a=sigma b+7d.
```

Then

```text
(c_2-sigma c_1)/7=13h d,                             (40r)
```

with the complete address table

| `(a,b)` | `sigma` | `d` |
|---|---:|---:|
| `(1,1)` | `-1` | `2` |
| `(1,2)` | `3` | `1` |
| `(2,1)` | `-2` | `4` |
| `(1,3)` | `2` | `1` |
| `(3,1)` | `-3` | `6` |
| `(4,1)` | `3` | `7` |
| `(2,3)` | `-3` | `5` |
| `(3,2)` | `2` | `5` |
| `(3,4)` | `1` | `5` |
| `(4,3)` | `1` | `7` |

Thus

```text
d in {1,2,4,5,6,7};                                  (40s)
```

the rows `(1,2),(1,3)` have descendant `13hd=C_2`, and the two rows with
first coordinate four climb one septimal digit.  This is an exact
THM-2377/Bockstein address sidecar; it does not by itself align the
descendant with the top-labelled target word.

Thus each middle-depth-two profile either has positive clean-hole mass,
or its correlated four-speed cage lies in (40o), with the sharper depth
split (40q).  If its low blockers have different `7`-adic valuations
before invoking THM-2391, the stronger uniform positive alternative
(40h) applies.

Together, (37) and (40e) give an unconditional positive toothpick cell in
`150` of the `165` profile rows: the `15` repeated-first rows and `135`
strict rows with `b>=3`.

## 7. What this still does not prove

The new word is exact but its type is limited.

- The low-blocker label at the parent can be simultaneous rather than an
  exclusive THM-2305 owner.
- The double pair lives on a thirteen-root predecessor, while the blocker
  owner which pays the hole lives one forward quotient level away.
- The selected ordinary edge is a lawful coefficient/root word, but no
  canonical expiration clock, endpoint-matched physical reference, or
  THM-2365 target covector has been identified.
- THM-2390/2391's one-double septimal word meets every generic clean
  parent by (22b), but on the transverse seven-sibling fibre rather than
  the thirteen-predecessor fibre.  No lawful cross-fibre phase or target
  intertwiner has been identified.
- For the last middle-depth-two profiles, (40o)--(40q) are finite
  compatible address banks, not an exclusion of either branch.

The first implication which fails in the unqualified THM-2388 finish is

```text
mu(Z intersection X)>=36/343
  does not imply
mu((Z intersection X) minus B)>0.                     (41)
```

Sections 5--6 repair (41) uniformly for `150` profile rows.  Even there,
positive coefficient-level target spectrum is not yet canonical target
landing.  No thirteen-adic row is removed, the scalar ledger remains
`165`, and LRC(14) remains open.

## 8. Exact companion

The dependency-free companion:

- verifies every rational identity in (25)--(40s);
- checks the exact same-line `1/91` interval geometry;
- exhausts the modulo-fourteen proof of the different-`7`-valuation
  `1/49` law;
- enumerates all fifteen double-pair types and confirms that an ordinary
  label outside the pair always exists;
- verifies the two root-count profiles in (13), the complete-word count
  `648648000`, and the sharper top-labelled `52/78` cell ledgers;
- checks the normalized and unnormalized target Parseval/fourth-power
  sums on every singleton and all thirteen adjacent edges;
- verifies the unique septimal excess word, its thirteen exact
  simultaneous/exclusive owner categories, the `338` joint-cell ledger,
  and the repeated/strict tensor floors in (22e), (38a), and (40f');
- verifies the full one-parameter product trade and the
  `theta=3/4459`, `delta=0` specializations;
- reconstructs the repeated-first and all `150` strict cage caps, checks
  that precisely the `135` rows with `b>=3` pass uniformly, and identifies
  the `b=4` tie family;
- independently enumerates the `52` preliminary middle-depth-two ratios,
  computes the exact compatible four-comb union for both orientations,
  recovers the ten-oriented/six-unordered bank in (40o), and checks the
  `M=2` singleton and empty `M>=3` congruence filters in (40q); and
- checks every balanced first-digit/Bockstein address in (40r)--(40s).

Run

```bash
python3 04-computation/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.py
python3 -O 04-computation/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.out
```

after LF normalization.  Every executable check raises explicitly under
optimized Python.
