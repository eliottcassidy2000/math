---
id: THM-2392
title: "Clean toothpick or bounded cross-ancestor cage"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. Assuming
  THM-2388's last-lane hypotheses and exact 36/343 hole/excess ledger, let
  delta be the mass of unit holes which are outside the quotient-blocker
  union. Every such parent has exactly one double-covered inverse root,
  guard root count four, and five ordinary root counts two. A deterministic
  ordinary label outside the double pair gives a literal adjacent two-root
  word whose normalized nonzero target spectrum is everywhere nonzero,
  with exact square/fourth-power sums 22/169 and 62/28561. A labelled-owner,
  pair, and translate cell has mass at least delta/390. Conversely, the
  blocker cage has mass at least 36/343-delta. Different 7-adic valuations
  force exact danger overlap 1/49; the last-lane hypotheses apply this to
  both high cross-ancestor pieces. Thus if delta<6/4459, one of the two
  remaining low cross pairs has reduced product at most
  1/[2(6/4459-delta)]. Every repeated-first profile forces
  delta>=1/26754. Among the 150 strict profiles, all 135 rows with b>=3
  force delta>=6042/9796423 and a fixed charged cell of mass at least
  1007/636767495. The fifteen b=2 rows either have positive clean mass or
  their one unresolved reduced cross pair lies in an explicit 124-ratio
  bank with product at most 345. This is a coefficient-level root word,
  not a canonical expiration target, branch exclusion, row decrement, or
  proof of LRC(14).
source: codex-2026-07-26-clean-toothpick-cross-ancestor-cage
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-sidecar
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
script: 04-computation/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.py
output: 05-knowledge/results/lrc14_clean_toothpick_cross_ancestor_cage_thm2392.out
script_sha256: 1cb8732b49d8b5541db41dd0d5dd5926b6bdc96298626b04d42d59eb4892d24f
output_sha256: 2b325668958a7329b96c9916ed91d5818fb6836fef7462d4f3fca01d119fb170
hash_basis: working-tree bytes (LF)
---

# THM-2392 -- a clean toothpick or a bounded blocker cage

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2388 finds exact signed mass in the final septimal
`k=2,(t,b)=(1,0)` lane, but it does not put that mass outside the
quotient blockers.  This theorem keeps that missing bit rather than silently
discarding it.  The result is the exact alternative

```text
positive clean-hole mass
  -> one double root + a labelled adjacent two-root charged word;

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

## 3. A literal two-root charged word

Fix a deterministic order on the five ordinary labels.  The double pair
contains at most two of them, so choose the first ordinary label `q`
outside that pair.  Its two roots in (12) are both singleton roots.
On the `q`-scaled thirteen-grid they are adjacent.  Hence their indicator
has the form

```text
A(r)=1_({r_0,r_0+eta})(r),             eta in F_13^*,

W(r)=1-A(r).                                             (15)
```

There are thirteen possible translates of this adjacent edge.  Partition
`S` by:

```text
chosen low-blocker owner:       2 possibilities,

double-pair label:             15 possibilities,

selected-q edge translate:     13 possibilities.       (16)
```

One literal charged cell `Y` therefore has

```text
rho:=mu(Y)>=delta/(2*15*13)=delta/390.                 (17)
```

If the complete low-blocker word rather than one deterministic owner is
required, the corresponding bound is `delta/585`.

Use the normalized target DFT

```text
a_k=(1/13)sum_(r in F_13)A(r)zeta^(-kr),

w_k=(1/13)sum_(r in F_13)W(r)zeta^(-kr),

zeta=exp(2 pi i/13).                                   (18)
```

For every `k!=0`,

```text
w_k=-a_k,

|a_k|^2
 =(2+2cos(2 pi k eta/13))/169
 >0.                                                   (19)
```

The strict inequality uses oddness of thirteen: `-1` is not a thirteenth
root of unity.  The exact all-colour ledgers are

```text
sum_(k!=0)|a_k|^2=22/169,

sum_(k!=0)|a_k|^4=62/28561.                            (20)
```

Consequently, on the fixed cell `Y`,

```text
integral_Y a_k conjugate(w_k)
 =-rho |a_k|^2<0,

|integral_Y a_k|
 >=(2rho/13)sin(pi/26),                               (21)

integral_Y |a_k|^2
 >=(4rho/169)sin^2(pi/26),                            (22)
```

and the sums in (20) acquire the factor `rho`.

Normalization is load-bearing.  With the unnormalized DFT
`A_k^sharp=sum_r A(r)zeta^(-kr)`, the two sums are `22` and `62`, not
`22/169` and `62/28561`.

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
  a fixed labelled charged cell has mass at least theta/390;

delta<theta:
  one of the two low cross pairs has
    rho>1/49+3/4459-theta/2
  and
    ab<[2(6/4459-theta)]^(-1).                        (30)
```

For example `theta=3/4459` gives

```text
charged-cell mass>=1/579670,

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
rho>=1/10434060.                                      (38)
```

On that one fixed cell, every nonzero target colour has the coefficient
and energy floors

```text
|integral_Y a_k|
 >=sin(pi/26)/67821390,

integral_Y |a_k|^2
 >=sin^2(pi/26)/440839035,                            (39)
```

while the exact summed floors are

```text
sum_(k!=0)integral_Y |a_k|^2
 >=11/881678070,

sum_(k!=0)integral_Y |a_k|^4
 >=31/149003593830.                                   (40)
```

Thus every repeated-first last-lane packet has a positive literal
two-root charged word.  This is not merely the formal alternative (30).

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
 >=(6042/9796423)/390
 =1007/636767495.                                     (40f)
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

rho(labelled charged cell)>=6/3767855.                (40h)
```

It remains to record the exact boundary when the two valuations in (40g)
are equal.  Divide the unresolved pair by its gcd and write it, up to
interchange, as

```text
(C_2,c_1)=g(a,b),

a<=b, gcd(a,b)=1, gcd(ab,91)=1.                       (40i)
```

The other five cage pieces have total upper bound

```text
2/91+25/1183+2/49=695/8281.                           (40j)
```

Therefore the extreme no-clean-hole case `delta=0` forces

```text
rho(a,b)>=36/343-695/8281
          =1219/57967,

Delta(a,b)/(ab)>=144/1183.                            (40k)
```

The coarse THM-1166 endpoint bound first gives `ab<=402`.  Exact
enumeration of the finite universe

```text
a<=b, gcd(a,b)=1, gcd(ab,91)=1, ab<=402,

[F((a+b) mod 14)-F((b-a) mod 14)]/(ab)>=144/1183     (40l)
```

leaves exactly `124` unordered pairs.  Across this bank

```text
ab<=345,                 min(a,b)<=11,

max(a,b)<=197.                                        (40m)
```

The product maximum occurs at `(3,115)` and the coordinate maximum at
`(1,197)`.  Thus each middle-depth-two profile either has positive
clean-hole mass, or its one unresolved cross ancestor lies in this
explicit finite bank.  If its low blockers have different `7`-adic
valuations, the stronger uniform positive alternative (40h) applies.

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
- THM-2390 forces a separate labelled seven-root partition or one-double
  word.  It supplies no positive intersection with `S` and no alignment of
  its septimal root with the thirteen-root toothpick.
- For the last middle-depth-two profiles, (40l) is a finite reduced-ratio
  bank, not an exclusion of either branch.

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

- verifies every rational identity in (25)--(40m);
- checks the exact same-line `1/91` interval geometry;
- exhausts the modulo-fourteen proof of the different-`7`-valuation
  `1/49` law;
- enumerates all fifteen double-pair types and confirms that an ordinary
  label outside the pair always exists;
- verifies the two root-count profiles in (13), the complete-word count
  `648648000`, and the `390/585` cell ledgers;
- checks the normalized and unnormalized target Parseval/fourth-power
  sums on all thirteen adjacent edges;
- verifies the full one-parameter product trade and the
  `theta=3/4459`, `delta=0` specializations; and
- reconstructs the repeated-first and all `150` strict cage caps, checks
  that precisely the `135` rows with `b>=3` pass uniformly, and identifies
  the `b=4` tie family; and
- independently enumerates the `124` middle-depth-two residual ratios and
  checks the sharp product/coordinate witnesses in (40m).

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
