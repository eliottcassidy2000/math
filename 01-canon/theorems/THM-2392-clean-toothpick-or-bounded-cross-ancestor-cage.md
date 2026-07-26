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
  blocker cage has mass at least 36/343-delta. Its two same-line pieces have
  exact mass 1/91, so a cross-ancestor pair has an explicit excess overlap;
  if delta<6/4459 its reduced product is at most
  1/(6/4459-delta). In every repeated-first profile (1,1,c), 5<=c<=19,
  THM-2263 bounds the whole cage by 17981/171366 and therefore forces
  delta>=1693/58778538 and a fixed charged cell of mass at least
  1693/22923629820. This is a coefficient-level root word, not a canonical
  expiration target, branch exclusion, row decrement, or proof of LRC(14).
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
script_sha256: PENDING
output_sha256: PENDING
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

The repeated-first thirteen-adic profiles fall strictly on the first side.
The strict profiles retain the quantitative alternative.

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

Subtract (24) from (8) and use the union bound on the four cross pairs
`(i,j)` with `j!=i`.  Some cross pair satisfies

```text
rho(c_i,C_j)
 >=(36/343-delta-2/91)/4

 =1/49+3/8918-delta/4.                                (25)
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
ab<=1/(6/4459-delta).                                 (28)
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
  some cross pair has
    rho>1/49+3/8918-theta/4
  and
    ab<(6/4459-theta)^(-1).                           (30)
```

For example `theta=3/4459` gives

```text
charged-cell mass>=1/579670,

or

cross reduced product ab<=1486.                       (31)
```

At the extreme `delta=0`, equation (25) gives

```text
rho(c_i,C_j)>=185/8918=1/49+3/8918,

ab<=743.                                              (32)
```

The linear trade (25) is the exact output of the cage mass, the two
same-line identities, and the four-pair union bound.  No independence of
the six intersections is assumed.

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
3. The two high cross pairs `(c_i,C_3)` have gap

   ```text
   d=(c-1)-1=c-2 in {3,...,17}.
   ```

   The largest THM-2263 upper endpoint on this range occurs at the even
   gap `d=4` and is

   ```text
   1/49+6/(49*13^4)=583/28561.                        (35)
   ```

Therefore the whole cage obeys

```text
mu(Gamma)
 <=2/91+2*(23/1092)+2*(583/28561)
 =17981/171366.                                       (36)
```

Combining (5)--(8) and (36) gives the unconditional clean-hole floor
within the THM-2388 candidate:

```text
delta
 >=36/343-17981/171366
 =1693/58778538
 >0.                                                   (37)
```

The labelled-owner charged cell from (17) consequently satisfies

```text
rho>=1693/22923629820.                                (38)
```

On that one fixed cell, every nonzero target colour has the coefficient
and energy floors

```text
|integral_Y a_k|
 >=(1693/149003593830)sin(pi/26),

integral_Y |a_k|^2
 >=(1693/968523359895)sin^2(pi/26),                   (39)
```

while the exact summed floors are

```text
sum_(k!=0)integral_Y |a_k|^2
 >=18623/1937046719790,

sum_(k!=0)integral_Y |a_k|^4
 >=52483/327360895644510.                             (40)
```

Thus every repeated-first last-lane packet has a positive literal
two-root charged word.  This is not merely the formal alternative (30).

## 6. What this still does not prove

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
- For strict profiles, (30) is a finite reduced-ratio alternative, not an
  exclusion of either branch.

The first implication which fails in the unqualified THM-2388 finish is

```text
mu(Z intersection X)>=36/343
  does not imply
mu((Z intersection X) minus B)>0.                     (41)
```

Section 5 repairs (41) only for the repeated-first profiles.  Even there,
positive coefficient-level target spectrum is not yet canonical target
landing.  No thirteen-adic row is removed, the scalar ledger remains
`165`, and LRC(14) remains open.

## 7. Exact companion

The dependency-free companion:

- verifies every rational identity in (25)--(40);
- checks the exact same-line `1/91` interval geometry;
- enumerates all fifteen double-pair types and confirms that an ordinary
  label outside the pair always exists;
- verifies the two root-count profiles in (13), the complete-word count
  `648648000`, and the `390/585` cell ledgers;
- checks the normalized and unnormalized target Parseval/fourth-power
  sums on all thirteen adjacent edges;
- verifies the full one-parameter product trade and the
  `theta=3/4459`, `delta=0` specializations; and
- reconstructs the THM-2263 repeated-first cage cap and positive margin,
  including the odd-gap-three versus even-gap-four boundary.

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
