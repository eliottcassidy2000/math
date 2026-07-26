---
id: THM-2440
title: "Sharp two-comb centred-window radii"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. If two
  radius-1/14 integer danger combs cover a pullback window
  {||ny||<rho} with rho>1/14, then both speeds are multiples of n.
  For almost-everywhere coverage rho<=15/182, with equality exactly
  for {n,13n}. For literal pointwise coverage rho<=15/196, with
  equality exactly for {n,14n}. Equivalently, after gcd
  normalization, the open almost-everywhere radius equals the closed
  radius and is at most 15/182 at {1,13}, while the literal connected
  radius of the open union is at most 15/196 at {1,14}.
  Thus the formerly reserved strict 15/182 statement is false:
  {1,13} misses the handoff x=1/14 in both open combs and has literal
  radius only 1/14.  The theorem can obstruct an a.e. two-comb cover
  of a normalized radius-16/182 window, but it does not itself
  produce such a cover, remove a scalar row, close the noncirculant
  graft, or prove LRC(14).
source: codex-2026-07-26-two-comb-centred-radius
depends_on: []
related:
  - THM-1094-exact-two-comb-component-theorem
  - THM-1147-exact-two-comb-gap-law
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - MISTAKE-273
  - MISTAKE-274
script:
  - 04-computation/lrc14_sharp_two_comb_centred_radius_thm2440.py
  - 04-computation/lrc14_sharp_two_comb_strict_radius_referee_thm2440.py
output:
  - 05-knowledge/results/lrc14_sharp_two_comb_centred_radius_thm2440.out
  - 05-knowledge/results/lrc14_sharp_two_comb_strict_radius_referee_thm2440.out
script_sha256:
  - ff6045d45842d2daa5ffb7c97481be7b93495bb303e33b0cc377167c2d5f54fb
  - aa6a27fd5e9bb67cc9c63142b0531ee0d2ef6551ac8ec4ce5d5ac2923fab9bea
output_sha256:
  - 774b94744cc4260be6295ac623d59eda7fdf682c70fcf3428f4997705049085f
  - 06de603f79561ac7c4be1eb8387a0254e3dc5571dc03792840821f34b83ca454
hash_basis: working-tree bytes (LF)
---

# THM-2440 -- the two sharp centred radii are different

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

For a positive integer `m`, put

```text
D_m^o  = {x in R : ||mx|| <  1/14},
D_m^cl = {x in R : ||mx|| <= 1/14}.                         (1)
```

For positive speeds `a,b`, define

```text
r_o(a,b)  = sup{R : (-R,R) is contained in D_a^o union D_b^o},
r_cl(a,b) = sup{R : [-R,R] is contained in D_a^cl union D_b^cl},
r_ae(a,b) = sup{R : (-R,R) is covered a.e. by D_a^o union D_b^o}.
                                                                    (2)
```

Let

```text
g=gcd(a,b),        p=a/g,        q=b/g,        p<=q.          (3)
```

The normalized radii are `rho_*=g r_*(a,b)`.  Equivalently, after the
pullback coordinate `y=gx`, they are the radii of the coprime pair
`{p,q}`.

## Theorem

For every positive pair,

```text
rho_ae = rho_cl <= 15/182,                                   (4)
```

and equality holds exactly when

```text
{p,q}={1,13},        equivalently {a,b}={g,13g}.              (5)
```

For the literal open union,

```text
rho_o <= 15/196,                                             (6)
```

and equality holds exactly when

```text
{p,q}={1,14},        equivalently {a,b}={g,14g}.              (7)
```

In unnormalized coordinates, the equality radii are respectively
`15/(182g)` and `15/(196g)`.

There is a stronger pullback form.  Let `n,A,B` be positive integers,
let `rho>1/14`, and write

```text
W_n(rho)={y in R/Z: ||ny||<rho}.                            (8)
```

Then

```text
W_n(rho) subset_ae D_A^o union D_B^o
  -> n|A, n|B, rho<=15/182,                                (9)

W_n(rho) subset D_A^o union D_B^o
  -> n|A, n|B, rho<=15/196.                               (10)
```

Equality in (9) occurs exactly for `{A,B}={n,13n}` and equality in
(10) exactly for `{A,B}={n,14n}`.  Both equality pairs really cover:
the first almost everywhere and the second pointwise.

## 1. A covered pullback forces the common divisor

Every centre `j/n` of `W_n(rho)` belongs to

```text
D_A^cl union D_B^cl.                                      (11)
```

Otherwise an open neighborhood of that centre would have positive
uncovered measure.  For one speed `u`, put

```text
d=gcd(u,n),                  N=n/d.                        (12)
```

As `j` runs modulo `n`, the values `uj/n` run through the `N`th roots,
each `d` times.  Hence the exact number of centres in `D_u^cl` is

```text
d(2 floor(N/14)+1).                                       (13)
```

If `N>=3`, (13) is at most `n/3`.  If `N=2`, it is `n/2`, and
the hit set is the unique index-two subgroup of the cyclic centre
set.  Two proper masks cannot cover all `n` centres: their cardinality
bound is at most `2n/3`, at most `5n/6`, or exactly the same `n/2`
subgroup according as zero, one, or two masks have `N=2`.  Therefore
at least one speed is divisible by `n`; exchange labels and write

```text
A=an.                                                      (14)
```

Almost-everywhere coverage and the union-measure bound first give
`rho<=1/7`.  For every `t` in the nonempty interval

```text
1/(14a)<t<min(rho,13/(14a)),                               (15)
```

all `n` points

```text
y_j=(j+t)/n,                   j in Z/nZ,                  (16)
```

lie in `W_n(rho)` and avoid `D_A^o`.  If `n` did not divide `B`, the
values `B y_j` would contain a full coset of `N>=2` equally spaced
circle points.  Such a coset cannot lie in the single arc
`||x||<1/14`.  Thus for every `t` in (15), at least one label `j`
also avoids `D_B^o`.  The finite union over labels covers the
positive-length `t`-interval, so some fixed label fails on a set of
positive measure, contradicting (9).  Hence

```text
n|A,                         n|B.                          (17)
```

Literal coverage implies a.e. coverage, so the same divisibility
conclusion applies to (10).  Pullback by `x=ny` now reduces both
statements exactly to the normalized pair problem below.

## 2. Why the a.e. and closed radii agree

On every bounded interval, the difference between
`D_a^cl union D_b^cl` and `D_a^o union D_b^o` is a finite set of tooth
endpoints.  Closed coverage therefore implies open coverage almost
everywhere.

Conversely, if a point in an open centred window lies outside the
closed union, an open neighborhood of that point lies outside it as
well.  That neighborhood has positive measure, contradicting a.e.
coverage by the open union.  Shrinking the window away from its two
outer endpoints gives

```text
r_ae(a,b)=r_cl(a,b).                                        (18)
```

This equality is only an equality of suprema.  It does **not** say that
the open union contains its internal handoff points.

## 3. The fragmentation bound

For every real interval `I` of length `L`,

```text
|I intersect D_m^o| <= L/7 + 6/(49m).                       (19)
```

Indeed, scale `I` by `m` and write `mL=k+s`, where `k` is an integer
and `0<=s<1`.  The `k` full unit periods contribute `k/7`.  The
remaining interval meets the circle danger arc, of length `1/7`, in
at most `min(s,1/7)`.  Its excess above the mean `s/7` is at most

```text
max_(0<=s<=1) (min(s,1/7)-s/7)=6/49.                        (20)
```

Division by `m` proves (19).  This is the sharp discrepancy estimate
also used in THM-1094, but its proof here is self-contained.

Suppose `[0,R]` is covered a.e. by `D_p^o union D_q^o`.  Subadditivity
and (19) give

```text
R <= 2R/7 + 6/(49p) + 6/(49q),

1/p+1/q >= 35R/6.                                          (21)
```

We use (21) only to leave a few elementary boundary branches; no
bounded search supplies an infinite quantifier.

## 4. Sharp closed/a.e. radius

Set

```text
R_ae=15/182,       (35/6)R_ae=175/364.                      (22)
```

Assume `rho_ae>=R_ae`.  Since `p<=q` and `(p,q)=1`, (21) gives:

- `p>=5` is impossible because `1/p+1/q<=2/5`;
- `p=4` is impossible because `q>=5` and
  `1/4+1/5<175/364`;
- for `p=3`, every `q>=7` is impossible because
  `1/3+1/7<175/364`.

The remaining `p=3` pairs are `{3,4}` and `{3,5}`.  Their centred
closed component ends at the central `3`-tooth endpoint `1/42`, well
before `R_ae`.

For `p=2`, its central tooth ends at `1/28` and its next positive
tooth starts at `13/28>R_ae`.  Hence the whole interval

```text
J=(1/28,R_ae),        |J|=17/364.                           (23)
```

would have to lie, up to endpoints, in one `q`-tooth.  Distinct
`q`-teeth have positive gaps, so a connected interval covered a.e.
by `D_q` cannot cross between them.  A tooth has length `1/(7q)`;
thus (23) forces `q<=3`.  The only coprime possibility `q=3` has
central endpoint `1/42<1/28` and next positive tooth starting at
`13/42>R_ae`, so it fails.

It remains to take `p=1`.  The central `1`-tooth ends at `1/14`, so

```text
J=(1/14,R_ae),         |J|=1/91                             (24)
```

must lie in one `q`-tooth.  Hence `q<=13`.  For `q<=12`, the first
positive tooth begins at

```text
13/(14q)>1/14,                                             (25)
```

leaving a gap immediately after the central tooth.  For `q=13`, its
first positive closed tooth is exactly

```text
[13/(14*13),15/(14*13)]=[1/14,15/182].                     (26)
```

Together with the central `1`-tooth, (26) covers `[0,R_ae]` in the
closed convention and covers it a.e. in the open convention.  The
open gap immediately to the right of `15/182` forbids any extension.
Reflection supplies the negative half.  This proves (4)--(5).

## 5. Sharp literal-open radius

Set

```text
R_o=15/196,        (35/6)R_o=175/392.                       (27)
```

If `rho_o>=R_o`, the open window is also covered a.e., so (21)
applies.  It leaves only:

```text
p=4: q=5,
p=3: q in {4,5,7,8},
p in {1,2}.                                                 (28)
```

The pair `{4,5}` stops at `1/56`, and every displayed `p=3` pair stops
at `1/42`: none of their first noncentral teeth begins before `R_o`.

For `p=2`, the interval after its central tooth has length

```text
R_o-1/28=2/49.                                              (29)
```

One `q`-tooth can cover (29) only when `q<=3`; the sole coprime
candidate `{2,3}` stops at `1/28`.

For `p=1`, the required tail has length

```text
R_o-1/14=1/196,                                             (30)
```

so one-tooth coverage forces `q<=28`.  Unlike the a.e. problem, the
internal point `x=1/14` must itself be covered.  The `1`-comb misses
it, while

```text
x=1/14 lies in D_q^o
  iff ||q/14||<1/14
  iff 14 divides q.                                         (31)
```

The only possibilities up to `28` are `q=14,28`.  For `q=14`, the
first positive tooth is

```text
(13/196,15/196).                                            (32)
```

It strictly overlaps the central `1`-tooth and gives literal coverage
of `(-R_o,R_o)`.  For `q=28`, the tooth crossing the seam ends at

```text
29/392 < 30/392=15/196.                                    (33)
```

so it fails.  The gap immediately after (32) prevents extension.
This proves (6)--(7).

## 6. The reserved statement and the LRC scope

The formerly reserved claim combined the first numerical constant
with the second endpoint convention.  Its minimal witness is the
claimed extremizer itself:

```text
{p,q}={1,13},       x=1/14,
||x||=||13x||=1/14.                                        (34)
```

Both strict combs miss (34), so their literal connected radius is only
`1/14`, not `15/182`.  MISTAKE-274 records the correction lineage.

The valid a.e./closed theorem does show that two normalized combs
cannot cover a centred window of radius `16/182`, because

```text
16/182 > 15/182.                                           (35)
```

Using (35) in LRC(14) still requires a lawful reduction that produces
an a.e. two-comb cover of that entire normalized pullback window while
retaining its labels and scale.  THM-2440 does not supply that
reduction, does not decrement any of the `165` scalar rows, and does
not close the surviving `c_3<=M` noncirculant-graft branch.

## 7. Exact referee

The dependency-free `Fraction` referee:

1. reconstructs closed and literal-open centred components directly
   from tooth endpoints;
2. verifies every arithmetic branch in Sections 3--5;
3. reproduces the strict seam witness (34);
4. scans all `9,276` coprime pairs with
   `1<=p<=40` and `p<=q<=400`, finding unique maxima `{1,13}` and
   `{1,14}` in the two conventions; and
5. checks the scaled equality families through `g=50`.

Normal and optimized Python runs execute the same explicit checks and
byte-match the stored output.  The finite scan is a hostile control;
the proof of the universal statement is Sections 1--5.
