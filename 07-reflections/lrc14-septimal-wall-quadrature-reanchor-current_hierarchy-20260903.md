---
title: "LRC14 septimal wall quadrature, arbitrary CRT intersections, and valuation re-anchoring"
status: >
  PROVED ELEMENTARY + VERIFIED-EXACT, promoted in-session to
  THM-4370; LRC(14) remains open. A consistent probability measure on a
  7-adic subtree of the signed u-walls makes every odd mask have mass exactly
  zero or 1/7 according to valuation, gives an exact even-anchor trichotomy,
  and yields a sharp six-versus-seven cover theorem and a direct
  six-of-twelve anchor-valuation restriction, without requiring the
  auxiliary wall scale to occur in the row. Pairwise-distinct lower
  valuations have exact product intersections; equality at seven forces one
  common valuation. Common-height owners are exact cyclic APs, with a
  complete first-shell classification. A THM-4335 coarse-defect lift closes
  the smallest literal tiler but not every residue-equivalent tiler.
  Arbitrary wall-mask intersections have an exact generalized-CRT formula.
  None of this retains physical tooth scale or completed renewal traces.
source: current_hierarchy / LRC14 continuation session, 2026-09-03
artifacts:
  - 04-computation/lrc14_wall_septimal_quadrature_reanchor_current_hierarchy_20260903.py
  - 05-knowledge/results/lrc14_wall_septimal_quadrature_reanchor_current_hierarchy_20260903.out
related:
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
  - THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law
  - THM-4367-lrc14-active-first-exit-scale-collision-classification
  - THM-4372-lrc14-sharp-current-layer-cake-and-critical-fibre-transport
---

# Septimal wall quadrature and valuation re-anchoring

## 1. Inheritance and the live board

The closest proved mechanism is THM-4348's signed wall set and exact
gcd-quotient capacity formula. Its corrected near miss is equally important:
the wall certificate stopped at three residual masks and declared `7u|h` the
untreated resonance. The canonical hostile is THM-4348's quotient-three
anchor-plus-two partition. The least-used sidecar is not another scalar
capacity, but the natural address projection from the signed `7u`-walls to
the signed `u`-walls.

The live board was

```text
signed wall address    <-> generalized CRT congruence class
7-adic address tower   <-> consistent probability measure
strict endpoint        <-> exact zero-versus-1/7 mass
distinct valuations    <-> exact 7^(-r) intersections
cover equality         <-> common valuation / cyclic AP owner partition
even anchor            <-> valuation trichotomy
resonant old anchor    <-> valuation-matched tail re-anchor
first-shell owner word <-> coarse defect escape
wall mask              <-> physical tooth scale/address warning
```

THM-4363/4365/4367 supply the guardrail on the last line: equal masks or
finite quotient states do not preserve tooth diameters, partial traces,
first exits, or physical binders. The present lane uses a wall only when the
wall itself is a genuine safe physical point. It makes no renewal-chain
inference from a mask identity.

## 2. Exact arbitrary-fold wall-mask calculus

For positive odd `u`, identify the `2u` points of THM-4348's wall set with
addresses

```text
(n,sigma), n in Z/uZ, sigma in {-1,+1},
x_(n,sigma)=(14n+sigma)/(14u) mod 1.                  (1)
```

For a positive speed `s`, put `d=gcd(u,s)`, `q=u/d`, and `S=s/d`. On the
fixed sign sheet `sigma`, strict danger is equivalent to the existence of an
integer

```text
y=14Sn+S sigma-14qm,             -q<y<q.              (2)
```

Thus `y=S sigma mod 14`, and each admissible `y` determines exactly one class

```text
a(y)=S^(-1)(y-S sigma)/14 mod q.                       (3)
```

Let `A_(s,sigma)` be the resulting set of classes in `Z/qZ`. For speeds
`s_1,...,s_k`, a choice `a_i in A_(s_i,sigma)` has a common wall address if
and only if

```text
a_i=a_j mod gcd(q_i,q_j)             for all i,j.     (4)
```

Every compatible tuple has exactly `u/lcm(q_1,...,q_k)` lifts modulo `u`.
Consequently the exact `k`-fold intersection size is

```text
sum_(sigma=+-1) sum_(compatible tuples a)
    u/lcm(q_1,...,q_k).                                (5)
```

Inclusion-exclusion applied to `(5)` is an exact union formula for any
finite residual bank. If the `q_i` are pairwise coprime, compatibility is
automatic and the fixed-sign count factors as

```text
u/(product_i q_i) product_i |A_(s_i,sigma)|.          (6)
```

This answers the formal three-or-more-intersection problem left by THM-4348.
It does not by itself force overlap: generalized CRT can make different
class choices incompatible, and the sharp seven-tail family below has
pairwise disjoint masks on the selected support.

## 3. The 7-adic wall measure

Write

```text
u=7^a v,                    7 does not divide v.       (7)
```

For each sign, let `n_sigma mod v` be the unique solution of

```text
14 n_sigma+sigma=0 mod v.                              (8)
```

Define the support

```text
Y_u={(n_sigma+jv,sigma): sigma=+-1, 0<=j<7^a}.        (9)
```

It has exactly `2*7^a` distinct walls. Give every support wall mass
`1/(2*7^a)` and every other wall mass zero; call this probability measure
`mu_u`.

The construction is a literal filtration. On addressed walls define

```text
pi_(7w,w):X_(7w)->X_w,       (n,sigma)|->(n mod w,sigma). (10)
```

Then

```text
Y_(7w)=pi_(7w,w)^(-1)(Y_w),   (pi_(7w,w))_*mu_(7w)=mu_w, (11)
```

and every supported parent has exactly seven supported children, each with
conditional weight `1/7`. This is a finite reverse-martingale/address-tree
interpretation, not Haar measure on all `X_u`.

For every positive odd speed `s`, the exact law is

```text
mu_u(D_s)= 1/7,  if nu_7(s)<a,
             0,  if nu_7(s)>=a,                       (12)
D_s={x:||sx||<1/14}.
```

To see this, a factor `7` in the speed makes the phase constant on every
seven-child fibre and descends the mask one level. Once the speed is not
divisible by `7`, its phases on the seven children are a full `1/7` grid. An
open interval of length `1/7` contains one grid point unless two consecutive
points are its endpoints. Endpoint equality would give

```text
s(14n+sigma)=+-u mod 14u.                              (13)
```

At a nontrivial septimal level, the right side is zero modulo `7`, whereas
the left side is `s sigma`, nonzero modulo `7`. Hence exactly one child is
strictly dangerous. If every `7`-factor in `u` is first removed from the
speed, the two base walls reduce to odd multiples of `1/14`, which are
equality-safe or farther away. This proves `(12)` including strictness.

## 4. Exact even-anchor trichotomy

Put `b=nu_7(h)`. On the support, `(8)` writes `14n+sigma=vp` with `p` odd
and not divisible by `7`, so

```text
2h x = hp/7^(a+1).                                     (14)
```

Therefore

```text
mu_u(D_(2h)) = 1,    if b>a,
                  0,    if b=a,
                  1/7,  if b<a.                       (15)
```

For `b>a`, `(14)` is integral. For `b=a`, it is a nonzero seventh, hence has
distance at least `1/7`, with equality still safe. For `b<a`, descend the
`b` constant fibres and inspect the next seven-grid fibre. It has exactly one
strict hit. Here endpoint equality is impossible already by parity: the
centered congruence would equate the even integer `2c(14n+sigma)` with the
odd integer `+-u'` modulo the even modulus `14u'`.

Equations `(12)` and `(15)` are stronger than an upper bound. The later cover
conclusions use only the safe inequality `<=1/7`; strict noncoverage comes
from a total mass at most `6/7<1`, never from pretending that `<=1/7` is
itself strict.

## 5. Six versus seven, and the anchor-valuation restriction

No collection of at most six positive odd speeds covers `Y_u`, since its
union has `mu_u`-mass at most `6/7<1`. For `a=0`, `(12)` is stronger: the two
base support walls are simultaneously safe for every odd speed.

The number seven is sharp whenever `7|u`. The seven distinct positive odd
nonmultiples

```text
s_j=1+2ju,                   j=0,...,6,                (16)
```

partition `Y_u` by strict danger. Indeed, on a sign-`sigma` wall their phases
are the seven translates `x+j sigma/7`. No endpoint occurs because `u=0 mod
7`, while `s_j(14n+sigma)=sigma mod 7`.

For an anchored row, the decisive auxiliary scale need not be a row speed.
Given any positive `h`, set

```text
a=nu_7(h),                     u=7^a.                 (17)
```

Then the anchor has mass zero by `(15)`. Every odd tail at valuation at least
`a` also has mass zero by `(12)`, and every lower-valuation tail has mass
exactly `1/7`. Therefore an anchor with any finite odd-tail bank is safe on
`Y_(7^a)=X_(7^a)` whenever at most six tails have valuation below
`nu_7(h)`. In particular, a strict failure with twelve odd tails forces at
least seven lower-valuation tails, or equivalently at most five tails in the
upper cone

```text
{s:nu_7(s)>=nu_7(h)}.                                 (18)
```

This is the headline six-of-twelve restriction. It does not assume that a
tail has valuation exactly `nu_7(h)`. When `7` does not divide `h`, `a=0`
and there are no lower-valuation odd tails, so the fixed walls `+-1/14`
already give safety.

For comparison, if one deliberately anchors the quadrature to an odd row
speed `u`, one gets the more general local alternatives:

- if `nu_7(h)=a`, the anchor has mass zero, so at most six lower-valuation
  other tails cannot cover;
- if `nu_7(h)<a`, the anchor costs exactly `1/7`, so at most five
  lower-valuation other tails cannot cover;
- if `nu_7(h)>a`, this particular support is wholly anchor-dangerous and
  gives no certificate.

There is also a tail-based re-anchor specialization. Suppose a `2+12` row has
twelve distinct positive odd tails `W` and contains `r` with

```text
nu_7(r)=nu_7(h).                                       (19)
```

If at most six members of `W\{r}` have lower `7`-valuation, the row is safe
on `Y_r`. Thus strict failure forces at least seven lower-valuation tails for
every `r` satisfying `(19)`. Equivalently, if one matching tail exists, a
strictly bad row can have at most five of its twelve tails at valuation at
least `nu_7(h)`. Six of twelve in that upper valuation cone force safety.

There is a purely divisibility-based specialization, useful when entering
from THM-4348's nesting. If at most six other tails are nonmultiples of `r`,
then every remaining tail has the form `cr`. Since both it and `r` are odd,
the integer multiplier `c` is odd. Hence on every `r`-wall

```text
cr(14n+sigma)/(14r)=cn+c sigma/14                     (19a)
```

is strictly safe. The anchor is safe by `(19)`, and the at most six
nonmultiples have total mass at most `6/7<1`. Therefore the row is safe. In a
strictly bad `2+12` row, every valuation-matched tail has at least seven
nonmultiples among the other eleven; in particular, six tails divisible by
one such `r` are impossible.

This applies even if a previously selected base `u` lies in THM-4348's total
resonance `7u|h`: re-anchor at a different tail `r` satisfying `(19)` and use
the exact physical `r`-walls.

## 6. Mixed-height independence and equality collapse

The filtration contains more than the one-mask mass law. If positive odd
speeds `s_1,...,s_r` have pairwise distinct valuations below
`a=nu_7(u)`, then on each sign sheet of `Y_u`

```text
|Y_u^sigma intersection intersection_i D_(s_i)|=7^(a-r),
mu_u(intersection_i D_(s_i))=7^(-r).                 (19b)
```

The mechanism is conditional, not global probabilistic independence. Sort
the valuations. At the first nonconstant child fibre, the least-valuation
speed takes exactly one child while every higher-valuation mask is constant.
Remove that child layer, divide the higher speeds by `7`, and iterate. The
corrected endpoint congruence from the theorem,

```text
s(14(n+jw)+sigma)=+-7w mod 98w,                     (19c)
```

is essential: modulo `7` it is impossible for the `7`-unit speed, so the
open mask has one child rather than zero.

If exactly seven lower masks cover the support, their seven masses sum to
one, so equality in the union bound makes them pairwise disjoint. Two
different valuations would instead meet in mass `1/49` by `(19b)`.
Therefore every exact seven-cover has one common lower valuation `e`. In a
primitive critical `2+12` row with exactly seven lower tails, the other five
tails and anchor all have valuation at least `a>e`; hence the common divisor
`7^e` forces `e=0`. This argument intentionally stops when there are eight
or more lower tails, where equality need not hold.

## 7. The common-height owner partition

Write a common-height tail as `s=7^e r` and put `U=7^(a-e)`. With
`delta_r` the unique member of `{-5,-3,-1,1,3,5}` congruent to `r mod 14`,
the exact positive-sheet owner block is

```text
J_r={ell:-U<delta_r+14ell<U},
A_r^+=r^(-1)(J_r+(delta_r-r)/14) mod U,
A_r^-=-A_r^+.                                       (19d)
```

The set `J_r` contains `U/7` consecutive integers, so each block is a cyclic
AP of length `U/7` and unit step `r^(-1)`. Seven masks partition both sheets
if and only if the seven positive blocks partition `Z/U`.

At the first shell `U=7`, the block is one point

```text
kappa(r)=r^(-1)(delta_r-r)/14 mod 7.                 (19e)
```

The complete odd-unit residue atlas modulo `98` is

```text
kappa  residue classes
0      +-1, +-3, +-5
1      +-13, +-33, +-39
2      +-17, +-27, +-37
3      +-9, +-25, +-41
4      +-19, +-31, +-43
5      +-11, +-29, +-47
6      +-15, +-23, +-45.                            (19f)
```

Thus a literal seven-bank tiles the wall exactly when it selects one class
of each colour. The smallest positive choice is

```text
L_0=(1,9,11,13,15,17,19),
kappa(L_0)=(0,3,5,1,6,2,4).                         (19g)
```

The deeper computation is intentionally scoped as **FINITE-EXACT / OPEN**.
At `U=49`, all `147` distinct blocks have exactly `21` seven-partitions; at
`U=343`, all `1,029` blocks have exactly `147`. Every one is a parallel lift
with signed representatives forming `R+2U Z/7` modulo `14U`. A targeted
`U=2401` search among `7,203` blocks finds `46` disjoint partners for the
speed-one block and its unique exact completion

```text
(1,4801,4803,9603,9605,14405,14407).                (19h)
```

This strongly suggests a cyclic-AP tiling rigidity theorem, but does not
prove one. Pairwise disjointness has many more solutions than complete
seven-way tiling, so no all-depth statement was promoted.

## 8. Coarse defect: the physical coordinate after wall equality

At `a=1`, write every physical time above `y` as `t_j=(y+j)/7`. For a lower
unit `s`, let

```text
O_s(y)={j:||s(y+j)/7||<1/14}.                       (19i)
```

This has at most one element; a half-integer equality gets no owner. Let
`P_L` be the locus where seven lower tails own a permutation of the seven
lifts, and let `E_L` be its complement. THM-4335's owner-permutation theorem
specializes exactly as follows. For `h=7H` and upper tails
`7q_1,...,7q_m`, `m<=5`, strict failure is equivalent to

```text
E_L subset D_(2H) union D_(q_1) union ... union D_(q_m). (19j)
```

The reason is literal fibre geometry: on `P_L` the lower tails cover every
lift; on `E_L` there is a free lift, while the anchor and upper tails are
constant across the fibre with phases `2Hy,q_1y,...,q_my`. Thus failure
requires

```text
lambda(P_L)>=(6-m)/7.                                (19k)
```

For `L_0`, an exact sweep of its `77` distinct internal half-integer
breakpoints (`78` open cells) gives only two permutation intervals:

```text
P_(L_0)=(1/18,3/38) union (35/38,17/18),
lambda(P_(L_0))=8/171,
lambda(E_(L_0))=163/171>6/7.                         (19l)
```

The owner vectors on the two cells are `(0,3,5,1,6,2,4)` and
`(6,3,1,5,0,4,2)`. Hence this literal lower family can never be completed to
a bad row by an anchor plus at most five upper-cone tails.

It is nevertheless an exact hostile to the compressed THM-4370/4372 data.
For anchor `14`, lower tails `L_0`, and upper tails
`(7,21,35,63,77)`, every `X_7` wall is owned by exactly one lower tail while
the anchor and upper tails are safe. In THM-4372 notation, `n_a=5` and
`Z=0`, so the retained current transport is zero. The restored coarse defect
finds the escape; explicitly

```text
t=15/182,       clearance=1/14,       unique binder=13. (19m)
```

The physical warning is sharp. The residue-equivalent parallel tiler
`(1,15,29,43,57,71,85)` has permutation mass
`7122/46835>1/7`, so the scalar defect measure no longer closes it. The
mod-98 wall atlas does not determine the coarse owner evolution: adding `98`
preserves a wall mask but changes physical tooth scale and address.

## 9. Physical controls and stopping boundary

First, a control shows that the headline result genuinely needs no
exact-height tail. Take `h=420` and

```text
W_0=(1,3,5,9,11,13,49,147,245,343,441,539).          (20)
```

Exactly six tails have valuation below `nu_7(h)=1`, the other six have
valuation strictly above it, and none has valuation exactly one. At the
auxiliary wall `u=7`,

```text
t=13/98,        min_(s in {840} union W_0)||st||=13/98>1/14. (21)
```

For `h=420`, take the twelve odd tails

```text
W=(1,7,21,35,63,77,91,105,119,133,147,161).          (22)
```

The original tail `u=1` lies in `7u|h`. Re-anchor at `r=7`. Of the other
eleven tails, ten are multiples of `7` and one is a nonmultiple. Their
integer quotients by `7` are odd. The first audited safe wall is

```text
(n,sigma)=(1,-1),       t=13/98,
min_(s in {840} union W)||st||=1/14,                  (23)
```

with physical binding speeds exactly `7,91,105`. Thus the theorem does not
merely find a quotient state; it returns a real safe time with its binders.

Three hostile controls delimit the claim.

1. For `7|u`, `(16)` is an exact seven-mask cover of the support. The
   six-tail threshold cannot be increased by this measure.
2. At `u=5,h=7`, the anchor and four odd residual speeds
   `(71,81,83,87)` cover all ten walls `X_5`. This lies in the `b>a` branch
   where `(15)` correctly offers no help.
3. At `u=11,h=420`, the anchor and five lifted odd residual speeds
   `(159,331,523,665,835)` cover all twenty-two walls `X_11`; the complete
   displayed speed set including `u=11` is primitive. This prevents a
   blanket anchor-plus-five improvement without a valuation hypothesis.

Finally,

```text
M_u(s+14uL)=M_u(s)                                     (24)
```

for every integer `L`, so the hostiles can be moved to larger physical
speeds without changing their wall masks. But the danger-tooth diameter is
`1/(7s)`, and its address and reach also change. In the language of
THM-4363/4365/4367, `(24)` preserves the finite mask while destroying the
physical scale and potentially the partial renewal trace. The theorem stops
at exact wall safety and exact mask intersection; it does not promote mask
equivalence to component coverage, greedy-chain equality, or first-exit
equality.

## 10. Exact audit

The script checks:

- `35,000` direct masks against `(2)--(3)`, including `26,424` strict
  endpoint occurrences;
- `4,065` arbitrary intersections and `1,500` inclusion-exclusion unions,
  with folds through six, against direct wall enumeration;
- `107,143` odd-speed valuation laws and `192,000` even-anchor trichotomy
  cases for every odd `u<400` in the declared bounded banks, including
  `10,800` anchor rechecks through an independent Fraction-circle path;
- all address projections and seven-child fibre sizes in that universe;
- `57` sharp seven-tail partitions for `7<=u<800`, also checked pointwise
  through the independent Fraction-circle path;
- `3,104` pairwise-distinct-valuation intersections through fold four on
  composite support scales, including `152` independent Fraction-circle
  rechecks and exact counts on each sign sheet;
- the complete `42`-residue mod-`98` owner atlas, `2,450` formula/direct
  owner-block comparisons, exhaustive `U=49,343` cover enumeration, and the
  targeted speed-one completion at `U=2401`;
- the exact `8/171` and `7122/46835` permutation measures by two independent
  owner predicates, the retained-invariant hostile and its physical witness,
  and `26,440` exact fibre checks of the coarse-defect equivalence;
- the no-exact-height control, physical tail re-anchor witness, and both
  hostile cover families.

Normal, optimized, and hash-seeded replays match the frozen output
byte-for-byte.
