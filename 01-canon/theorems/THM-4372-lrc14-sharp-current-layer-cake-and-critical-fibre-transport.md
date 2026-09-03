---
id: THM-4372
title: "LRC(14) sharp-current layer cake and critical septimal-fibre transport"
status: >
  PROVED ANALYTIC + VERIFIED-EXACT; LRC(14) OPEN. The sharp absolute-current
  rebate is exactly one half of a threshold-three event plus one half of a
  third binomial tail moment. The 7-adic reverse-martingale filtration is
  monotone for this convex cost. At an even anchor's critical septimal
  height, every runner becomes one directed unit edge on a seven-point fibre
  and the anchor deletes exactly one vertex. The resulting edge-transport
  budget gives an explicit anchored nonlinear-current lower bound. Coupled
  with THM-4370, every strict counterexample lies in the five-tail upper cone
  and retains at least 4/7 of its deeper-shell rebate. This coefficient is
  sharp for the retained upper-cone fibre data even after restoring the
  residue-labelled edge types. For
  twelve tails, high-shell current at least four forces a positive core
  rebate, while level-three high current can be hidden at all six safe
  vertices of a physical critical fibre. A
  four-tail physical row is the minimal exact obstruction to the stronger
  naive 6/7-survival claim. The theorem supplies a new q-cubic certificate,
  not a closure of the arbitrary one-shell or anchored LRC(14) branch.
source: anchor_current / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope
  - THM-4346-lrc14-halfturn-current-brownian-kernel-and-cubic-moment-boundary
  - THM-4370-lrc14-septimal-wall-quadrature-and-valuation-reanchor
related:
  - THM-2216-residual-capacity-hinge-gram-law
  - HYP-3223-lrc14-green-current-lorentzian-exchange-angles
  - THM-4338-lrc14-all-pair-rank-four-cubic-majorant-and-uniform-one-appender
artifacts:
  - 04-computation/lrc14_anchor_current_tail_filtration_anchor_current_20260903.py
  - 05-knowledge/results/lrc14_anchor_current_tail_filtration_anchor_current_20260903.out
  - 07-reflections/lrc14-anchor-current-tail-filtration-anchor_current_20260903.md
script_sha256: 5a79bcff44bfde4e9d655f1640fe3ee2078abd4e9a5ef65f3f5f201a44dc414f
output_sha256: 6276649efa7fcde418387cd397bae72a2e42e3d8818b29aa0f8076fb685f5cd2
hash_basis: raw LF bytes
audit: >
  PASS. An independent integer dynamic program checks all 169 current-level
  and transport-budget entries through twelve against the closed water-fill
  formula and verifies convexity of every resulting envelope. Exact rational
  wall sweeps check the theorem on 90 structured or deterministic-pseudorandom
  profiles, the minimal four-tail hostile, a twelve-tail hostile, the exact
  six-edge level-three cancellation, the first positive locally sharp fibre,
  the exact global `4/7` upper-cone equality, and its primitive twelve-tail
  `7/2/3` local embedding. Common-dilation tests check every current threshold
  1,...,12 against the exact residual-wall formula at reduced remainders
  0,1,6. Normal, optimized, and nondefault-hash-seed runs reproduce the
  frozen output byte-for-byte. The analytic proof below is not finite-range.
---

# THM-4372 -- sharp-current tails and the critical septimal fibre

**PROVED ANALYTIC + VERIFIED-EXACT. LRC(14) REMAINS OPEN.**

## 1. Definitions and inheritance

Put, up to null endpoints,

```text
A(x)=1_{||x||<1/14},
sigma_w(t)=A(wt)-A(w(t+1/2))                         (1)
```

for a positive odd speed `w`. For a finite family `W`, define the two sheet
depths and their signed current by

```text
alpha(t)=sum_(w in W) A(wt),
beta(t) =sum_(w in W) A(w(t+1/2)),
C(t)=alpha(t)-beta(t)=sum_(w in W) sigma_w(t).        (2)
```

All integrals below use ordinary `dt` on the half-base `[0,1/2)`. Every
display is insensitive to endpoint conventions.

The closest proved mechanisms are THM-4345's exact half-periodic anchor-strip
transducer and THM-4346's divisor-square/reverse-martingale current
filtration. The canonical hostile is the one-shell geometric chain
`{1,9,...,9^11}` from THM-4346. The corrected near miss is to retain only
quadratic current energy: even its exact anchor restriction can miss the
positive nonlinear certificate. The least-used sidecar is the Green-current
view of exact-shell current as vertex divergence. THM-2216 supplies the
closest prior integer layer-cake mechanism; no novelty is claimed for layer
cake in isolation. The new content is its synthesis with the critical
septimal anchor fibre and the resulting transport envelope.

## 2. Exact binomial and tail representation

The sharp current-only rebate from THM-4345 is

```text
g(0)=g(1)=g(2)=0,       g(3)=1/2,
g(d)=d(d^2-6d+11)/12                         (d>=4). (3)
```

For every integer `d>=0`, it has the exact binomial form

```text
g(d)=1/2 [1_{d>=3}+binom(d-1,3)],                     (4)
```

where the binomial term is declared zero below its natural range. Its
successive differences are therefore

```text
Delta g(1)=Delta g(2)=0,
Delta g(3)=1/2,
Delta g(j)=1/2 binom(j-2,2)                  (j>=4). (5)
```

For a measurable region `E`, define the current-tail masses

```text
T_j(E;C)=integral_E 1_{|C|>=j}dt.                    (6)
```

If `|W|=12`, integer layer cake gives the exact identity

```text
integral_E g(|C|)dt
 =1/2 T_3(E;C)
  +1/2 sum_(j=4)^12 binom(j-2,2)T_j(E;C).            (7)
```

For a general finite family the upper limit is `|W|`. Thus the nonlinear
rebate is precisely a threshold-three atom plus a positive third-binomial
tail moment, rather than an unspecified higher-moment correction.

### Proof

For `d>=3`, expansion of `(4)` gives `(3)`; the values `d=0,1,2` are direct.
Alternatively, `(5)` follows by subtraction, since for `j>=4`

```text
g(j)-g(j-1)=(j-2)(j-3)/4=binom(j-2,2)/2.             (8)
```

Now sum `g(d)=sum_(j>=1) Delta g(j)1_{d>=j}` pointwise and integrate. This
proves `(7)`.

## 3. Exact residual-wall transport of every tail

Let `B` be a finite odd family and let `C_B` be its current. Fix an even
anchor `2h` and an odd common multiplier `m`. In THM-4345's notation put

```text
c=gcd(m,h)=gcd(m,2h),
D=m/c=7q+r,       0<=r<=6,
H=2h/c,
R_(H,q,r)={x: ||Hx+q/2||<r/14}.                      (9)
```

The indicator `1_{|C_B|>=j}` is half-periodic, because
`C_B(t+1/2)=-C_B(t)`. Applying THM-4345 separately at every level gives

```text
T_j(anchor strip;C_(mB))
 =[q T_j([0,1/2);C_B)+T_j(R_(H,q,r);C_B)]/D.        (10)
```

Summing `(10)` with the weights `(5)` gives the identical formula for the
whole rebate `integral g(|C|)`. At `r=0`, every tail and the rebate itself
split exactly in the ratio `1:6` between anchor strip and core.

Equation `(10)` is a specialization of the general observable transducer in
THM-4345, not a new common-dilation LRC family. Its role here is to sharpen
the missing sidecar: one needs the nested current-tail vector on the single
residual wall, not merely its quadratic energy.

## 4. Convex order and anchor-measurable martingale steps

Let `Phi` be the even piecewise-linear interpolation of the values
`g(|n|)`, `n in Z`. The nonnegative-side slopes are `(5)` and are
nondecreasing. Reflection across zero shows that `Phi` is convex on `R`.

Split the current into exact septimal speed shells:

```text
W_e={w in W: nu_7(w)=e},
C_e=sum_(w in W_e) sigma_w,
M_e=sum_(j>=e) C_j.                                  (11)
```

By THM-4346,

```text
M_e=Q_e C,
(Q_e f)(t)=7^(-e)sum_(k=0)^(7^e-1)f(t+k/7^e).       (12)
```

The `Q_e` are reverse-martingale conditional expectations. Jensen gives

```text
integral Phi(C)>=integral Phi(M_1)>=integral Phi(M_2)>=... . (13)
```

This is convex order. It does not imply monotonicity of each individual tail
mass in `(6)`.

Now fix the anchor-safe indicator and its septimal height:

```text
G_h(t)=1-A(2ht),              a=nu_7(h).              (14)
```

For `e+1<=a`, `G_h` is invariant under translation by `1/7^(e+1)` and hence
is measurable for `Q_(e+1)`. Conditional Jensen and self-adjointness give

```text
integral G_h Phi(M_e)
 >=integral G_h Phi(M_(e+1)).                         (15)
```

Indeed,

```text
Phi(M_(e+1))=Phi(Q_(e+1)M_e)<=Q_(e+1)Phi(M_e),
```

and multiplying by `G_h` preserves the integral because
`Q_(e+1)G_h=G_h`. Iteration yields the weighted-Jensen reduction

```text
integral G_h g(|C|)>=integral G_h g(|M_a|).          (16)
```

The next step is the first at which the anchor is not measurable. The
critical fibre below resolves part of that boundary exactly.

## 5. Critical septimal-fibre transport theorem

Put

```text
n_a=|W_a|,
Z=M_(a+1)=sum_(nu_7(w)>a) sigma_w.                   (17)
```

For integers `s,n>=0`, define

```text
u=min(n,6(s-2)_+),       u=6q+r,       0<=r<=5,
L_n(s)=(6-r)g(s-q)+r g(s-q-1),                        (18)
```

where the second term is omitted when `r=0`.
Here and below, when a translated argument is nonpositive, extend the integer
cost by the convention `g(k)=0` for every integer `k<=2`. This agrees with
`(3)` on its original domain.

> **Critical-fibre transport bound.** For every finite family of positive
> odd speeds and every positive integer `h`, one has
>
> ```text
> integral_[0,1/2) G_h(t)g(|C(t)|)dt
>  >=(1/7)integral_[0,1/2)L_(n_a)(|Z(t)|)dt.         (19)
> ```

### Proof: one exact-shell runner is one unit edge

By `(16)`, it is enough to bound the core rebate of

```text
M_a=C_a+Z.                                           (20)
```

Both `M_a` and `G_h` are invariant under translation by `1/7^a`. Thus they
descend to that quotient circle, where translation by `1/7^(a+1)` has order
seven. Consider, away from the null wall set, one resulting seven-point fibre

```text
t+j/7^(a+1),                j=0,...,6.               (21)
```

Every speed of valuation greater than `a` changes its phase by an integer on
this fibre, so `Z` is constant. Write `h=7^a h_0`, with `7` not dividing
`h_0`. The seven anchor phases differ by `2h_0/7`; hence exactly one fibre
vertex is anchor-dangerous and the other six are safe.

If `w=7^a v` has exact valuation `a`, then `v` is a unit modulo seven and the
seven phases of `wt` differ by `v/7`. The lower-sheet danger arc has length
`1/7`, so it selects exactly one vertex. The translated upper-sheet danger
arc also selects exactly one, distinct vertex. Thus the vector of `sigma_w`
on the fibre has one `+1`, one `-1`, and five zeros. Equivalently, it is the
divergence of one directed unit edge. Therefore `C_a` is the vertex divergence
of a directed multigraph with `n_a` unit edges.

There is a residue-labelled refinement. If `p` is the positive vertex and
`q` the negative vertex for `w=7^a v`, then

```text
v(q-p)=3 or 4                           (mod 7),      (22)
```

according to which side of the positive tooth centre contains `t`. The
bound below forgets this edge type, but `(22)` is a faithful sidecar for any
future sharpening.

Let the common high current be `z`; reverse all edge orientations if needed
so that `z>=0`. Write the divergences on the six safe vertices as `d_i` and
put `r_i=(-d_i)_+`. Each unit edge supplies only one negative endpoint, so

```text
sum_(i=1)^6 r_i<=n_a.                                (23)
```

Since `g` is even after composition with absolute value and nondecreasing on
the nonnegative integers,

```text
g(|z+d_i|)>=g(max(z-r_i,0)).                         (24)
```

The saving obtained by reducing a current level `k` by one is
`g(k)-g(k-1)`. These savings are nondecreasing in `k` by `(5)`. Hence a
fixed reduction budget is optimally allocated by repeatedly lowering a
currently highest one of the six equal stacks, and no budget is useful after
a stack reaches the zero plateau `0,1,2`. Thus at most `6(z-2)_+` units are
spent; Euclidean division of the spent budget gives exactly `(18)`. The sum
of the six safe costs is at least `L_(n_a)(z)`.

Average this inequality over the fibres `(21)`. Their normalized seven-point
average supplies the factor `1/7`. Every integrand is half-periodic, so
passing from the full circle to the ordinary half-base changes both sides by
the same factor two. This proves `(19)`.

**QED.**

### Equality boundary and general convex costs

The relaxed water-fill value is attained when all useful negative divergence
is placed as evenly as possible on the six safe vertices and the compensating
positive divergence is hidden at the danger vertex. Sections 8.2--8.3 give
physical local attainments at the first zero and positive boundaries.

The fibre proof works for any even convex lattice cost that is nondecreasing
on the nonnegative integers. The six-stack water-fill replaces `(18)`. If
the cost has a zero plateau through level `k`, the useful budget stops at
`6(s-k)_+`. Formula `(18)` is the sharp-current case `k=2`.

## 6. Consequences for deeper shells and twelve tails

### 6.1 Sparse critical shell

For `n<=6`, formula `(18)` (with the just-stated convention at `s=0`) gives

```text
L_n(s)=(6-n)g(s)+n g(s-1)>=(6-n)g(s).                (25)
```

For `n>=6` the corresponding right side is zero. Hence

```text
integral G_h g(|C|)
 >=(6-n_a)_+/7 integral g(|Z|).                      (26)
```

Thus a nonzero high-shell rebate always leaves a positive certified fraction
when the critical shell contains at most five tails.

### 6.2 Deeper reverse-martingale tails

Write `n_a=6q+r`. Extend the integer values of `g` piecewise linearly to a
function on `R`, with value zero on the whole half-line `x<=2`. Formula `(18)`
may then be written without the budget truncation as

```text
L_(n_a)(s)=(6-r)g(s-q)+r g(s-q-1).                  (27)
```

It is a nonnegative combination of translates of a convex nondecreasing
function. Consequently the even extension `L_(n_a)(|x|)` is convex. For
every `b>=a+1`, unweighted conditional Jensen yields

```text
integral G_h g(|C|)
 >=(1/7)integral L_(n_a)(|M_b|).                     (28)
```

The choice `b=a+1` is strongest in convex order; a deeper tail may be easier
to evaluate arithmetically.

### 6.3 A twelve-tail high-current bound

Now assume `|W|=12`. Let `N_>` be the number of tails above height `a`.
Pointwise,

```text
|Z|<=N_>,       n_a+N_><=12,
```

so `n_a<=12-|Z|`. Because `L_n(s)` decreases with `n`, put

```text
Lambda(s)=L_(12-s)(s),       s=0,...,12.             (29)
```

The exact table is

```text
s:       0 1 2 3  4    5  6  7  8     9   10  11  12
Lambda:  0 0 0 0  2  11/2 15 38 78 279/2 227 345 498. (30)
```

Direct comparison with `(3)` gives

```text
Lambda(s)>=2g(s)1_{s>=4}.                            (31)
```

Therefore

```text
integral G_h g(|C|)
 >=(1/7)integral Lambda(|Z|)
 >=(2/7)integral g(|Z|)1_{|Z|>=4}.                  (32)
```

The first bound in `(32)` is itself an exact positive tail functional. The
successive weights of `Lambda` at thresholds `4,...,12` are

```text
lambda_j:
2, 7/2, 19/2, 23, 40, 123/2, 175/2, 118, 153.       (33)
```

Thus

```text
(1/7)integral Lambda(|Z|)
 =(1/7)sum_(j=4)^12 lambda_j T_j([0,1/2);Z).         (34)
```

The threshold boundary is sharp at the fibre level: level three can be
hidden, whereas every level at least four forces a positive rebate under the
twelve-tail budget.

### 6.4 Counterexample-local `4/7` survival

The wall quadrature and the current filtration meet at exactly the same
septimal height. Suppose `|W|=12` and the anchored row `{2h} union W` is
strictly bad. THM-4370 gives

```text
n_a+N_> = #{w in W:nu_7(w)>=a} <= 5.                (34a)
```

If `n_a>=3`, then `N_> <=2`, so `|Z|<=2` pointwise and `g(|Z|)=0`. If
`n_a<=2`, the coefficient in `(26)` is at least `4/7`. In both cases,

```text
integral G_h g(|C|)>=(4/7)integral g(|Z|).           (34b)
```

More precisely, a positive right side is possible only when
`n_a<=2`, `N_> >=3`, and `n_a+N_> <=5`; then `(26)` retains respectively at
least `6/7`, `5/7`, or `4/7` for `n_a=0,1,2`. Thus the six-edge
level-three cancellation below is sharp for the unrestricted transport
theorem, but it cannot occur at the critical height of a strict LRC(14)
counterexample. The lower septimal shells, which THM-4370 forces to contain
at least seven tails, remain outside this rebate.

The coefficient `4/7` is sharp for exactly the upper-cone information retained
here, even after restoring the edge types `(22)`. Take

```text
h=1,        W_a={1,3},        W_>={7,21,35}.         (34c)
```

For every `0<t<1/490`, on the fibre `t+j/7`, with the danger vertex listed
first,

```text
Z=(3,3,3,3,3,3,3),
C_a=(2,-1,0,-1,0,0,0),
C_a+Z=(5,2,3,2,3,3,3).                              (34d)
```

The speed-one edge has `(p,q)=(0,3)` and the speed-three edge has
`(p,q)=(0,1)`. Hence their typed displacements are respectively
`1(3-0)=3` and `3(1-0)=3 mod 7`: both positive endpoints are hidden at the
danger vertex and their negative endpoints are distinct safe vertices. The
safe cost is exactly

```text
2g(2)+4g(3)=2=4g(3).                                (34e)
```

An exact rational wall sweep upgrades this local attainment to the global
five-tail equality

```text
integral g(|Z|)=1/70,
integral G_1 g(|C|)=2/245
                    =(4/7)(1/70)
                    =(1/7)integral L_2(|Z|).        (34f)
```

The valuation profile required by THM-4370 is also physically compatible
with the same local upper-cone fibre. At `h=7`, take

```text
W_<={5,9,11,15,17,29,31},
W_a={7,21},
W_>={49,147,245}.                                   (34g)
```

These are twelve distinct odd speeds and `gcd({14} union W)=1`, with counts
`7/2/3`. At `t=1/1000000`, the upper projection `M_a` is
`(5,2,3,2,3,3,3)`, while the lower-shell current is
`(7,0,-1,0,-1,-1,-1)` and the raw total is
`(12,2,2,2,2,2,2)`, whose six safe rebates all vanish. This does not
contradict `(34b)`, because the lower-shell Jensen step is integrated, not
fibrewise. Neither row is claimed to be a strict counterexample. Therefore a
stronger coefficient conditional on strict badness would have to use global
time coherence or one-sheet information; it cannot arise from a pointwise
critical-fibre refinement based only on the typed edges and lower-shell count.

## 7. Proof-facing q-cubic certificate

Let

```text
p(n)=1-n+binom(n,2)-binom(n,3).                      (35)
```

THM-4345 proves the sharp pointwise inequality

```text
p(min(alpha,beta))
 >=[p(alpha)+p(beta)]/2+g(|C|).                      (36)
```

On the half-base, the measure-preserving reflection `t -> 1/2-t` preserves
`G_h` and swaps `alpha` with `beta`: for odd `w`, periodicity and evenness of
`A` give

```text
A(w(1/2-t))=A(w(t+1/2)),
A(w(1-t))=A(wt).                                    (37a)
```

Define

```text
B_3(G_h)=integral_[0,1/2)G_h(t)p(alpha(t))dt
        =integral_[0,1/2)G_h(t)p(beta(t))dt.         (37)
```

Combining `(19)` and `(36)` gives the anchored certificate

```text
integral_[0,1/2)G_h p(min(alpha,beta))
 >=B_3(G_h)+(1/7)integral_[0,1/2)L_(n_a)(|Z|).       (38)
```

Equations `(26)`, `(28)`, and `(32)` may replace the last term by their
simpler lower bounds. On the strict-counterexample locus, `(34b)` instead
gives

```text
integral G_h p(min(alpha,beta))
 >=B_3(G_h)+(4/7)integral g(|Z|).                    (38a)
```

Thus the critical-fibre transport enters the actual q-cubic route, not
merely an auxiliary current norm.

## 8. Sharp controls and hostile examples

### 8.1 Minimal failure of naive `6/7` survival

The tempting strengthening

```text
integral G_h g(|C|)>=(6/7)integral g(|Z|)            (39)
```

is false for a physical integer-speed row. Take `h=1`, so the anchor is
speed two, and take

```text
W={1,7,21,35}.                                       (40)
```

An exact rational wall sweep gives

```text
high-shell rebate                         =1/70,
core rebate of {7,21,35} alone            =3/245=(6/7)(1/70),
core rebate after adjoining speed 1       =1/98,
slack in (39)                              =-1/490.  (41)
```

The corrected bound `(19)` is equality on `(40)`: `n_a=1` and its right side
is `1/98`. This witness is minimal in tail count. A positive right side in
`(39)` needs at least three high runners because `g` vanishes through current
two. With only those high runners, the reduced degree is seven and `(10)`
gives exact `6/7` core survival. One critical-height runner is therefore the
smallest possible perturbation. For distinct positive odd speeds at `h=1`,
`7,21,35` are the three smallest high runners and `1` is the smallest
critical runner.

A twelve-tail physical control with the smallest possible maximum speed for a
three-runner high shell is

```text
W={1,5,7,9,11,13,17,19,21,29,31,35}.                (42)
```

Its exact values are

```text
integral G_1 g(|C|)=23423/8369690,
slack in (39)=-79063/8369690.                        (43)
```

Thus padding to the physical twelve-tail count does not repair `(39)`.

### 8.2 Complete loss at level three

At `t=1/10000` on the anchor's seven-fibre, the nine-tail family

```text
{1,3,5,7,9,11,13,21,35}
```

has exactly six critical-height runners and constant high current `3`. In the
order danger vertex followed by the six safe vertices, its exact-shell and
total currents are

```text
C_a=(6;-1,-1,-1,-1,-1,-1),
C_a+Z=(9;2,2,2,2,2,2).                              (43a)
```

Every safe rebate is zero. Thus six edges really suffice to erase all six
copies of the level-three high rebate; no positive universal coefficient of
that tail can follow from only the data retained in `(19)` once `n_a>=6`.

The twelve-tail padded control

```text
{1,3,5,7,9,11,13,15,17,19,21,35}                   (44)
```

has the same constant high current and total currents

```text
(12;1,1,1,2,2,2).                                   (45)
```

It is included as a physical twelve-tail control, but the nine-tail row above
is the sharp edge-budget witness.

### 8.3 Sharp first positive level

For

```text
W={1,3,5,7,9,11,13,17,19,21,35,49}                 (46)
```

on the same fibre and at the same time, the constant high current is `4` and
the total fibre currents are

```text
(12;2,2,3,3,3,3).                                   (47)
```

The safe cost is exactly

```text
4g(3)+2g(2)=2=Lambda(4).                            (48)
```

Thus the first positive entry of `(30)`, and the coefficient `2/7` in the
second bound `(32)`, are sharp at the physical fibre level.

### 8.4 Common-dilation and canonical one-shell controls

For the THM-4345 base

```text
B=(1,3,5,7,9,11,13,15,17,19,21,45),       h=420,   (49)
```

the exact program verifies `(10)` simultaneously at every threshold
`1,...,12` for multipliers `1,13,49,127`, covering reduced remainders
`1,6,0,1`. At multiplier `49`, the reduced degree is seven and

```text
strip rebate=7133471/79214135,
core rebate =42800826/79214135=6 times strip rebate. (50)
```

This multiplier-`49` row is deliberately nonprimitive and serves only as the
exact-independence normalization control.

Conversely, the THM-4346 geometric chain

```text
W_9={1,9,...,9^11}                                   (51)
```

occupies one septimal shell. For `h=1`, the high projection `Z` in `(17)` is
zero. If `7|h`, then `a>=1` and the weighted-Jensen projection `M_a` is
already zero. Hence `(19)` is exactly silent on this canonical hostile. The
new bound is a mixed-height certificate, not a solution of the arbitrary
one-shell branch.

## 9. Scope, connection contract, and reproduction

The theorem gives a positive anchored rebate whenever the critical shell has
at most five tails and the high-shell rebate is nonzero, or whenever the
high-shell current reaches absolute level at least four. It does **not**
prove LRC(14):

1. the one-sheet term `B_3(G_h)` in `(38)` can still be negative;
2. the high projection can vanish, as in `(51)`;
3. level-three current can be transported entirely into the deleted vertex
   for the unrestricted theorem, although THM-4370 excludes the required
   six critical tails on the strict-counterexample locus;
4. the relaxed bound forgets the residue-labelled endpoints `(22)` and their
   coherence as `t` moves.

The strongest survivor of the false assertion `(39)` is `(19)`. The missing
coordinate is not another scalar energy but the critical-shell unit-edge
budget. In the Green-current language, the relevant object is integer vertex
divergence with a marked deleted vertex, rather than effective resistance of
an unmarked positive-conductance graph.

```text
source:       sharp current cost plus the 7-adic reverse martingale
target:       nonlinear rebate on the even-anchor-safe core
map:          critical seven-fibre; each exact-shell runner is a unit edge
preserved:    integer current, septimal height, tails, deleted anchor address
destroyed:    residue-labelled endpoints and time-coherent edge motion
sidecar:      n_a and, for a refinement, the edge types in (22)
hostile:      {1,7,21,35} refutes naive 6/7 survival by exactly 1/490
decisive test: exact rational wall sweep and the D=7 common-dilate control.
```

Reproduce the exact audit from the repository root with

```bash
python3 -B 04-computation/lrc14_anchor_current_tail_filtration_anchor_current_20260903.py
```

The frozen output is
`05-knowledge/results/lrc14_anchor_current_tail_filtration_anchor_current_20260903.out`.
