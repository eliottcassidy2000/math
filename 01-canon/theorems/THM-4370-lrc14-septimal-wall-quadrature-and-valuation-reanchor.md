---
id: THM-4370
title: "LRC(14) septimal wall quadrature, arbitrary CRT calculus, and valuation re-anchor"
status: >
  PROVED ELEMENTARY + VERIFIED-EXACT; LRC(14) OPEN. A canonical uniform
  measure on 2*7^nu7(u) signed u-walls gives every positive odd speed
  strict-danger mass exactly 1/7 below the u-valuation and zero at or above
  it. The even anchor has exact masses 1/7,0,1 as its valuation lies below,
  at, or above the u-valuation. Taking the auxiliary scale
  u=7^nu7(h), whether or not it is a row speed, proves that strict failure in
  the anchored 2+12 branch requires at least seven tails below the anchor's
  7-valuation, so at most five tails lie in the upper valuation cone. Six odd
  masks never cover the support, while an explicit seven-tail family
  partitions it. Pairwise-distinct lower valuations have exact product
  intersections, forcing every exact seven-cover to one common valuation.
  Common-height masks are cyclic owner APs, completely classified at the
  first shell. A coarse owner-defect criterion then closes the smallest
  literal first-shell tiler. Arbitrary finite wall-mask intersections and
  unions have an exact generalized-CRT formula. These are pointwise wall and
  owner certificates, not physical renewal-chain or first-exit
  classifications.
source: current_hierarchy / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
related:
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
  - THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law
  - THM-4367-lrc14-active-first-exit-scale-collision-classification
  - THM-4372-lrc14-sharp-current-layer-cake-and-critical-fibre-transport
script: 04-computation/lrc14_wall_septimal_quadrature_reanchor_current_hierarchy_20260903.py
output: 05-knowledge/results/lrc14_wall_septimal_quadrature_reanchor_current_hierarchy_20260903.out
reflection: 07-reflections/lrc14-septimal-wall-quadrature-reanchor-current_hierarchy-20260903.md
script_sha256: f021693acabeebf178ae7190fbad63967fd081b9f38187dab5d6ebc4152e9220
output_sha256: e27bff14b1a556fb6bcfa0ca8e81bef089089cbcdf57643e1f88f30f324d91b4
hash_basis: raw LF bytes
audit: >
  PASS. Integer/Fraction-exact controls compare 35,000 masks with the
  quotient formula, including 26,424 strict endpoint occurrences; compare
  4,065 generalized-CRT intersections and 1,500 inclusion-exclusion unions
  through fold six with direct enumeration; verify 107,143 odd-speed mass
  laws, 192,000 even-anchor trichotomy cases, 10,800 of those by an
  independent Fraction-circle path, and every support projection for odd
  u<400 in the declared banks; verify 57 sharp seven-tail partitions by both
  centered-residue and independent Fraction paths; verify 3,104
  distinct-valuation intersections through fold four, including 152
  independent Fraction checks; classify all 42 first-shell unit residues,
  compare 2,450 owner blocks by independent paths, exhaust the U=49 and 343
  partitions, and target the speed-one U=2401 completion; verify both exact
  coarse permutation measures, 26,440 fibre implications, the retained-data
  hostile and its physical witness; and reconstruct the earlier controls.
  Normal, optimized, and hash-seeded runs reproduce the frozen output
  byte-for-byte.
---

# THM-4370 -- septimal wall quadrature and valuation re-anchor

**PROVED ELEMENTARY + VERIFIED-EXACT. LRC(14) REMAINS OPEN.**

## 1. Signed walls and arbitrary finite intersections

For a positive odd integer `u`, use THM-4348's `2u` distinct physical walls

```text
X_u={x_(n,sigma)=(14n+sigma)/(14u) mod 1:
     n in Z/uZ, sigma in {-1,+1}}.                    (1)
```

For a positive integer speed `s`, its strict wall mask is

```text
M_u(s)={x in X_u:||sx||<1/14}.                        (2)
```

Put

```text
d_s=gcd(u,s),             q_s=u/d_s,       S_s=s/d_s. (3)
```

For each sign `sigma`, define `A_(s,sigma) subset Z/q_s Z` as follows. For
every integer `y` satisfying

```text
-q_s<y<q_s,               y=S_s sigma mod 14,         (4)
```

include the unique residue `a mod q_s` satisfying

```text
S_s a=(y-S_s sigma)/14 mod q_s.                       (5)
```

The right side is integral by `(4)`, and the residue is unique because
`gcd(S_s,q_s)=1`. The definition also has its unique meaning in the trivial
ring `Z/Z` when `q_s=1`.

> **Arbitrary-fold CRT wall formula.** Let `s_1,...,s_k` be any finite list
> of positive integer speeds. For fixed `sigma`, choose
> `a_i in A_(s_i,sigma)`. Call the tuple compatible when
>
> ```text
> a_i=a_j mod gcd(q_(s_i),q_(s_j))       for every i,j. (6)
> ```
>
> Then
>
> ```text
> |intersection_i M_u(s_i)|
>  =sum_(sigma=+-1) sum_(compatible (a_1,...,a_k))
>      u/lcm(q_(s_1),...,q_(s_k)).                    (7)
> ```
>
> Inclusion-exclusion using `(7)` gives the exact union size for every
> finite mask bank. If the `q_(s_i)` are pairwise coprime, compatibility is
> automatic and the contribution of sign `sigma` is
>
> ```text
> u/(product_i q_(s_i)) product_i |A_(s_i,sigma)|.    (8)
> ```

### Proof

At `x_(n,sigma)`, strict danger is equivalent to an integer `m` satisfying

```text
|s(14n+sigma)-14um|<u.                                (9)
```

Divide by `d_s` and set

```text
y=14S_s n+S_s sigma-14q_s m.                         (10)
```

Then `(9)` is exactly `(4)`, while division of
`y-S_s sigma` by `14` gives `(5)`. Distinct admissible `y` give distinct
classes modulo `q_s`: equality of their classes would make their difference
a nonzero multiple of `14q_s`, impossible inside the interval of length
`2q_s` in `(4)`.

The generalized Chinese remainder theorem says that the simultaneous
classes `n=a_i mod q_(s_i)` are soluble exactly under `(6)`. Every `q_(s_i)`
divides `u`, so each compatible tuple has one solution modulo their lcm and
therefore exactly `u/lcm_i(q_(s_i))` lifts modulo `u`. Different class tuples
have disjoint solution sets. Sum over the two disjoint sign sheets to obtain
`(7)`. Inclusion-exclusion and the pairwise-coprime specialization are
immediate. **QED.**

## 2. Canonical septimal wall quadrature

Write

```text
u=7^a v,                    a=nu_7(u),       7 does not divide v. (11)
```

For each `sigma in {-1,+1}`, let `n_sigma mod v` be the unique solution of

```text
14n_sigma+sigma=0 mod v.                              (12)
```

For `v=1`, this is the unique class `0`. Define

```text
Y_u={x_(n_sigma+jv,sigma):
     sigma=+-1, 0<=j<7^a}.                            (13)
```

These are exactly `2*7^a` distinct points of `X_u`. Define the probability
measure

```text
mu_u=(1/(2*7^a)) sum_(x in Y_u) delta_x.              (14)
```

Thus every support point has weight exactly `1/(2*7^a)`; all other points of
`X_u` have weight zero.

This measure is consistent along the address tower. For positive odd `w`,
define

```text
pi_(7w,w):X_(7w)->X_w,
pi_(7w,w)(x_(n,sigma))=x_(n mod w,sigma).             (15)
```

Then

```text
Y_(7w)=pi_(7w,w)^(-1)(Y_w),
(pi_(7w,w))_*mu_(7w)=mu_w.                            (16)
```

Every point of `Y_w` has the seven children with addresses
`(n+jw,sigma)`, `0<=j<7`, each of conditional weight `1/7`.

> **Odd-speed quadrature law.** For every positive odd integer `s`,
>
> ```text
> mu_u(M_u(s)) = 1/7,       if nu_7(s)<a,
>                    0,       if nu_7(s)>=a.           (17)
> ```

### Proof

On the seven children of `x_(n,sigma) in X_w`, a speed not divisible by `7`
has phases

```text
s x_(n+jw,sigma)=s x_(n,sigma)/7+j s/7 mod 1.        (18)
```

Thus the child phases form a full seven-point grid. An open circular interval
of length `1/7` contains exactly one grid point unless two adjacent grid
points are its two endpoints. Endpoint equality at child `n+jw` would imply

```text
s(14(n+jw)+sigma)=+-7w mod 98w.                       (19)
```

The right side is zero modulo `7`, while the left side is `s sigma`, nonzero
modulo `7`. Hence every supported parent has exactly one strictly dangerous
child.

If `s=7s'`, its phase is constant on a child fibre and equals the `s'` phase
at the projected parent:

```text
7s' x_(n+jw,sigma)=s' x_(n mod w,sigma) mod 1.        (20)
```

Therefore each common factor `7` in `s` and `u` descends both one level with
mass unchanged. If `nu_7(s)<a`, descend until `(18)--(19)` apply, giving
mass exactly `1/7`. If `nu_7(s)>=a`, descend all `a` levels. At the two base
support points, `(12)` writes `x=p/14` with `p` odd. Multiplication by the
remaining odd speed gives an odd multiple of `1/14`, whose circle distance
is at least `1/14`. Strict danger is absent. This proves `(17)`. **QED.**

The tower `(15)--(16)` may equivalently be viewed as a finite filtration:
conditional on a supported parent, a first lower-valuation odd mask occupies
exactly one of seven children. This interpretation does not assert
independence between different speeds.

## 3. Exact even-anchor trichotomy

Let `h` be positive and put `b=nu_7(h)`. Then

```text
mu_u(M_u(2h)) = 1,       if b>a,
                       0,       if b=a,
                       1/7,     if b<a.               (21)
```

### Proof

At a supported address, `(12)` writes

```text
14n+sigma=vp,                 p odd, 7 does not divide p. (22)
```

Consequently

```text
2h x_(n,sigma)=hp/7^(a+1).                            (23)
```

If `b>a`, this is integral, so every support point is strictly dangerous. If
`b=a`, it is a nonzero seventh modulo one and its distance is at least
`1/7`; all support points are safe.

If `b<a`, descend the `b` constant fibres as in `(20)`. The remaining even
speed is not divisible by `7`, so every next seven-child fibre has exactly
one hit unless an endpoint occurs. Such an endpoint would require an even
integer of the form `2c(14n+sigma)` to be congruent to the odd integer
`+-u'` modulo the even modulus `14u'`, impossible by parity. Thus exactly one
of seven children is strictly dangerous and the mass is `1/7`. **QED.**

The endpoint semantics are load-bearing: `(17)` and `(21)` are exact
equalities. In the cover arguments below, each relevant mask has mass
`<=1/7`, while strict noncoverage comes from the separate strict inequality
`6/7<1`.

## 4. Sharp six-versus-seven wall theorem

No collection of at most six positive odd speeds covers `Y_u`. Indeed, the
union bound and `(17)` give mass at most `6/7<1`. If `7` does not divide `u`,
then `a=0` and the stronger statement holds: both points of `Y_u` are
simultaneously safe for every positive odd speed.

The number seven is sharp for this support whenever `7|u`. Put

```text
s_j=1+2ju,                    j=0,...,6.               (24)
```

These are seven distinct positive odd speeds, and none is divisible by `u`.
Their strict masks partition `Y_u` into seven equal pieces, each of
`mu_u`-mass `1/7`.

### Proof of sharpness

For `x=x_(n,sigma)`,

```text
s_j x=x+j sigma/7 mod 1.                              (25)
```

Thus the seven phases form a full `1/7` grid. Boundary equality would imply
`s_j(14n+sigma)=+-u mod 14u`. Modulo `7`, its right side is zero, while its
left side is `sigma` because `s_j=1 mod 7`. Hence no endpoint occurs, exactly
one `j` is strictly dangerous at each support wall, and the seven masks are
disjoint and exhaustive. **QED.**

## 5. Distinct-valuation intersections and equality rigidity

The filtration gives an exact intersection law that does not require the
quotients in `(3)` to be coprime.

> **Distinct-valuation intersection law.** Let `s_1,...,s_r` be positive odd
> speeds whose valuations `e_i=nu_7(s_i)` are pairwise distinct and satisfy
> `e_i<a`. Then
>
> ```text
> |Y_u^sigma intersection intersection_i M_u(s_i)|=7^(a-r)
>                                                    for sigma=+-1,
> mu_u(intersection_(i=1)^r M_u(s_i))=7^(-r).        (25a)
> ```

### Proof

Order the valuations increasingly. Descend the first `e_1` constant levels
using `(20)`. At the next seven-child fibre, the reduced `s_1` is a
`7`-unit, so `(18)--(19)` put it on exactly one child. Every remaining speed
is still divisible by `7`, hence its mask is constant on that fibre and is
the pullback of its divided speed on the parent. Conditional on any parent
belonging to all remaining masks, exactly one of its seven equally weighted
children belongs to the whole intersection. This contributes a factor
`1/7`. Divide the remaining speeds by `7` and repeat. The residual
valuations stay pairwise distinct and below the residual support height, so
induction gives `r` factors of `1/7`. **QED.**

This pins down equality in the sharp cover theorem. Suppose seven lower-
valuation odd masks cover `Y_u`. Each has mass `1/7`, so equality holds in
the union bound and the seven restricted masks are pairwise disjoint. If two
speeds had different valuations, `(25a)` at fold two would give their
intersection mass `1/49`, a contradiction. Therefore

```text
seven lower masks cover Y_u  =>  all seven have one valuation e<a. (25b)
```

In particular, consider a primitive anchored `2+12` row at the critical
scale `u=7^a`, `a=nu_7(h)`, with exactly seven lower-valuation tails and five
upper-cone tails. If the seven lower masks cover `X_(7^a)`, then `(25b)`
gives one common value `e<a`. The anchor and the five other tails have
valuation at least `a`, so every row speed is divisible by `7^e`. Primitivity
forces

```text
e=0.                                                        (25c)
```

The assertion is deliberately about an exact seven-mask equality case. With
eight or more lower masks, union-bound equality and hence the disjointness
step need not hold.

## 6. Common-height owner blocks

The remaining equality case has an exact finite description. Suppose a
lower speed has the common form

```text
s=7^e r,       7 does not divide r,       U=7^(a-e). (25d)
```

Its mask on `X_(7^a)` depends only on `n mod U`; each owner address below has
`7^e` lifts on each sign sheet. Let

```text
E={-5,-3,-1,1,3,5},
delta_r=the unique element of E congruent to r mod 14,
J_r={ell in Z:-U<delta_r+14ell<U}.                   (25e)
```

Then `J_r` consists of exactly `U/7` consecutive integers. On the positive
sign sheet the owner block is the cyclic arithmetic progression

```text
A_r^+=r^(-1)(J_r+(delta_r-r)/14) mod U,              (25f)
A_r^-=-A_r^+ mod U.                                  (25g)
```

Indeed, apply `(4)--(5)` with `q_s=U`, `S_s=r`, and write the admissible
centered residue as `y=delta_r+14ell`. Strict endpoint equality cannot occur:
it would make the `7`-unit `delta_r` equal to `+-U` modulo `7`. Formula
`(25f)` follows, and reflection `(n,-1)<->(-n,+1)` gives `(25g)`. Thus seven
common-height masks partition `X_(7^a)` exactly when their seven positive
blocks partition `Z/U`; the negative-sheet partition is then automatic.

### Complete first-shell classification

When `U=7`, one has `J_r={0}`. Define

```text
kappa(r)=r^(-1)(delta_r-r)/14 mod 7.                 (25h)
```

The positive block is the singleton `{kappa(r)}`. Consequently seven masks
cover if and only if their `kappa` values exhaust `Z/7`. The complete fibres
among the `42` odd `7`-units modulo `98` are

```text
kappa   r mod 98
  0     +-1, +-3, +-5
  1     +-13, +-33, +-39
  2     +-17, +-27, +-37
  3     +-9, +-25, +-41
  4     +-19, +-31, +-43
  5     +-11, +-29, +-47
  6     +-15, +-23, +-45.                            (25i)
```

For example, the smallest positive choice

```text
L_0=(1,9,11,13,15,17,19)                            (25j)
```

has positive owners `(0,3,5,1,6,2,4)` and partitions both sign sheets.

> **FINITE-EXACT / OPEN boundary.** Exhaustive exact enumeration finds `147`
> distinct owner blocks and `21` exact seven-block partitions at `U=49`, and
> `1,029` blocks and `147` partitions at `U=343`. Every enumerated partition
> is a parallel lift: after choosing signs `epsilon_i in {+-1}`, there is a
> unit `R mod 2U` such that
>
> ```text
> epsilon_i r_i=R mod 2U,
> {epsilon_i r_i mod 14U}={R+2Uj:j in Z/7}.          (25k)
> ```
>
> At `U=2401`, a targeted search among all `7,203` blocks finds `46` blocks
> disjoint from the speed-`1` block and the unique completion
>
> ```text
> (1,4801,4803,9603,9605,14405,14407).               (25l)
> ```
>
> These are finite-exact data, not a proved all-depth rigidity theorem.
> Proving that every cyclic-AP partition `(25f)` has form `(25k)` remains
> open; pairwise disjointness alone does not supply it.

## 7. First-shell coarse-defect escape

There is a physical-scale refinement when `a=1` and exactly seven lower
tails remain. It is the `d=7` specialization of THM-4335's owner-permutation
characterization, written here on the critical wall fibre. For seven
positive odd `7`-units `L`, put, on the coarse circle `R/Z`,

```text
O_s(y)={j in Z/7:||s(y+j)/7||<1/14},
P_L={y:the seven singleton owner labels O_s(y), s in L,
       form a permutation of Z/7},
E_L=(R/Z)\P_L.                                       (25m)
```

Each `O_s(y)` has at most one element; away from finitely many equality
points it has exactly one. Thus `E_L` is exactly the locus on which at least
one of the seven lifts `(y+j)/7` is not killed by the lower tails.

Let `h=7H`, with `7` not dividing `H`, and let the remaining `m<=5` odd tails
be `7q_1,...,7q_m`. Then the full row is strictly bad if and only if

```text
E_L subset D_(2H) union D_(q_1) union ... union D_(q_m),
D_c={y:||cy||<1/14}.                                 (25n)
```

For the forward direction, take `y in E_L` outside the right side and choose
an unowned lift `t=(y+j)/7`. Every lower tail is safe there, while

```text
14H t=2Hy mod 1,             7q_i t=q_i y mod 1.     (25o)
```

Hence the anchor and every upper tail are also safe, contradicting strict
failure. Conversely, on `P_L` the lower tails kill all seven lifts, while on
`E_L` the inclusion says that one coarse speed kills the entire fibre. This
proves the equivalence. Since every `D_c` has Lebesgue measure `1/7`, strict
failure requires

```text
lambda(E_L)<=(m+1)/7,       lambda(P_L)>=(6-m)/7.    (25p)
```

In particular, for every `m<=5`, the uniform sufficient escape condition is
`lambda(E_L)>6/7`, equivalently `lambda(P_L)<1/7`.

For the wall tiler `L_0` in `(25j)`, exact rational endpoint sweep gives

```text
P_(L_0)=(1/18,3/38) union (35/38,17/18),
lambda(P_(L_0))=8/171,
lambda(E_(L_0))=163/171>6/7.                         (25q)
```

Off those equality points, the owner is equivalently

```text
kappa_s(y)=-s^(-1) floor(sy+1/2) mod 7.             (25r)
```

For completeness, each owner changes only at a point
`(2m+1)/(2s)`. Sorting the common endpoints gives `78` open cells. The only
permutation cells are the two intervals in `(25q)`, with owner tuples in the
order `(1,9,11,13,15,17,19)` equal respectively to

```text
(0,3,5,1,6,2,4),       (6,3,1,5,0,4,2).             (25s)
```

Thus `(25p)` proves a safe physical lift for every choice of `H` and the at
most five upper tails attached to this literal lower family.

The same family is a sharp hostile to the information retained by the wall
mass and the THM-4372 critical transport alone. Take

```text
h=7,   anchor=14,
L=L_0,
upper tails=(7,21,35,63,77).                         (25t)
```

At every point of `X_7`, exactly one lower mask is dangerous, while the
anchor and all five upper masks are safe. In THM-4372's notation `a=1`,
`n_a=5`, and `Z=0`, so its retained transport right side is identically
zero. There is no forced escape from those two compressed invariants.
Nevertheless `(25q)` closes the literal row, and an explicit witness is

```text
t=15/182,       min ||st||=1/14,       unique binder s=13. (25u)
```

Here the minimum ranges over all thirteen displayed speeds. This is not a
counterexample: it isolates the lower owner-defect/physical-scale coordinate
discarded by the mass and deeper-current summaries.

The closure in `(25q)` is literal-speed, not merely residue-class data. The
parallel wall tiler

```text
(1,15,29,43,57,71,85)                               (25v)
```

has the same exact one-owner-per-class wall property, but its physical
permutation locus has measure `7122/46835>1/7`. Thus the scalar defect-
measure test does not close every first-shell tiler. Adding `98` to one speed
preserves its `X_7` wall mask but can change the coarse owner evolution,
exactly as the physical-scale warning in Section 11 requires.

## 8. Direct anchor-valuation and six-of-twelve theorem

Let `W` be any finite set of positive odd speeds. Given the positive anchor
parameter `h`, put

```text
a=nu_7(h),                    u=7^a.                  (26)
```

The auxiliary `u` need not occur in `W`. Here `v=1`, so `(13)` gives
`Y_u=X_u`, with exactly `2*7^a` equally weighted physical walls. By `(21)`,
the anchor `2h` has mass zero. By `(17)`, every tail at valuation at least
`a` has mass zero and every lower-valuation tail has mass exactly `1/7`.

> **Direct anchor-valuation theorem.** If
>
> ```text
> #{s in W:nu_7(s)<nu_7(h)}<=6,                       (27)
> ```
>
> then some `x in X_(7^nu_7(h))` satisfies
>
> ```text
> ||2hx||>=1/14,              ||sx||>=1/14 for every s in W. (28)
> ```

Indeed, the strict-danger union has `mu_u`-mass at most `6/7<1`. The point
outside it is the claimed genuine physical wall.

For twelve distinct positive odd tails, strict failure therefore forces

```text
#{s in W:nu_7(s)<nu_7(h)}>=7,
#{s in W:nu_7(s)>=nu_7(h)}<=5.                        (29)
```

Thus six of the twelve tails in the upper valuation cone force safety. No
tail at exact anchor height is required. In particular, if `7` does not
divide `h`, the lower set in `(27)` is empty and the fixed walls `+-1/14`
prove safety for every odd-tail bank.

## 9. Tail re-anchor and divisibility specialization

The more general quadrature permits re-anchoring at an odd row speed. Suppose
`W` contains `r` with

```text
nu_7(r)=nu_7(h).                                      (30)
```

If at most six members of `W\{r}` have valuation below this common value,
then `{2h} union W` is safe at a point of `Y_r`, by `(17)`, `(21)`, and the
union bound.

There is a purely divisibility-based specialization. If at most six members
of `W\{r}` are not divisible by `r`, then the row is safe. Every divisible
speed is `cr` with integer `c`; since both speeds are odd, `c` is exactly odd,
and at every `r`-wall

```text
cr(14n+sigma)/(14r)=cn+c sigma/14                    (31)
```

has circle distance at least `1/14`. The anchor is safe on `Y_r` by `(30)`,
and the at most six nonmultiples have total danger mass at most `6/7<1`.

Consequently, in a strictly bad anchored `2+12` row, every tail satisfying
`(30)` has at least seven nonmultiples among the other eleven. In particular,
strict failure cannot have six of its twelve tails divisible by one such `r`,
counting `r` itself. This remains valid even if a previously chosen base
speed `u_0` lies in THM-4348's formerly untreated resonance `7u_0|h`: use
the exact physical walls `Y_r subset X_r` instead.

## 10. Exact physical controls

The headline theorem genuinely needs no exact-height tail. Take

```text
h=420,
W_0={1,3,5,9,11,13,49,147,245,343,441,539}.          (32)
```

Exactly six tails have valuation below `nu_7(h)=1`; the other six have
valuation strictly above it, and none has valuation exactly one. At the
auxiliary wall

```text
x=13/98 in X_7,
min_(s in {840} union W_0)||sx||=13/98>1/14.          (33)
```

For a physical re-anchor inside the old resonance, take

```text
h=420,
W={1,7,21,35,63,77,91,105,119,133,147,161}.          (34)
```

The original choice `u_0=1` has `7u_0|h`. Re-anchor at `r=7`, which satisfies
`nu_7(r)=nu_7(h)=1`. Among the other eleven tails, exactly ten are multiples
of `7` and exactly one, the tail `1`, is a nonmultiple. Each quotient of an
odd multiple by `7` is odd. The supported wall

```text
(n,sigma)=(1,-1),               x=13/98              (35)
```

satisfies

```text
min_(s in {840} union W)||sx||=1/14,                  (36)
```

with physical binding speeds exactly `7,91,105`. Both controls are primitive
thirteen-speed relative rows because they contain speed `1`; they are safe
controls, not LRC(14) counterexamples.

## 11. Sharp and hostile boundaries; physical scale retained

Three exact controls delimit the theorem.

1. The seven speeds `(24)` strictly cover `Y_u` for every `7|u`, so the
   six-mask conclusion cannot be extended to seven for this measure.
2. For `u=5,h=7`, the even anchor `14` and the four positive odd residuals
   `(71,81,83,87)` cover all ten points of `X_5`. This is the `b>a` branch in
   which `(21)` correctly gives no certificate.
3. For `u=11,h=420`, the anchor `840` and the five positive odd residuals
   `(159,331,523,665,835)` cover all twenty-two points of `X_11`. Including
   the base speed `11`, the displayed speed set is primitive. Thus there is
   no unconditional anchor-plus-five wall escape at an arbitrarily chosen
   row-speed scale without the valuation hypothesis.

For every integer `L` for which the speed remains positive,

```text
M_u(s+14uL)=M_u(s).                                   (37)
```

This permits independent lifting of a hostile wall mask to larger physical
speeds, as in items 2--3. It does **not** preserve the physical danger-tooth
diameter `1/(7s)`, tooth address, farthest reach, partial renewal trace, or
first exit. Accordingly, THM-4363/4365/4367 prevent upgrading `(37)` to a
renewal-chain equivalence. THM-4370 proves exact mask calculus and genuine
pointwise wall witnesses only; it does not erase physical scale.

## 12. Exact audit scope

The primary script uses integer centered residues with convention

```text
-7u<e<=7u,
```

and strict danger `|e|<u`. It checks:

- `35,000` individual masks against `(4)--(5)`, including `26,424` actual
  endpoint equalities;
- `4,065` intersection instances and `1,500` union instances through fold
  six against direct enumeration;
- `107,143` odd-speed cases of `(17)`, `192,000` anchor cases of `(21)`,
  including `10,800` rechecks through an independent Fraction-circle path,
  and every seven-child support projection in the declared odd-`u<400`
  banks;
- `57` instances of the sharp partition `(24)`, each also checked pointwise
  through the independent Fraction-circle path;
- `3,104` distinct-valuation intersections through fold four on composite
  support scales, with exact counts on each sign sheet and `152` independent
  Fraction-circle rechecks;
- all `42` first-shell unit residues, `2,450` formula/direct owner-block
  comparisons, exhaustive exact partitions at `U=49,343`, and the targeted
  speed-one completion at `U=2401`;
- both exact permutation measures in Section 7, the retained-invariant row
  and its physical witness, and `26,440` pointwise coarse-fibre equivalence
  checks;
- the no-exact-height witness `(32)--(33)`, the physical re-anchor witness
  `(34)--(36)`, and both hostile full-wall covers.

Normal Python, optimized Python, and a nondefault hash-seeded Python replay
produce the frozen output byte-for-byte.
