---
id: THM-3416
title: "Zero-mode-cochain global rank-six support"
status: >
  PROVISIONAL PROOF CANDIDATE / AUDIT REQUIRED + VERIFIED-EXACT companion.
  The proposed all-q theorem says that the global zero-mode-cochain rank is
  six exactly when one of 11,15,23,25 divides q and none of 8,9,10,12 does.
  The converse is analytic: two exact complement classifiers remove orders
  three and five, and a reflection-orbit/CRT argument removes every mixture
  of orders 17 and 29 with arbitrary low-density companions.  The independent
  primitive cap-six census is FINITE-EXACT only through Q=200.  No arbitrary
  physical-time or LRC(14) conclusion is claimed.
source: q8-multiblocker-2026-08-15
audit: pending independent proof/type/replay audit
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3414-fixed-zero-six-owner-base-classification
  - THM-3415-zero-mode-cochain-global-rank-five-support
related:
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor
  - THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs
script: 04-computation/lrc_zero_mode_cochain_rank6_support_thm3416.py
output: 05-knowledge/results/lrc_zero_mode_cochain_rank6_support_thm3416.out
script_sha256: c946cf9a66daa13da790cc9b1129993f9b5c1a8a2bdc6dec1bcc07b644989122
output_sha256: 733a3ed02910348b87b561df7f1c79eef2dc431a8d702f2c0f807db33c1298cc
semantic_sha256: 99892baf39b3d2b1b6a802bf21e0fe4164f155030d8ad051bc7ae26513b01ca3
hash_basis: LF-normalized bytes
---

# THM-3416 -- zero-mode-cochain global rank-six support

**PROVISIONAL PROOF CANDIDATE / AUDIT REQUIRED + VERIFIED-EXACT companion.**

## 1. Exact statement

Retain THM-3415's notation.  For `q>=2`, `rho_ZMC(q)` is the least number of
distinct positive transverse owners whose strict danger blocks cover all
`q` cyclic sheets at one common THM-3398 mode centre and whose complete mode
cochain vanishes.  Put `rho_ZMC(q)=infinity` if there is no such family.

Define

```text
L={8,9,10,12},       N={11,15,23,25},
S=L union N={8,9,10,11,12,15,23,25}.                 (1)
```

Then the proposed theorem is

```text
rho_ZMC(q)<=6  iff  s|q for some s in S,              (2)

rho_ZMC(q)=6   iff  n|q for some n in N
                       and l does not divide q for every l in L.   (3)
```

The lower layers from THM-3415 are

```text
rho_ZMC(q)=4 iff 8|q or 9|q,
rho_ZMC(q)=5 iff (10|q or 12|q) and 8 does not divide q
                                      and 9 does not divide q.     (4)
```

Thus `(3)` is an exact rank statement, not merely a six-owner construction.
It concerns the zero-complete-cochain/common-mode-centre category.  It is not
a classification of arbitrary common physical times or positive cochains,
and it gives no LRC(14) ledger decrement.

## 2. Inheritance and connection contract

THM-3405 reduces every zero-cochain family, after a common sheet relabelling,
to a fixed-zero layer or a Boolean half-twist layer.  THM-3415 then supplies
the exact primitive divisor minimum.  For an active speed family,

```text
U=dV,  gcd(V)=1,  g=gcd(q,d),  Q=q/g,                 (5)
```

the cover is the `g`-fold pullback of a primitive `Q`-sheet cover at
`epsilon=0` or `1`.  Its normalized residues satisfy

```text
gcd(M_epsilon,r_1,...,r_k)=1,
M_0=Q, M_1=2Q.                                        (6)
```

Only the necessity of this reduction is used in the converse.

| field | exact connection |
|---|---|
| source | primitive THM-3405 Boolean half-twist covers by at most six blocks |
| target | quotient-order cylinders on odd reflection orbits |
| map | send residue `r` to `m=Q/gcd(Q,r)`, then reduce its mask to `Z/mZ` |
| preserved | strict danger predicate, block size, cyclic pullback, reflection, and cover |
| destroyed | the location of individual points inside a reflection pair and the residue lift |
| required sidecar | the common fixed point of `ell -> -1-ell` and CRT independence of the 17/29 cylinders |
| cheapest decisive tests | the sharp low block at `Q=43` and the mixed hostile `Q=17*29=493` |

The proof board is therefore: exact block density; order-three complement;
order-five complement; reflection pairs; CRT cylinder product; positive
atoms.  The complement classifiers remove the two tempting small orders,
while the reflection sidecar closes the mixed large-order case that density
alone leaves open.

## 3. Exact half-twist block counts

On the primitive half layer put

```text
H_(Q,r)={ell in Z/QZ:
         ||r(2ell+1)/(2Q)||<1/14},
m=Q/gcd(Q,r).                                         (7)
```

The mask in `(7)` is the pullback of a mask on `Z/mZ`.  Let `z(m)` be the
zero-like count and `a(m)` the odd-coset count.  Direct reduction modulo
`2m` gives

```text
z(m)=1+2 floor((m-1)/14),
a(m)=2 ceil(floor((m-1)/7)/2),

h(m)=a(m)                     if m is even,
h(m)=max(a(m),z(m))           if m is odd.             (8)
```

Indeed, when the reduced numerator is odd it permutes the odd residue coset
modulo `2m`; an even reduced numerator is possible only for odd `m` and gives
the zero-like coset.  Strict endpoints are already built into `(8)`.

Both counts satisfy

```text
h(m)<=(m+6)/7.                                        (9)
```

For `z`, multiply by seven and use the defining floor.  For `a`, write
`k=floor((m-1)/7)`: if `k` is even then `a=k`, and if `k` is odd then
`a=k+1`; both give `7a<=m+6`.

If six blocks cover, at least one has density at least `1/6`.  The tail
`m>=37` is excluded by `(9)`, and direct evaluation for `m<=36` gives the
complete all-order list

```text
6h(m)>=m iff
m in H6={3,5,8,9,10,11,12,15,17,22,23,24,29,36}.     (10)
```

Call `Q` **S-free** when no member of `S` divides it.  Since every selected
order divides `Q`, `(10)` says that an S-free six-cover contains an owner of
order

```text
3, 5, 17, or 29.                                      (11)
```

The next two sections show that none is possible.

## 4. Exact complement classifiers eliminate orders three and five

### 4.1 An order-three anchor

Every nonempty order-three half block is the same pullback singleton.  Its
complement has density `2/3`.  If five remaining blocks cover that
complement, one of them contributes at least `2/15` of all sheets.

For a companion of order `m`, reduce the intersection exactly to
`lcm(3,m)`.  If `gcd(3,m)=1`, work on `3m` sheets.  Each covered quotient
point has three lifts and exactly one belongs to the order-three class, so
the maximum complement contribution is `2h(m)`.  The threshold is therefore

```text
5h(m)>=m.                                             (12)
```

By `(9)`, `(12)` forces `m<=15`; direct evaluation of the coprime cases
leaves `m in {5,8,10}`.

If `3|m`, multiplication by a unit makes the removed order-three class one
colour of the centred interval.  Hence the complement contribution `c`
satisfies

```text
c<=ceil(2h(m)/3).                                     (13)
```

The necessary inequality `c/m>=2/15`, together with `(9)` and
`ceil(x)<=x+2/3`, forces `m<=32`.  Exact evaluation of the multiples of
three through 32 leaves `m in {9,12,15}`.  Thus the all-order classifier is

```text
an order-three anchor forces a companion order in
{5,8,9,10,12,15}.                                     (14)
```

For reference, the qualifying rows `(m,lcm(3,m),complement size,maximum
contribution)` are

```text
(5,15,10,2), (8,24,16,4), (9,9,6,2),
(10,30,20,4), (12,12,8,2), (15,15,10,2).              (15)
```

### 4.2 An order-five anchor

An order-five half block is likewise one pullback singleton; its complement
has density `4/5`.  One of five companions must contribute at least `4/25`.
For `gcd(5,m)=1`, reduction to `5m` gives complement contribution `4h(m)`,
and again the condition is `(12)`.  The coprime survivors are `{3,8,9}`.

For `5|m`, the balanced-colour bound is

```text
c<=ceil(4h(m)/5).                                     (16)
```

Using `ceil(x)<=x+4/5`, the threshold `c/m>=4/25` and `(9)` again force
`m<=32`.  The divisible survivors are `10,25`.  Therefore

```text
an order-five anchor forces a companion order in
{3,8,9,10,25}.                                        (17)
```

The qualifying rows are

```text
(3,15,12,4), (8,40,32,8), (9,45,36,8),
(10,10,8,2), (25,25,20,4).                            (18)
```

Now suppose that `Q` is S-free.  If order three occurs, `(14)` leaves only
order five; then `15=lcm(3,5)|Q`, contradiction.  Hence order three is
absent.  If order five occurs, `(17)` leaves either order three or an order
forcing a member of `S` to divide `Q`; both are impossible.  Thus neither
order in `{3,5}` can occur in an S-free six-cover.

## 5. Reflection and CRT eliminate every 17/29 mixture

By `(10)--(11)` an S-free cover must now contain an order-17 or order-29
owner.  Every remaining non-17/non-29 selected order has block density at most

```text
7/43.                                                 (19)
```

For `m>=43`, this follows from `(9)` because
`(m+6)/(7m)<=7/43`.  Direct evaluation below 43, after removing the orders
already excluded, has maximum `6/37<7/43`.  The constant is sharp at the
hostile order `m=43`, where `h(43)=7`; nevertheless six such blocks have
total capacity only `42/43`.

The half-twist predicate is invariant under

```text
sigma(ell)=-1-ell.                                     (20)
```

On `Z/17Z` and `Z/29Z`, `sigma` has one common fixed point and all other
orbits are pairs.  An order-17 block has size at most three, so `a` such
blocks cover at most `1+2a` quotient points.  An order-29 block has size at
most five, so `b` such blocks cover at most `1+4b` quotient points.  Hence
the fractions missed by all high blocks are at least

```text
A_a=(16-2a)/17 if a>0, and A_0=1,
B_b=(28-4b)/29 if b>0, and B_0=1.                     (21)
```

These are cylinder bounds, not an assumption that the masks coincide.  If
both orders occur, then `493|Q`, and the equal-fibre map

```text
Z/QZ -> Z/17Z x Z/29Z                                (22)
```

makes the set missed by every high block have density at least `A_a B_b`.
This exact product is the sidecar lost by a scalar density estimate.

Pad a cover by empty slots to six owners and put `c=6-a-b`.  The remaining
blocks cover at most `7c/43`.  The three possible cases give

```text
b=0, a>=1:
  A_a-7(6-a)/43=(33a-26)/731 >= 7/731;

a=0, b>=1:
  B_b-7(6-b)/43=(31b-14)/1247 >= 17/1247;

a,b>=1:
  A_a B_b-7(6-a-b)/43
   =(344ab+1043a+699b-1442)/21199
   >=644/21199.                                       (23)
```

Every right-hand side is positive.  Thus the low blocks cannot fill the
reflection/CRT core missed by the high blocks.  This contradiction proves:

```text
every primitive half-twist cover by at most six owners has Q divisible
by a member of S.                                      (24)
```

The mixed hostile `Q=493` is load-bearing: pairwise density alone does not
remember that order-17 and order-29 complements occupy independent CRT
coordinates, whereas `(22)` does.

## 6. Fixed zero and the all-q converse

THM-3414 proves that a fixed-zero cover by at most six owners exists exactly
when its sheet number is divisible by one of

```text
15,16,18,20,24.                                       (25)
```

Each number in `(25)` has a divisor in `S`: respectively `15,8,9,10,12`.
Therefore every primitive zero-layer cover also forces an `S` divisor.
Together with `(24)` and the primitive divisor reduction `(5)--(6)`, every
global cover by at most six owners forces `s|q` for some `s in S`.

Conversely, the four new atoms have the following exact covers.  Here
`epsilon=0` denotes fixed zero and `epsilon=1` the half twist.

```text
Q=11, epsilon=1: (1,2,3,5,7,9),
  block sizes (2,1,2,2,2,2), a partition;

Q=15, epsilon=0: (1,2,3,4,5,7),
  multiplicity one on 14 sheets and six at sheet 0;

Q=15, epsilon=1: (1,4,6,7,8,10),
  multiplicity one on 14 sheets and four at sheet 7;

Q=23, epsilon=1: (1,4,5,7,9,11),
  block sizes (4,3,4,4,4,4), a partition;

Q=25, epsilon=1: (1,9,10,11,19,21),
  block sizes (4,4,5,4,4,4), a partition.             (26)
```

The lower atoms `8,9,10,12` are inherited from THM-3415.  THM-3405 dilation

```text
(Q,u,c) -> (lambda Q,lambda u,c/lambda)               (27)
```

pulls every atom back to each positive multiple, preserving the cover and
zero cochain.  This proves `(2)`.  Subtracting the exact lower support `(4)`
gives `(3)`.

## 7. Periodic and arithmetic corollaries

The exact rank-six predicate is periodic modulo

```text
lcm(S)=455400.                                        (28)
```

Inclusion-exclusion, independently matched by a direct residue census, gives
exactly `55000` accepting residues.  Consequently

```text
#{q<=X:rho_ZMC(q)=6}=(25/207)X+O(1),
sum_(q<=X,rho_ZMC(q)=6) 1/q=(25/207)log X+O(1).        (29)
```

The rank-four, rank-five, and rank-six layers are disjoint and have densities
`2/9`, `4/45`, and `25/207`.  Hence

```text
#{q<=X:rho_ZMC(q)<=6}=(149/345)X+O(1),
sum_(q<=X,rho_ZMC(q)<=6) 1/q=(149/345)log X+O(1).      (30)
```

Two exact transport sidecars follow.

For Fibonacci numbers `F_n` with `n>=3`, the ranks of apparition at the bases
in `(1)` are

```text
z(8,9,10,12,11,15,23,25)=(6,12,15,12,10,20,24,25).   (31)
```

Using `gcd(F_m,F_n)=F_gcd(m,n)`, divisibility by each base occurs exactly on
multiples of its rank.  Thus

```text
rho_ZMC(F_n)=6 iff
(10|n or 24|n or 25|n) and 6 does not divide n
                                  and 15 does not divide n.        (32)
```

This indicator has period 600, with 48 accepting residues and density
`2/25`.

On the Berggren `U`-spine `q(n)=4n^2+12n+11` for `n>=0`, the polynomial is
odd and has no root modulo five; its discriminant is a nonsquare modulo 23.
Hence among the new bases only 11 can divide it, exactly for `n=0,8 mod 11`.
Among the lower bases only 9 can divide it, exactly for `n=1,5 mod 9`.
Therefore

```text
rho_ZMC(4n^2+12n+11)=6 iff
n=0 or 8 (mod 11), and n is not 1 or 5 (mod 9).        (33)
```

There are 14 accepting classes modulo 99.  Rank five is absent on this
spine.

## 8. Exact companion and evidence boundary

The standard-library companion freezes:

1. the exact formula `(8)`, bound `(9)`, and list `(10)`;
2. both all-order complement cutoffs, including every boundary tie;
3. all order-17/order-29 mask sizes and exact union maxima through six;
4. the `Q=493` cylinder product and every rational deficit in `(23)`;
5. all five displayed positive controls in `(26)`;
6. an independent augmented-prime-breaker primitive cap-six census for both
   twists through `Q=200`; and
7. the period, Fibonacci, and Berggren arithmetic in `(28)--(33)`.

The finite primitive census is explicitly **FINITE-EXACT only through
`Q=200`**.  It is a hostile audit of the analytic proof, not an extrapolation
and not a dependency of `(2)--(3)`.  No primitive pattern beyond that cutoff
is claimed.  The theorem has no tournament: the decisive observable is a
reflection-orbit cylinder with a CRT sidecar, not a pairwise orientation.
