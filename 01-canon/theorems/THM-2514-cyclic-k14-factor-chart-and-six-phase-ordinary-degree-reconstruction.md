---
id: THM-2514
title: "Cyclic K14 factor chart and six-phase ordinary-degree reconstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  At the critical
  dimensions 13 x 7, the ordered representative cut is a transversal of
  the nonzero antipodal classes of F_13.  Resolving its thirteen diagonal
  pairs through one marked vertex Omega identifies every affine-cut chart
  bijectively with E(K_14).  The thirteen horizontal rows are the cyclic
  one-factorization; the seven cut columns are one hub star and six
  Hamilton 13-cycles.  For a rational row-zero defect, take ordinary
  weighted degrees on K_14 in the seven translated cut charts.  Their sum
  is zero, and any six of the seven degree vectors reconstruct all 78
  defect coordinates exactly.  Fourier reconstruction has multiplier
  K_(lambda,beta)+K_(-lambda,beta)-1, nonzero precisely on the six
  nontrivial cut characters; the marked vertex recovers the alpha=0
  septimal modes.  The full affine CRT gauge acts by vertex relabelling
  and cut permutation.  One antipodal pair of Radon marginals still has
  rank only 24 and kernel dimension 54.  Septimal translation exchanges
  the unique star column with a cycle column, which no K_14 vertex
  automorphism can realize, so the cut phase remains a necessary sidecar.
  This is a signed static root chart, not runner-labelled Boolean
  owner/time/deep ancestry, and it excludes no live row or LRC(14) case.
source: codex-2026-07-27-k14-degree-reconstruction
depends_on:
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2509-antipodal-radon-cospan-and-lossless-septimal-chart
related:
  - THM-128-z13-circulant-maximizer-structure
  - THM-913-parallel-class-book-drawing
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2510-truncated-radon-tight-frame-and-integral-dipole-boundary
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
script: 04-computation/lrc14_k14_degree_reconstruction_thm2514_referee.cpp
output: 05-knowledge/results/lrc14_k14_degree_reconstruction_thm2514_referee.out
script_sha256: 034b5c300575155303ee4ca18a5e948c6cb08d70fb1367d3ea0d1aa4fe942056
output_sha256: 2197609b9877239180f0e06924654f0ca9c49bdf48d7a4210c4e49ae4f0b6d10
hash_basis: working-tree bytes (LF)
---

# THM-2514 -- the septimal strip is a cyclic `K_14` factor chart

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2509 keeps the joint two-leg Radon chart because its two marginals lose
`54` dimensions.  At `13 x 7`, the joint chart has a more intrinsic
description.  Its `91` cells are exactly the `91` edges of `K_14`; the
septimal rows label one edge in each round of the standard cyclic
one-factorization.

This is not only a relabelling theorem.  If the diagonal pairs are resolved
through the fourteenth vertex and one takes ordinary weighted degrees, then

```text
six translated cut phases x thirteen independent vertex degrees
  = 78 defect coordinates.                                      (1)
```

The equality in (1) is an exact reconstruction theorem.  It exposes a
converse-invariant root-current carrier while also proving why the cut phase
cannot be quotiented away.

## 1. The critical antipodal chart

Write

```text
I={0,1,...,6} subset F_13,
rho:F_7 -> I                                                   (2)
```

for the ordered representative section.  Fix

```text
tau in F_13^*,       a in F_7^*,       c in F_7,              (3)
```

and put `s=ar+c`.  The ordered cospan map is

```text
Pi_(tau,a,c)(h,r)
 =(h+tau rho(s),h-tau rho(s)).                                (4)
```

For the baseline cut `(a,c)=(1,0)`, the ordered image is

```text
{(x,y):(x-y)/(2tau) lies in I}.                               (5)
```

Its inverse on that image is

```text
h=(x+y)/2,
r=rho^(-1)((x-y)/(2tau)).                                     (6)
```

Thus every ordered fibre has size one on (5) and zero off it.  More is true:

```text
I intersection (-I)={0},
I union (-I)=F_13.                                            (7)
```

Consequently (5) consists of the diagonal together with exactly one
orientation of every off-diagonal unordered pair.  After quotienting by leg
swap, (4) is a bijection

```text
F_13 x F_7  ->  Sym^2(F_13).                                  (8)
```

At `tau=0`, by contrast, all seven cells `(h,r)` collapse to `(h,h)` and
their signed pushforward is zero by the row-sum law.  Nonzero slope is a
sharp hypothesis.

## 2. Resolving the diagonal gives `E(K_14)`

Add one marked vertex `Omega`.  Replace a diagonal pair by a spoke and keep
every off-diagonal pair as an ordinary edge:

```text
Phi_(tau,a,c)(h,r)
 = {Omega,h},                              if ar+c=0 in F_7,
 = {h+tau rho(ar+c),h-tau rho(ar+c)},      otherwise.          (9)
```

Equation (8) makes (9) a bijection

```text
F_13 x F_7  ->  E(K_14).                                      (10)
```

The inverse is explicit.  A spoke `{Omega,h}` has `s=0`.  For a finite edge
`{x,y}`, its midpoint is `h=(x+y)/2`; its half-difference determines one
antipodal class, and (7) gives the unique `s in {1,...,6}` representing that
class after division by `tau`.  Finally

```text
r=a^(-1)(s-c).                                                 (11)
```

For fixed `h`, the seven edges form

```text
M_h={{Omega,h}}
    union {{h-u,h+u}:u in F_13^*/{+-1}}.                       (12)
```

This is a perfect matching of the fourteen vertices.  Every finite edge has
one midpoint, so the thirteen matchings `M_h` partition `E(K_14)`.  Their
finite parts are exactly THM-913's parallel classes

```text
x+y=2h                       mod 13.                            (13)
```

The orthogonal column view is equally rigid.  For one fixed cut value `s`:

```text
s=0:       the thirteen edges form the hub star K_(1,13);
s!=0:      the thirteen edges form the Hamilton cycle
           h -> {h-tau s,h+tau s}.                             (14)
```

The cycle is connected because its step `2tau s` is nonzero in the prime
field.  Hence

```text
E(K_14)=one 13-edge star disjoint-union six Hamilton C_13 cycles. (15)
```

The ordered finite arcs in (4) have difference

```text
y-x in -2tau{1,...,6}
       =tau{1,3,5,7,9,11}.                                    (16)
```

Thus they are precisely the multiplier orbit of THM-128's odd-step
`C_13` circulant tournament.  Replacing `tau` by `-tau` takes the global
converse while leaving the unoriented `K_14` edge chart unchanged.

## 3. Row-zero defects are factor-balanced edge weights

Let

```text
d:F_13 x F_7 -> Q,
sum_r d(h,r)=0                         for every h.             (17)
```

Push `d` through (9), without summing fibres, to obtain an edge weighting

```text
W_(tau,a,c)^d:E(K_14)->Q.                                    (18)
```

Because (9) is bijective, (18) is linearly lossless.  Under the chart, the
row-zero law is exactly

```text
sum_(e in M_h) W_(tau,a,c)^d(e)=0          for every h.         (19)
```

Conversely, pulling back any edge weighting satisfying (19) gives a unique
row-zero array.  The thirteen equations have disjoint edge supports and are
independent, so this image has dimension

```text
91-13=78.                                                       (20)
```

## 4. Ordinary weighted degrees and the seven-phase relation

Let

```text
G_(tau,a,c)^d(v)=sum_(e incident to v)W_(tau,a,c)^d(e),
                         v in F_13 union {Omega}.               (21)
```

These are ordinary undirected weighted degrees: every spoke and every finite
edge contributes once at each of its two endpoints.  Put

```text
r_c=-a^(-1)c                    in F_7.                         (22)
```

For finite `x`, direct incidence counting gives

```text
G_(tau,a,c)^d(x)
 =R_(tau,a,c)^d(x)+R_(-tau,a,c)^d(x)-d(x,r_c),                 (23)

G_(tau,a,c)^d(Omega)=sum_h d(h,r_c).                           (24)
```

The correction in (23) is load-bearing.  The two Radon marginals count the
zero-cut diagonal twice at `x`; the resolved spoke `{Omega,x}` contributes
only once.

Now fix one cell `(h,r)` and vary `c`.  The value `s=ar+c` runs through all
of `F_7`.  Its seven chart edges meet every vertex of `K_14` exactly once:
the zero value supplies `{Omega,h}`, while the six nonzero values pair `h`
with all twelve other finite vertices.  Therefore, coefficientwise in `d`,

```text
sum_(c in F_7)G_(tau,a,c)^d(v)
 =sum_(h,r)d(h,r)=0                    for every v.             (25)
```

Each individual degree vector also has total vertex sum zero, because its
total is twice the total edge weight and (17) makes that total zero.  Hence
one phase has at most thirteen independent coordinates.

## 5. Fourier diagonalization of the degree bank

Fix primitive roots `zeta=zeta_13` and `xi=zeta_7`.  Use the unnormalized
defect transform

```text
dtilde(alpha,beta)
 =sum_(h,r)d(h,r)zeta^(-alpha h)xi^(-beta r).                  (26)
```

For finite vertices, put

```text
Ghat_(tau,a,c)(alpha)
 =sum_(x in F_13)G_(tau,a,c)(x)zeta^(-alpha x).                (27)
```

Retain THM-2508's geometric kernel

```text
K_(lambda,beta)=sum_(s=0)^6 zeta^(-lambda s)xi^(-beta s)       (28)
```

and define its converse-symmetrized, spoke-corrected form

```text
L_(lambda,beta)
 =K_(lambda,beta)+K_(-lambda,beta)-1.                          (29)
```

Fourier transform (23) in `c`.  The substitution `c=s-ar` gives, for
`alpha!=0`,

```text
sum_c Ghat_(tau,a,c)(alpha)xi^(-beta c)
 =L_(alpha tau,beta)dtilde(alpha,-beta a).                     (30)
```

At the marked vertex, (24) gives the missing horizontal-zero slice:

```text
sum_c G_(tau,a,c)(Omega)xi^(-beta c)
 =dtilde(0,-beta a).                                          (31)
```

### The multiplier never vanishes on a nontrivial cut character

Let `lambda,beta!=0` and set `z=zeta^(-lambda)`.  Equation (29) is

```text
L_(lambda,beta)
 =1+sum_(s=1)^6 xi^(-beta s)(z^s+z^(-s)).                     (32)
```

Replace `z^(-s)` by `z^(13-s)`.  The right side is a polynomial of degree at
most twelve over `Q(xi)`.  The coprime cyclotomic fields have intersection
`Q`, so `Phi_13` remains irreducible over `Q(xi)`.  If (32) vanished, that
polynomial would be a scalar multiple of `Phi_13`.  Its constant coefficient
is one, so every coefficient would have to equal one.  In particular

```text
xi^(-beta s)=1                 for every s=1,...,6,             (33)
```

forcing `beta=0`, a contradiction.  Thus

```text
L_(lambda,beta)!=0             when lambda,beta!=0.             (34)
```

At `beta=0`, the twelve nonzero powers in (32) occur once each, so

```text
L_(lambda,0)=Phi_13(z)=0.                                     (35)
```

This is the exact boundary: the degree bank sees all six nontrivial cut
characters, while the row-zero law already kills
`dtilde(alpha,0)` for every `alpha`.

## 6. Any six cut phases reconstruct the defect

Given the degree vectors for any six values of `c`, equation (25) recovers
the omitted seventh vector.  Fourier transform the complete phase list.

- Equations (30) and (34) recover all `72` values
  `dtilde(alpha,beta)` with `alpha,beta!=0`.
- Equation (31) recovers the remaining six values with
  `alpha=0,beta!=0`.
- The seven values with `beta=0` are zero by (17).

Fourier inversion now reconstructs all `78` independent coordinates of
`d`.  Therefore

```text
d -> (G_(tau,a,c)^d)_(c in C)             |C|=6               (36)
```

is injective for every fixed `tau,a!=0` and every six-element
`C subset F_7`.  Its exact rank is `78`.

This is dimension-minimal among banks of ordinary degree vectors: each
`G_c` lies in the thirteen-dimensional zero-total hyperplane of `Q^14`, so
fewer than six phases have rank at most `65`.  In fact the exact referee
finds rank `13` for every single phase and rank `78` for all `504`
six-phase banks obtained by choosing `tau`, `a`, and one omitted phase.

There is a useful intermediate form.  For the six slope classes
`tau modulo {+-1}`, the converse-invariant looped degrees

```text
R_tau+R_(-tau)                                                (37)
```

have joint rank `72` and exact kernel

```text
d(h,r)=b(r),                 sum_r b(r)=0.                     (38)
```

Indeed, on every nonzero horizontal frequency they evaluate the even
Laurent polynomial `P(z)+P(z^(-1))` on all six inversion classes.  If all
six values vanish, then

```text
z^6(P(z)+P(z^(-1)))
```

vanishes at all twelve nontrivial thirteenth roots.  It also vanishes at
`z=1` because `P(1)=0`, so this degree-at-most-twelve polynomial is zero
(equivalently, it is `c Phi_13` with `c=0`).  Laurent-coefficient comparison
forces `P=0`.  Horizontal frequency zero is precisely (38).  Resolving the
diagonal through `Omega` and translating the cut restores those last six
modes.

In particular, if a row-zero array also has zero column sums, as does a
doubly centred ANOVA interaction, then (38) intersects it trivially and the
six converse-invariant sums already reconstruct it in one chosen cut chart.
This conditional algebraic corollary is not a physical-current theorem.

## 7. Exact affine covariance and the carry obstruction

Let

```text
g=(A,H;B,C) in AGL_1(F_13) x AGL_1(F_7)
```

act by

```text
d^g(h,r)=d(A^(-1)(h-H),B^(-1)(r-C)).                          (39)
```

On `K_14`, let

```text
phi_(A,H)(Omega)=Omega,
phi_(A,H)(x)=Ax+H                    for x in F_13.             (40)
```

Substitute `h=Ah_0+H`, `r=Br_0+C` into (9).  One obtains the exact edge and
degree covariance laws

```text
W^(d^g)_(tau,a,c)(e)
 =W^d_(A^(-1)tau,aB,aC+c)(phi_(A,H)^(-1)e),                   (41)

G^(d^g)_(tau,a,c)(v)
 =G^d_(A^(-1)tau,aB,aC+c)(phi_(A,H)^(-1)v).                   (42)
```

Thus the thirteenth-root carry is a geometric vertex relabelling, while the
septimal carry is a finite permutation of the cut chart.  There is no
representative-wrap error once the chart bundle is kept.

But the cut phase cannot be discarded.  In one chart, the semantic column

```text
r=-a^(-1)c                                                     (43)
```

is the hub star; every other column is a Hamilton cycle.  A nonzero
septimal translation changes which semantic column has type star.  No
`K_14` vertex automorphism can realize that exchange, because the two colour
classes have distinct degree multisets:

```text
star:                 (13,1,1,...,1),
cycle plus Omega:     (2,2,...,2,0).                           (44)
```

Equation (44) is the graph-theoretic form of THM-2507's cut seam.  The full
bundle is covariant; one unmarked quotient is not.

## 8. The two marginals still lose fifty-four dimensions

For a fixed nonzero `tau`, the paired marginal map on the row-zero space is

```text
d -> (R_tau d,R_(-tau)d).                                     (45)
```

At horizontal frequency zero, (45) kills the six-dimensional vertical
space.  At each of the twelve nonzero frequencies, the allowed polynomials
have degree at most six and vanish at one, hence dimension six.  Evaluation
at `zeta^(-alpha tau)` and `zeta^(alpha tau)` gives two independent
functionals.  Therefore

```text
rank(45)=12*2=24,
dim ker(45)=78-24=54.                                         (46)
```

The boundary has an explicit integral hostile for every `tau`.  Put

```text
T(k)=12, if k=0 mod 13, and T(k)=-1 otherwise,

d_tau(h,0)=-T(h),
d_tau(h,1)= T(h)+T(h+tau)+T(h-tau),
d_tau(h,2)=-T(h)-T(h+tau)-T(h-tau),
d_tau(h,3)= T(h),
d_tau(h,r)=0 for r>=4.                                       (47)
```

Then `d_tau` is integral, row-zero, nonvertical, has `L1=168`, and satisfies

```text
R_tau d_tau=R_(-tau)d_tau=0.                                  (48)
```

The joint edge chart retains (47), while both marginals erase it.  The
ordinary-degree reconstruction theorem succeeds only because it retains the
marked diagonal resolution and sweeps the cut phase.

## 9. General antipodal boundary

The same cardinality mechanism holds beyond `13 x 7`.  Let `p` be odd prime,
`2<=q<p`, and use the consecutive section

```text
I_q={0,...,q-1} subset F_p.                                   (49)
```

The unordered map

```text
(h,r) -> {h+tau r,h-tau r}                                   (50)
```

is injective exactly when `I_q` contains no two distinct nonzero antipodes,
equivalently

```text
q<=(p+1)/2.                                                    (51)
```

At equality it is bijective onto `Sym^2(F_p)` and, after diagonal
resolution, onto `E(K_(p+1))`.

THM-2507 gives at least `p-q+1` good nonzero slopes.  Among the
`(p-1)/2` antipodal slope pairs, this forces at least

```text
max(0,(p-2q+3)/2)                                             (52)
```

completely good pairs.  The bound is sharp in its abstract trace class.
Choose a set `B` of `q-2` bad slopes meeting as many antipodal pairs as
possible and form

```text
P_B(X)=(X-1) product_(b in B)(X-zeta_p^b).
```

Writing `P_B=sum_r D_rX^r`, the THM-2507 trace array

```text
d(h,r)=Tr_(Q(zeta_p)/Q)(D_r zeta_p^h)
```

is integral and row-zero, and its Radon output vanishes exactly at slope
zero and at the chosen slopes `B`.  Thus it realizes the extremal bad-pair
mask, not merely its combinatorial count.
Thus the threshold for a forced complete pair is

```text
p>=2q-1,                                                       (53)
```

the same as (51).  The case `(p,q)=(13,7)` is the critical equality: the
unordered chart just fills pair space and the abstract Radon count forces
exactly one pair in the worst case.

This section is a scoped corollary about the chart and the THM-2507 root
count.  The six-phase degree reconstruction above is stated and proved only
at `13 x 7`.

## 10. Exact gain and physical nonclaims

The proved object can be read in three compatible ways:

```text
13 x 7 row-zero strip
  = factor-balanced weights on the cyclic one-factorization of K_14
  = six ordinary-degree root currents plus a retained cut chart.          (54)
```

This is a genuine relational reframe.  It connects the Radon strip to the
parallel-class and circulant-tournament objects already present in the repo,
and it converts a lossless `91`-entry joint chart into a lossless bank of
one-dimensional vertex currents.

It does not close a physical LRC(14) step:

- `d` and `W` are signed; a degree is not the measure of one Boolean event;
- the thirteen finite vertices are root labels and `Omega` is the resolved
  diagonal mark, not proved identities of the fourteen runners;
- changing the cut phase swaps the star and cycle semantic columns, so the
  chart remains a required sidecar;
- the theorem supplies no source/arrival type, owner, clock, terminal word,
  old address, or deep ancestry sheet;
- pointwise nonzero degree vectors can cancel after integration over a
  moving parent; and
- the original THM-2436 application is in an already-empty branch.  A lawful
  signed transplant such as THM-2512 still does not force its non-replica
  alternative or turn a diagonal contraction into a Boolean current.

Thus no live row is excluded and LRC(14) remains open.  The next decisive
test is whether one factor-balanced edge weight in (54) can be realized on
the same owner/time/deep ancestry cospan as THM-2471/2478 without erasing
the marked cut phase.

## 11. Exact dependency-free referee

Compile and run

```bash
c++ -std=c++20 -O2 \
  04-computation/lrc14_k14_degree_reconstruction_thm2514_referee.cpp \
  -o /tmp/thm2514 && /tmp/thm2514

c++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_k14_degree_reconstruction_thm2514_referee.cpp \
  -o /tmp/thm2514-opt && /tmp/thm2514-opt
```

Both executions reproduce

```text
05-knowledge/results/lrc14_k14_degree_reconstruction_thm2514_referee.out
```

byte-for-byte.  The dependency-free referee checks:

- all `504` affine-cut charts as bijections with `E(K_14)`, all `6,552`
  matching rows, and all `3,528` star/cycle columns;
- `183,456` edge-covariance rows for four generators of the full product
  affine gauge;
- `550,368` ordinary-degree/Radon identities and `78,624` seven-phase
  zero-sum identities;
- exact rank `78` over `F_547` for all `504` six-phase banks, with every
  single phase having rank `13`;
- all `72` nontrivial multipliers in (34), all twelve zero-character
  vanishings, `404,352` coefficientwise instances of (30), and `33,696`
  marked-vertex instances of (31);
- rank `24` and kernel dimension `54` for all twelve antipodal marginal
  maps, together with all twelve integral hostiles (47);
- rank `72` for the six converse-invariant pair sums; and
- the injectivity and sharp good-pair boundaries through every odd prime
  `p<=31`.

The integer `547` is prime (trial division through `23` suffices), and
`547-1=2*3*7*13`.  The referee explicitly finds an element whose powers at
the four prime cofactor exponents are nontrivial, hence the field contains
primitive seventh and thirteenth roots.  A nonzero reduced multiplier is
also an exact characteristic-zero nonvanishing certificate.  Normal and
optimized transcripts and the source
and output hashes agree with the frontmatter. **QED.**
