---
id: THM-2260
title: "Positive finite-fibre capacity and thin-predicate boundaries"
status: >
  PROVED. A constant-Jacobian finite-fibre pushforward
  has image measure at least source mass times the Jacobian divided by the
  maximal fibre occupancy; complex weights inherit the support conclusion
  under a uniform open-half-plane sidecar on each fibre. This
  abstracts THM-2257 and gives two concrete positive transfers: every
  nonempty saving ideal in a fixed finite knot-prime alphabet has asymptotic
  box density one, and a monic planar Keller map of generic degree N obeys
  lambda_4(F(A)) >= |Jac F|^2 lambda_4(A)/N. The same audit proves sharp
  stopping boundaries. Every monic Keller image has lambda_4-null
  complement, so coarse image support misses the resultant nonproper curve;
  the DC2 grade-six response has infinite fibres and is not
  continuation-constant; the Mahler ordinary-integer itineraries project
  onto every finite carry cylinder while forming a countable null set;
  coefficient-uniform GMC2 support transfer fails already in two-element
  dihedral fibres; and tournament forced-cell support does not determine its
  multiplicity energy. No new JC2, DC2, Mahler-Z,
  positive-knot-catalyst, or LRC closure is claimed.
source: codex-2026-07-25-image-sieve-cross-frontier
depends_on:
  - THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation
  - THM-2111-effective-compound-root-bound-for-one-variable-constant-terms
  - THM-2172-frobenius-collapse-of-mahler-and-twisted-cyclic-packets
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-2240-dc2-grade-response-gauge-is-not-a-continuation-state
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2249-directed-triangle-forced-quotient-frustration
  - THM-2251-two-level-principal-ray-threshold-and-exposed-geodesic-face
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends
---

# THM-2260 -- when an image sieve transfers, and when it cannot

THM-2257 succeeds by replacing the mass of

```text
A=E_H intersection D_q
```

with the support of its image under the `169`-fold circle map. Its two
load-bearing properties are:

```text
the source mass is nonnegative;
every target fibre contains at most K_r source points.              (1)
```

This theorem extracts that mechanism, transfers it where (1) is genuinely
present, and gives exact hostile boundaries where one of its two clauses
fails or where the target predicate lives on a null sidecar.

## 1. The positive finite-fibre capacity lemma

Let `(X,mu)` and `(Y,nu)` be measure spaces, let `pi:X->Y` be measurable,
and let `A subset X` have finite measure and measurable image `pi(A)`. Put

```text
n_A(y)=#(A intersection pi^(-1)(y)).                           (2)
```

Assume that `n_A` is measurable and that a constant `J>0` gives the
disintegration identity

```text
integral_Y n_A(y) dnu(y)=J mu(A).                            (3)
```

If

```text
n_A(y)<=K 1_(pi(A))(y)                                      (4)
```

almost everywhere for an integer `K>=1`, then

```text
nu(pi(A))>=J mu(A)/K.                                       (5)
```

Indeed, integrate (4) and use (3). Equality in (5) requires the occupied
fibres to have multiplicity `K` almost everywhere on the image.

For the degree-`d` circle map `x->dx`, normalized Haar disintegration gives
`J=d`. Thus (5) is exactly THM-2257's inequality

```text
measure(169A)>=169 measure(A)/K_r.                           (6)
```

The sets in THM-2257 are finite unions of intervals up to endpoints, and
the circle covering is open, so the image-measurability hypothesis is
automatic there.

The lemma is deliberately a positive-mass theorem. Cardinality is its
finite version: for a map of finite sets and nonempty `A`,

```text
|pi(A)|>=|A|/max_y |A intersection pi^(-1)(y)|.             (7)
```

No topology, tournament orientation, or independence assumption is hidden
in either statement.

## 2. A sector is a sufficient coefficient sidecar; the orbit sum is exact

There is a weighted version, but it requires more than bounded fibre size.
Let every occupied fibre be finite, assign a nonzero complex weight `w_x`
to each `x in A`, and define

```text
W(y)=sum_(x in A intersection pi^(-1)(y)) w_x.              (8)
```

Suppose that for almost every `y` there is a unit complex number `zeta_y`
and one fixed `gamma>0` such that

```text
Re(conjugate(zeta_y) w_x)>=gamma |w_x|                      (9)
```

for every `x` in that fibre. Then

```text
|W(y)|
 >=Re(conjugate(zeta_y)W(y))
 >=gamma sum_(x in A intersection pi^(-1)(y))|w_x|.         (10)
```

Consequently

```text
{y:W(y)!=0}=pi(A),                                         (11)
```

up to null sets. This nonzero set is the coefficient-support meaning of
`support(W)` used below, and every capacity conclusion from Section 1
transfers.
A common sector of angular width `2 theta<pi` gives (9) with
`gamma=cos(theta)`.

The full orbit sum (8) is the exact coefficient sidecar; (9) is a convenient
uniform sufficient condition for its nonvanishing, not a necessary condition
for a particular fibre. Without any phase restriction, no support theorem
follows from `K` alone: in one two-point fibre, weights `1` and `-1` have
`K=2` and nonempty image support before weighting, but push forward to zero.
Bounded multiplicity controls how much positive mass can pile up; it does
not control coefficient phase.

## 3. NC2/GMC2: Frobenius transfers a packet, not positive image mass

THM-2070 supplies an internal hostile witness, not merely the abstract pair
`1,-1`. Put

```text
f(u)=u^2+u+u^(-1)-u^(-2).                                (12)
```

Its support has a balanced word at every length `m>=2`, but

```text
f(-u^(-1))=-f(u),                                        (13)
```

so `CT(f^m)=0` for every odd `m`.

For odd `m`, the ordered coefficient words of total charge zero are paired
by applying the involution in (13) to every letter. If `c_q` is the
coefficient at charge `q`, comparison of coefficients in (13) gives

```text
c_(-q)=-(-1)^q c_q.                                      (13a)
```

Hence a balanced word `(q_1,...,q_m)` and its image have weight ratio

```text
(-1)^m (-1)^(q_1+...+q_m)=(-1)^m.                        (13b)
```

Because the support has no charge zero, no ordered word is fixed, so its
involution orbit has exactly two elements. If ordered words are first
grouped into multinomial multiplicity vectors, a fixed vector would use
equal counts at charges `q` and `-q`, hence would have even total length;
there is again no fixed vector for odd `m`. In either representation the
two orbit weights are opposite. This is an exact `K=2` failure of (11),
inside an actual horizontal Gaussian face by THM-2070.

Thus THM-2257's positive image argument does not handle arbitrary radial
polynomial coefficients. THM-2022 succeeds for a different reason: after
Kummer removes every off-layer channel, Lucas and Frobenius preserve the
complete algebraic face sum as

```text
Q -> Q^p,                                                (14)
```

including all collisions. It does not replace `Q` by its support.

A second support-only boundary blocks a direct sharpening of THM-2111 by
subset-product compression. Let `1<=a<d` and choose

```text
alpha_i=2^(2^i),              0<=i<d.                    (15)
```

The products over the `a`-subsets are all distinct, because their base-two
exponent sums uniquely recover the subsets. Hence the compound map

```text
S -> product_(i in S) alpha_i                             (16)
```

has fibre cap one and image size exactly

```text
binom(d,a).                                               (17)
```

This example occurs on a literal fibre of THM-2111's root family: put

```text
P(X)=product_i(X-alpha_i),
R(X)=X^a-P(X),
Phi(X,t)=X^a-tR(X).                                      (18)
```

Then `Phi(X,1)=P(X)`, while `R` has nonzero constant and leading
coefficients. Therefore a support/multiplicity argument on the compound
roots alone recovers the current binomial degree and cannot reduce it to
`d`. Any sharp `M+N=d` improvement must use extra algebraic contact
structure, not finite-image compression.

## 4. JC2 and DC2: null codimension versus infinite memory

### 4.1. A planar Keller map satisfies the capacity inequality

Let

```text
F=(P,Q):C^2->C^2,              Jac_C(P,Q)=kappa!=0,       (19)
```

with `P` monic in `y`, and let `N` be the generic degree from THM-2241.
For a Borel set `A subset C^2` of finite four-dimensional Lebesgue measure,
the continuous image `F(A)` is analytic and hence Lebesgue measurable. The
real area formula gives

```text
integral_(C^2) # (A intersection F^(-1)(z)) d lambda_4(z)
 =|kappa|^2 lambda_4(A).                                (20)
```

The real Jacobian determinant is `|kappa|^2`. One may prove (20) by taking
a countable injective local-chart cover for the etale map, refining it to
disjoint Borel pieces contained in those charts, and applying ordinary
change of variables on each piece.

THM-2241 proves

```text
#F^(-1)(z)<=N                                           (21)
```

at every target point: it equals `N` off the resultant leading-coefficient
zero set and drops on that set. Sections 1 and (20)--(21) therefore give

```text
lambda_4(F(A))>=|kappa|^2 lambda_4(A)/N.                (22)
```

This is a genuine transfer of the finite-fibre capacity lemma. On the
automorphism branch `N=1`, every occupied fibre has multiplicity one and
the area formula makes (22) an equality.

### 4.2. Why it cannot see the planar Jacobian obstruction

Let `c(U,V)` be THM-2241's leading resultant coefficient. That theorem gives

```text
C^2\V(c) subset F(C^2),                                 (23)
```

because every point off `V(c)` has exactly `N` preimages. The zero set of a
nonzero polynomial in two complex variables has four-dimensional Lebesgue
measure zero. An elementary proof writes the polynomial in `V`: outside
the measure-zero set of `U` where all its coefficients vanish, each vertical
line has only finitely many roots, and Fubini applies.

It follows that every monic Keller map in this scope has

```text
lambda_4(C^2\F(C^2))=0.                                (24)
```

If `c` is constant, THM-2241 proves that `F` is an automorphism. If a
nonautomorphic planar Keller map existed, its nonproper set would be the
nonconstant curve `V(c)`, but (24) would still report full-measure image.
Thus four-dimensional image support forgets the codimension-one resultant
sidecar. The capacity constant (22) retains the generic degree `N`; proving
`N=1`, not proving full-measure image, is the unresolved algebraic
obligation. This is only a stopping boundary for the coarse support
observable. It is not a converse, an equivalence, or a claim that every
refined fibre-capacity argument is incapable of detecting `N` or `c`.

### 4.3. The DC2 response does not have finite fibres

For THM-2240's grade-six response map `d_6`, the hypothesis (4) fails before
any measure issue arises. A fixed response fibre contains the arbitrary
polynomial axes

```text
J(q,u),       a(u),       b(u),                              (25)
```

so it is infinite-dimensional. Moreover the exact next residual is not
constant on that fibre: the homogeneous syzygy `C=1` changes it by

```text
8q(4u^2-13u+13).                                        (26)
```

There is therefore no finite `K` to divide by, and quotient image support
cannot be a continuation state. A representative-sensitive next-rung
sidecar is mandatory. This is an infinite-memory failure, distinct from
JC2's null-codimension failure.

## 5. Mahler: exact support transfer and a thin inverse-limit target

There are two different Mahler objects in the current canon.

First, THM-2172's good-prime twisted cyclic quotient has

```text
r -> pr mod n                                           (27)
```

as a permutation when `p` does not divide `n` and the twist is nonzero.
Its coefficient weights are also nonzero. Thus `K=1` and Frobenius
preserves labelled coefficient support exactly. This is a positive,
lossless image-support transfer. At either bad-prime boundary THM-2172 gives
a nonzero element whose Frobenius image vanishes, so nonzero weights are
load-bearing.

Second, THM-2228's length-`m` carry map

```text
{0,1}^m -> Z/2^m Z,             prefix -> r_m             (28)
```

is a bijection. It also has `K=1`. Nevertheless the infinite words arising
from positive ordinary integers form a countable, dense, and co-dense set
`I_+` in the binary shift. Hence:

```text
every finite cylinder meets I_+;
every finite cylinder meets its complement;
Bernoulli measure(I_+)=0.                                (29)
```

In particular, the projection of `I_+` onto every finite image (28) is the
entire image even though `I_+` is null. The `K=1` prefix-support data in
(28), by itself, cannot detect eventual stabilization to an ordinary
integer.

The two exact hostile controls from THM-2228 separate the missing sidecars:

```text
(100)^infinity: all real suffix tails are safe,
                but the 2-adic state is -9/19 and never stabilizes;

the itinerary of A=1: ordinary stabilization holds,
                      but the real tail inequality fails.            (30)
```

Thus the Mahler `Z`-number predicate is an intersection of a real tail
condition with a thin termination condition. The image sieve sees every
finite carry address but neither their infinite conjunction nor ordinary
termination.

## 6. Knots: a nonempty finite-alphabet saving ideal has density one

Let `I subset N_0^r` be any nonempty additive upper ideal, and choose

```text
a=(a_1,...,a_r) in I.                                   (31)
```

For

```text
B_n={0,...,n}^r
```

and `n>=max_i a_i`, upper closure gives

```text
{b in B_n:b_i>=a_i for every i} subset I intersection B_n.
```

The translation map from the smaller orthant box into this set is injective,
so it is the `K=1` case of (7). Therefore

```text
|I intersection B_n|
 >=product_(i=1)^r (n-a_i+1),                           (32)

|B_n\I|
 <=sum_(i=1)^r a_i (n+1)^(r-1),                         (33)

lim_(n->infinity) |I intersection B_n|/|B_n|=1.          (34)
```

Equation (33) follows by covering the complement of the displayed orthant
by the `r` slabs `b_i<a_i`.

Apply this to THM-2191. Fix oriented prime-knot types
`P_1,...,P_r` and identify their connected-sum monoid with `N_0^r`.
For every knot `K` and saving threshold `s`, if

```text
I_s(K) intersection <P_1,...,P_r>                       (35)
```

is nonempty, then its exponent vectors have asymptotic box density one.
The same holds for optimal contexts whenever the chosen alphabet contains
one. Thus one catalyst in a finite alphabet makes almost every sufficiently
large box context at least as effective.

This does not produce the first catalyst. THM-2251's exact word-metric
models realize an arbitrary principal-ray threshold

```text
N in N_0 union {infinity}.                              (36)
```

For finite `N`, the saving support is a tail and hence has density one; for
`N=infinity`, it is empty. The abstract image-capacity theorem is therefore
sharp: it amplifies a nonempty support but cannot prove that the support is
nonempty or bound its first element.

## 7. Tournaments: support detects the zero layer, not the scale

For THM-2249's directed-triangle quotient, put

```text
Sigma(X)={(i,j):X_ij X_(i+1,j-1)>0}.                     (37)
```

Here both indices are modulo three. Because every term is nonnegative,

```text
F(X)=sum_(i,j)X_ij X_(i+1,j-1)=0
 iff Sigma(X)=empty.                                    (38)
```

THM-2249 then proves that the empty support is exactly a whole-block cyclic
rotation. This is a valid discrete image-support transfer.

Concretely, take as source atoms the forced reversed vertex pairs counted
by THM-2249 and send a pair from transport cells `(i,j)` and
`(i+1,j-1)` to the label `(i,j)`. Its fibre multiplicity is exactly
`X_ij X_(i+1,j-1)`, and its image is `Sigma(X)`. Thus (38) is an actual
finite pushforward, not a tournament metaphor.

Binary support does not retain the positive frustration scale. At margin
`N=2`, the two transport matrices

```text
X =[[0,0,2],          X'=[[0,0,2],
    [1,1,0],              [0,2,0],
    [1,1,0]],             [2,0,0]].                    (39)
```

both have row and column sums two and exactly the same active cells

```text
Sigma(X)=Sigma(X')={(0,2),(1,1),(2,0)}.                 (40)
```

But their three active products are respectively

```text
(2,1,2) and (4,4,4),
```

so

```text
F(X)=5,                 F(X')=12.                       (41)
```

Thus forced-cell support preserves the automorphism zero layer but destroys
the multiplicity energy needed for THM-2249's sharp floor and pinned
transport certificates. The full transport matrix `X`, the rotation
assignment costs, and the residual `G-F` remain the required sidecars.

## 8. Transfer ledger and scope

The exact source-to-target contracts are:

| Source | Target and map | Preserved | Lost coordinate / needed sidecar | Cheapest hostile test |
|---|---|---|---|---|
| THM-2257 LRC overlap | image support under `x->169x` | positive mass after a fibre cap | root labels and endpoint state; retain `K_r` and the resonant strip parameter | the four `K_r=25` classes |
| NC2/GMC face channels | channel-orbit support or `p`-dilation | balanced feasibility; the complete sum only under Frobenius | arbitrary coefficient phase; retain the sector or full orbit sum `Q` | THM-2070's odd `K=2` dihedral pairs |
| monic Keller map | measurable image under `F:C^2->C^2` | area capacity and the explicit degree parameter `N` | the null resultant curve; retain `c(U,V)` and fibre degree | `C^2\V(c)` already has full measure |
| DC2 grade response | current-response quotient | grade-six cancellation | infinite next-residual memory; retain a representative or next-rung residual | `C=1` splitter (26) |
| Mahler carry words | finite prefix under `c->r_m` | every finite address, with `K=1` | real tail and ordinary stabilization; retain both | `(100)^infinity` and `A=1` |
| knot contexts | saving support under addition by a fixed context | every already-realized saving | existence and first threshold; retain a minimal context antichain | arbitrary `N`, including infinity |
| tournament transport | active forced-cell support `Sigma(X)` | the automorphism zero layer | cell multiplicities and internal kernel; retain `X,A_a,D,G-F` | matrices (39) |

The reusable rule is therefore:

```text
positive mass + bounded fibres
  -> quantitative image capacity;

complex weights
  -> first retain a sector or the complete orbit sum;

thin/null target predicate
  -> retain its algebraic, termination, or continuation sidecar.     (42)
```

This theorem adds the finite-alphabet density-one knot corollary and the
Keller area-capacity formula. It does not improve THM-2111's first-return
bound, prove a positive knot catalyst, decide a Mahler `Z`-number, detect
the JC2 nonproper curve from measure, terminate the DC2 lift, strengthen
THM-2249's frustration floor, or alter LRC(14)'s 165-row frontier. QED.
