---
id: THM-2516
title: "Anchored moment simplex disintegration and owner-star recovery"
status: >
  PROVED + INDEPENDENTLY AUDITED.  On a finite cyclic probability space, an
  unanchored m-th moment
  is the exact average of fibres with m-1 difference coordinates, reflecting
  C_D^m/C_D,diag.  Multiplying by one positive marked anchor adds one vertex
  and one difference: h A_j^m is the average of the nonnegative owner-supported
  m-arm fibres E_x H(x) product_q F_j(x+u_q).  Hence any nonzero linear
  coefficient of an m-th moment survives on one owner-supported rational
  simplex fibre, with magnitude at least h times the original coefficient.
  For a C_7-by-C_13 response bank, one primitive mixed coefficient then
  Galois-saturates all 72 colours and all 5,184 THM-2508 cut coefficients.
  For Boolean densities, coincident arms collapse exactly to the previous
  moment order, giving a literal toothpick self-similarity.  At m=2 the
  carrier is an owner-rooted triangle/two-arm star.  Its arms remain unordered
  and need not be lawful phase-bank, clock, collision, or deep-ancestry edges;
  no one-point Boolean current, row exclusion, or LRC(14) follows.
source: codex-2026-07-27-anchored-moment-simplex
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2515-haar-self-correlation-disintegration-and-rational-shift-recovery
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
---

# THM-2516 -- moments are difference simplices; an owner adds one arm

**PROVED + INDEPENDENTLY AUDITED.**

THM-2515 disintegrates an entrywise square into one-difference
autocorrelation fibres.  Its sharp hostile shows why that operation can lose
a separately marked owner anchor.  The repair is not to force the owner onto
the same edge.  Keep it as a distinguished vertex and disintegrate the other
moment legs relative to it.

The resulting coordinate count is structural:

```text
unanchored m-th moment:      m points modulo diagonal translation
                            = m-1 differences;

owner-anchored m-th moment:  one marked point plus m moment points
                            modulo diagonal translation
                            = m differences.                   (1)
```

For the square, (1) is an owner-rooted triangle.  This is the smallest
canonical difference chart produced by the product geometry, not a claim of
absolute minimality among every possible nonlinear observer.

## 1. The unanchored moment quotient

Let `J` be a finite label set, let `D>=1`, and let

```text
F_j:C_D -> Q,                           j in J,

A_j=1/D sum_x F_j(x).                                          (2)
```

Fix an integer `m>=1`.  For

```text
t=(t_2,...,t_m) in C_D^(m-1)                                 (3)
```

(with no coordinates when `m=1`), define

```text
U_t(j)=1/D sum_x F_j(x) product_(q=2)^m F_j(x+t_q).           (4)
```

Then

```text
1/D^(m-1) sum_t U_t(j)=A_j^m.                                (5)
```

Indeed,

```text
(x,t_2,...,t_m)
 -> (x_1=x,x_2=x+t_2,...,x_m=x+t_m)                           (6)
```

is a bijection `C_D^m -> C_D^m`.  The average in (5) is therefore the
product of `m` independent copies of the mean in (2).

Equivalently, simultaneous translation is the diagonal action on `C_D^m`,
and (3) is a global difference chart on its quotient.  At `m=2`, equations
(4)--(5) are exactly THM-2515's autocorrelation disintegration.

## 2. A marked anchor adds one difference coordinate

Let

```text
H:C_D -> Q,
h=1/D sum_x H(x),                                              (7)
```

and for

```text
u=(u_1,...,u_m) in C_D^m                                     (8)
```

define the anchored moment fibre

```text
T_u(j)=1/D sum_x H(x) product_(q=1)^m F_j(x+u_q).             (9)
```

Then

```text
1/D^m sum_u T_u(j)=h A_j^m.                                 (10)
```

For fixed `x`, each sum over `u_q` is `D A_j`.  Hence the left side of
(10) is

```text
1/D^(m+1) sum_x H(x) product_q sum_(u_q)F_j(x+u_q)
 =1/D^(m+1) (D h)(D A_j)^m
 =h A_j^m.                                                    (11)
```

Alternatively, `(x,u_1,...,u_m)` parametrizes the full `m+1`-point
configuration

```text
x, x+u_1, ..., x+u_m                                        (12)
```

while `u` alone parametrizes the quotient by common translation and `x`
parametrizes the integrated diagonal fibre.  The function `H` marks the
first point, so there is no legal gauge that sets one of the `u_q` to zero
while leaving an arbitrary anchor density unchanged.

Let `K` be a subfield of `C` and let

```text
Lambda:Q^J -> K                                               (13)
```

be linear.  Equations (10)--(11) imply

```text
h Lambda(A^(circ m))!=0

 => there is u in C_D^m with Lambda(T_u)!=0,                  (14)

max_u |Lambda(T_u)|>=h|Lambda(A^(circ m))|                   (15)
```

when `h>0`.  If `H` and all `F_j` are nonnegative rational functions, every
selected table `T_u` is nonnegative and rational.  More importantly, every
one of its entries is represented by an integrand containing the same
nonnegative anchor factor `H(x)` of positive mean.  This is owner support at
the marked base point; it is not a claim that the selected table's own
untwisted anchor cell has positive triple overlap.

## 3. Circle and base-`13` versions

For a finite bank of rational step functions on `T=R/Z`, take a common grid
denominator `D` and identify cell values as in THM-2515.  Then (9) is

```text
T_u(j)=integral_T H(x) product_(q=1)^m F_j(x+u_q/D) dx,       (16)
```

and all conclusions are exact.

The Haar version holds for bounded functions on any compact abelian
probability group:

```text
integral_(u_1,...,u_m) T_u(j) du_1...du_m
 =h A_j^m.                                                    (17)
```

For bounded step data, `u->Lambda(T_u)` is continuous: telescope the finite
product and use translation continuity with all other factors bounded.
Arbitrary `L^2` inputs alone are not claimed at higher moment order.  If the
integral in (17) is nonzero, its nonzero locus contains an open subset of
`T^m`.  Consequently, for arbitrarily large `L`, all arms can be chosen on
one base-`13` grid:

```text
u_q=k_q/13^L,                  q=1,...,m.                      (18)
```

At rational shifts the table remains rational.  As in THM-2515, the
common-`D` selection gives the exact magnitude invoice (15); the dense
base-`13` selection guarantees nonvanishing but not that same invoice.

## 4. Primitive cut spectrum on an owner-supported fibre

Now take

```text
J=F_7 x F_13                                                   (19)
```

and suppose, for one `kappa,b!=0`,

```text
Lambda_(kappa,b)(A^(circ m))
 =1/91 sum_(ell,s)A_(ell,s)^m
       xi_7^(kappa ell)zeta_13^(b s)
 !=0.                                                         (20)
```

Let `H` be any positive rational step anchor with `h>0`.  Equations
(14)--(15) choose an owner-supported rational fibre `T_u` with this mixed
coefficient nonzero.  Rational Galois transitivity then gives

```text
Lambda_(kappa',b')(T_u)!=0

for every kappa' in F_7^*, b' in F_13^*.                     (21)
```

ANOVA-centre `T_u`, transpose the interaction to a row-zero
`F_13 x F_7` array, and apply THM-2508.  Exactly as in THM-2515,

```text
Psi^u_(tau,a)(alpha,beta)
 =91 K_(alpha tau,beta)Lambda_(beta a,-alpha)(T_u)
 !=0                                                          (22)
```

for every nonzero `tau,a,alpha,beta`.  Thus all `5,184` primitive cut
coefficients survive on one fibre whose every entry retains the common
anchor factor `H(x)`.

At a fixed lawful THM-2449 clock in its baseline owner chart, take the genuine
Boolean owner/word event before the deep-label sum,

```text
H_own(x)=U_(0,0)(x)d_(j,0)(x)Q^epsilon_0(Rx).                 (22a)
```

If `N(x)=sum_r Delta_r(x)`, then `N in {1,2}` almost everywhere and

```text
F_(0,0)=H_own N,
a=A_(0,0)=integral H_own N.                                  (22b)
```

Thus, with `h=integral H_own`,

```text
a/2<=h<=a.                                                    (22c)
```

Use `H=H_own` in (9) and the full response densities `F_(ell,s)` on its
arms.  The selected fibre is then literally supported on the marked Boolean
owner event, and its common-grid coefficient satisfies the strengthened
invoice

```text
|Lambda(T_u)|>=(a/2)|Lambda(A^(circ m))|.                     (22d)
```

This still does not assert that the selected table's own `(0,0)`
`(m+1)`-fold overlap entry is positive.

The theorem is conditional only on the moment coefficient in (20).  It does
not claim that every response table has such a coefficient at every order.
When a first-or-second-moment dichotomy supplies the square channel, `m=2`
is the relevant application.

## 5. Toothpick self-similarity and the triangle carrier

The symmetric group `S_m` permutes the arms of (9), and `T_u` is invariant
under this permutation.  If every `F_j` is Boolean, then on a diagonal face

```text
u_r=u_s,
```

the two repeated factors satisfy `F_j^2=F_j`.  Deleting either repeated arm
therefore gives exactly the anchored `(m-1)`-st moment fibre.  Thus the
simplex bank has literal recursive faces:

```text
m-arm carrier
  -- collide two arms --> (m-1)-arm carrier.                   (23)
```

For non-Boolean weights, the diagonal retains `F_j^2`; Boolean idempotence is
load-bearing.  In particular, the live THM-2449 density after its deep-label
sum can take the values one and two.  The recursive face law is an abstract
Boolean boundary, not an assertion about that summed response bank.  Only an
arm--arm collision `u_r=u_s` has this recursion.  An arm at the owner point
`u_q=0` does not collapse unless the separate anchor `H` idempotently absorbs
that particular `F_j`.

At `m=2`, the three points are

```text
x                    marked owner vertex,
x+u, x+v             two moment legs.                          (24)
```

They form an owner-rooted triangle, equivalently a two-arm toothpick star;
the third edge has difference `v-u`.  The two arms are unordered because
the square is symmetric.  THM-2515's equal-cell hostile proves that
collapsing this data to the single unanchored difference `v-u` can separate
the spectral signal from the owner anchor.  Retaining `(u,v)` is the
canonical repair supplied by (10).

## 6. Exact remaining physical boundary

For a phase bank `{d(cx-r/p):r in C_p}`, an arm `u_q` preserves the bank as a
label permutation only if

```text
p c u_q in Z.                                                  (25)
```

This condition is also sufficient for bank-set preservation; keeping every
label fixed requires the stronger `c u_q in Z`.  For a delayed word, define

```text
Stab_T(Q)={t in T:Q(x+t)=Q(x) almost everywhere}.              (26)
```

Then `Q(13^K x)` is literally unchanged under the arm exactly when

```text
13^K u_q mod 1 lies in Stab_T(Q).                              (27)
```

The simpler condition `13^K u_q in Z` is sufficient, and is necessary when
`Q` has trivial translation stabilizer.  Every arm must satisfy these
conditions simultaneously.  The open spectral
locus in Section 3 need not meet this finite lawful product grid; a retained
septimal factor with `13`-unit speed again forces each base-`13` arm to be
zero.  Thus owner support has been restored algebraically, but lawful
ancestry has not.

The carrier (24) also does not identify its two arms with THM-2471 collision
edges, THM-2478's old deep sheet, or an address/gain edge.  In the THM-2449
application every vertex retains the same fixed terminal-word label and
clock, but one common word phase is not proved; agreement is governed by
the exact stabilizer test (27), with integrality only a sufficient shortcut.
The carrier is a nonnegative
three-point cospan observable, not a one-point Boolean current.  The next
decisive test is now precise: realize both differences
from one marked owner root inside the existing collision/ancestry graph, or
prove that the surviving semantic rows forbid such a triangle.  No live row
is removed and LRC(14) remains open.

Two independent audits rederived the all-`m` normalization (including exact
rational controls through `m=4`), quotient coordinate count, magnitude
selection, base-`13` density argument, Galois/cut signs, and Boolean
arm-collision recursion.  They also repaired three scope points before final
promotion: `x` parametrizes the diagonal fibre rather than the quotient;
the live summed response is not Boolean; and word invariance is governed by
the stabilizer (27), not by integrality as a necessary condition. **QED.**
