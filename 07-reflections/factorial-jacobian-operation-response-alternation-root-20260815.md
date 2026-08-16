# Factorial--Jacobian alternation: characters, currents, and response length

**Status:** CURRENT RESEARCH SYNTHESIS, 2026-08-15.  The theorem statements
below are routed to THM-3465/3466 and the strengthened THM-3383.  `FC(3)`,
`HFC(3)`, `JC(2)`, and `DC(2)` remain open.

## 1. Portfolio and inheritance

The session alternated rather than merging the frontiers.

| lane | inherited mechanism | canonical hostile / near miss | least-used sidecar |
|---|---|---|---|
| Anchor: factorial exact frontier | THM-3152/3201 first Euclidean triple and multi-place Newton barcodes | THM-3210's invisible/invisible/visible cancellation ray; a two-step multiplication kernel can exit at the third power | the full progressive degree-set trace, not the final survivor count |
| Niche: constant-J HFC | THM-3303 boundary collision and THM-3328 cone/anti-tangent passport | the slit-annulus seam; scalar moments and first tangent data descend | Hermitian dagger and a pointed boundary-current primitive |
| Wildcard: polynomial effectivity | THM-3383's oriented Laurent intersection and THM-3397's two-boundary filtration | killing the cyclic torsor and class group still leaves poles | the operation induced by the surviving polynomial coordinate on the additive quotient |

The incoming branch-archaeology emphasis on typed source currents suggested
retaining the triangle edge labels.  The incoming THM-3397 effectivity message
then supplied the hostile that forced a uniform-response audit.  Neither was
used as historical truth: the exact current was rederived by Stokes, and the
Laurent claim was checked directly against THM-3383's membership criterion.

## 2. Live concept board

| concept | source -> target map | status after the session |
|---|---|---|
| pure nonreal cyclic character | `g -> (g,g^dagger) -> {g,g^dagger}` | PROVED all-degree Keller rigidity; triangle HFC intersection empty |
| mixed-character top layer | ordinary leading form plus dagger | THM-3470 closes one reflection-even triangular quartic slice at moment 5; general cell OPEN |
| factorial Boolean face descent | `prod(1-partial_i)` from bulk to coordinate faces | PROVED identity; external access boundary survives |
| triangle current | `dbar(g)` under multiplication by `g` and boundary integration | PROVED exact length-two Krylov block under HFC+constant J |
| pointed clutch | `d kappa=g^2 dbar(g)` on labelled source edges | PROVED existence/degree; universal sheet separation REFUTED |
| Laurent effectivity response | multiplication by the polynomial target on `B/(A intersect B)` | PROVED locally nilpotent with unbounded exact response lengths |
| prime-power carry face | first Euclidean triple at `d=a p^k+1` | PROVED complete singleton ledger for every `2<=a<p` |
| binary-submask compiler | `NP_2(A_N^(N+1))` intersect odd reset | PROVED all `2^s p^k` windows with `2<=2^s<p` |
| all-divisor digit compiler | exact positive `F` suffix intersect THM-3161 `G` hull | PROVED pair barcode for every `p|d-1`; 33/38 residuals close |
| seven-exit factorial boundary | first Euclidean triple at `d=2501,2502` | both closed; the wider `2501..2600` audit is intentionally not yet a theorem here |

## 3. First alternation: factorial lowering becomes Keller flux

For the factorial functional

```text
L_n(h)=integral_(R_+^n) h exp(-sum x_i),
```

integration by parts gives

```text
B_I(h)=L_n(prod_(i in I)(1-partial_i)h).                 (1)
```

At a null power,

```text
m L_n(f^(m-1)partial_i f)=-B_i(f^m).                    (2)
```

Thus derivative multipliers are exactly face debt; they are not consequences
of the scalar moment orbit.  Divisibility by `x_i` is the sharp hostile: it
kills the entire corresponding face response while leaving the interior
uncontrolled.

For a planar Keller pair `{P,Q}=c`, Hamiltonian integration by parts makes a
related multiplier endogenous, but only with both exponential drift and
labelled face flux retained.  This is the lawful FC--JC interface.  Dropping
the drift or merging the two faces fails already on the triangular
automorphism `(P,Q)=(x,y+x^2)`.

## 4. Second alternation: the Hermitian bracket closes every pure character

Let an orientation-preserving finite Euclidean rotation act in conjugate
coordinates by

```text
rho(z,w)=(xi z,xi^(-1)w).
```

If `g o rho=eta g` with `eta` nonreal and `g^dagger` is the coefficient-
conjugate `z<->w` mate, then a constant nonzero real Jacobian is equivalent to
`{g,g^dagger}` being a nonzero constant.  If the top degree is `d>1`, its top
bracket vanishes.  Equal-degree binary forms with zero bracket are
proportional:

```text
f=w^dF(z/w), h=w^dH(z/w)
=> {f,h}=d w^(2d-2)(F'H-FH').                            (3)
```

But the two top forms have the distinct characters `eta,eta^(-1)`.  Hence
`d=1`; moreover `eta` must be one of the two source linear characters.

On the triangle, the only survivors are `Az` and `Aw`, and

```text
<z^3>=<w^3>=1/10.                                      (4)
```

So no pure nontrivial `C3` eigencell, in any degree, meets the constant-J HFC
sector.  This does not close the cyclic support-five moment cell: it says its
Jacobian is zero or nonconstant.

The sidecars are sharp.  Without dagger,

```text
(g,h)=(w+Bz^2,z),             {g,h}=-1,                 (5)
```

and for the real order-two character the odd shear `(x,y+x^3)` is nonlinear
with Jacobian one.  Finite character syntax alone is therefore not the
theorem.

The first mixed cell behaves differently but is still finite.  On the
normalized triangle, with `U=x-y` and `M=x+y-2/3`, the real family

```text
g=aM+b(U^2-1/6)+c(U^4-1/15)+iU                         (5a)
```

has `Jac(Re g,Im g)=-2a` and an explicit triangular polynomial inverse when
`a!=0`.  Exact elimination shows that moments two through four have `24`
reduced complex solutions, exactly four real; all four are tame automorphisms
with nonzero support in every `C3` character and vanish through moment four.
Moment five excludes each one, both by a unit Groebner ideal and by a coprime
lex remainder.  Thus small HFC prefixes can occur even on automorphisms; the
fifth-moment gate is universal only in this three-parameter shear slice.

## 5. Third alternation: HFC moments become a boundary Krylov block

For `g=p+iq` on an oriented triangle and `J={p,q}`,

```text
integral_boundary g^n dbar(g)
  =-2ni integral_Delta g^(n-1)J dA.                      (6)
```

When `J=c!=0` is constant and all positive HFC moments vanish, the response of
`alpha=dbar(g)` under `T(alpha)=g alpha` is

```text
C(alpha)=0,
C(T alpha)=-2ic Area(Delta)!=0,
C(T^j alpha)=0 for j>=2.                                 (7)
```

Modulo the all-future invisible kernel, `[alpha],[T alpha]` form exactly one
nilpotent Jordan string of length two.  The all-order HFC hypothesis proves
the terminal absorption; two sampled zeros do not.

The exterior Cauchy transform is the pure dipole

```text
integral_boundary dbar(g)/(zeta-g)
  =-2ic Area(Delta)/zeta^2                               (8)
```

on the unbounded complementary component.  Consequently zero lies on the
boundary image or in a bounded complementary cell.  Also

```text
d kappa=g^2 dbar(g)                                      (9)
```

has a periodic piecewise-polynomial primitive on the labelled source
boundary.  It is a lawful clutch coordinate, but not a universal separator.
At an anti-tangent collision its first target jet agrees automatically; the
first possible difference is a signed-curvature term in the second jet.  At
`g=0` that jet is blind, and on the slit-annulus seam the whole primitive
descends as `(g^3-rho^3)/3`.

## 6. Fourth alternation: classify response length before using a finite bank

The univariate factorial hostile

```text
f=x^2+(-4+2i)x+(2-2i)
```

satisfies

```text
L(f)=L(f^2)=0,       L(f^3)=32+80i,       L(f f^dagger)=8. (10)
```

Thus `ker L intersect M_f^(-1)ker L` is only a depth-two address.

There is also a moving-depth face hostile.  If `ord_(x_i)(f)=e_i>=1`, then
for every `beta_i<m e_i`, repeated weighted integration by parts encounters
no boundary jet and gives

```text
L_n(partial^beta f^m)=L_n(f^m).                           (10a)
```

Hence on a null orbit every fixed pure-derivative bank becomes blind for all
sufficiently large powers; reaching the boundary requires derivative order
growing linearly with `m`.  THM-2846's first-three-null polynomial, multiplied
by the remaining coordinates, realizes this exactly on a finite prefix.
This still does not endogenize polynomial multipliers.

At the opposite extreme, in THM-3383's positive Laurent orientation let
`E=A intersect B=C[u,v,y^e]`, `B=C[u,t]`, and `v=ut`.  Multiplication by `u`
acts on the additive quotient `B/E`, and for every `q>=1`,

```text
u^k t^q notin E for k<q,        u^q t^q=v^q in E.        (11)
```

The negative orientation swaps `u,t`.  Every class is eventually killed, but
the exact killing times are unbounded.  The finite residue field, torsor, and
class group therefore do not bound polynomiality debt.

Equations `(7)`, `(10)`, and `(11)` give three genuinely different regimes:

```text
uniform length two / finite prefix with later exit /
objectwise finite but uniformly unbounded.                         (12)
```

This is the session's reusable cross-frontier result.  Before turning any
sidecar into a finite certificate, classify which regime its native operation
occupies.

## 7. Prime-power carry face and exact factorial boundary

The exact `d=2501` signal suggested looking at its predecessor

```text
d-1=2500=4*5^4.
```

That pattern survives proof rather than only interpolation.  For `p>=5`,
`H=p^k`, `2<=a<p`, and `d=aH+1`, the resonant rows `F=A_(aH-1)`,
`G=A_(aH)` and their first full Euclidean row share the reduced slope

```text
sigma=2(H-1)/((p-1)H),                                  (13)
```

with exact common capacity

```text
min(a-1,(p-1)/2)H.                                      (14)
```

Lucas and Kummer expose the coefficient valuations; two base-`p` digit-sum
supporting inequalities isolate all equality anchors.  The extra work is the
Euclidean projection: its shifted `F` term has a unique minimum because the
three scalar coefficients reduce to `-3,-6,2`, all units for `p>=5`.  This is
the missing sidecar in the older verified large-divisor face probe.  At
`a=1` only a zero-capacity supporting contact remains.  At `p=3` the
Euclidean-row shift changes, and allowing `a>=p` introduces a second carry
scale; neither boundary is cosmetic.

The local address in fact closes for every `2<=a<p`.  When
`a<=(p-1)/2`, the displayed `G` face spans its full degree.  Above that
threshold, digit-sum supporting lines give complete `G` slopes
`sigma,tau,rho_G`, while `F` has only `sigma` in common because its terminal
slope `rho_F` is strictly larger.  Thus

```text
D_p={0,p^k,2p^k,...,T p^k},   T=min(a-1,(p-1)/2).        (14a)
```

At the former first unaudited row `d=2501,r=2499`, both THM-3201 engines give
the identical progressive factor-degree trace

```text
p=2: {256,2048,2304},
p=3: {256},
p=5: empty.                                             (15)
```

The common semantic digest is
`365533925519a4d8d44db78394f0785e87be5f4cc03e0a98d759f93609fb09ee`.
The `p=5` polygon has the single relevant slope `312/625`, cap `1250`, and
denominator `625=5^4`; because `a=4=p-1`, its complete local degree set
`{625,1250}` is structural and incompatible with degree `256`.  This closes the row
FINITE-EXACT and extends the exact quadratic window boundary to `r=2499`.

The binary place then changes the quantifier.  THM-3161's proved digital
skeleton specializes to a complete `2`-polygon with one edge for each set bit
of even `N`; every factor degree is a binary submask of `N`.  Intersecting
this with (14a) proves all windows

```text
d=2^s p^k+1,                 2<=2^s<p,                   (15a)
```

for every `k>=1`.  Thus `d=2501` sits in the all-height family
`d=4*5^k+1`, rather than being only an isolated computation.

The next row is a sharper adaptive-place control.  At `d=2502`, the
`61`-adic singleton reset forces a common-factor degree to be divisible by
`61`.  THM-3161's exact `41`-adic polygon of the second resonant row has
block denominators `41` and `41^2`, so it independently forces divisibility
by `41`.  A common factor also divides the degree-`2499` Euclidean row, but
`41*61=2501`; hence none exists.  The old reset-plus-`2,3` trace leaves the
hostile address `{61}`, while the adaptive `61,41` trace is empty.  This
extends the contiguous exact quadratic boundary to `r=2500` and identifies
the right design principle: compile divisor places from `d-1`, retaining
their full capacities, instead of asking a fixed small-prime bank to guess
the useful scale.

The adaptive observer can be made exact for every divisor prime, not only a
one-digit prime power.  For `p|N=d-1`, the first resonant row has a strict
two-plateau digit recursion (an odd-anchor recursion at `p=2`); every possible
positive-hull vertex is an actual coefficient unit.  Its unresolved Kurepa
prefix has nonpositive slope, while THM-3161's second-row polygon has only
positive slopes.  Their full common ledger is therefore computed from digit
weights alone, with no large coefficients.  Two implementations apply this
to the `38` seven-exit residuals through `d=2600`: `33` close, and only

```text
d=2516,2564,2571,2576,2586
```

retain explicit degree packets.  The first is
`{503,1006,1509,2012}`, so the contiguous quadratic boundary reaches
`r=2513`.  These survivors identify the exact boundary of the divisor-only
compiler; a nondivisor place or a new sidecar is genuinely required there.

## 8. Typed bridge and next work

The honest joint carrier is

```text
(moment multiplication orbit, Hermitian Poisson defect,
 labelled face/current response, effectivity quotient).           (16)
```

- **Preserved:** the exact factorial/HFC moments, real-star structure,
  constant-J predicate, source orientation, and native multiplication.
- **Destroyed by the tempting scalarizations:** mixed-character ancestry,
  edge ownership, target-cell topology, boundary valuations, and uniform
  response length.
- **Cheapest next tests:** solve the degree-at-most-four mixed-character
  Keller/HFC joint ideal; classify anti-tangent polynomial edge pairs on which
  both `kappa` and its curvature jet descend; and compile adaptive divisor
  places for the non-power-of-two residuals rather than extending a fixed
  prime bank by pattern alone.

The main reframing is deliberately modest: the frontiers communicate through
paired native operations and their response modules, not through a common
scalar discriminant.  That produced one all-degree restricted theorem, one
new current/topology passport, one unbounded effectivity obstruction, and a
new exact factorial row without converting any of them into a full-conjecture
claim.
