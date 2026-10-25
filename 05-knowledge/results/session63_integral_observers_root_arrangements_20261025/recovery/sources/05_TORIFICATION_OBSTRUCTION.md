# The birational threshold: five roots can be toric; six already cannot be torified

**Status: PROVED in characteristic zero.** The complete-intersection and
surface-theoretic mechanisms are classical. The uniform coordinate-projection
translation and its explicit hourglass application are new deductions in this
packet, not claims of historical novelty. The full nine-root surface's
256 nodes, general type, p_g=111 and K^2=576 are recovered known results;
see the primary-source note below.

The relative branch obstruction in `RESULT.md` is not by itself a global
torification obstruction. This note supplies the missing birational argument.
It does **not** rule out blue presentations or every framework called F1.

## 1. A local-to-global square-root-cover theorem

Let k be algebraically closed of characteristic zero. Let L_0,...,L_(n-1)
be pairwise nonproportional linear forms on k^3 spanning the dual space,
with n>=4. Assume that at most three of the lines `{L_j=0}` meet at any
point of P^2. Write t_3 for the number of triple points.

The forms embed P^2 as a linear plane Lambda in P^(n-1). Define X_L as the
inverse image of Lambda under the coordinate-squaring map

    [x_0:...:x_(n-1)] -> [x_0^2:...:x_(n-1)^2].

Equivalently X_L is cut out by the n-3 independent linear relations among
the L_j, with each L_j replaced by x_j^2. Let H denote its hyperplane class.

### Theorem T1

1. X_L is a geometrically integral normal surface and a complete intersection
   of n-3 quadrics in P^(n-1). It has degree 2^(n-3).
2. Its map to P^2 is finite locally free of degree 2^(n-1).
3. Its only singularities are `t_3*2^(n-4)` ordinary double points (A1 nodes),
   with exactly 2^(n-4) nodes above every triple point.
4. If pi:Y_L->X_L is the minimal resolution, then

       K_(Y_L)=pi^*((n-6)H),
       K_(Y_L)^2=(n-6)^2*2^(n-3).

   In particular n=6 gives a K3 surface, and every n>=7 gives a surface of
   general type. The geometric genus is

       p_g(Y_L)=[z^(n-6)] (1-z^2)^(n-3)/(1-z)^n,

   with p_g=0 when n<6. Its first values for n=5,...,9 are 0,1,7,31,111.

**Proof of the algebraic assertions.** Choose three independent L_j. Each
remaining squared coordinate is a linear combination of these three squared
coordinates. Successively adjoining the other root variables by those monic
quadratics gives a free module of rank 2^(n-3) over the polynomial ring in
the first three roots. These monic equations are a regular sequence, so the
projective surface is a complete intersection with the stated degree.

The coordinate-squaring map on P^(n-1) is finite locally free of rank
2^(n-1): on the chart where target coordinate j is nonzero its preimage has
x_j!=0, and the other normalized root coordinates satisfy n-1 independent
monic square equations. Pulling back along Lambda preserves that property.

In K=k(P^2), the n-1 square classes `L_j/L_0`, j!=0, are independent. A
nonempty product has odd valuation on any selected line L_j=0, since the
lines are distinct. The generic field obtained by adjoining their square
roots therefore has degree 2^(n-1). The finite flat algebra on every
affine base chart injects into this generic field, hence is a domain. This
proves integrality, and the argument works after algebraic closure.

**Proof of the singularity assertion.** At a base point choose one
nonvanishing L_j and normalize its root to one. Each other root whose
radicand is a unit gives an etale local extension, so it can be separated
into its two local branches when classifying geometric completed local
rings. At an ordinary point or a single line, the root chart is smooth.
At a double crossing, use the two vanishing line ratios as parameters
s,t; adjoining u^2=s and v^2=t still gives the smooth local ring k[[u,v]].

At a triple crossing the three vanishing line ratios satisfy
`L_3/L_j=alpha (L_1/L_j)+beta (L_2/L_j)` with alpha,beta nonzero. The
completed local equation is

    k[[u,v,w]]/(w^2-alpha u^2-beta v^2),

an ordinary double point. There are n-3 nonzero root coordinates at the
point. Their sign choices, modulo simultaneous projective sign, number
`2^(n-3)/2=2^(n-4)`. These give precisely that many distinct nodes over
the triple point. There are no other base-point types by hypothesis.
The displayed smooth or A1 local rings also prove normality.

**Proof of the canonical and genus assertions.** Complete-intersection
adjunction gives

    omega_(X_L)=O_(X_L)(2(n-3)-n)=O_(X_L)(n-6).

Resolving each ordinary double point is crepant. This can be checked
directly: blow up its ambient smooth threefold at the origin. The ambient
canonical divisor acquires 2E, while the strict transform of the quadratic
surface subtracts 2E; adjunction cancels these contributions. Thus the
canonical divisor of the resolution is the stated pullback.

Normality and proper birationality give pi_*O_Y=O_X, so the spaces of
sections of the canonical bundle are equal by the projection formula.
The complete-intersection Koszul resolution gives the homogeneous Hilbert
series `(1-z^2)^(n-3)/(1-z)^n`, and the intermediate cohomology of line
bundles on projective space gives the needed equality between its degree
n-6 component and H^0(X_L,O(n-6)). For n<6 this negative line bundle has
no sections. The degree and pullback intersection formula give K^2.

For n>=7 the canonical divisor is the pullback of a positive multiple of
an ample divisor, so its sections grow quadratically: it is big and the
resolution is of general type. At n=6 the canonical bundle is trivial.
Also H^1(X_L,O)=0 by the same Koszul resolution. A1 singularities do not
change H^1(O) under resolution: their exceptional P^1 has normal bundle
O(-2), so the successive exceptional-neighborhood quotients are O(2m)
with vanishing H^1; the theorem on formal functions gives R^1 pi_*O=0.
Thus H^1(Y_L,O)=0, establishing the usual algebraic definition of a K3
surface. QED.

The no-fourfold-intersection hypothesis is load-bearing. With four or more
vanishing lines, the completed local cover has more than one quadratic
relation and the A1/crepant calculation above does not apply automatically.
Distinct lines and characteristic zero are also retained, not silently
replaced by a characteristic-dependent incidence matroid.

## 2. Apply the theorem to the actual magic labels

In homogeneous parameter coordinates [E:x:z], the nine lines are

    E,
    E+x, E-x, E+z, E-z,
    E+x+z, E-x-z, E+x-z, E-x+z.

Their vanishing lines have exactly eight triple points and no point on four
lines. Four triple points are at infinity E=0, one for each opposite pair.
The other four, in the chart E=1, are `(x,z)=(1,0),(-1,0),(0,1),(0,-1)`.
All remaining intersections are double. This can be checked directly from
the displayed linear forms; `check_torification.py` also independently lists
every projective intersection and every incident-line set over Q.

Every subset of at least four of these lines remains without fourfold
intersections. Every subset of at least four also spans the three-dimensional
dual space, since otherwise all its lines would pass through one projective
point. Thus Theorem T1 applies to every n-root coordinate projection for
4<=n<=9.

| Root labels retained | n | Triple points | Nodes | omega_X | K_Y^2 | p_g(Y) |
|---|---:|---:|---:|---|---:|---:|
| E,E±x,E±z | 5 | 2 | 4 | O(-1) | 4 | 0 |
| E,E±x,E±z and either one of E±(x+z) | 6 | 2 | 8 | O(0) | 0 | 1 |
| E,E±x,E±z,E±(x+z) (one hourglass) | 7 | 3 | 24 | O(1) | 16 | 7 |
| All nine magic labels | 9 | 8 | 256 | O(3) | 576 | 111 |

The six-root row is one specified intermediate projection, not all six-root
projections' node number. Every six-root projection is K3, but its number
of nodes depends on which labels were retained. Similarly every seven-root
projection has p_g=7 and K^2=16, but need not have 24 nodes.

In the hourglass case the three triple points are precisely the intersections
of E=0 with x=0, z=0, and x+z=0. In E=1, an intersection of x=±1 and z=±1
has x+z in {-2,0,2}, not ±1, so there is no finite triple point.

For completeness, the bounded exact census of all labelled subsets is:

| n | Number of subsets by triple-point number t_3 |
|---:|---|
| 5 | 0:28, 1:76, 2:22 |
| 6 | 0:2, 1:20, 2:46, 3:16 |
| 7 | 2:4, 3:18, 4:12, 5:2 |
| 8 | 4:1, 5:4, 6:4 |
| 9 | 8:1 |

These counts are on the binomial(9,n) labelled coordinate subsets, NOT on
isomorphism types or D4-orbits. They serve as FINITE-EXACT controls for the
theorem's different node counts; classification of those isomorphism types
is not asserted.

## 3. The actual finite-torification obstruction

### Corollary T2

For every coordinate projection retaining at least six of the nine root
coordinates, neither the resulting surface nor any of its nonempty Zariski
open subsets admits a finite torification by split tori over C. In particular
this holds for H+, H− and the full nine-root surface M.

**Proof.** A finite torification of an irreducible surface is a finite
disjoint decomposition by locally closed split tori whose field points
cover those of the surface. At least one stratum must have dimension two:
finitely many lower-dimensional locally closed sets cannot cover an
irreducible surface. Such a torus stratum is dense and locally closed,
hence open. It identifies the function field with C(u,v), so the surface
is rational. The same implication holds for any nonempty open subset.

But Theorem T1 gives p_g>=1 for a smooth projective resolution of every
n>=6-root surface. Geometric genus is a birational invariant of smooth
projective surfaces, while a rational surface has p_g=0. This is a
contradiction. QED.

This is a stronger statement than the branch obstruction in `RESULT.md`,
but has a different proof and scope. The branch point blocks the specified
finite toric morphism over the given base torus; here a birational invariant
blocks **every finite split-torus decomposition** of the source surface,
even after deleting the branch divisor or other proper closed subsets.

It follows in particular that these surfaces are not split toric varieties.
Consequently the ordinary normal toric/monoidal F1 route cannot describe
them after complex base extension. This is not a theorem against blueprints,
which explicitly retain additive relations, or against every other use of
the phrase F1 geometry. The blue cone pushout in `RESULT.md` continues to
exist exactly as proved.

## 4. Inheritance, connection map, and literature guard

- Closest inherited mechanism: Session58 proves these are complete
  intersections and computes their dualizing sheaves, but does not classify
  their resolutions. This note supplies the missing local singularity check
  before using canonical positivity as a birational obstruction.
- Hostile near miss: a positive/ample dualizing sheaf on a badly singular
  projective surface alone cannot justify the same conclusion for a
  resolution. Here the A1 computation and explicit discrepancy cancellation
  are the decisive repairs.
- Source -> target: labelled root coordinates -> branched cover of the
  parameter plane -> local A1 models -> crepant resolution.
- Preserved: the root-cover function field and labelled equations.
- Changed: nodes become (-2)-curves; deleting arithmetic boundary divisors
  does not change that function field.
- Sidecar: the entire incidence matroid of the branch lines over the chosen
  characteristic, not just the number of roots.
- Cheapest hostile test: look for four concurrent branch lines; then audit
  the crepant step rather than applying the node formula.

**Primary-source recovery.** Auel and Singer, *The algebraic geometry of
3-by-3 magic squares of squares*, arXiv:2609.09351v1, Sections 1.1--1.3,
Lemma 1.1 and Theorem 2, already give the full surface's 256 ordinary double
points, general type and invariants p_g=111, K^2=576:
https://arxiv.org/html/2609.09351v1 . Those full-surface claims are recovered
controls here, not discoveries. Our use does not import their Picard-rank
computation, rational-point conjectures, or claim to solve distinct positive
integer-square existence. The exact hourglass and all-coordinate-subset
application follows from the self-contained theorem above.

The standard torification definition is that of Lopez Pena and Lorscheid,
*Torified varieties and their geometries over F1*,
https://arxiv.org/abs/0903.2173 . The implication from finite torification
to rationality for an irreducible variety is proved above rather than
assumed from a numerical q=1 heuristic.

## 5. Independent topological check: the node count cancels after resolution

Let t_2 be the number of double points of the line arrangement. Then
`t_2+3t_3=binom(n,2)`. Its complement, open line pieces, double points and
triple points have Euler characteristics respectively

    3-2n+t_2+2t_3,
    2n-2t_2-3t_3,
    t_2,
    t_3.

The **reduced inverse images** give etale covers of these four strata of
degrees `2^(n-1),2^(n-2),2^(n-3),2^(n-4)`, respectively. Scheme-theoretic
fibers over branch strata are not being called reduced or etale here.
Euler characteristic is additive, ignores nilpotents, and multiplies under
the indicated finite covering maps. Replacing each
node by P^1 adds one further unit per node. Consequently

    e(Y_L)=2^(n-4)(24-8n+2t_2+6t_3)
          =2^(n-4)(n^2-9n+24).

For n=5,...,9 this gives 8,24,80,256,768, independent of how many triple
points the particular subset has. It agrees with the independent Noether
identity `e(Y)=12(1+p_g)-K^2`, using H^1(O_Y)=0. In particular the
hourglass resolution has Euler characteristic 80, whereas the singular
hourglass has Euler characteristic `80-24=56`. No finite-field point-count
polynomial is asserted by this complex topological calculation.

## 6. Reproduction

Run `python3 check_torification.py`. It enumerates every intersection of
the nine rational lines, checks the exact subset census and the inherited
5/6/7/9-root rows, verifies all complete-intersection Hilbert coefficients,
and checks local ordinary-double-point ranks at every triple intersection.
The output is `torification_results.json`. These exact bounded checks do
not replace the characteristic-zero proof or assert a reduction theorem at
the primes where the line-incidence matroid changes.
