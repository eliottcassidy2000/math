# Independent audit of the three-dimensional collision-period quotient

**Status: INDEPENDENT ANALYTIC/SOURCE AUDIT PASS + EXACT REPLAY.**
This audit accepts the theorems and scoped application in
[the primary note](planar_jc48_sep06_collision.md). The full proof and
standalone source were read, the named inherited tangent and seminormal
claims were reopened, and normal and optimized executions were compared
with the frozen transcript. No mathematical correction remains.

## 1. Exact quotient and first-order collision compatibility

Let the nonzero labelled `t_i` span a three-dimensional vector space `V`.
The common-motion map to `N=direct_sum V/(k t_i)` is injective: its kernel
would be the intersection of all the tangent lines, which is zero under
the spanning hypothesis. Tangential reparametrization leaves every
`t_i wedge n_i` unchanged, and a common motion is killed by every scalar
tangent relation. Thus the stated period map and its inclusion of common
motion are well defined.

Independently, regard the zero-period condition as factorization through
the surjection `k^m -> V`, sending the ith basis vector to `t_i`. It says
precisely that the map assigning `t_i wedge n_i` to that basis vector
factors uniquely as `A:V -> wedge^2 V`. A nonzero volume form identifies
this map with the bilinear form `b(u,v)=vol(u wedge A(v))`. Its diagonal
quadratic form vanishes on every `t_i`.

If that quadratic form is zero, polarization makes `b` alternating.
The map `w -> [(u,v) -> vol(u wedge v wedge w)]` is an isomorphism
from `V` onto its alternating bilinear forms. Hence `A(v)=v wedge w`,
and each normal class is the class of the same `w`. Conversely every
quadric vanishing at the `t_i` has a symmetric polarization. Its associated
`A` satisfies `t_i wedge A(t_i)=0`, so `A(t_i)` lies in the two-dimensional
subspace `t_i wedge V`, supplying the required normal class. This proves
both exactness and surjectivity in

```
0 -> V -> ker(period) -> I_2(t) -> 0.
```

The resulting dimension is `6-rank(E_2)`. It counts all two-form slots;
there is no unnoticed choice of a single area form. Repeated or rescaled
tangent directions do not increase the quadratic evaluation rank. At least
six projectively distinct directions are therefore necessary for
completeness, and the coordinate vectors together with their three pairwise
sums give the asserted sufficient six-direction example.

## 2. All-size conic hostile

For `t_i=(1,a_i,a_i^2)`, `n_i=(0,1,2a_i)`, the literal cross product is
`(a_i^2,-2a_i,1)`, the value of one fixed linear map on `t_i`. All scalar
relations consequently annihilate every wedge coordinate. Three distinct
parameters already span `V`, by the Vandermonde determinant.

The independent common-motion contradiction is exact: the first
coordinate of `n_i+v_i t_i=w` gives `v_i=w_1`; the second at two distinct
parameters gives `w_1=0,w_2=1`; the third would then require every `2a_i`
to be the same number. Thus no such motion exists. The quadratic form is
`2XZ-2Y^2`. For five or more distinct parameters, its uniqueness up to
scale follows from the degree-at-most-four substitution polynomial.
This establishes the arbitrary-size statement analytically; the finite
three-, five-, and nine-branch controls are not used to extrapolate it.

## 3. Actual suspension and formal graph criterion

The inherited statements were checked in **THM-4381**,
[exceptional-quartic seminormalization](../../01-canon/theorems/THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification.md),
**THM-4404**,
[descended two-form cokernel](../../01-canon/theorems/THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel.md),
**THM-4411**,
[surface collision transgression](../../01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md),
and **THM-4412**,
[seminormal suspension](../../01-canon/theorems/THM-4412-exceptional-quartic-seminormal-suspension-compensator-firewall.md).
They supply `r(-1)=r(0)=r(1)=0`, the derivative rows
`C'=(3,3,3)`, `E'=(-9,4,9)`, and the nonzero retained derivative period.

The determinant's cofactors are independently `(15,-54,39)`, so it is
exactly `3(5r'(-1)-18r'(0)+13r'(1))`. Therefore the three augmented
tangents are independent, the scalar relation space is zero, and the
first-order quotient has dimension three. The source's three derivative
triples are explicitly synthetic controls; they are not substituted for
the inherited nonzero period of the actual degree-175 generator.

A prose ambiguity was clarified before acceptance. The smooth target
means the ambient Russell surface `C^2E=B(B+4)` times an affine line,
locally with coordinates `(C,E,D)`. It does **not** mean the singular
embedded `Spec(S[D])`. At the retained point the derivative of the surface
equation with respect to `B` is `-4`, justifying that local chart. The
current primary note makes this distinction explicitly.

For the actual graph `(C(x),E(x),r(x)+H(x,s))`, the first two target
coordinates are independent of `s`. Suppose labelled sections have the
same image and all earlier section coefficients are zero. At the next
order their common `(C,E)` coefficient must lie in each of the three
distinct tangent lines `(3,-9),(3,4),(3,9)`, so it vanishes. Since the
tangent vectors are nonzero, every section coefficient vanishes too.
Induction fixes all three sections and their common first two coordinates.
The last coordinate is common exactly when
`H(-1,s)=H(0,s)=H(1,s)`. Monic division, coefficient by coefficient in `s`,
is equivalent to `H=a(s)+x(x^2-1)J(x,s)`.

Thus the all-order iff retains moving labelled sections and common target
motion; it does not assume their constancy at the outset. The polynomial
graph controls `H=sx` and `H=sx(x^2-1)` have the asserted opposite outcomes.
This vertical graph operation differs from THM-4411's deformation of the
planar compiler polynomial, so its criterion does not contradict the
older surface-period or ninth-source results. No descended primitive,
elimination of the extra coordinate, or Keller conclusion is inferred.

## 4. Source and replay

The Fraction implementation constructs the relation-period matrix, the
common/tangent-motion matrix, and the quadratic evaluation matrix by
separate formulas. The quotient rank formula, gauge invariance, invertible
coordinate change, tangent rescaling, synthetic vertical cubes and monic
division controls all match their declared universes. No producer is
imported and all checks remain active under optimization.

```
python3 -B 04-computation/planar_jc48_sep06_collision.py > /tmp/planar_jc48_sep06_collision_audit_normal.out
python3 -B -O 04-computation/planar_jc48_sep06_collision.py > /tmp/planar_jc48_sep06_collision_audit_optimized.out
```

Both replays are byte-identical to
[the frozen output](planar_jc48_sep06_collision.out), with **257 exact
gates**. Freshly recomputed SHA256 values agree with the primary note:

```
source c6363ae79380172c134fcf73b9b8d16b80520a4d6b4c32d65e2046335e4cf5ba
output 0d45a4b8e9a1769584b4e2ab426e2d1829c3514627bfc418ceba30966365d8f2
semantic 487244421195a2ad8153c51d057fa5f54695cb39395dcb6954735f32fbb6fef9
```

The primary source and output were not edited. This audit accepts the
stated first-order linear theorem and actual formal graph criterion, not
a planar Jacobian or global polynomial-termination result.
