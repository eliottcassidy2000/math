# Quadratic sidecars for collision periods after seminormal suspension

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The theorem below identifies the exact information lost by scalar
tangent-relation periods in a three-dimensional smooth target. It gives an
arbitrarily large explicit family of period-invisible collision failures
and an exact graph criterion for the exceptional seminormal suspension.
It does not concern a planar Keller map or settle JC(2).

## 1. Inheritance and the precise question

The closest mechanism is
[THM-4411 / first-order-collision-transgression-seminormal-tradeoff](../../01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md):
in a surface, all scalar tangent-relation periods detect first-order
collision compatibility. Its hostile `h=x` pays the missing period by
splitting the collision. The
[ninth-source theorem](overnight8_20260906_jc_collision_period.md) retains
that obstruction for the specified collision-preserving late source jets.
Those statements are inherited, not reproved here as new stability results.

[THM-4381 / exceptional-quartic-seminormalization-and-conductor-fibre-classification](../../01-canon/theorems/THM-4381-exceptional-quartic-seminormalization-and-conductor-fibre-classification.md)
and
[THM-4412 / exceptional-quartic-seminormal-suspension-compensator-firewall](../../01-canon/theorems/THM-4412-exceptional-quartic-seminormal-suspension-compensator-firewall.md)
give `S^sn=S[r]`, where `r` vanishes at every conductor branch and its
derivative supplies the one missing triple jet. The added coordinate makes
the three retained tangent vectors independent. The corrected near miss is
to read the disappearance of their old linear relation as completeness of
the same period test after the target dimension changes. The least-used
sidecar is the quadratic evaluation map on the tangent directions.

The live board has five entries: surface-period completeness; the augmented
triple's independent tangents; all two-form slots versus common motion;
quadrics through projective tangent directions; and actual source graphs
in the embedded suspension. The operation is stabilization followed by a
fresh kernel computation, as required by the repair-quotient and controlled
forgetting cards in `00-navigation/META-PATTERNS.md`. Targeted searches of the
named collision/transgression canon and JC result files found no existing
version of the exact quadratic quotient below. This is a repository
inheritance check, not a claim of literature priority.

## 2. Exact first-order quotient in a three-dimensional target

Let `k` have characteristic zero, let `V` be three-dimensional, and let
nonzero labelled vectors `t_1,...,t_m` span `V`. Set

```
N = direct_sum_i V/(k*t_i),
R = {ell in k^m : sum_i ell_i*t_i=0}.
```

A tuple `n in N` represents branch velocities modulo tangential
reparametrization. It preserves the collision to first order precisely
when it lies in the diagonal common-motion image

```
delta:V -> N,                  w |-> ([w],...,[w]).
```

This map is injective. Define the full scalar-relation wedge-period map

```
P:N -> Hom(R,wedge^2 V),
P(n)(ell)=sum_i ell_i*(t_i wedge n_i).                 (1)
```

Testing (1) against every target two-form is equivalent to `P(n)=0`.
It retains all two-form slots, not merely one selected area form.
Tangential choices do not affect (1), and `P delta=0`.

Fix a nonzero volume form on `V`, and write `I_2(t)` for the vector space
of homogeneous quadratic forms vanishing at every `t_i`. Then

```
0 -> V --delta--> ker P -> I_2(t) -> 0                (2)
```

is exact. In particular, if `E_2` has rows

```
(a_i^2,b_i^2,c_i^2,a_i*b_i,a_i*c_i,b_i*c_i),
                    t_i=(a_i,b_i,c_i),
```

then

```
dim(ker P / delta V)=6-rank(E_2).                      (3)
```

Consequently all these scalar periods are complete for first-order
collision compatibility **if and only if** no nonzero quadric vanishes
on all projective tangent directions. At least six distinct directions
are necessary. Six suffice: take the three coordinate vectors and their
three pairwise sums. Six points lying on a conic need not suffice, and
arbitrarily many such points need not suffice either.

**Proof.** For `P(n)=0`, the assignments
`t_i |-> t_i wedge n_i` respect every linear relation. Since the `t_i`
span `V`, they extend uniquely to a linear map
`A:V -> wedge^2 V`. The volume form gives a bilinear form

```
b(u,v)=vol(u wedge A(v)).
```

Its diagonal quadratic form `q(v)=b(v,v)` vanishes on every `t_i`.
If `q=0`, polarization makes `b` alternating. In dimension three every
alternating bilinear form is uniquely `b(u,v)=vol(u wedge v wedge w)`
for some `w`. Hence `A(v)=v wedge w`, and each `n_i-w` is tangent to
`t_i`: this is exactly the common-motion image.

Conversely, take any `q in I_2(t)` and polarize it to the symmetric
bilinear form with diagonal `q`. The volume identification gives a linear
`A`. Since `vol(t_i wedge A(t_i))=q(t_i)=0`, one has
`A(t_i) in t_i wedge V`. Choose `n_i` with
`A(t_i)=t_i wedge n_i`. These choices produce an element of `ker P`
mapping to `q`, proving surjectivity and (2). A different volume form
rescales the quotient identification by a nonzero scalar. Completeness
and the dimension formula are independent of that choice. This proves
the theorem.

The two-dimensional theorem is consistent with this calculation: there,
every linear map `V -> wedge^2 V` is already `v |-> v wedge w`.
There is no extra symmetric coordinate. Passing to dimension three is
the first place that this argument changes.

## 3. An arbitrary-size hostile with all periods retained

Choose any `m>=3` distinct scalars `a_i` and put

```
t_i=(1,a_i,a_i^2),               n_i=(0,1,2a_i).       (4)
```

The tangents span `V`. In the usual cross-product coordinates for
`wedge^2 V`,

```
t_i cross n_i=(a_i^2,-2a_i,1)=M*t_i,
M=[[0,0,1],[0,-2,0],[1,0,0]].                          (5)
```

Thus every scalar relation among the `t_i` kills their full wedge
responses. Nevertheless no common `w` and scalars `v_i` can solve
`n_i+v_i t_i=w`. The first coordinate makes every `v_i=w_1`;
the second, at two distinct `a_i`, forces `w_1=0` and `w_2=1`;
the third would then require `w_3=2a_i` for all `i`, impossible.
This is an exact failure of sufficiency, not a failure to test enough
relations or two-forms. The missing quadratic form is

```
q(X,Y,Z)=2XZ-2Y^2.
```

For at least five distinct parameters this is the whole one-dimensional
miss: substituting `(1,a,a^2)` into a general quadric gives a polynomial
of degree at most four, so five zeros force it to be a multiple of
`XZ-Y^2`. Three independent directions already give the smallest possible
spanning hostile, but the all-size conic construction shows that branch
count alone cannot repair the test.

These are actual first-order map germs, for example the disjoint labelled
branches `f_i(u,s)=u*t_i+s*n_i` into affine three-space. They are not
asserted to be branches of a polynomial Keller map.

## 4. The actual seminormal suspension and an earlier source graph

Use the exact exceptional curve and field of THM-4381. The ambient target
is the smooth Russell surface `C^2 E=B(B+4)` times an affine line. Near
the retained point it has smooth coordinates `(C,E,D)`, with `D=z+r(x)`.
This ambient target is distinct from the singular embedded suspension
`Spec(S[D])`. On the source section `z=0`, the tangent columns at
`x=-1,0,1` are

```
T_-=(3,-9,r'(-1)), T_0=(3,4,r'(0)), T_+=(3,9,r'(1)).
det[T_- T_0 T_+]=3*(5r'(-1)-18r'(0)+13r'(1)) !=0.     (6)
```

The nonzero period in (6) is inherited from
[THM-4404 / exceptional-quartic-descended-two-form-seminormal-cokernel](../../01-canon/theorems/THM-4404-exceptional-quartic-descended-two-form-seminormal-cokernel.md).
There are no nonzero scalar tangent relations. Equation (3) therefore
gives **three** undetected first-order collision coordinates for arbitrary
ambient normal velocities. Adding the seminormal jet removed the old
scalar relation, but did not make common-motion compatibility automatic.

There is also a complete actual graph criterion, without restricting to
first order. Let `H(x,s) in s K[x][[s]]` and take the source graph

```
z=H(x,s),
Gamma_s(x)=(C(x),E(x),r(x)+H(x,s)).                    (7)
```

Then the three labelled branches have one common formal target image
if and only if

```
H(-1,s)=H(0,s)=H(1,s),
equivalently H(x,s)=a(s)+x(x^2-1)J(x,s).               (8)
```

For necessity, the first two coordinates in (7) do not vary with `s`.
At the first nonzero order of any proposed section or common-image
motion, its coefficient must belong to all three distinct planar tangent
lines `(3,-9),(3,4),(3,9)`. Their intersection is zero. Induction fixes
every section at its base point and fixes the common `(C,E)` image.
Since `r` vanishes at the retained points, equality of the last coordinate
is exactly (8). Sufficiency uses these fixed sections. Monic division
by `x(x^2-1)` proves the equivalent polynomial remainder statement
coefficient by coefficient. For polynomial `H`, both `a,J` are polynomial
in `s` as well.

In particular the lawful polynomial graph `z=s*x` splits the triple at
first order while all scalar tangent-relation periods are vacuous. The
positive graphs `z=s` and `z=s*x(x^2-1)` retain the collision exactly.
The latter is a nonconstant graph, so constancy as a polynomial in `x`
is not necessary. These graph changes are realized by the polynomial
source shear `z |-> z+H(x,s)` when `H` is polynomial; this does not make
their coefficient functions descend through `S`.

## 5. Connection contract, loss, and stopping boundary

The source is the normal-velocity space of one labelled collision in a
smooth three-dimensional target. The target is the scalar-relation
two-form period packet (1), and the map is differentiation followed by
wedge contraction and summation. It preserves every necessary scalar
period identity and forgets exactly the quadratic quotient (2), after
tangent and common-motion gauges. The needed sidecar is the symmetric
part of `A`, or equivalently its six quadratic coefficients modulo the
evaluation constraints. The cheapest decisive test is the three-branch
instance of (4); the six-direction nonconic control proves sharpness of
the repair.

For the exceptional suspension the map from the old surface into its
augmented tangent space preserves the point collision and adds the missing
seminormal jet, but changes the dimension and the tangent-relation kernel.
The all-order graph criterion (8) retains the source section explicitly.
It is not the moving graph/edge-primitive complex of THM-4067, a new
descended differential on the planar core, or an elimination of `z`.
It neither contradicts the fixed ninth-source theorem nor excludes the
unrestricted earlier source-jet space. No Keller or polynomial-termination
conclusion follows.

The new question is whether a proposed stabilization or earlier-jet
operation carries these quadratic collision coordinates together with
its descended differential. Repeating the scalar period test alone cannot
answer that question, even with arbitrarily many branches.

## 6. Verification manifest

The [standalone standard-library source](../../04-computation/planar_jc48_sep06_collision.py)
and [frozen output](planar_jc48_sep06_collision.out) pass **257 always-active
exact gates**, with byte-identical normal and optimized output. The source
imports no inherited producer. It constructs the collision and period
matrices separately from the quadratic evaluation matrix, then checks
the quotient dimension identity on eight named configurations. Their
unseen dimensions are `3,2,1,0,0,0,1,1`, in the declared order:

```
three independent; four general; five coordinate controls;
six complete; seven complete; eight with repeated directions;
five conic; nine conic.
```

The three conic hostiles use exactly `3,5,9` parameter values. Each is
checked before and after tangent/common-motion gauge changes, an
invertible target coordinate change, and independent tangent rescalings.
The suspension determinant is checked through its symbolic cofactors.
Three explicitly synthetic derivative triples, each with nonzero retained
period, supply all `3*27=81` retained-value cube controls for the vertical
graph criterion. They are not presented as numerical evaluations of the
actual degree-175 seminormal generator. The eight named graph polynomials
check the monic-divisor criterion, including both the actual `h=x` hostile
and nonconstant `h=x(x^2-1)` positive control. The exceptional application
uses the proved nonzero period of the actual generator, not a synthetic
substitution for that dependency.

Reproduce from the repository root:

```bash
python3 -B 04-computation/planar_jc48_sep06_collision.py > /tmp/planar_jc48_collision_normal.out
python3 -B -O 04-computation/planar_jc48_sep06_collision.py > /tmp/planar_jc48_collision_optimized.out
cmp /tmp/planar_jc48_collision_normal.out /tmp/planar_jc48_collision_optimized.out
cmp /tmp/planar_jc48_collision_normal.out 05-knowledge/results/planar_jc48_sep06_collision.out
```

```
source SHA256
c6363ae79380172c134fcf73b9b8d16b80520a4d6b4c32d65e2046335e4cf5ba
output SHA256
0d45a4b8e9a1769584b4e2ab426e2d1829c3514627bfc418ceba30966365d8f2
semantic SHA256
487244421195a2ad8153c51d057fa5f54695cb39395dcb6954735f32fbb6fef9
```

The source and output are frozen. Universal claims rest on the displayed
proofs, not this finite control bank. The
[independent analytic/source audit](planar_jc48_sep06_collision_audit.md)
passes, including normal/optimized/frozen replays and the distinction
between the smooth ambient target and the singular embedded suspension.
