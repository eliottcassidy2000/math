# Session 144: equivariant composition and disjoint orbit-collapse divisors

**Status: PROVED by the elementary arguments below; FINITE-EXACT companion replay passed normally and under `-O`. The proof was independently audited twice; see `independent/INDEPENDENT_AUDIT.md` for the second audit and its reproducible certificates.** No complete classification of Keller maps is asserted. Source snapshot: `origin/main = 5e1cedf452c7cc1ec447205939a9d6cf4c0e8743`.

## Inheritance, scope, and the classification question

The closest proved mechanisms are the Keller composition monoid in `01-canon/theorems/THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule.md`, the all-degree weighted lifts and primitive symmetric monodromy in `THM-3438-weighted-lift-keller-degree-spectrum.md`, and the mixed-weight quotient formula in `THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md`.

THM-4727 correctly identifies the seed tangent-sweep family's affine orbit scale. Its literal suggestion that existence of higher-degree orbit scales remains open overlooks composition. This note closes that **existence** question and records a gauge-invariant distinction in the equivariant category. It does not classify primitive solutions, prove new inequivalence under arbitrary changes of coordinates from orbit degrees alone, or provide a new primitive counterexample.

Live concepts: composition, quotient loss, geometric degree, primitive monodromy, and orbit-collapse components. The hostile is the square of the already verified untwisted map, not an untested deformation.

Work over an algebraically closed field of characteristic zero. The exact formulas have rational coefficients. The torus action is

`lambda.(x,y,z)=(lambda x,lambda^(-1)y,lambda^(-2)z)`.

Let `F` be the verified Alpoge map of THM-1300, and reorder its coordinates:

`G=(F_3,F_2,F_1)`.

This has the same source and target weights, determinant `2`, and generic geometric degree `3`. Put

```
v=xy,  t=x^2 z,  l=2-3v-t,
K=t(1+v)^2+v^2(4+3v),
A=v+3K,  B=(1+v)K.
```

Then `G=(xl,A/x,B/x^2)`, with all apparent divisions polynomial after substitution. Its quotient is

`H(v,t)=(l A,l^2 B)`, and `Jac(H)=2l^2`.

The source quotient map forgets the orbit coordinate `x`. It therefore forgets which orbit points survive over the target coordinate hyperplane. The orbit scale, polynomial lifting conditions, and reconstruction on `x l !=0` are necessary sidecars; quotient ramification is not ramification of the lifted Keller map.

## 1. Composition theorem

Define, for integers `j>=0`,

`f_j=l composed with H^j`, and, for `m>=1`, `l_m=product_(0<=j<m) f_j`.

Then:

1. `G^m` is equivariant, polynomial, noninjective, with determinant `2^m` and geometric degree `3^m`.
2. Its orbit scale is exactly `l_m`, and its quotient is `H^m`.
3. `Jac(H^m)=2^m l_m^2`.

**Proof.** If an equivariant map has orbit scale `a` and quotient `P`, and a second one has scale `b` and quotient `Q`, their composite has scale `a(b composed with P)` and quotient `Q composed with P`. Induction gives the formula. The chain rule gives both determinant statements; multiplicativity of finite function-field degrees gives `3^m`. Noninjectivity follows directly by applying `G^(m-1)` to a collision of `G`. This is also a direct witness for every ambient dimension at least three after identity padding. ∎

In particular,

`l_2=l(2-3lA-l^2B)`.

Its total degree is `7`, its value at zero is `4`, and its leading homogeneous part is

`v^3(t+3v)^4`.

The complete exact factor expression is recorded in `verification.json`. It is unnecessary to expand the three-dimensional composite to verify this theorem.

## 2. Exact degree growth for every iterate

For `j>=1`,

`(deg H_1^j,deg H_2^j)=(4*5^(j-1),6*5^(j-1))`,

where `H_i^j` means the ith coordinate of the iterate, not the power of one polynomial. Consequently

`deg(l_m)=(3*5^(m-1)-1)/2`,

beginning `1,7,37,187,937,...`.

**Proof.** The leading homogeneous coordinate terms of `H` are

`-3v^2(t+3v)^2` and `v^3(t+3v)^3`,

so their degrees are `4,6`. Suppose polynomial inputs `V,T` have degrees `a<b`. In `l(V,T)` the unique leading term is `-T`. In `A(V,T)` the unique leading term is `3TV^2`; in `B(V,T)` it is `TV^3`. Therefore the two output leading terms are `-3T^2V^2` and `T^3V^3`, which cannot cancel, and the degree pair is `(2(a+b),3(a+b))`. This recurrence preserves strict inequality and multiplies the sum by five. For `j>=1`, `deg(f_j)=deg(H_2^j)=6*5^(j-1)`, because its degree exceeds that of `H_1^j`; `deg(f_0)=1`. Degrees add under products, giving the formula. ∎

This is a structural induction, not extrapolation from the finite checks at `m<=3`.

## 3. Pairwise disjoint smooth boundary levels

For every `i<j`,

`f_j-2 belongs to the ideal (f_i)`.

In particular `(f_i,f_j)=(1)`. Each curve `V(f_j)` is nonempty and smooth, each polynomial `f_j` is squarefree, and the curves for different `j` are pairwise disjoint. Hence the reduced divisor `V(l_m)` is smooth and has at least `m` distinct irreducible components.

**Proof of disjointness.** Both coordinates of `H` are divisible by `l` and `H(0)=0`. Both coordinates of `H^(i+1)` are therefore divisible by `f_i`. Applying further iterates keeps them divisible, since every coordinate has zero constant term. Evaluating `l=2-3v-t` gives `f_j=2 modulo f_i`.

**Proof of smoothness.** At a point of `V(f_j)`, no earlier `f_i` vanishes. The composition determinant is

`Jac(H^j)=2^j product_(i<j) f_i^2`,

so it is nonzero there. The differential of `f_j=l composed with H^j` is the nonzero row `(-3,-1)` times an invertible matrix. It never vanishes on `V(f_j)`. This proves smoothness and squarefreeness. Each `f_j` is nonconstant by Section 2, hence has a nonempty zero set over an algebraically closed field. Pairwise disjoint nonempty curves supply at least one distinct irreducible component per level. ∎

We do **not** claim that every `f_j` is irreducible; the lower bound is sufficient and proved.

## 4. Complete zero-fiber set of the quotient tower

As sets of geometric points,

`(H^m)^(-1)(0)={0} union V(l_m)`.

The point `0` lies on none of the curves because every `f_j(0)=2`. Every positive-dimensional component is one of the smooth disjoint components above, and every component collapses to the quotient origin. This is a statement about the quotient; `G^m` has no finite ramification and no positive-dimensional fibers.

**Proof.** At `l=0`, both coordinates of `H` vanish. At `l!=0`, they vanish exactly when `A=B=0`. The identity

`3B-(1+v)A=-v(1+v)`

then forces `v=0` or `v=-1`. At `v=-1`, `A=2`, impossible. At `v=0`, `B=t`, so `t=0`. Thus `H^(-1)(0)=V(l) union {0}`. Taking successive inverse images proves the formula. ∎

The result is set-theoretic. It does not assert a reduced scheme structure for the whole fiber or describe its embedded multiplicities.

## 5. What equivalence does and does not preserve

Fix this torus action. An equivariant polynomial automorphism `T` has first coordinate `x a(v,t)` and its inverse has first coordinate `X b(V,T)`. Composing implies

`a(v,t) b(H_T(v,t))=1`.

Thus `a` is a nonzero constant. The induced quotient `H_T` is a polynomial automorphism because the inverse descends too.

If `G` is replaced by `T_left composed with G composed with T_right`, its orbit scale becomes a nonzero constant times `l composed with H_(T_right)`. Consequently the number of irreducible components of its orbit-zero divisor is unchanged. Total polynomial degree need not be unchanged.

Every nonconstant affine orbit scale has exactly one irreducible zero component. For `m>=2`, `l_m` has at least two. Thus `G^m` is **not equivariantly two-sided equivalent to any map with nonconstant affine orbit scale**. This answers the coordinate-artifact objection in the appropriate category.

Under unrestricted polynomial source/target automorphisms, the chosen torus action is not fixed. The orbit divisor need not be an invariant of that broader equivalence. A separate invariant does separate `G^2` from the primitive degree-nine weighted lift of THM-3438: `G^2` has a proper degree-three intermediate function field; the weighted degree-nine map has geometric monodromy `S_9` with maximal point stabilizer `S_8`, hence no proper intermediate field. That is already implicit in THM-3438's atom/composite distinction, not a new classification theorem.

## 6. Boundaries and next question

The seed has affine scale, its composites have unbounded scale degree and unbounded smooth boundary-component count, and quotient collapse is compatible with an everywhere-etale lift. Thus the right unresolved problem is not whether nonlinear orbit scales exist. It is:

- Which nonlinear-scale maps are composition atoms?
- Which are equivalent to affine-scale atoms under a specified equivalence group?
- Can arbitrary equivariant Keller maps be factored into a controlled menu of primitive weighted lifts?
- Which additional boundary and monodromy data distinguish atoms of equal geometric degree?

None is answered by this tower. Neither dimension two nor general non-equivariant Keller maps is classified here.

## Exact companion and audits

Run in this directory:

```
python -B verify_equivariant_composition.py > verification.json
python -B -O verify_equivariant_composition.py > verification.optimized.json
cmp verification.json verification.optimized.json
```

The exact companion directly checks the original 3-by-3 determinant and rational collisions; the quotient identities and determinant; the direct expanded determinant of `H^2`; degree pairs `(4,6)` and `(20,30)`; scale degrees `1,7,37`; comaximality and squarefreeness for `f_0,f_1,f_2`; monomial weight-module controls; the corrected power law; and positive and hostile controls. There are no `assert` gates. The universal statements rely on the proofs above. No formal Isabelle/Lean build was executed in this worker.
