# Universal boundary-tower corollary in the one-positive equivariant category

**Status: PROVED by the argument below; independent audit passed; see the separately delivered `UNIVERSAL_COROLLARY_AUDIT.md`.** This is an additive generalization of the S144 composition mechanism, not part of the fixed-map exact-degree calculation and not a classification of arbitrary Keller maps.

Let `k` be algebraically closed of characteristic zero, let `w_i>=1`, and let `F:A^n->A^n` be a polynomial Keller map equivariant for the **same** diagonal action with weights `(1,-w_2,...,-w_n)` on source and target. If source and target begin in different orders, invertible linearization forces the weight multisets to agree, and a target permutation puts them in this situation.

Write `x=x_1`, `v_i=x^(w_i)x_i`, and `W=sum_(i>=2) w_i`. The invariant ring is `k[v_2,...,v_n]`. There are polynomials `l,A_2,...,A_n` such that

```
F_1=x l(v),
F_i=x^(-w_i) A_i(v),
H_i(v)=l(v)^(w_i) A_i(v),
```

where the apparent divisions in the displayed `F_i` are interpreted through their polynomial lifts. `H` is the induced quotient endomorphism. If `det JF=c!=0`, then

`Jac(H)=c l^(W-1)`.

## The universal claims

1. `a:=l(0)!=0`, `H(0)=0`, and `H(V(l))={0}` when the divisor is nonempty.
2. `V(l)` is smooth (possibly disconnected), hence `l` is squarefree when nonconstant.
3. Let `f_j=l composed with H^j`, `j>=0`. For `i<j`, `f_j-a` is divisible by `f_i`. Hence the curves/hypersurfaces `V(f_j)` are pairwise disjoint.
4. Every `V(f_j)` is smooth. If `l` is nonconstant, all are nonempty, and the orbit-zero divisor of `F^m`, defined by `product_(j<m) f_j`, has at least `m` distinct irreducible components.
5. The iterated quotient satisfies

`Jac(H^m)=c^m product_(j<m) f_j^(W-1)`.

6. Every equivariant polynomial automorphism has constant orbit scale. Therefore a map with nonconstant `l` is a Keller nonunit, and component counts of the orbit-zero divisor are preserved under equivariant source/target automorphisms.

## Proof

The origin is the only point fixed by the entire torus. Equivariance gives `F(0)=0`, and the derivative at zero is an invertible intertwiner. Since the positive weight `1` has multiplicity one, its corresponding one-by-one block is `l(0)`, so `a!=0`. Each `A_i(0)=0`: evaluating `F_i` on `(x,0,...,0)` gives zero, because a negative-weight polynomial contains no pure nonnegative power of `x`. Thus `H(0)=0`. Its displayed formula shows that every coordinate of `H` is divisible by `l`, proving the collapse assertion.

On the source chart `x!=0`, the change to `(x,v)` is invertible. The first output coordinate is `x l(v)`. At `l(v)=0`, its differential is `x dl`; this cannot vanish because the full Jacobian is invertible. Hence `dl` is nonzero at every point of `V(l)`, proving smoothness.

Every coordinate of `H^(i+1)` is divisible by `f_i`. Further composition preserves divisibility because `H(0)=0`. Therefore `l(H^j)-l(0)` is divisible by `f_i`, proving `(f_i,f_j)=(1)` and disjointness.

The chain rule gives the iterated determinant formula. On `V(f_j)`, no earlier `f_i` vanishes. Thus `JH^j` is invertible there (also when `W=1`), and `d(l composed with H^j)` is nonzero because `dl` is nonzero on `V(l)`. This proves smoothness of every level. `H` is dominant, since its determinant is not identically zero. Thus each pullback of a nonconstant `l` is nonconstant. Over an algebraically closed field it has a nonempty zero set, and pairwise disjoint nonempty divisors supply at least one distinct irreducible component per level.

Finally, if an equivariant map and its polynomial inverse have orbit scales `l` and `r`, their composition has scale `l(r composed with H)=1`. Polynomial units are constants, so `l` is constant. Such an automorphism induces a quotient automorphism; left/right composition changes a general scale only by a nonzero constant and pullback by that quotient automorphism. This preserves its component count. ∎

## Scope

No exact degree-growth law is asserted here. The closed formula `deg(l_m)=(3*5^(m-1)-1)/2` belongs to the particular Alpoge tower proved separately. Likewise, the exact zero-fiber formula `{0} union V(l_m)` requires the extra seed assertion `H^(-1)(0)={0} union V(l)`, which was checked for that seed but is not assumed universally.

These are equivariant-category statements. Arbitrary polynomial coordinate changes can change the chosen torus action, so orbit-zero component counts are not asserted to be invariants of unrestricted Aut×Aut equivalence. No global topological covering or surjectivity claim is made: Keller maps here are etale polynomial endomorphisms, and finite unramified coverings arise only after restricting the target to a suitable open set away from nonproper values.

## Positive-coordinate fibers and one non-slice direction

Assume `l` is nonconstant. Then the positive-weight component has the exact geometric fibers

```
F_1^(-1)(0) ~= A^(n-1) disjoint-union (G_m times V(l)),
F_1^(-1)(c) ~= A^(n-1) minus V(l)       for every c!=0.
```

Indeed `F_1=x l(v)`. At `x=0`, all invariants vanish and `l(v)=l(0)!=0`, so the full coordinate hyperplane `A^(n-1)` is one zero-fiber component and meets none of `l=0`. On `x!=0`, the orbit coordinates are an isomorphism to `G_m times A^(n-1)`, giving the second zero-fiber part. For `c!=0`, the equation determines `x=c/l(v)` uniquely on the complement of `V(l)`, with reconstruction `x_i=v_i(l(v)/c)^(w_i)`. These maps are regular on the displayed open varieties and inverse to the invariant-coordinate map.

Consequently `F_1` is neither a coordinate nor a stable coordinate: its zero fiber is reducible, and remains reducible after a product with any affine space, unlike an affine-space coordinate fiber. Its Keller dual derivation `D_1`, characterized by `D_1 F_i=delta_(1i)`, is not locally nilpotent either. If it were, `F_1` would be a slice and the slice theorem would identify all fibers with `Spec(ker D_1)`. But the zero fiber is reducible and every nonzero fiber is a nonempty irreducible open subset of affine space, a contradiction.

This is a claim about the **positive-weight component and its dual derivation only**, not all coordinate components or arbitrary linear combinations. It uses no finite-field point count. The fixed seed's three-direction non-slice statement is a stronger, separate result of THM-4562. Root independently audited this corollary; the separate exact reconstruction companion gives finite controls.

## Finite exact higher-dimensional controls

The separate `verify_padded_boundary.py` tests the untwisted three-variable seed times identity in four and five dimensions, with seven choices of extra negative weights. Every direct quotient determinant equals `2 l^(W-1)`. Normal and `-O` output files are byte-identical. These check signs and exponents on explicit witnesses; the all-map corollary is proved above, not inferred from these examples.
