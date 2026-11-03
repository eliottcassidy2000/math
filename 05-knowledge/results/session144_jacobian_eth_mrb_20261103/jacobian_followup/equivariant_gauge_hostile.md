# Arbitrarily high orbit-scale degree is already a primitive gauge artifact

**Status: PROVED, root independently checked the construction. FINITE-EXACT companion uses no assert gates and passes identically normally and under `-O`.** This is a correction to the *formulation* of the remaining classification problem, not a new inequivalent primitive counterexample.

For every integer `k>=2`, the elementary shear

`T_k(x,y,z)=(x, y+x^(2k-1) z^k, z)`

is a polynomial automorphism with determinant one and inverse obtained by subtracting the added monomial. That monomial has weight `(2k-1)-2k=-1`, so `T_k` is equivariant for weights `(1,-1,-2)`. Its quotient is the polynomial automorphism

`Q_k(v,t)=(v+t^k,t)`.

Let `G=(F_3,F_2,F_1)` be the degree-three untwisted witness, whose orbit scale is `l=2-3v-t`. Then `G composed with T_k` has orbit scale

`l composed with Q_k=2-3v-t-3t^k`.

Its total degree is exactly `k`, but its zero divisor is just the graph

`v=(2-t-3t^k)/3`,

isomorphic to an affine line. The new map is equivariantly two-sided equivalent to `G`; it has the same determinant `2`, geometric degree `3`, and primitive function-field extension. Thus **primitive nonlinear orbit scales of every polynomial degree already exist** even without composition of nonunits.

This is the cheapest hostile to using total orbit-scale degree as a moduli invariant. The composition tower offers a genuinely different phenomenon: its orbit-zero divisor has at least two irreducible components from the second iterate onward, so an equivariant polynomial coordinate change cannot turn it into an affine line. Those examples, however, are composite.

The meaningful unresolved question is therefore:

> Are there primitive equivariant Keller maps which are **not equivalent under equivariant polynomial source and target automorphisms** to the known affine-scale seed families?

Raw scale degree does not answer it. Neither does the existence of a primitive map written with nonlinear scale. Component count, the isomorphism type of the scale-zero divisor, geometric monodromy, and compatible boundary data are possible invariant sidecars. Their sufficiency for classification remains unproved.

Run `verify_equivariant_gauge.py` normally and with `-O`. It checks both inverse compositions, determinant, weighted equivariance, quotient formula, scale degree, and the affine-line parametrization for `k=2,...,6`. The arbitrary-k statement follows from the displayed formulas.
