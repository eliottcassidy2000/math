# Final bounded audit: field completion and tournament rank obstruction

## Singer completion

**PASS** `supersingular/SINGER_CLOCK_FIELD_COMPLETION.md`. Given a linear
Singer permutation sigma and a nonzero identity mark e, every nonzero vector
is sigma^i(e). Multiplication by that vector is exactly the linear operator
sigma^i, proving distributivity. The remaining nonzero multiplicative axioms
are the cyclic group law. Uniqueness correctly requires the prescribed
realization of sigma as multiplication by sigma(e), not merely linearity.

Changing e to sigma^t(e) changes the operation by sigma^-t and produces a
clock-preserving isomorphic, but not identical, named field. Replacing sigma
by a primitive power preserves the product and changes the marked one-step
clock. Neither periodic exponent residues nor a named field encode unbounded
integer winding. Source characteristic-five branch Frobenius is explicitly
not the auxiliary characteristic-two field's squaring map.

`check_singer_completion_independent.py` independently uses polynomial long
division in the six previously audited primitive charts. It verifies all
products and distributive triples, all31 identity phases in each chart,
and the nonlinear-clock failure(1,0,1), whose two outputs are16 and6.
The six source additions are distinct. Normal/-O outputs agree. The prior
incidence audits supply the independently checked map into these polynomial
charts; the new checker does not reimport producer modules.

## Raw odd-entry skew matrices cannot model larger coranks

**PROVED and audited for all orders.** Let M be an integral skew matrix with
every off-diagonal entry odd. Modulo2, M=I+J. If n is even then J²=nJ=0,
so(I+J)²=I and det(M) is odd. Thus M is invertible overQ.
If n is odd, every even principal(n-1)-minor is invertible by that argument,
while skew-symmetry gives det(M)=-det(M)=0 overQ. Therefore

    corank_Q(M)=n modulo2.

The n=1 case is the zero matrix with corank1. This includes every raw
unit-weight tournament dominance matrix, extending the Session96 even-order
observation. It cannot supply alternating matrices of corank>=2 for the
PPVW rank model. Changing the matrix, weights or observable would require
an explicit different construction; a shared21 supplies no such map.

## Same transitive orientation, every positive even rank

**PROVED and audited.** For2r<=n and r>=1 set

    M_ij=(j-i)*sum_(k=0)^(r-1)(ij)^(2k),  i<j,

and extend skew-symmetrically. Let V have columns1,i,...,i^(2r-1), and let
J_r be the direct sum of r blocks[[0,1],[-1,0]]. Each column pair contributes
(ij)^(2k)(j-i), so M=V J_r V^T. Distinct positive row labels make V a
full-column-rank Vandermonde matrix. The map V^T is surjective, J_r is
invertible and V is injective; hence rank(M)=2r.

All upper-triangular entries are strictly positive integers. Thus the same
fixed transitive sign orientation permits every positive even rank up to n,
while forgetting weights loses these rank distinctions. Rank0 is not claimed
when there are nonzero edges. This positive counterpart is an exact additive-
weight information-loss mechanism, not a realization of elliptic curves,
Sha groups or the PPVW probability distribution.
