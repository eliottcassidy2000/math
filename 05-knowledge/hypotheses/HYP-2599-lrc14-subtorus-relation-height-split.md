# HYP-2599: LRC(14) Subtorus Relation-Height Split

**Status:** computationally supported proof program; not a proof of LRC(14).

Script: `04-computation/lrc14_subtorus_relation_lattice_codex.py`  
Output: `05-knowledge/results/lrc14_subtorus_relation_lattice_codex.out`

## Claim

The pure LRC14 S3 density floor should be organized by the affine subtorus
relation lattice

`Lambda_aff(E) = {n in Z^k : sum_i n_i = 0 and sum_i n_i e_i = 0}`,

not by raw spread alone.

For the smooth empty-arc minorant

`G(v)=sum_gaps (gap-2/7)_+`,

the Fourier identity is

`int_0^1 G(E*x) dx = (5/7)^k + sum_{n in Lambda_aff(E), n != 0} prod_i psi_hat(n_i)`.

Thus a proof of the pure floor can plausibly split into:

1. **Low relation height:** primitive affine relations with small coefficient
   product occur.  These are finite relation patterns and should be handled by
   exact enumeration / interval arithmetic.
2. **High relation height:** every nonzero relation vector has large
   coefficient product.  Then the explicit product kernel gives a uniform
   lattice-tail lower bound for `int G(E*x) dx`, hence for `mu(E)`.

The point is that spread is a proxy for this relation height.  A large-spread
shape such as `{0,1,500,501}` has a short support-4 affine relation, while a
large-spread dissociated shape has no short relations and should sit close to
the iid ceiling.

## Evidence

The exact bounded bank scanned all primitive shapes with `0 in E`, `k=5..8`,
and `max(E)<=14`, totaling `9373` exact rows.

The low-`mu` rows are relation-rich:

- `k=7` minimum in the bank is `E=(0,2,3,4,5,6,8)` with
  `mu=13/35`, `25` short support-3 relations and `100` short support-4
  relations at coefficient bound `4`.
- `k=8` minima begin with five-point run shapes such as
  `E=(0,2,3,4,5,6,8,11)`, with `mu=71/220`, `29` short support-3
  relations, and `164` short support-4 relations.
- A relation-sparse comparison `E=(0,1,4,9,13,14)` has
  `mu=25576/28665`, far closer to the iid ceiling, with much smaller
  triple-decay and short-relation counts.

The proof-carrier tournament ranked relation observables by their ability to
predict the negative correction `F(k)-mu(E)`.  The order was transitive:

`triple_decay > small_triples > additive_quad > inverse_spread > run_mass > max_run > spread`.

Fingerprints:

- score histogram `0,1,2,3,4,5,6`;
- directed 3-cycles `0`;
- SCC sizes all `1`;
- Hamiltonian path count `1`.

The important correction is that every triple of points on a line has an
affine relation.  What matters is not the existence of a support-3 relation
but its coefficient product.  This is exactly the product-kernel decay visible
in the `G` Fourier identity.

## Relation to Existing LRC14 Lanes

HYP-2597 and HYP-2598 solve the fixed universal-center skeleton when the
small part is all-odd or 3-free.  HYP-2599 is aimed at the complementary
large-spread recurrence: the non-universal intervals are controlled by the
short vectors of `Lambda_aff(E)`.

HYP-2593 through HYP-2595 then remain the finite-placement side: after a
continuous reservoir is found, the 14-colored CRT grid must hit it with a
bounded discrepancy loss.

## Next Lemmas

1. Define a relation height
   `H(E)=min prod_i max(1,|n_i|)` over primitive nonzero
   `n in Lambda_aff(E)` with support at least `3`.
2. Prove an explicit tail inequality for the `G` kernel when `H(E)>H0`.
3. Enumerate all support-3 and support-4 low-height affine patterns for
   `k<=13` and exact-check the corresponding `G` or `mu` floors.
4. Carry the surviving reservoir through HYP-2593/HYP-2595's colored placement
   machinery.

This would replace the informal spread bound by an intrinsic subtorus
relation-lattice theorem.
