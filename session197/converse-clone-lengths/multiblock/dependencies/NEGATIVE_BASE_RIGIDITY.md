# Ternary rigidity of marked self-converse negative-unit clone clocks

## Theorem

Let (T,r) admit an anti-involution fixing r. Replace r by a transitive m-vertex module, retaining its exterior incidences. If, for some real A>0 and real b>1, the resulting Hamilton-path count satisfies

    H_m = A b^m - 1

for all sufficiently large positive integers m, then b=3 and A is an integer of the form2q^2 with q odd, with A>=242. The identity then holds for every positive m. No attainment at242 is asserted.

Thus ternary is the only nonconstant exponential base available to a marked self-converse negative-unit transitive-cloning clock. This is a restriction on that precise mechanism, not on unrestricted Hamilton-path counts or on asymmetric clone sources.

## Proof

Write a_k for the number of words that use each nonroot vertex once and k indistinguishable root symbols, with no adjacent root symbols and with every consecutive pair a tournament arc. A word with k roots has k-1 nonempty internal root-to-root fragments. Permuting those internal fragments while leaving the outer prefix and suffix in place is a valid free action of the symmetric group S_(k-1). It is free because every internal fragment contains at least one labelled exterior vertex, and distinct fragments have disjoint vertex sets. Therefore

    (k-1)! divides a_k.

This argument includes k=1, where the group is trivial. It does not require self-converseness.

The exact transitive-clone expansion is

    H_m = sum_(k>=1) a_k sum_(j=1)^k (-1)^(k-j) binom(k,j) j^m,

for positive m. There are finitely many nonzero a_k. Let K be the largest index with a_K>0. Its exponential coefficient is a_K, and all other exponential bases are strictly smaller integers. If H_m=A b^m-1 eventually, comparison of leading growth forces b=K and A=a_K. Independence of the remaining integer exponential bases forces the equality for every positive m. In particular b and A are positive integers and (b-1)! divides A.

An anti-involution of a tournament has at most one fixed vertex. Since r is fixed, the source has odd order. Extend the anti-involution over the clone chain by reversing its linear order. For every even m the clone tournament has even order and remains self-converse. The audited even-order congruence gives H_m=1 modulo4.

If b is even, take any even m>=2. Then A b^m-1=-1 modulo4, a contradiction.

If b is odd, then b^m=1 modulo4 for even m, so A-1=1 modulo4 and hence A=2 modulo4. But an odd b>=5 has4 dividing (b-1)!, contradicting its divisibility into A. The only remaining b>1 is3.

The audited square theorem for a tau-fixed root with cross-matching number two gives A=2q^2 or4q^2. Since A=2 modulo4, only2q^2 with q odd remains. Its complete order-eleven exclusion, combined with the strong-core bound, gives A>=242. This completes the proof.

## Dependencies and finite controls

The root-run expansion and the maximal-run/matching equality are proved in the frozen `MINIMUM_TERNARY.md`. The marked square theorem and amplitude242 bound are `SYMMETRIC_SQUARE.md` (independently reviewed PASS). The congruence is in the peer's `session197-orbit-wildcard/EVEN_CONVERSE.md` (independently reviewed PASS).

No finite computation is required for the new base exclusion itself. The small standalone `verify_negative_base_rigidity.py` checks the free-fragment divisibility on every labelled rooted tournament through order4 and verifies the arithmetic exclusions over a broad finite control range. Those controls are not used as an all-order proof.
