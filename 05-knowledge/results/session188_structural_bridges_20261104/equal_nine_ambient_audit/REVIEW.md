# Independent review of the equal-profile cyclic-child ambient theorem through order nine

**PASS.** This review concerns the frozen `equal_nine_ambient` packet, proof SHA256 `f5c32f6d93743d31083532ba9c3884685b96ef288524004b657ab4e466b0bbed`. Its new conclusion is common child orders eight and nine; together with the earlier certified range it gives orders one through nine. It does not make an inference about larger orders.

## Mathematical review

1. The adjacent-position-pair group acts freely on all permutations, with orbit size `2^floor(n/2)`. A directed Hamilton word must orient every chosen pair correctly; at most one orbit member has even these necessary adjacencies. Thus `H <= n!/2^floor(n/2)`. This requires no assumption about edges joining pairs.
2. The classical odd-cycle-cover profile expansion has nonnegative same-parity lower weights and top weight one. Since each transitive profile has Hamilton coordinate one, the Hamilton bound gives precisely the stated total lower-weight budgets 2,519 and 22,679. Allowing all real nonnegative weights in that simplex enlarges the actual tournament domain and is safe.
3. Homogenization uses one slack coordinate and total simplex weight one. The homogeneous cubic coefficients are equivalently its symmetric trilinear functional applied to simplex vertex profiles, with the ordinary multinomial multiplicity. This independent equivalent computation verifies both the indexing and all diagonal/repeated-coordinate factors.
4. All higher ordered-cover coordinates are nonpositive, and the inherited finite-support sum-of-squares identity makes the first coordinate strictly negative. The degree/support tail vanishes, so there is no unchecked infinite family of coordinates.
5. The ports are physical names with `x -> y`. Their old-core masks are the first and second child sets. At every exterior vertex, both ports receive the original root's incidence. Therefore the port-augmented core is a common module. Its exact Hamilton difference is a nonnegative combination of its ordered-cover differences, and the first coefficient is `H(Q)>0`. This proves strict separation in every fixed ambient, not just in a transitive exterior. Setting all exterior masks to zero would not justify this lift.
6. The retained positive order-ten coefficient concerns a larger artificial continuous simplex. It is a boundary of this sufficient coefficient test, not an actual tournament counterexample or an exact response collision.

## Independent executable evidence

`audit.py` imports no producer implementation. It enumerates complete old-colour run words with both physical ports by an endpoint recursion, symmetrizes the literal source-minus-target counts, and independently recovers `K/(1-XYZ)^2` throughout the exponent box through ten. It enumerates ordinary cyclic run words separately, multiplies their series, and contracts with simplex vertex profiles computed by inclusion-exclusion for ordered surjections.

Every retained coefficient is compared, not just its sign or a sample: 2,373 total coefficients, comprising the earlier 893 and the new 1,480, including all zeros. The latter span 53 higher-cover coordinates. Complete source coefficient lists remain in the producer receipts; the audit records the count, minimum, maximum and a hash of each independently reconstructed complete list. It also independently recovers the reported positive order-ten coefficient `20098555126088950347840000` on that packet's larger domain.

The orbit argument has literal permutation controls for all 1,099 labelled tournaments through order five. Both producer modes are replayed and equal the frozen receipts. Independent audit normal and optimized outputs are byte-identical. All checks use explicit exceptions rather than optimization-sensitive assertions.

## Provenance check

The cited [Grinberg–Stanley detailed paper](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf), Theorem 1.39, supplies the classical odd-cycle-cover expansion. The inspected [Alon paper](https://www.tau.ac.il/~nogaa/PDFS/hamilton.pdf) and [Divoux–Wang manuscript](https://arxiv.org/html/2610.23755v1) provide historical and asymptotic context only. No asymptotic constant from them is substituted into the finite coefficient certificate. The adjacent-swap bound and nonpositive-coefficient method are correctly presented without a novelty claim.

The result remains an equal-complete-profile cyclic-quotient theorem on fixed named child sets and a fixed exterior. It does not assert arbitrary tournament reconstruction, unequal-profile separation, or arbitrary changing exterior contexts.
