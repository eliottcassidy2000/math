# Prime and seven-prime states do not close under bounded Collatz returns

**Proved.** For every finite depth L≥1, and each c∈{1,7}, infinitely many odd primes q have this property: every positive-depth vertex in the full mixed-sign odd Collatz tree rooted at cq, up to depth L, has at least two distinct prime divisors ≥11. This is true with the *actual*, unrestricted valuations. One can impose any fixed reduced progression on q, and even arrange the property for both roots q and 7q simultaneously. It follows that a bounded reverse-block descent cannot work on only prime and seven-prime states, no matter how the signs or valuations are selected.

## Exact cylinders followed by a hostile sieve

Write D_s(x)=(3x+s)/2^v for odd x>0, s=±1, v=v₂(3x+s). Fix a prescribed reduced progression q≡r mod m, and choose any odd q₀>1 in it (if m is even its reduced residue is necessarily odd). Evaluate both sign trees on the starts q₀ and 7q₀, using their actual valuations. Each path of length k≥1 has an affine expression

    D_{s_k}⋯D_{s_1}(cq)=(c·3^k q+B_w)/2^A,

where A is the sum of the path's valuations, B_w=s₁3^(k−1)+∑_{i>1}sᵢ3^(k−i)2^(a₁+⋯+a_(i−1)) is odd, hence nonzero. Take K larger than every A in the two trees, and at least v₂(m). On q≡q₀ mod 2^K these affine formulas hold on the entire trees: at a prefix of total valuation A, the output differs from the old output by c·3^k(q−q₀)/2^A, an even integer. Equivalently the numerator retains its exact valuation at every prefix. Each positive q in this class has positive actual outputs since both odd maps take positive odd inputs to positive odd outputs.

For each nonempty word and each start multiplier c assign two distinct fresh primes ℓ≥11, avoiding all primes dividing m·c·3·B_w and avoiding those previously used. Set

    q ≡ −B_w (c·3^k)^(-1) mod ℓ.

This residue is nonzero. The odd part of q≡r mod m and the prescribed dyadic class are compatible, and all fresh prime constraints are coprime to their moduli. The Chinese remainder theorem gives a single progression q≡R mod M with gcd(R,M)=1. [Dirichlet's theorem](https://arxiv.org/abs/0808.1408) gives infinitely many primes q. Each endpoint is divisible by its two assigned primes: the numerator vanishes at those primes and the denominator is invertible there.

In fact assign any fixed number J of distinct fresh primes to each map. Then each endpoint has at least J distinct prime factors, all larger than any prescribed bound. This excludes any readiness set requiring bounded prime-factor support, even if fixed multiplicative tags other than 7 are added. It does **not** imply that prime targets lack decreasing Collatz blocks. It says all bounded routes for some prime targets must end at a composite with as many distinct prime factors as specified. Composite-state readiness, an unbounded block-depth mechanism, or arithmetic outside these affine maps is essential. As q grows within the progression the finitely many endpoints all tend to infinity, so no finite extra seed list eliminates the obstruction.

## Connections, and what does not transfer

- [Session 305's exact cylinder compiler](../session305_collatz_core_changes_20261113/CAP4.md) uses exact valuation cylinders to guarantee descent to *some* smaller allowed odd number. This theorem uses the same exact domain, but takes the whole sign tree and makes its endpoint image arithmetically hostile. Restricting the induction to primes and seven-primes was a new restriction not guaranteed by the original compiler.
- The old [Collatz route memory packet](../collatz_atom_route_memory_20261003.md) distinguishes the actual signed affine map from merely knowing that two vertices share a structural shape. Here the previously easy-to-drop affine constant B_w is what allows the congruence to remain in a reduced prime class.
- Older prime-gap and exponent-sieve work exploits small-prime divisibility to organize families; this argument runs the sieve adversarially, reserving a different prime pair for every sign word. There is no appeal to finding simultaneous primes for several forms: only q is required prime. The other forms are deliberately forced composite.

## A direct determinant-to-clock incompatibility

Use the [polarized determinant source](../target_collatz_transitions_20261113/POLARIZED.md) with module counts h,j and two-color responses p,q. For e=6(pj−hq)+hj=±1, even p,q give e≡hj mod12. The repaired target satisfies Y≡−hj≡−e mod12. Thus the determinant ±1 mechanism alone misses prime targets in classes ±5 mod12 and seven-prime targets with prime cofactor in classes ±1 mod12. These missed classes are infinite by Dirichlet.

On the [Π clock chart](../target_collatz_transitions_20261113/collatz/CLOCK_BOUNDARY.md), H=60P_M(3)P_N(2)+H(M)H(N), so H≡H(M)H(N) mod12. Consequently a determinant-repaired target cannot have a Π return preserving the old product hj, even if *both* original modules are replaced by nonisomorphic modules with the same product. A stronger congruence is H≡H(M)H(N) mod360: the two colors can be swapped freely on a nonempty graph, while the three colors can be rotated freely, giving divisibility of P_N(2) by 2 and P_M(3) by 3. The permissible new module product is therefore pinned to Y mod360.

There is a direct integer-size obstruction to literally keeping the two original modules as the two slots, in either order. Nonnegative falling-factorial coefficients give P_M(3)≥(3/2)p and P_N(3)≥(3/2)q, so either Π slot assignment has count at least 90pq+hj. But p≥2h and q≥2j imply Y=36pq+6hq+18pj−hj≤48pq−hj. Ordinal padding in either slot only increases the clock count; ordinally merging the two old modules into one slot, with any nonempty module in the other, also increases beyond Y by multiplicativity and the same inequalities. This excludes specified cross-chart reuses, not all ways to embed old modules into other skeletons. The repair must change more than slot names or transitive padding.
