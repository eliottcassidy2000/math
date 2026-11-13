# Fixed component banks cannot renew determinant repairs or exact-target clocks

**Scope.** These are finiteness theorems for a fixed finite repertoire of strong components under arbitrary ordinal composition. They do not exclude new strong components, interleaving joins, changing exterior clocks, or a full Chan–Pak proof by another route. The two parts extend the inherited fixed-core padding and fixed-amplitude obstructions in [the target-transition packet](../target_collatz_transitions_20261113/README.md).

Let H(A) count Hamilton paths and P_A(k) count ordered k-slot path covers including empty slots. For ordinal join every global slot is uniquely the concatenation of its restrictions, including empty slots, so both H and every P(k) are multiplicative. Fix a finite alphabet of actual nonempty components A, and allow a module to be any nonempty word in that alphabet under ordinal join. This includes arbitrarily many repetitions, rearrangements, mixtures, or transitive padding if the alphabet includes TT1. Let S include the prime divisors of H(A), P_A(2), P_A(3) across this alphabet, and 2,3,5. All such observables are positive S-units, and H is odd, P(2) is even, P(3) is odd and divisible by 3: a transposition of slots 2 and 3 leaves only covers using slot 1 alone fixed, while cycling the three slots fixes no cover of a nonempty module.

We use the [Evertse–Schlickewei–Schmidt theorem](https://annals.math.princeton.edu/2002/155-3/p04): a linear equation in a finite-rank multiplicative group has finitely many solutions with no vanishing nonempty proper subsum on the left. We claim finiteness, not an effective enumeration of the resulting witness set.

## 1. The polarized determinant

For modules M,N, use (h,p)=(H(M),P_M(2)), (j,q)=(H(N),P_N(2)). The inherited polarized repair has output Y=3A+2e, with e=6pj+hj−6hq. At e=±1, apply the theorem to the three positive units pj,hj,hq with signed coefficients. The only possible vanishing proper subsum is 6pj−6hq: the other mixed-sign pair hj−6hq cannot vanish since q≥2j. In that degenerate case hj=e, hence e=1, h=j=1 and p=q. A tournament with one Hamilton path is transitive (by the odd-cycle formula, or injective extension of paths in a directed triangle), so these are precisely the transitive diagonal responses.

In the nondegenerate case there are finitely many term triples. Each fixes hj and hence finitely many h,j; then pj and hq determine p,q. Thus only finitely many nontransitive compensated response pairs occur. In fact each response has finitely many physical module words: P_A(2)≥2 for every letter, so a fixed P(2) bounds word length. This is stronger than fixed-core padding: even a recursive grammar with arbitrarily many alternating SCCs chosen from the finite alphabet cannot keep supplying determinant compensation.

## 2. An exact-target clock refill

In the inherited clock frame write H=60X+h, where X=P_M(3)P_N(2) and h=H(M)H(N). Call v2(P_N(2)) its dyadic exponent (for TT_k this is k). Suppose another presentation within the same chart has H=60Y+j, with positive Y,j formed the same way. The fixed-intercept case h=j forces X=Y and thus equality of dyadic exponents, since each P_M(3) is odd. So holding the product h fixed does not refill the clock even if both modules change without limit. If a proposed presentation increases the exponent to r, then its product satisfies the necessary exact congruence h'≡H mod (45·2^(r+2)), from 3|P_M(3) and 60=15·4.

Now assume both presentations use the fixed finite alphabet and H is prime or 7 times a prime. For h≠j rewrite the collision as

    X/Y + (h/Y)/60 − (j/Y)/60 = 1.

Every ratio is an S-unit. No proper nonempty subsum vanishes: h−j≠0; 60X−j≠0 since 60X is even and j odd; the other pair has two positive terms. There are finitely many such ratio triples. For each triple choose its unique primitive positive integer tuple (X0,h0,Y0,j0), i.e. gcd of its four entries 1. Every integer tuple on that ray is an integral multiple d of this tuple. The target H=d(60X0+h0) with 60X0+h0≥61. For a prime target d=1. For a target 7 times a prime, d can only be 1 or 7 since the other factors left by choosing d equal to the cofactor or the whole target are at most 7. Hence there are finitely many exact prime or 7p clock values admitting a changed-intercept collision from the fixed alphabet. Degenerate collisions cannot increase the exponent, as proved above. This rules out unbounded target-neutral refill within that finite-component clock grammar, not refill after an actual chart change.

## 3. Why determinant outputs cannot carry the old intercept into the clock

Since p,q are even, e≡hj mod12 and Y≡−hj mod12. A clock presentation of the same Y needs its new intercept product h'≡Y mod360, since P_M(3) is divisible by 3 and P_N(2) is even. Thus an exact e=±1 output cannot keep the old product hj when moving into the clock: the sign of that product modulo 12 must reverse. The same necessary condition shows that determinant ±1 cannot be the universal last move: its outputs lie only in classes ±1 mod12, omitting prime classes ±5 mod12 and seven-prime targets 7p with p≡±1 mod12. This is a chart obstruction, not nonattainment of those targets.

The same-sector refill analogy in the separate hypergraph note is now literal: one must return to the right intercept sector before replenishment counts. Finite palettes cannot supply a universal source even with arbitrary reorderings. The proof target becomes a move that constructs new strong components and reestablishes the correct intercept and port guards simultaneously.
