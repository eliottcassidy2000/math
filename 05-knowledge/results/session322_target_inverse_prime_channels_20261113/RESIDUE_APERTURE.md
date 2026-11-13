# Every allowed residue in both selective-ray quotients

**Theorem.** For either of the two marked selective quotients from [Session 321](../session321_selective_returns_and_fano_wiring_20261113/SELECTIVE_RAYS.md), any modulus N and any residue a modulo N compatible with oddness, there are literal nonempty tournaments U,V whose strong inflation has count congruent to a modulo N. If N is odd every residue is allowed. This is simultaneous over all prime powers of N, not a claim about equality to the prescribed absolute integer.

Use the inherited bivariate layout formula. In terms of the ordered-slot polynomial Q, the two quotients evaluate to

    H_out = A(U) Q_V(2) - H(U)H(V),
    A(U) = -3H(U)+6Q_U(2)+c Q_U(3),   c=2 or 6.

The idea is to annihilate the positive term at all requested primes, using V at primes other than 3, and the independent U coefficient at 3. The negative term then performs the exact CRT readout.

## One-apex pattern formula

Let P_d(n) be the increasing chain of dn vertices and one apex, whose out-cut on that chain is the repeated pattern (10^(d-1))^n. Use d=2 or 3. It is strong for n≥1 (the apex beats the first chain vertex and loses to the last). Set

    Phi_n(r)=sum_{k=0}^{n-1}(n-k) r^k
            =n/(1-r)-r(1-r^n)/(1-r)^2  (r≠1).

At r=1 the finite sum is n(n+1)/2. Then

    H(P_d(n))=1+n+(2^(d-1)-1) Phi_n(2^d).

Also, for odd moduli, its normalized two-slot response is

    Q_{P_d(n)}(2)/2^(dn) = 2 + A_d Phi_n((3/2)^d) + B_d Phi_n(2^(-d)),
    A_d=sum_{j=1}^{d-1} 3^(j-1)/2^j,     B_d=sum_{j=1}^{d-1} 1/2^j.

To see both identities, any directed cycle uses the apex. It goes from the apex to a selected 1, through an increasing subset, and back from a selected 0. Their distance is dk+j, with 1≤j<d and n−k choices of the pair. Odd total cycle order selects an even number of intervening vertices. The odd-cycle packing formula (at most one cycle here) gives, at general z, the base term z^(dn+1) and the contribution

    (n-k) z^(dn-dk-j) [(z+1)^(dk+j-1)+(z-1)^(dk+j-1)]

for each pair of distances. Evaluating at 1 and 2 gives the formulas. The replay uses a division-free matrix evaluation for arbitrary rings and checks it by literal subset DP at small orders.

## V: annihilate every non-three prime power while keeping its count a unit

Write N=2^k 3^a M with M coprime to 6. For a prime power p^e dividing M, choose a pattern P_d(n), d=2 for p≥7 and d=3 for p=5. Initially require n a multiple of p−1. The exponential ratios in its normalized two-slot expression are then 1 to the n-th power modulo p, and both denominators 1−r are units. For d=2 the normalized response reduces modulo p to 2+4n/15, and H reduces to 1+2n/3. For d=3,p=5 the normalized response is 2−2n, and H is 1+4n/7. Thus the two-slot roots give H congruent respectively to −4 modulo p and 3 modulo 5, both units.

Lift the two-slot root to p^e along n in a fixed progression of step D divisible by p−1 and coprime to p. Indeed, shifting n by j D p^t changes the linear term of the normalized response by a fixed unit times j p^t modulo p^(t+1), while the exponential terms do not change there: each base r has r^D≡1 modulo p. This gives a unique digit at each level.

Make this component's H a unit at every other odd prime of N as follows. For 3 and the primes q<p in M, require at the outset that n be a multiple of q and q−1. These constraints only enlarge D by factors coprime to p. The formula for H gives H≡1 modulo q, including when 2^d≡1 (then the triangular sum n(n+1)/2 vanishes). For q>p, process in increasing order after the p-lift. Enlarge the step by q−1 to fix the exponential phase at q, without changing the previously chosen representative. The enlarged step is coprime to q, since all old step factors are smaller primes. Along the progression H modulo q is affine in n with nonzero slope 2^(d−1)/(2^d−1) if 2^d is not 1; in the exceptional case it is a nonconstant quadratic (only d=3,q=7 occurs here). Choose n modulo q outside its roots and retain this constraint. These are finite arithmetic selections, not a random search for tournaments.

Ordinally join the resulting components over the primes of M, then add a transitive singleton tail of length k (or a singleton if this would otherwise be empty). Slot polynomials and Hamilton counts multiply under ordinal joins. This constructs a nonempty V with Q_V(2) divisible by 2^k M and H(V) coprime to N. The later outer inflation is strong even though V itself may not be.

## U: prescribe its count residue and kill the coefficient at 3

First prescribe H(U) modulo N to be −a_out/H(V), where a_out is the requested output residue. Ordinally join b directed triangles to the construction below, where b is the requested 3-valuation truncated at a, so these triangles supply 3^b. If a>b, the remaining requested unit at 3 can be made 1 modulo 3: when necessary multiply by a fixed physical count 2^L+1 with L even and 2^L≡1 modulo M. A cyclic wrap of TT_L supplies this factor, which is a unit at every prime of N and is 2 modulo 3. All these fixed factors can be divided out of the relevant congruences, after reducing the modulus by 3^b.

For the remaining source use P_2(n), with exact count

    h(n) = (4^(n+1)+6n+5)/9.

At 2^k, take n large and choose n modulo 2^(k−1); the exponential term vanishes and the linear term runs through exactly the odd residues. At 3, choose n≡0 modulo 3 when a>b; that phase has h≡1. Its higher digits can be chosen freely: for t≥1,

    h(n+j·3^t)−h(n) ≡ j(n+1)3^t (mod 3^(t+1)).

Indeed (4^(3^t)−1)/3^(t+1)≡7 modulo 9, starting with t=1 and preserved by cubing. Substitution in h and expansion gives the displayed congruence. The derivative n+1 is a unit on the chosen phase, so this lifts along any existing step coprime to 3. At each p≥5 in increasing order, enlarge the current step by p−1 (coprime to p and compatible with all earlier phases), and lift along that progression. The exponential term is fixed modulo p at the first digit, and fixed modulo p^(t+1) under shifts by the step times p^t. The linear term 6n/9 has a unit slope there. This solves all count congruences simultaneously. Increasing n by the final period makes it positive and as large as needed for the 2-primary argument.

Now append a further transitive tail TT_J to U. It preserves H(U) and scales Q_U(2),Q_U(3) by 2^J,3^J. The existing triangle factors ensure both Q_U(2) and H(U) have exact 3-valuation b if b<a: this uses Q(2)≡±H modulo 3 for the remaining unit-count tournament, by the parity of the slot polynomial. Therefore, when a>b+1, choose J to solve

    2^(J+1) [Q_U(2)/3^b] ≡ H(U)/3^b (mod 3^(a-b-1)).

Powers of 2 exhaust the units modulo every power of 3 (LTE gives order 2·3^(t−1)), so this is possible with J arbitrarily large in its congruence class. Choose J large enough also to annihilate c·3^J Q_U(3) modulo 3^a. If a≤b+1 the first two terms are already divisible by 3^a, so only the latter condition remains. Consequently A(U) is zero modulo 3^a for either value of c.

The resulting literal inflated tournament is strong because its fixed outer quotient is strong. Its first term vanishes modulo N, and its second gives exactly the requested residue. This settles the entire congruence feasibility of the two selective geometries; an exact ending theorem still needs absolute size, a target-decreasing return, or other constructive control not implied by these congruences.

## The same family is extremely sparse in absolute size

The construction above uses U of the form P_2(n) with an ordinal join of b triangles, at most one cyclic wrap of a transitive run, and a transitive tail. It uses V as an ordinal join of any multiset of the patterns P_2(n), P_3(n), and transitive padding. Even if all parameters are made completely free, this family's number of distinct absolute output counts up to X is at most

    exp(O(sqrt(log X))).

Indeed Q_U(2)≥2H(U), and Q_V(2)≥2H(V) by allowing either slot to be empty. The negative term is therefore more than absorbed by the 6Q_U(2)−3H(U) contribution. The odd-cycle packing term with no selected cycles gives Q_U(3)≥3^|U| and Q_V(2)≥2^|V|. Thus H_out≥2·3^|U|·2^|V|, forcing total input order O(log X).

At bounded order the number of U parameter tuples is polynomial in the order. V has at most two possible pattern types of each part size, plus padding; their order in the ordinal join does not change either H or Q(2), the only coordinates read by the output. The distinct V profiles are bounded by three-colored integer partitions. For size at most L, their count is exp(O(sqrt L)): for x=exp(−1/sqrt L), the logarithm of the product generating function is bounded by a constant times

    sum_{j≥1} [j(exp(j/sqrt L)−1)]^(-1) = O(sqrt L),

and extracting terms up to L costs at most x^(−L). Polynomial factors are absorbed by the same bound. This explicit family is therefore a physical example of complete simultaneous congruence coverage coexisting with subpolynomial absolute coverage. Improving only the local moduli of this construction cannot complete the Chan–Pak inverse.
