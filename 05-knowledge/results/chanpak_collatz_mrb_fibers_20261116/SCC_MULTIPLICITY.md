# Exact Hamiltonian-count fibers: an inverse-multiplicity analogue
     
Fix a positive odd integer h. Write F_h(z) for the ordinary order generating series of tournament isomorphism classes of Hamiltonian-path count h, including the empty tournament at h=1 (order zero). Write A_h(z) for the strong-class generating series, excluding the empty tournament but allowing the singleton. No all-order enumeration is needed for the following law:
     
`F_h(z) = sum_{k>=0} [sum_{d_1...d_k=h, d_i>1} A_{d_1}(z)...A_{d_k}(z)]/(1-z)^{k+1}.`
     
All nontrivial SCCs of a tournament are linearly ordered by their intercomponent arcs. A Hamiltonian path visits them in that order, contiguously, and arbitrary choices of Hamiltonian paths inside the SCCs concatenate. This gives the product rule for counts and an isomorphism classification by the ordered word of SCC classes. Deleting all singleton SCCs leaves an ordered list of k nontrivial components, with k+1 independent gaps for singleton insertions. This proves the identity without any converse quotient or ordering overcount.
     
A strong tournament on m>=3 has a directed Hamiltonian cycle (Camion). The m breaks of this cycle give distinct Hamiltonian paths, so m<=h. For completeness: a longest directed cycle exists in a nontrivial strong tournament. A vertex outside it with mixed adjacency could be inserted along a cycle edge, so all outsiders point uniformly either into or away from the cycle. If there are any outsiders, strong connectivity guarantees both types and an arc from the type dominated by the cycle to the type dominating the cycle; the pair can then be inserted along any cycle edge, another contradiction.
     
Consequently A_1=z and A_h is a polynomial of degree at most h. The total order s contributed by any nontrivial SCC list is at most sum d_i <= prod d_i = h. Its individual contribution to [z^n]F_h is exactly binom(n-s+k,k) for n>=s and zero otherwise. Thus for n>=h, [z^n]F_h is a polynomial in n of degree equal to the maximum attainable k. There is no cancellation: its leading coefficient is
     
`(1/k!) sum_{d_1...d_k=h, d_i>1} A_{d_1}(1)...A_{d_k}(1)`.
     
For h=1 this gives identically one; for an omitted h it gives zero, without assigning a degree to the zero polynomial. The order h is a crude uniform threshold; an exact useful threshold is the maximum total order of any compatible nontrivial SCC list. The minimum order, by contrast, follows the same product formula in the min-plus semiring, taking the least achievable total order over compatible lists. These are three projections of a single SCC language: minimization (birth), summation (multiplicity), and feasibility (spectrum).
     
Examples. If p is an attained odd prime other than 3, the count-p fiber is exactly eventually linear, with slope the number of strong p-classes. For p*3^j its eventual polynomial has degree j+1 and leading coefficient A_p(1)/j!: the top-degree lists are j directed triangles and a single p-core in one of j+1 positions. The only strong 3-core is the directed triangle, by the order bound. Along 7*3^j, using the exclusions 7,21 and a 63-witness, the first 7-bearing factor must consume at least two factors of 3. Exactly one can carry the 7. For j>=2 the maximum is k=j-1, achieved exactly by one H=63 core and j-2 triangles. The resulting leading coefficient is A_63(1)/(j-2)!. The exceptional chain is not just missing two initial values: it has an exact loss of two polynomial degrees relative to any other attained prime chain with prime other than 3.
     
## Converse involution
     
Let G_h(z) count self-converse tournament classes by order and C_e(z) count self-converse strong classes by order (C_1=z). Then
     
`G_h(z) = sum_{d^2 e=h} F_d(z^2) (1_{e=1}+C_e(z)).`
     
A self-converse class has its SCC word equal to its reversal with every SCC replaced by its converse. The first half is free, the middle is empty or a self-converse strong SCC, and the second half is uniquely determined. This gives an eventual quasipolynomial of modulus at most two and exactly determines the number of two-element converse orbits as (F_h-G_h)/2 coefficientwise. For prime h, each eventual order parity has one fixed class per self-converse strong h-core of that parity. This identity is also unrestricted in h: summing first over h yields the usual ordered-SCC generating series and its converse-palindrome version.
     
## Actual target selection is equivalent to strong selection on explicit atoms
     
The odd integers except 7,21 are a multiplicative monoid. Its unique set of irreducibles is:
     
- all odd primes except 7;
- 7p for odd primes p other than 3,7;
- 49,63,343.
     
To generate an input divisible by 7, pair off 7's if their valuation is even. For valuation odd, absorb one 7 into a 7p if such a prime divisor is available, or into 63 if at least two 3's are available. If neither holds, valuation one gives precisely the two exclusions; otherwise use a single 343 and pair off the rest. The displayed irreducibility checks follow directly from prime factorization and the exclusions. Other inputs just factor into odd primes. Thus after 49,63,343 are realized, Chan–Pak is equivalent to physical strong witnesses for the odd primes other than 7 and physical strong witnesses for all 7p with p not 3 or 7. Any witness at one of these atoms must have just one nontrivial strong component: all other components are singleton padding.
     
This is analogous to the distinction in Carmichael's totient conjecture between an image problem and finite inverse multiplicity; the inverse question here must specify strong cores or a fixed order. At unrestricted order any existing witness has infinitely many singleton-padded realizations, so an unqualified 'no singleton fibers' claim is vacuous. A finite census also cannot establish a universal absence of singleton strong fibers. In this setting the exact slope of an atomic fiber isolates *true strong multiplicity* from padding. The unbounded missing step remains a strong selector on the two prime lanes; none of the SCC, padding or converse identities generates a new atomic target.
     
## Repo lineage to retain
     
The historical affine Carmichael transfer is not a theorem that transfers totient multiplicities to tournament multiplicities. See `05-knowledge/results/carmichael_affine_transfer_20261019/`, the THM-4884/HYP-9477(A) finite-multiplicity corrections recorded in the Mistakes Log, and Session224/225 for the distinction between minimum-order class-universality and existence of a strong first witness. Ordered SCC products and converse palindromes are inherited. The claims here are the joint fixed-count eventual-polynomial formula, its leading coefficient and the explicit two-degree comparison of prime chains; no novelty relative to uninspected history is asserted.

## Which leading degree detects the 7p lane?

Write h=7^r·3^a·b with gcd(b,21)=1, and let Ω count prime factors with multiplicity. Suppose each prime dividing b is attained. If r=0, the degree is Ω(h). For even r≥2, it is Ω(h)−r/2 (49 absorbs exactly two sevens per component, and splitting them would require an extra non-seven prime for each singly carried seven). For odd r≥3, it is Ω(h)−(r+1)/2 (use one 343, then 49 blocks; all lists have at least this loss). A uniform proof of the lower bound on loss is to sum over 7-bearing SCCs: if an SCC contains s factors of 7, its loss Ω(d)−1 is at least s/2, since s=1 requires a companion prime and s≥2 gives s−1≥s/2. Round up the sum.

For r=1 and b>1, the degree equals Ω(h)−1 **if and only if** some prime p|b has an actual strong 7p-core. Indeed that maximum leaves only one loss, so the 7-bearing component must have exactly two prime factors. The partner cannot be 3. If no such core exists and a≥2, the degree is exactly Ω(h)−2 by using one 63 and all other factors separately. If no such core exists and a≤1, a larger 7-bearing composite may still realize h; the degree is then at most Ω(h)−2, and may not exist. This specifies where all counting improvements relying only on the SCC product do and do not see the missing prime lane.

`python check_scc_multiplicity.py` cross-checks the ordered-SCC and converse formulas against all tournament classes through order five using literal permutation path tests; it is an implementation control, not the evidence for the identities.
