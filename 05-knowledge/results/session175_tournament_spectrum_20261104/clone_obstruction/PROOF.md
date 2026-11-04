# Why transitive clone inflation cannot complete the prime lane

Session 175. **Proved universal obstruction to a scalable constructor class.** This does not obstruct arbitrary nonuniform growing interfaces, and does not claim any particular uncovered integer is globally missing. The earlier complete-prefix statement strictly below80557 is historical evidence, not the current repository frontier or nonattainment at80557. [Session134's two-apex provenance](../../session134_constructive_spectrum_frontier_20261101/two_apex/README.md) already credits inherited coverage through99,999 except7 and21; the present clone obstruction does not extend or retract that finite coverage.

## 1. Constructor class

Fix a finite library Q_1,...,Q_K of tournaments of order at most d. For each quotient vertex i replace it by a nonempty transitive tournament of order r_i, orienting every edge between modules according to Q. Allow arbitrary positive integers r_i, reversal, and ordered sums of any number of resulting tournaments. This includes every transitive twin inflation of the 1536 seven-vertex interface states in the sibling `interface` packet, together with arbitrary SCC products.

**Theorem.** Let L=floor(log_2(N−1)), for N≥3. The number of distinct *prime* Hamilton-path counts at most N supplied by this entire constructor class is at most

    K d L^d.                                           (1)

Consequently the proportion of primes at most N supplied by any fixed finite library tends to zero. In particular such a library cannot complete the prime lane, nor therefore prove completeness of all permissible odd counts. The same obstruction applies if the library is first enlarged by any finite number of additional fixed cores.

## 2. Three elementary lemmas

**Insertion monotonicity.** H(T)≤H(T+x). For every Hamilton path of T there is an available position to insert x: an endpoint works, or along the path the relation to x changes from old vertex beating x to x beating the next vertex. Choose the first available position by a fixed rule. Deleting x recovers the old path, so the resulting map is injective. Iteration proves H(T)≥H(T[S]) for every induced subtournament.

**Triangle through every strong vertex.** If Q is strong of order at least 3 and v is a vertex, write O for its out-neighbors and I for its in-neighbors. Both are nonempty. If no triangle v→a→b→v exists, all edges between I and O point I→O. Then no vertex in {v}∪O can reach I, contradicting strongness. Thus each quotient vertex belongs to a directed triangle.

**Inflated triangle count.** If the three triangle modules have orders r,1,1 and the order-r module is transitive, then

    H = 2^r+1.                                         (2)

Write the cyclic orientation X→a→b→X. Paths with consecutive a,b correspond to assigning each vertex of X to the increasing segment before a or the increasing segment after b, giving 2^r choices. The only path without consecutive a,b is b,X,a. Indeed a has only b as an outgoing neighbor and b only a as an incoming neighbor; nonconsecutiveness forces a last and b first. This proves (2) exactly, including r=1.

## 3. Proof of the bound

First suppose the quotient Q is strong of order at least 3. A singleton strong quotient is excluded from this case: its transitive inflation has H=1, independent of its size. For every module i choose the two other quotient vertices in a triangle through i, and retain one vertex of each of their modules. This induces the tournament of (2). By monotonicity,

    H(Q[T_(r_1),...,T_(r_q)]) >= 2^(r_i)+1

for every i. Therefore H≤N implies every r_i≤L. There are at most L^q≤L^d parameter tuples, so at most that many distinct values, regardless of collisions. Notice that this bounds **all** small values of a strong-template family, not only primes.

For a general quotient Q, its strongly connected components are uniquely totally ordered, and the inflated tournament is their ordered sum. A nontrivial quotient component C is strong; its inflation remains strong, as cyclic routes through other modules connect any desired internal vertices. A singleton quotient component becomes a transitive module and contributes H=1. Counts multiply across ordered sums because a Hamilton path can never return to an earlier component.

Each nontrivial component has H≥3 by (2). Therefore a prime total can come only from exactly one nontrivial component; any other components must contribute 1. Its module parameters satisfy the same bound as above. There are at most d components per template and at most K templates. Summing gives (1). This overcounts some templates, which is harmless. Arbitrary further ordered sums cannot introduce a new prime value: a product of positive integers is prime only when one factor is that prime and all others are 1. Reversing all arcs preserves H, so reversal introduces no values. QED.

## 4. Why polylogarithmic coverage misses infinitely many primes

No analytic number theory assumption is needed. For m≥1,

    binom(2m,m) >= 4^m/(2m+1).

In the prime factorization of this binomial coefficient, the p-adic exponent is
sum_j(floor(2m/p^j)−2floor(m/p^j)); each term is 0 or 1, and there are at most floor(log_p(2m)) terms. Thus each prime contributes a prime power at most 2m, giving

    binom(2m,m) <= (2m)^pi(2m).

Combining,

    pi(2m) >= [2m log2 − log(2m+1)] / log(2m).

Hence there are at least a constant times N/logN primes up to N, whereas (1) is O((logN)^d). Their ratio tends to zero. Removing the exceptional prime 7 makes no difference. This proves the claimed limitation with an elementary argument.

## 5. Consequence for the current nonuniform interface

The sibling interface packet exhibits actual strong witnesses 77,91,133 and unit-response components, so its nonuniform edges genuinely escape the mirrored-wrapper obstruction. The present theorem says that merely replacing its finitely many vertices by increasingly long transitive twin blocks cannot turn those successes into universal prime coverage. A scalable positive program must let the number or structure of distinct interface roles grow, or use internally nontransitive modules whose arithmetic cannot be reduced to this finite clone library.

The theorem does not forbid infinitely many prime values in one template, prove any fixed prime unattainable globally, or show that nonuniform growing interfaces fail. It rules out a specific attempted completeness strategy. It also does not assume 80557 is a missing value.

## 6. Checks and attribution

`check.py` independently verifies insertion monotonicity on all labelled tournaments through order 5, the triangle-through-vertex condition on all nontrivial strong labelled tournaments of orders 3 through 5, the exact inflated-triangle formula through r=9, and the parameter bound for all blowups with module sizes 1,2,3 of all nontrivial strong labelled quotients of orders 3 and 4. These finite checks supplement, not replace, the proofs. Run normal and `python3 -O` versions and compare receipts.

The S174 synthesis already uses insertion monotonicity and the special count 2^r+1; SCC multiplicativity and the interface response construction are inherited. The finite-template obstruction and its application to clone-based prime-lane coverage are derived here, with no historical novelty claim.
