# Independent audit: multiblock symmetry completion and orbit series

**PASS, no correction.** The finite union of weight equalities, exact mutable-block completion minimum, ordinary and twisted cycle indices, labelled mass and asymptotic centralizer factors all hold at the stated fixed canonical quotient scope. No Hamilton-count-preserving or unrestricted spectrum conclusion follows.

## Frozen source

- PROOF.md: `934f76db619e39187d26c83f8c944694ab47371971077d32f3cdf81fdc8a2be5`.
- verify.cpp: `82bb34b479fe21a16519a566cee9ab97f06f83c79cb2f6f0d85ddb62f4841936`.
- Supplied optimized/debug receipts and both copied source replays: `b2a796bd10f56e8c0c17553d2b9c22ea7d3fcd774b406a8e9d241c8d69219969`.

The source directory is session197-orbit-wildcard/multiblock-orbits/. The mathematical dependencies are the previously reviewed canonical-block proof c3516df7… and converse-mass proof07f4b20d…. The copied source is standalone. Original files were not changed.

## Canonical language and physical completion

The maximal transitive-module partition is canonical. For a pair-free quotient, a larger transitive module in a positive blow-up would, by repeated union with the blocks it meets, induce a larger transitive module in the quotient. Therefore the original chains remain exactly the canonical blocks at every positive weight vector. This argument does not require strongness or full module primeness.

Any isomorphism must permute the canonical blocks and lift an automorphism of the quotient; order inside each transitive block fixes its lift uniquely. An anti-isomorphism reverses each block order and lifts a quotient anti-automorphism. Thus exact map sets, not just their sizes, are controlled by weight-preserving quotient maps. Repeated or simultaneous vertex cloning changes only the indicated block sums, proving the affine length equations at all lengths.

A weight-preserving anti-automorphism has order twice an odd number. Its appropriate odd power is a weight-preserving anti-involution. Thus all successful descendants can be tested against the finite anti-involution set, and the fixed-block/pair constraints suffice. The singleton quotient uses the identity and is treated separately.

For each involution pair, any descendant must give the two blocks a common weight at least their initial maximum. An immutable member can therefore be feasible only if it already has that maximum. When feasible, raising the deficient mutable member to that maximum attains cost |w_i−w_j|. Distinct involution pairs are disjoint, so these costs add with no interaction. Taking the minimum over involutions proves both the lower bound against arbitrary clone sequences and an explicit attaining construction. A fixed block needs no change. If no involution is feasible, changing roots or using more operations cannot repair the obstruction. This also correctly handles empty mutable sets and non-self-converse quotients.

## Orbit fibres, labelled mass and asymptotics

For ordinary classes, the standard orbit stabilizer identity gives the automorphism cycle index. For self-converse classes, each invariant weight orbit contributes exactly |Aut(Q)| to the sum of fixed pairs (weight, anti-automorphism); each non-invariant orbit contributes zero. Dividing by |Aut(Q)| is therefore the correct twisted trace, with no extra factor2 or centralizer divisor at this stage.

The unique chain lifts identify the full graph stabilizer with the weight stabilizer. Consequently the ordinary automorphism mass is the number of positive weight vectors divided by |Aut(Q)|, proving the stated labelled formula n! binom(n−1,q−1)/|Aut(Q)|. These are unmarked graph classes, not clone recipes or rooted classes. For q1 all three series are z/(1−z), and q2 cannot be canonical.

The anti-involution conjugacy theorem gives exactly |G|/C leading anti terms, where C is the centralizer of one involution. An anti-automorphism has at most one fixed vertex. Every longer cycle is even with length2 modulo4: an odd cycle of length>1 would yield an odd anti power fixing a nontrivial induced tournament pointwise. Non-involutory anti maps contain a cycle of length at least6, and hence have strictly fewer cycles than an involution. Their rational generating functions have lower pole order at every root of unity. The ordinary identity term likewise has strictly greater pole order than all nonidentity automorphisms.

For odd q=2r+1, the leading series is z^q/[C(1−z)(1−z²)^r]. Only z1 has pole order r+1, yielding n^r/(C2^r r!). The z−1 pole has order r and cannot change the leading constant. For even q=2r, all anti cycles are even, so odd coefficients vanish exactly. On the even subsequence, the leading coefficient is C^−1 binom(n/2−1,r−1), giving the stated factor2^(r−1), not2^r. Combining with the ordinary identity term proves both ratios and their exponents. The singleton and empty anti-coset exceptions are explicit.

Independent exhaustive quotient permutations give the Paley group order21, seven anti-involutions and centralizer3. Its other14 anti maps have cycle type1*6. These yield exactly the stated rational functions and the ratio constant105. The class-sampling interpretation is correctly distinguished from Hamilton weighting or all-tournament sampling.

## Independent controls

The standalone checker imports no author code. It exhaustively permutes each of the11 named quotient seeds to obtain its groups. On actual blow-ups it uses a fixed degree-ordered backtracking search, rather than the source's dynamically selected vertex search. Maximal transitive modules are independently tested by all subsets and absence of directed triangles. Ordinary weight orbits are formed by disjoint-set union rather than lexicographic canonical keys.

Most importantly, the completion check is a different optimization algorithm. For each mutable-block mask, it starts a multi-source reverse breadth-first search at every self-converse weight vector in the entire box{1,2,3}^q. An allowed edge decrements one mutable coordinate by one. Its exact distances to every starting vector are compared with the proposed formula. This checks unreachable cases as well as all optimal finite costs, rather than repeating the source's direct minimization over dominating vectors.

Optimized and assertion-enabled unoptimized receipts agree on:

- 11 canonical quotient seeds, including the singleton, weak/non-self-converse examples, the six-vertex ribbon and Paley7;
- 2,140 actual blow-ups and4,280 exact map-set comparisons;
- 520 full maximal-block checks;
- 82 ordinary/twisted orbit coefficient checks;
- 368,286 reverse-BFS completion comparisons, covering every mutable mask;
- 25 anti-involution centralizer/conjugacy checks.

The supplied source's complete larger control receipt also replays exactly in both modes. These tests corroborate the implementation; the all-order completion lower bound and pole-order analysis above establish the unrestricted lengths and asymptotics within the fixed quotient language.

## Reproduction and hashes

Compile `check.cpp` as C++17 and run it. It has no file dependencies and writes only its receipt to stdout. The source copy can be compiled separately. No audit step changes the source packet.

- check.cpp: `6adcf4a72af64a28265b6d9605887995a4717f02852ebbaf393c4a504de48d16`.
- independent-normal.json: `58c45059befcf95f538b66d343a5f9c0ae9c7b5f0104162af8bd6686ae2ba6f7`.
- independent-debug.json: `58c45059befcf95f538b66d343a5f9c0ae9c7b5f0104162af8bd6686ae2ba6f7`.
