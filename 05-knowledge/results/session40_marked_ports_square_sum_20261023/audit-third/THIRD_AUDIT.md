# Third independent audit: prime-cell carriers and color-word limits

## Verdict

PASS within the explicit scopes below. This packet extends, without modifying, the first and second immutable audits. It does not certify a universal tournament-metagraph classification.

## Five-cell repair

`check_five_cell_orbits.py` independently enumerates all 120 matching permutations, their rotation double orbits, converse orbits, and Hamiltonian counts by vertex DFS. There are eight ordinary types, six self-converse types and seven merged types. In lexicographic representative order the Hamiltonian counts are 8265,8245,8221,8201,8193,8175,8175,8245; J is 10,4,−2,−4,2,0,0,−10; K is 5,3,1,−3,−1,5,−5,−5.

Exact J and J modulo16 classify the seven merged types. J modulo8 alone has four buckets. The pair (J modulo8,K modulo4) distinguishes all eight ordinary types. Since K is odd, the refinement adds one bit. Testing all64 linear combinations agrees with the proof that aJ+bK modulo8 has fixed parity and at most four labels. This is a restriction on linear compression, not arbitrary encoding. K is not globally converse-odd. The intrinsic short arcs are independently verified by common-outneighbor counts.

## Prime matching formulas and constructive realizability

`check_prime_matching_counts.py` independently enumerates every permutation at p=3,5,7 and quotients by actual rotations and reflected-inverse converse. Ordinary/self-converse/merged counts are (2,2,2), (8,6,7), (108,34,71). The p=11 row is a formula evaluation, not exhaustive enumeration.

The all-prime proof was checked directly: the ordinary rotation fixed sets give N=((p−1)!+(p−1)^2)/p. In the converse coset, the equation becomes (q t_a r)^2=t_(b−a). For a=b there are I_p involutions. Otherwise any square root commutes with the nontrivial p-cycle, and its centralizer is cyclic, yielding exactly one square root. Thus F=(I_p+p−1)/p and M=(N+F)/2. The positive involution EGF bound establishes the stated asymptotics. The scalar-J obstruction counts possible values, not numerical coincidence.

`check_displacement_necklaces.py` separately tests all nonzero displacement words: 8 at p=3,1024 at p=5,279936 at p=7. Exactly2,24,720 satisfy the total-sum and distinct-proper-partial-sum criterion, giving2,8,108 rotation necklaces. It also checks necklace classification against permutation double orbits. Converse must reconstruct, invert and reflect the actual permutation; it is not generally reversal of the necklace. The constant p=5 word22222 becomes33333 under converse, a direct hostile control.

The automorphism and path-fiber claims in the prime note retain the actual directed path action, which is free: fixing a full ordered path fixes every vertex. They do not count all tournament classes.

## Square-sum color extension

`check_color_extension.py` verifies all56 explicit path witnesses at birth orders25 through80 and both distinct monochromatic order80 paths. It checks every path and ordered color label in the producer's already audited census25..32, then compiles an independent C++ plain path DFS for33..36. The C++ source is included; no previously computed producer path list is used in that new census.

| order | paths up to reversal | ordered-color labels with endpoints | histogram labels | largest ordered-word fiber |
|---|---:|---:|---:|---:|
|33|669|669|537|1|
|34|4041|4041|2158|1|
|35|17175|17175|5827|1|
|36|12960|12957|4947|2|

The independent searches visit479643168 partial paths. The explicit order36 collision has endpoints6 and32 and is checked vertex-by-vertex and edge-by-edge. Combined with endpoint uniqueness through24 and the independently checked25..32 census, this certifies first failure at36 within the complete finite range1..36. It is not an infinite existence or color reconstruction theorem.

The order80 examples remain a different, stronger monochromatic mechanism, not a claimed smallest collision. The three colors are birth categories, not a proper graph coloring or a usage partition. Actual shell sums plus a starting owner reconstruct through x_next=s−x, with domain and no-repeat checks; coarse color words cannot replace those owners.

## Final proof and synthesis review

The final prime-cell, five-cell, endpoint-hostile and path-scatter notes were checked with their exact source hashes recorded in MANIFEST.json. The final session synthesis and index accurately distinguish target-class decoding, actual path transport, and literal partner ownership; they scope regular-block counts to two cyclic prime cells rather than A000568. The separate regular-block packet supplies full coefficient proofs, independent literal coefficient multiplication and a matching-union census. Its normal and optimized receipts agree with the root's separately implemented coefficient DP.

## Reproduction

Run each Python checker normally and with `python -O`. The color wrapper requires g++ and automatically compiles its adjacent C++ source into a temporary directory. Place this audit folder adjacent to `square-sum/` (published layout), or retain the development `session40-square-sum/` sibling. All assertions essential to verification use explicit exceptions and survive optimized Python. Checkers overwrite only their own deterministic receipts. Exact source hashes and commands are in MANIFEST.json.
