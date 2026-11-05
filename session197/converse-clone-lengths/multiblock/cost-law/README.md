# Exact cost of self-converse completion in cyclic three-block tournaments

This independently reviewed addendum studies ordinary unmarked classes C3[Tr_a,Tr_b,Tr_c], with a,b,c positive. The cost D is the minimum number of vertices that must be added by transitive cloning to make the graph self-converse.

- [Proof and exact bivariate generating function](C3_COMPLETION_LAW.md)
- [Independent audit](audit/REPORT.md)
- [Standalone source checker](verify_c3_cost.cpp)
- [Standalone independent checker](audit/check.py)

The exact generating function is

    F(z,t)=z^3(1+t*z^3)/[(1-z)(1-z^2)(1-t*z^3)].

Here z marks order and t marks the optimal added order. It follows from D=min(|a-b|,|b-c|,|c-a|), an explicit sorted-triple bijection, and the complete physical minimum-completion theorem.

For a uniformly chosen ordinary class of order n, D/n converges to a distribution on[0,1/3] with survival function(1-3x)^2. The mean added order divided by n tends to1/9. The probability of requiring no added vertices is asymptotic to3/n. The proof distinguishes ordinary class multiplicity from automorphism or Hamilton weights and covers all endpoint cases.

This is a physical symmetry-completion law, not a transformation that preserves H or establishes prime/7p coverage.

## Reproduction

Both checkers have no input-file dependencies.

    g++ -O3 -std=c++17 verify_c3_cost.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify_c3_cost.cpp -o verify-debug
    ./verify-debug
    python3 audit/check.py
    python3 -O audit/check.py

The source outputs match c3-cost.txt and c3-cost-debug.txt. Independent outputs match audit/normal.json and audit/optimized.json. The source enumerates183800 classes through150 and independently optimizes76 small actual graphs. The independent checker enumerates437844 classes through200 using sorted partitions, checks the formal bivariate series, and uses literal degree-compatible permutations on actual enlarged graphs.

## Included unchanged dependencies

The frozen proof retains its original source name PROOF.md for the complete multi-block theorem. Its release alias is dependencies/MULTIBLOCK.md. That theorem and this addendum now both have independent PASS reviews. The earlier canonical-block/group proofs are included as dependencies/CONVERSE_CLONE_LENGTHS.md and CONVERSE_MASS.md, along with the canonical-block packet's six unchanged arithmetic prerequisites.

The existing main multiblock checkpoint is not edited by this addendum. No literature-priority claim is made. RELEASE_FILES.sha256 binds every other included file.
