# Self-converse clone lengths: finite or universal

This reviewed Session197 packet proves a structural decision theorem for one fixed root replaced by a nonempty transitive chain. The canonical maximal transitive-module partition changes only by enlarging its root block.

- [Proof](PROOF.md)
- [Independent map-for-map audit](audit/REPORT.md)
- [Standalone source verifier](verify.cpp)
- [Standalone independent checker](audit/check.cpp)

The self-converse clone lengths are either all positive integers or a finite set contained in explicit block-size collision lengths. For an n-vertex source, every finite exception is at most n-2. The threshold is sharp: C3[K1,Tr_(n-2),K1] is self-converse precisely at lengths1 and n-2 for n>=4.

Automorphisms stabilize outside the same collision set. Every eventual unit exponential clone clock A*b^m±1 has trivial root stabilizer and fully rigid sufficiently long outputs, giving exact labelled-orbit consequences. Infinitely-often self-converse negative-unit clocks inherit the ternary/odd-square obstruction after the root block is contracted.

The persistent self-converse positive-binary unit clocks are completely classified: amplitude2^(l-1), with equal transitive padding around C3[Tr_l,K1,K1], and the root at any chain position. A root-fixed anti-involution exists only at the central position of an odd chain. The distinction matters: the four-vertex source[6,12,8,1] gives3*2^m-1 but is self-converse only at m=1,2.

No result here concerns arbitrary adaptive graph reconstruction, the unrestricted spectrum, or prime infinitude.

## Reproduce the controls

Both checkers are standalone C++17 and require no bank or external input.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug
    g++ -O3 -std=c++17 audit/check.cpp -o independent
    ./independent
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 audit/check.cpp -o independent-debug
    ./independent-debug

The first output matches receipt.txt and receipt-debug.txt. The second matches audit/normal.txt and audit/debug.txt. The independent audit checks32110 actual clone graphs, compares automorphism and converse maps exactly with their weighted quotient lifts, and verifies the classification's rooted and labelled factors. These finite controls support the separately proved all-order statements.

## Included reviewed mathematical dependencies

Frozen source and review retain their session-local source names. Their unchanged release aliases are:

- General unit-clock/factorial/SCC reduction: dependencies/UNIT_CLOCK_REDUCTION.md.
- Marked positive-binary classification: dependencies/MARKED_BINARY.md.
- Marked negative-base rigidity: dependencies/NEGATIVE_BASE_RIGIDITY.md.
- The latter's quantitative square/amplitude bound: dependencies/SYMMETRIC_SQUARE.md.
- Its even-order self-converse congruence: dependencies/EVEN_CONVERSE.md.
- Root-run and matching conventions: dependencies/MINIMUM_TERNARY.md.

The structural module and stabilization arguments are proved in this packet. The arithmetic corollaries retain the dependencies' precise hypotheses and do not assert attainment at amplitude242. No literature-priority claim is made. RELEASE_FILES.sha256 binds every other included file.
