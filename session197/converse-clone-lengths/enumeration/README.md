# Exact binary-family orbits and sampling bias

This reviewed Session197 addendum counts the complete persistent self-converse positive-binary unit-clock family, as ordinary classes, rooted classes, converse orbits, automorphism mass, and Hamilton-normalized presentations.

- [All-order proof](PROOF.md)
- [Independent audit](audit/REPORT.md)
- [Standalone source verifier](verify.cpp)
- [Standalone independent checker](audit/check.cpp)

The family is Tr_d -> C3[Tr_l,K1,K1] -> Tr_d, with l>=1 and d>=0, and the marked clone root at any chain position. The underlying graph has order l+2+2d and H=2^l+1. The included complete classification is separate from this enumeration addendum.

The normalized self-converse fibre of each class is exactly1 for odd l, and2^(l/2)+1 for even l. This gives a concrete asymptotic sampling discrepancy: at even order2s, conditioned on this family, the strong-class probability is1/(s-1) under uniform ordinary-class sampling, but(2^(s-1)+1)/(2^s+s-3) under uniform Hamilton-normalized self-converse sampling. The limits are0 and1/2, respectively. At odd order both probabilities are1/s.

All class sequences are zero below order three. Ordinary class counts, automorphism mass, labelled counts and normalization fibres are kept distinct. The statements do not imply unrestricted spectrum coverage or prime infinitude.

## Reproduce the controls

Both C++17 programs are standalone and have no input-file dependencies.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug
    g++ -O3 -std=c++17 audit/check.cpp -o independent
    ./independent
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 audit/check.cpp -o independent-debug
    ./independent-debug

The source outputs match receipt.txt and receipt-debug.txt. Independent outputs match audit/normal.txt and audit/debug.txt. The checkers construct42 actual graphs through order14, independently count10950 Hamilton words and282 fixed normalized words, and verify exact normalized adjacency sets, centralizer fibres and class counts. These finite controls support the explicit all-order derivations rather than substitute for them.

## Included unchanged dependencies

The proof and audit retain their frozen source names. Their release aliases are:

- Complete physical classification and stabilization: dependencies/CONVERSE_CLONE_LENGTHS.md.
- Exact normalization fibres and total normalized universe: dependencies/CONVERSE_MASS.md.
- The classification's arithmetic dependencies: dependencies/UNIT_CLOCK_REDUCTION.md, MARKED_BINARY.md, MINIMUM_TERNARY.md, NEGATIVE_BASE_RIGIDITY.md, SYMMETRIC_SQUARE.md and EVEN_CONVERSE.md.

The arithmetic dependency files retain their precise marked-root hypotheses and do not assert attainment at amplitude242. No literature-priority claim is made. RELEASE_FILES.sha256 binds every other file in this packet.
