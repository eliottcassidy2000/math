# The full binary-mirror bin operator

This reviewed Session197 packet proves an exact full-polynomial substitution operator, including its empty-input boundary. It also gives explicit positive difference coefficients and a physical pair of tournaments with the same B(2)/H ratio but different full profiles.

- [Proof](PROOF.md)
- [Independent review](audit/REPORT.md)
- [Standalone C++ verifier](verify.cpp)
- [Independent odd-cycle and literal-word checker](audit/check.py)

For a nonempty marked side S, its converse mirror I, and any input T including empty, put f=B_T and let E_S sum residual bin profiles over paths ending at the mark. Then

    B_(F_T)(x)=B_S(x)^2 f(x)+x E_S(x)^2 f(x+1)-x E_S(-x)^2 f(x-1).

The empty input retains a genuine base-zero correction. The ten-vertex graph Q3 has H=161 and B(2)=10304, sharing ratio64 with Tr6, but their B(3) values are not proportional. Pair-based ratio cancellation therefore cannot by itself certify primitive integer content. This is not a total-positivity theorem or an unrestricted spectrum classification.

## Reproduction

Both checkers are standalone and require no external data or research helper files.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug
    python3 audit/check.py
    python3 -O audit/check.py

The C++ output matches receipt.txt and receipt-debug.txt. The Python output matches audit/normal.json and audit/optimized.json. The source checks1505 full physical substitutions; the independent checker uses325 physical odd-cycle profiles and2395 literal crossing-cycle words. All-order validity comes from the proof, not extrapolation from these controls.

## Included reviewed mathematical dependencies

The frozen proofs and review retain their original session scratch source names for hash binding. Their local release aliases are:

- Original path-polynomial PROOF.md: [dependencies/PROFILE_BASIS.md](dependencies/PROFILE_BASIS.md), with the self-contained odd-cycle/bin derivation in Section2.
- Original UNIVERSAL_CONTEXTS.md: [dependencies/UNIVERSAL_CONTEXTS.md](dependencies/UNIVERSAL_CONTEXTS.md). Its companion PROOF.md means PROFILE_BASIS.md above.
- Original ENDPOINT_RECIPROCITY.md: [dependencies/ENDPOINT_RECIPROCITY.md](dependencies/ENDPOINT_RECIPROCITY.md).
- Original binary-boundary PROOF.md: [dependencies/BINARY_BOUNDARY.md](dependencies/BINARY_BOUNDARY.md).

These dependencies are copied unchanged. The earlier [binary-boundary release](https://github.com/eliottcassidy2000/math/blob/01177cf7bf2e95de57bc40f9f5b7686e7ee2331c/session197/binary-boundary/README.md) includes its independent review and standalone controls. Dependency verifiers are not needed to run either checker in this packet. No literature-priority claim is made.

RELEASE_FILES.sha256 binds every other file in this release.
