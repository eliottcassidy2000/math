# Complete fixed-core inverse algorithm

**PROVED and independently audited PASS.** For any fixed finite bank of nonempty strong tournament cores, the star-closure template has an exact all-length membership algorithm using O(sqrt(N) (log N)^2) arithmetic steps and O(log N) live node records. This is template membership, not a proof of all-odd Hamilton-path-count coverage. The three-core occupancy conjecture remains OPEN.

[Proof and implementation domain](SUFFIX_INVERSE.md). The C++ executable implements the Session216 bank of80 core-cut states, for integer1<=N<=10^15. The abstract theorem permits arbitrary-precision inputs and any fixed finite strong-core bank. Arithmetic-operation complexity is distinguished from bit complexity in the proof.

## Dependencies

- [Session216 three-core construction](../../session216_overnight_constructions_20261107/three-core/PROOF.md): literal cores, exact state coordinates, strongness condition, counting formula and finite inverse bounds.
- [Session216 numeration](../../session216_overnight_constructions_20261107/numeration/PROOF.md): prefix-parity recursion and the earlier balanced-band inverse. The present theorem adds the reversed-suffix search and complete square-root bound.
- [Original80-state table](../../session216_overnight_constructions_20261107/three-core/CORE_STATES.txt): CORE_STATES.txt is its byte-identical local copy.

## Reproduce

From this directory, using C++17 and Python3 standard libraries:

    g++ -std=c++17 -O3 inverse.cpp -o /tmp/session217_inverse
    /tmp/session217_inverse CORE_STATES.txt 1 7 21 23 47 63 65 265 287 1237 4295971 9948407 1000000007 1000000000039 > /tmp/session217_inverse_cases.jsonl
    cmp /tmp/session217_inverse_cases.jsonl INVERSE_CASES.jsonl
    python3 check_suffix.py > /tmp/session217_suffix_check.json
    cmp /tmp/session217_suffix_check.json SUFFIX_CHECK.json

The full reconstruction verifier currently expects the compiled executable to be named inverse in this directory:

    g++ -std=c++17 -O3 inverse.cpp -o inverse
    python3 verify_combined.py
    rm inverse

It regenerates COMBINED_VERIFY.json, LITERAL_GRAPHS.json and CASE_GRAPHS.json. Run in a scratch copy if preserving release files byte-for-byte is desired. Importing check_suffix also reruns and prints its synthetic suffix checks; that is expected.

## Retained checks

-4,400 synthetic suffix-inverse controls, exhaustively compared with suffix words through length10.
-Complete membership support agrees with independent all-length parameter enumeration for all targets1..3000:29,952 legal candidate recipes and1,490 represented counts.
-60 literal strong graphs counted independently by vertex-subset Hamilton-path DP.
-9 larger/example recipes reconstructed as literal outgoing masks and checked for strong connectivity; their counts use the inherited exact formula, not large-graph subset DP.
-Permanent template omissions47 and63 recovered by complete all-length inverse searches. Neither is an unrestricted tournament-spectrum omission.

For N=1000000000039 the returned prefix-only recipe has core-cut index0, P=16231501406 and empty suffix. It is a concrete improvement on that target's previously retained balanced-band witness, not a global coverage claim.

The independent adversarial review and final source freeze are supplied alongside this packet by the separate reviewer. PUBLICATION_FILES.txt lists this packet's owned deliverables; SHA256SUMS records their exact bytes. Built executables, caches, provisional notes and the separate entropy audit are excluded.

## Independent audit and portable replay

[Independent review](audit/REVIEW.md) checks the proof, entire recipe sets and final DFS executable. Historical review commands retain their original source-directory names. From this release directory use:

    python3 audit/audit.py
    g++ -std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=all inverse.cpp -o /tmp/session217_inverse_audit_ubsan
    python3 audit/check_final_cpp.py

The last check writes its declared DFS_RECEIPT.json into audit/. PORTABLE_PATH.diff records the sole checker change: the relocated core-table path. The algorithms are unchanged. The returned recipe for 10^12+39 has 40 vertices; the earlier retained balanced-band recipe has 41. No minimum-order claim is made.
