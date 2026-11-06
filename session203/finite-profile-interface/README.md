# Finite weighted presentation of bounded-start profiles

PROVED, independently audited. This packet gives a physical finite-interface construction for every
actual rooted profile (B,s) with 1<=s<=K. A prescribed marked set is measured
by its total Hamilton-path start count, which stays within budget when a
chosen root is replaced by its initial-SCC entry set. The resulting register
bound is

    W(K)=1+max_{0<=a<=K}(m(a)+m(K-a)),
    m(0)=0, m(t)=floor(sqrt(4t-1)) for t>=1.

For K=26, W=15. A finite nonnegative path-forest transfer computes exact B
and s from register words. Capped vectors recognize each bounded start
budget by finite control, while uncapped vectors retain the actual B.

The numerical image is not asserted to be regular in any positional base,
semilinear, cofinite, or efficiently computable. F(8,137) remains OPEN.
The reviewed proof and dependencies determine the precise scope.

## Self-contained author controls

Run from this packet directory:

    python3 check_forest.py > /tmp/profile-forest-normal.json
    python3 -O check_forest.py > /tmp/profile-forest-optimized.json
    cmp FOREST_NORMAL.json /tmp/profile-forest-normal.json
    cmp FOREST_OPTIMIZED.json /tmp/profile-forest-optimized.json
    g++ -O3 -std=c++17 check_register.cpp -o /tmp/profile-register
    /tmp/profile-register REGISTER_FIXTURES.tsv > /tmp/profile-register-normal.json
    cmp REGISTER_FIXTURES_NORMAL.json /tmp/profile-register-normal.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check_register.cpp -o /tmp/profile-register-debug
    /tmp/profile-register-debug REGISTER_FIXTURES.tsv > /tmp/profile-register-debug.json
    cmp REGISTER_FIXTURES_DEBUG.json /tmp/profile-register-debug.json
    g++ -O3 -std=c++17 check_marked.cpp -o /tmp/profile-marked
    /tmp/profile-marked > /tmp/profile-marked-normal.json
    cmp MARKED_NORMAL.json /tmp/profile-marked-normal.json

The forest checker compares h,s,L,B,C2 and the rooted three-path response
against literal subset counting on all 1,099 labeled inputs through order 5,
plus 90 register words of order 9 with forgetting and label reuse.

The register fixture has 281 actual sources through order 13; each is checked
for strongness and its given start count, compiled to a register word, and
compared literally with the original graph. Its B column is inherited
metadata and is not used or certified by this register check. The fixture
is a deterministic selection from a larger research bank; its provenance
file records the selection and source hash. No completeness of that bank
is asserted or needed.

The marked-set checker tests all 17,239 nonempty marked subsets of total
start budget at most 26 in all 571 strong labeled inputs through order 5.
It verifies the general marked-budget construction, not only singleton roots.
All runtime checks remain enabled under optimization.

The separate 65,125-source research replay is not needed to reproduce any
claim in this packet and is deliberately excluded from release receipts.

Independent checks, frozen original inputs, scope analysis, and reproduction
commands are in audit/REVIEW.md. The audit separately verifies complete
intermediate forest patterns, saturation, reusable register words, and
weak-tail preservation of both B and s. PROMOTION.json records the reviewed
status/convention clarifications and the explicit-return cleanup.
