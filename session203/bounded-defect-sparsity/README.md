# All bounded physical Z defects have sparse profile support

**PROVED, independently audited.** For any fixed positive integer K, the number of distinct B values up to X realized by nontrivial positive-start rooted tournaments with physical Z defect D=B(ZU)-2B<=K is

    O_K(log^(2K+1) X).

Thus every fixed absolute defect cap has zero density among odd B values. This covers all physical sources, including D>s, rather than just the exact equality face or a selected carrier menu. The hidden constant depends on K; no bound with K growing with X is asserted.

The proof removes at most two exceptional vertices to obtain a tournament with at most K Hamilton paths and at most K marked exceptions. Its nontrivial SCC cores are finite; only linearly ordered transitive singleton runs remain unbounded. An exact positive sum of onto-map counts describes every Hamilton path of the physical B readout. Any run too long for a count at most X is invisible to that count and can be shortened for the counting upper bound. No preservation of D under this normalization is required.

The rooted-edge boundary is handled separately, and the singleton source is excluded from the physical-defect statement. Fixed finite staggered shifts preserve zero density. Therefore any bounded-start staggered construction covering every sufficiently large odd base must allow unbounded absolute defect. This does not establish a spectrum gap, disprove F(8,137), or prove the eventual threshold 23; both remain OPEN.

Read PROOF.md and audit/REVIEW.md. The reviewed proof and status-only delta are preserved. Dependency metadata is updated only to the now-published audit URL and theorem status, with the original and exact metadata diff retained. Immutable source blobs are listed in dependencies/PUBLISHED_PROVENANCE.json.

Reproduce the producer controls:

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > replay.json
    cmp replay.json NORMAL.json

The independent audit checks all 164,030 positive marked inputs through order 6, extracts actual homogeneous-run patterns from all 164,028 nonedge sources, and independently verifies 636 chain patterns, 9,096 literal expansions and 93,768 normalization cases. Empty and invisible runs are included. These controls supplement the all-order proof.

Every frozen audit entry and its original manifest are preserved. This focused release changes no navigation or theorem-index files.
