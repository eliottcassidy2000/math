# Independent audit: complete universal-source obstruction through order ten

Status: APPROVED FINITE-EXACT result, using the separately approved rooted-reciprocity theorem. No scalar spectrum exclusion, all-order paired-source obstruction, order-eleven census, or eventual-coverage conclusion is approved.

## Approved statement

No two physical rooted tournaments BOTH of order at most10 satisfy the width-eight universal all-O/Z/J paired-source criterion. The complete necessary-modular inventories contain0 sources through order7,41 order8 presentations,301 order9 presentations, and1412 order10 presentations. These are cut presentations over canonical unrooted deletion cores, not rooted-isomorphism counts. Their1549 full-profile-plus-label keys contain no cross-residue match.

This conclusion relies on the already independently approved reciprocity theorem, whose exact frozen bytes match the dependency here. It forces opposite order parity; starts(4,2),(10,8),(16,14),(22,20); eight-coordinate divisibility by3; and L=6 mod12, C=b=0 mod12, T=0 mod6. The singleton is excluded by odd start. Thus restricting to these necessary conditions cannot discard an admissible pair.

## Canonical completeness

The producer's ordered equitable refinement is isomorphism invariant: count vectors into ordered cells determine invariant splits. Selecting the first smallest nonsingleton cell is likewise invariant. Every individualization branch is traversed, with the selected vertex placed first and the remainder after it. At a discrete leaf the code is an actual relabeling; hence equal codes imply isomorphism. Isomorphisms biject the search trees, hence give the same minimum code. The routine need not give the minimum over ALL labelings; it only needs a complete canonical invariant, which it does.

Every one-vertex extension of every prior representative is processed. Induction from the singleton therefore proves completeness at each order; matching known class counts is only a control. All codes through order9 fit within36 bits; incidence masks fit within9 bits.

independent_generation.cpp uses a separate canonicalizer: incoming counts rather than outgoing counts, an initial degree partition, largest-cell branching rather than smallest-cell branching, and distinguished singleton LAST rather than first. It checks canonical equivalence against all permutations on all1099 labeled tournaments through order5. It then confirms6880 distinct inherited order8 classes,191536 distinct supplied order9 classes, and checks ALL1,761,280 order9 augmentations of the inherited order8 representatives. Every augmentation belongs to the supplied class set and every supplied class occurs. This is an independent closure check, not merely sampling relabelings or reproducing the class count.

The producer generator was also rebuilt and replayed: BANK2 through BANK9 are byte-identical to the frozen banks. Its own complete small-label, inherited-order8,191536-idempotence, and2000-relabel controls replay exactly. BANK9.txt.gz decompresses exactly to the bound plaintext. The trivial supplied BANK1.txt is correct; generate.cpp starts from that singleton internally but does not itself write BANK1.txt, so reproduction should retain or explicitly create that trivial input.

## Root-incidence compiler and safe filters

Every rooted U arises from its unrooted deletion Q and one root incidence subset. Transporting the cut under a relabeling of Q shows that a complete unrooted deletion bank suffices; no rooted class enumeration is required.

The first filters d=H(Q) divisible by3 and b=B_Q(2) divisible by12 are inherited necessary conditions. The minimum-positive-full-start filter is safe: s is the sum of selected nonnegative start counts, so if every positive summand exceeds22, no permitted positive s can occur. All positive bounded-start possibilities are in{2,4,8,10,14,16,20,22}. Start26 cannot pair, because its residue-one companion would require28 and exceed width8.

The linear formulas for s,e,L,R,T and the quadratic h,C formulas follow by splitting the root-containing path before and after the root. For C, the remaining vertices form the other labeled bin and the factor2 selects the root's bin. The term L+R-d includes the root-only path exactly once; the extra convolution includes nonempty before/after portions. For T, the root-only first path contributes b. No empty case is omitted or double-counted.

Only Section5 of the copied SUFFIX_RIGIDITY.md is needed here. Its root-incidence decomposition was independently checked in this audit. The packet's unrelated suffix theorem and all-tail claims are NOT being imported or approved by this audit.

independent_search.cpp computes all path tables by forward extensions in the graph and its converse, and computes the quadratic cover expression by first convolving the suffix with the leftover ordinary bin. It deliberately DOES NOT use the producer's minimum-start filter. At deletion order9 it checks33,955,840 root cuts after only the d,b filters, versus the producer's3,258,880 linear checks, and finds exactly the same1412 complete profiles, byte-equivalent as indexed profiles. The independent search also reproduces all smaller-order inventories. Thus the minimum-start pruning received a direct hostile audit.

## Literal profiles and matching

literal_profiles.cpp independently reconstructs every one of the1412 order10 candidate graphs and recomputes all nine coordinates with direct path and cover tables. All agree in optimized and debug builds. The direct root graph uses a separate labeling from deletion-bank generation. Producer verify_candidates also replays exactly.

The final key is eight nonstart coordinates plus s+2r, r=s mod3. Between residue1 and residue2, equal last coordinates mean exactly s1=s2+2. The combined search writes a detected match before any representative overwrite, so no overwrite can hide one. Independent combination instead groups the eight coordinates alone and directly compares allowed start pairs; it finds none and independently recovers1549 label keys.

The older TARGETS list has381 entries including80 start26 presentations. Those extras do not create false absence: the complete conclusion comes from the combined exact inventories, not from matching only this older target list. The source-order9 comparison count301 is consistent after removing the80 unpairable presentations.

## Separate asynchronous depth-seven scan

The broader scan is independently scoped to1874 supplied literal sources and arbitrary independent O/Z/J words of length at most7. Its source JSON and text agree exactly; full profiles are deduplicated and all starts are permitted. There are593 inherited-bank sources and1281 new-bank sources after deduplication. All1874 literal source profiles were independently reconstructed and recounted.

The executable replay gives6,146,720 presentations,2,335,127 full-profile-plus-label keys, and0 hits, exactly as recorded. These are finite construction-menu negatives, not absence of arbitrary order11..17 sources or all-depth extensions. Only the separate BOTH-orders-at-most10 result is complete over its full order universe.

At order10,10!*binomial(12,2) bounds the relevant three-bin counts. The depth-seven sources reach at most17 vertices, where17!*binomial(19,2) still fits signed64-bit. Actual source starts remain bounded because the inherited operations preserve them. No arbitrary-order numerical safety is asserted.

## Reproduction and provenance

The copied reciprocity dependency exactly matches the prior approved freeze. The deletion bank matches the separately audited inherited bank. All files in frozen/ are bound by FROZEN_SHA256SUMS; the compressed and uncompressed order9 banks were cross-checked. Producer source files were not modified.

    g++ -O3 -std=c++17 independent_generation.cpp -o /tmp/session216-independent-generation
    /tmp/session216-independent-generation frozen/dependencies/BANK_THROUGH8.txt frozen/BANK9.txt
    g++ -O3 -std=c++17 independent_search.cpp -o /tmp/session216-independent-search
    /tmp/session216-independent-search frozen/BANK9.txt frozen/CANDIDATES10.txt
    g++ -O3 -std=c++17 literal_profiles.cpp -o /tmp/session216-literal-profiles
    /tmp/session216-literal-profiles frozen/BANK9.txt frozen/CANDIDATES10.txt
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

INDEPENDENT_GENERATION.json, INDEPENDENT_SEARCH1..9.json, INDEPENDENT_COMBINE.json, the LITERAL receipts, and REPLAY.json record the independent outcomes. Producer replay outputs remain under replay/. The original proof labels its reciprocity dependency conditional because that dependency was awaiting audit when written; the independent reciprocity audit has now discharged that named dependency. No publishing or stage edits were performed by this worker. Parent remains sole publisher.
