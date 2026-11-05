# Independent audit: full C2 bilinear normal form and bounded inverse

**PASS, no correction.** The two-source graft theorem covers the full strong C2 class except the stated regular-five graph. Its ordered graph presentations are bijective under the cited class-count theorems. The scalar oracle remains an exact finite physical inventory, and the positive product bounds make the numerical-target inverse complete. No cofinite scalar-spectrum conclusion is implied.

## Frozen binding

Source: session197-hybrid-constructions/bilinear-c2-inverse/.

- PROOF.md: `5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813`.
- bilinear.py: `33cfe6cea487b17da44e083d7e81bd860dfeddcab3096deb6d135191e5537ac0`.
- verify.py: `8735852f72520deed64d70fd916faf521ad671f4ef783a71932b35dd23419bbc`.
- core.py: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`.
- Both copied source replays and all supplied validation receipts: `f8907d7e56531c8d49cff29f531eda0d2740e2bca7430d9872a6feca3d9fb757`.
- SUPPORT.json: `427003e29b23e8045ae2c8c7efdf29ed6fb6378e7385ad21c00396563738b6cc`.

The exact dependency hashes are recorded in the source RELEASE_FILES.sha256. I read the named ROOTED, START_TWO, C2_COMPLETE and FIBONACCI_CLASS snapshots, including the rooted-source and nonregular unrooted-C2 generating functions used in the injectivity argument. Sources were copied before replay because the verifier writes SUPPORT and VALIDATION beside itself. Original files remain unchanged.

## Physical source oracle

A positive root-start count puts the root in the initial SCC. Count1 forces every later SCC to have H1 and therefore to be a singleton. This proves the strong-source plus unique weak-tail decomposition. Appending a sink gives exactly (D,0,2U,d): the sink can occupy the final position of either reserved path slot, including an initially empty first slot. Ordinary sewing gives the stated first three responses, and deleting the newborn gadget root leaves the old graph followed by a transitive chain, so its fourth response is exactly old D.

The inequalities used to bound the oracle are physical. Cutting a Hamilton path before the root injects it into the U covers. Inserting the root just before the first out-neighbor in any deletion path gives a non-ending path, proving D−E>=d for every nontrivial source. Such an out-neighbor exists because a root-start Hamilton path visits another vertex. The U>=order injection supplies at least one reserved cover for every possible length of its initial rooted path segment.

Both moves strictly increase U. Therefore every ancestor of a cap-B state stays within the cap, and L<=B/U is safe because a possible accepted successor satisfies U_new>=D_new>=LU. Deduplicating full quadruples preserves all future scalar responses; recording one physical recipe per quadruple suffices for scalar construction. This is not a claim that distinct rooted graphs have distinct quadruples, or that every nonnegative quadruple is realizable.

## Graft identity and pair forms

The arbitrary-G path identity has the correct correction term. A path using a reversed seam has unique exterior run pattern P,X,Q,Y, with X ending at r and Q starting at one selected entry. The only remaining forbidden join requires X to be the singleton r and P to end at the other entry; those cases are exactly dK. Both reversed seams cannot be used, because both leave the same r. For a cycle with |G|>=2, the G path starts away from r, so the incoming T-to-G join is automatically legal. This proves the stated cycle product and explains why a right singleton must be excluded.

The Hamilton-path injection into R−K is valid: cut before a unless its predecessor is b, then cut before b. Both bad predecessor relations cannot occur in one linear path. Concatenation recovers the original path. Together with L>=d+1 for a nontrivial endpoint-one G, this proves the two unsigned bounds without cancellation assumptions.

Both left forms check literally. In N, K_ac=U and K_ca=E can also be seen directly: the a-start path in the latter cover is forced to be singleton a, while deleting the last a in the former leaves an arbitrary c-reserved cover. In M, the only u-start path shapes avoiding c are u and u,b. They yield d+E, while the other cover orientation yields U−d. The fourth source coordinate is therefore essential and is not an independent free parameter.

The START_TWO block theorem places every marked C1 pair in exactly the required adjacent or first/last positions. Deleting its one or two degree-one vertices leaves an actual root-start-one source, including weak cases. The singleton-left M case is precisely the triangle redundancy and can be replaced by the nontrivial N source (1,0,2,1).

## Full-class reconstruction and graph injectivity

For every nonregular strong C2 graph, remove an outdegree-one vertex v and cut its two cycles after v. The remaining marked graph has start count2. Its initial SCC has start count2 and every later SCC is a singleton, because a nontrivial strong component has H>=3. Following the unique-entry START_TWO reductions must eventually reach its first two-entry reduction, with a strong C1 core and two start-one entries.

All earlier discarded roots point back only to their preceding root, which lies outside the final C1 core. All discarded SCC tails are dominated by that core, and v points only to the outer closing root. Therefore the only reverse cross edges into the core are those from the first pair-entry root r to its selected pair. The bundled right graph contains r and v, so it is nontrivial. The cycle product forces its ending count at r to be1. Its converse is thus an actual source in S*. This reconstructs the entire graph; no unproved suffix normalization or count-preserving rebasing is used.

The normal-form construction is well-defined on ordered rooted-source graph classes. It is surjective onto all nonregular C2 graph classes, and cannot produce the R5 exception: its sole order-five presentation has H13. The source class generating function W=B/(1−z) counts actual rooted classes with their unique weak tails. The presentation series (z+z²)(W−z)² equals the independently proved nonregular C2 series coefficient by coefficient. These are finite sets at every order, so the surjection is bijective. Numerical quadruple deduplication is not substituted for graph-class counting in this argument. Ordered sources remain essential even though both scalar forms are symmetric.

## Positive inverse bounds and scope

The coordinate change x=(D−E−d,U−D,E,d) is invertible and nonnegative on actual nontrivial sources, with last coordinate at least1 and coordinate sum U. Direct expansion gives exactly the two displayed symmetric matrices. Their nonnegative remainders after subtracting the stated all-ones/last-coordinate terms prove N>=3U1U2+1 and M>=6U1U2+U1+U2+1. Thus both source U values are at most floor((n−1)/6), since the other source has U>=2. Symmetric pair ordering and product pruning are safe.

The cap oracle plus the two forms and the explicit15 exception therefore gives a complete inverse for the full C2 scalar class. Every retained quadruple has a physical recipe and every accepted value has an actual graph. The polynomial claim is in the numerical target, not its bit length; source inventory, transition work and the B8 pair bound are all finite polynomial numerical work. Intermediate arithmetic has O(log n) bits. No graph-order census, Busch cutoff, saturation or formal inverse closure is hidden in the argument.

The graft root need not remain start-two: its start count is2d. The H41 example with start count6 is a valid and necessary scope hostile. The result does not prove that C1 union C2 is cofinite or has only the three stated exceptional scalar values.

## Independent controls

The standalone C++ checker imports no source code. It counts all small Hamilton paths, anchored cycles and reserved covers by literal vertex permutations and cuts. Its graph factorization is deliberately different: it exhausts every possible left vertex subset and right graft root, checks all cross arcs, tries both physical pair forms, and reconstructs every edge under an explicit map. It does not reuse the source SCC peeling algorithm.

Optimized and assertion-enabled unoptimized runs agree on:

- 3,744 labelled C2 graphs through order6;
- exactly one nontrivial ordered source presentation for each of the3,720 nonregular graphs, with24 labelled R5 graphs separately;
- 7,410 arbitrary right-gadget graft identities by literal path/cycle counting;
- 3,389 literal left-response controls from all start-one sources through order5;
- 103 literal full source quadruples from the cap40 oracle;
- an independent cap310 inventory of2,048 source quadruples;
- all4,192,256 unordered bilinear comparisons, without the source's product pruning;
- exactly the same884 target values through1861 as the frozen SUPPORT.json;
- all883 nonexceptional reconstructed witness graphs counted by a backward path/cycle dynamic program, maximum order14.

The independent source constant calculation scans the binary gadget word and uses direct reverse-complement bits, rather than the source's recursive L compiler. Its source recipes are physically reconstructed. The backwards large-graph counter differs from the source's forward extension counter; small graph identities and factorizations use literal permutations instead of either dynamic program.

A separate standalone algebra checker expands all32 matrix entries and their32 positive remainders. It independently generates rooted strong class coefficients from the canonical SCC recurrence, adds weak tails, convolves ordered source pairs and compares121 graph-class coefficients with the nonregular C2 rational-function recurrence. Initial nonregular counts are1,5,18,60,191,589,1775 at orders5–11.

The source normal/optimized replays reproduce the frozen receipt exactly. The finite controls corroborate different parts of the implementation; the cited complete classifications and the all-order reconstruction/bounds supply completeness beyond the tests.

## Reproduction and hashes

Compile check.cpp as C++17 and run it; run check_algebra.py with Python3 and its -O form. Neither independent checker reads external data or writes files. The independent JSON includes the full bounded support, so it can be compared directly with the source SUPPORT file. Source replay must use the copied source directory because it writes receipts beside itself.

- check.cpp: `c49b7effeb569f3155cb5cdc51c1c9d35abb7b154b6b7bd1ac764c38fe8026f8`.
- independent-normal.json: `2d1069599b6b968355aa3ee883f67c7c0fb07606befa34f144a34fee2197c7f8`.
- independent-debug.json: `2d1069599b6b968355aa3ee883f67c7c0fb07606befa34f144a34fee2197c7f8`.
- check_algebra.py: `516c0c1a534d656f0597a791360bcc8250631639ff74d5ae848765647ce721b2`.
- algebra-normal.json: `1089374786fe1709b6838aab7d3affd9cbcc0789f3bbcb77518c8f2b8dc8c579`.
- algebra-optimized.json: `1089374786fe1709b6838aab7d3affd9cbcc0789f3bbcb77518c8f2b8dc8c579`.
