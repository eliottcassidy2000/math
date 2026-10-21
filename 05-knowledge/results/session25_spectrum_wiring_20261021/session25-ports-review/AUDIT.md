# Independent audit: named boundary signatures

## Verdict and frozen files

APPROVED for the mathematical claims of sections 1–3 and the stated abstract-versus-geometric / counting-versus-support scope. This is not an audit of the separately cited two-port theorem, tournament packet, or A.9 recursion.

Frozen SHA-256:

- PORT_SIGNATURES.md: b5d64e5063c1f918ef6c18afd51524e0cc04268f6d88a9a6bc5df76303e8d686
- verify_ports.py: 3091bb365ceef16927c0618e0ff48f8553cfaa6823054ae7f9198a73c427f875

No mathematical revision is required. Two editorial improvements: explicitly specify finite loopless graphs (parallel edges may be suppressed), and change the original checker description “verifies their number” to “reports their number” unless including the independent Stirling recurrence check delivered here. Original checker correctly enumerates restricted-growth words and checks full rank, but its own code has no separate count oracle.

## Proof audit

A boundary assignment has an equality partition pi and exactly (k)_b assignments realize each b-block partition. Colour permutation acts transitively on those assignments, proving f is well-defined. A colouring of a glued union corresponds bijectively to two extensions with identical boundary assignment. The hypothesis that the interiors are disjoint is essential and is explicitly present; overlapping boundary edges are harmless because the same inequality can be tested twice without changing validity.

For necessity, equality of counts gives weighted orthogonality of the signature difference to every boundary-only graph indicator u_E. For each sigma, the product of 1-d_ij over same-block pairs and d_ij over cross-block pairs is exactly its indicator. Expanding gives an integral linear combination of u_E, so the sigma coordinate of the difference times (k)_b is zero. The factor is positive for every allowed partition, proving equality coordinatewise. This does not assert existence of a single graph implementing the signed sum.

The minimal-state result is about colour assignments under all fixed-name continuations, not arbitrary encodings of fragment vectors. Two different equality partitions disagree on at least one pair; the single edge on that pair distinguishes their feasibility as well as their fixed-assignment extension count. Equal partitions are related by a colour permutation. Hence the number of states is sum_{b=0}^{min(r,k)} S(r,b), with the usual S(0,0)=1. This is not a port-permutation quotient.

At three ports, exchanging ports 2 and 3 transforms edge12 into edge13, but keeping the edge12 continuation fixed gives one edge versus a two-edge star. Their counts are k^2(k-1) and k(k-1)^2; their difference is k(k-1). It is positive for k>=2, and zero at k=1 as stated by the example's scope. Permuting every operand together preserves gluing; independently permuting one requires matching interface information. For two ports all equality partitions are fixed by S2. At k=1 there is only one equality state for every r, so no S3 loss is claimed.

For r=0 the sole partition and sole empty assignment give f_G=Col_k(G), with (k)_0=1. The product formula is ordinary disjoint-union multiplicativity. For k=1 any graph containing an edge has zero extensions; edgeless fragments have one. If loops are included by convention, every looped fragment has zero signature and all arguments remain valid. Source looplessness is thus a convention issue, not a counterexample. Geometry compatibility and restricted realizable continuations are correctly excluded from the universal minimality claim. Boolean support suffices for feasibility, while actual counts are needed for counting.

## Primary source verification

Directly read https://arxiv.org/pdf/math/0408232v1 (12-page PDF, arXiv v1 18 August 2004).

Theorem 2.2 is on printed page 5 (PDF index 4), Lemma 2.4 on printed page 6 (PDF index 5), so the packet's pages 5–6 reference is correct. The theorem is for twin-free weighted targets, with positive vertex weights and real symmetric edge weights, and equates connection rank to the number of automorphism orbits on ordered labelled tuples. Source inputs are finite loopless multigraphs; gluing identifies equally named vertices. Its target convention represents absent edges and loops by weight zero.

Specialization: K_k has equal positive vertex weights, off-diagonal adjacency weight 1 and diagonal 0. Distinct vertices have different adjacency rows (including at their own positions); K_1 is vacuously twin-free. Automorphisms are the full colour symmetric group and tuple orbits are precisely equality partitions. Optional vertex-weight normalization scales individual graph homomorphism values without changing the rank or orbit conclusion. Lemma 2.4 identifies indistinguishable boundary maps up to target automorphism. The paper's source notation k counts boundary labels; the packet instead uses r for labels and k for colours. This is a valid relabelling.

The source is appropriately cited as prior general context, not as a novelty claim or substitute for the stronger elementary boundary-only spanning proof. Version-pinning the displayed link to /pdf/math/0408232v1 would improve reproducibility, but the inspected statements match.

## Independent verification

Run from the workspace root:

    python3 session25-ports-review/independent_check.py
    python3 -O session25-ports-review/independent_check.py

Both passed with byte-identical JSON receipts. The checker was written independently and uses union-find inclusion–exclusion for complete chromatic counts instead of the original checker’s full named-colouring enumeration. Restricted-growth partitions are generated recursively and counted independently against a Stirling recurrence.

Coverage:

- 36 partition-count controls, r=0..5, k=1..6.
- 11,033 direct signed expansion identities, including empty products.
- 10,704 pairwise distinct-state single-edge distinctions.
- 16,660 gluing pairs: all simple graphs on r named ports plus one private vertex for each operand, r=0..3, k=1..4. Interiors are mapped to different vertices. This includes shared boundary edges, impossible boundary assignments, empty boundary, one-colour cases, and disconnected fragments.
- Closed-form hostile examples k=1..8.
- Explicit loop-zero and repeated-edge invariance controls k=1..4.

The original checker’s exact rational elimination, falling factorial, private-vertex remapping, and RuntimeError checks were inspected and are correct. Its checks survive Python -O because they do not rely on assert. Its finite ranges do not prove the arbitrary-size theorem; the preceding proof does.

Original checker reproduction completed in both normal and -O modes: 12,322 checks PASS, outputs byte-identical to each other and the supplied receipt.json. Saved as original_receipt.json and original_receipt_optimized.json. The frozen source hashes were rechecked unchanged immediately before the optimized reproduction.

## Final editorial revision approval

APPROVED revised PORT_SIGNATURES.md SHA-256: 0d535d3d8500c498a18d9de136940abf1407603eb6265e98af0be03461d5d04e. Reversing exactly the two recommended wording substitutions reproduces the originally audited hash. Thus the only changes are “finite graph” to “finite loopless graph” and “verifies their number” to “reports their number”. Both recommendations are resolved. No mathematical claims or checker code changed; the earlier proof and receipts apply unchanged. The checker remains SHA-256 3091bb365ceef16927c0618e0ff48f8553cfaa6823054ae7f9198a73c427f875.
