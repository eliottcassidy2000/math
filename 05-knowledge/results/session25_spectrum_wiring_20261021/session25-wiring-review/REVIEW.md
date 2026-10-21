# Session25 wiring and labelled holonomy: independent audit

Verdict: APPROVED for the precise claims in the frozen producer packet. No required mathematical revisions. This review approves an exact finite wiring certificate, elementary all-order transfer/group statements, and the stated all-order labelled odd-girth family. It does not approve a merged-metagraph odd-girth theorem, all-order shape classification, connectivity theorem, or non-leaf-bridge exclusion.

## Frozen producer target

Producer directory: `/workspace/shared/session25-wiring`.

- RESULT.md: 15a2d1b6320acd163ac34081e5c3cf92a9a22b6aaef7ce06e2f3a17b753898a4
- check.py: c0d9978d55edbcc373e7c98ad26454d30c56c0f510a95af088bce5bd755c9423
- check.out: 52366b9c653da21e57f72fccf4e44d1fafd3d2cfeae8ac3d1b297d5b080afdea
- holonomy.py: db7ab82c361122368388ce7315e072977ac5f15c163474f8afb3f3541c304d49
- holonomy.out: 869c0838653e5f10be5bebacacf424852fef125f3dc9905a012091a95e439d9d

These five match the producer SHA256SUMS. Both scripts were rerun without modifying them; both outputs match the frozen outputs byte-for-byte. The repository was only read.

## Mathematical audit

1. Agreement iff is exact, not merely necessary. Pulling B back by a bijection separates exactly the retained path pairs from all flipped pairs. With n-1 arcs and both directed degrees at most one, acyclicity gives one path component; weak connectivity gives the same. The C3 plus isolate hostile invalidates the degree/count-only shortcut.

2. The automorphism divisor is exactly a_A a_B: each source path has a_B target identifications, then the source automorphism group acts freely on full ordered paths. It is not legitimate to divide by an unrelated module-action size. Blue means the actual path reversal is an anti-automorphism. In pinned coordinates complement preserves blue, is fixed-point-free for n>=3, and pairs diagonal incidences, giving the factor two on loops. The mirror action pairs black lines, including loops, after converse merging. Thus off-diagonal black division by two and diagonal black division by four are correct. The definitions and factors agree with THM-4729 section 0 and THM-4731 statement 3, inspected directly.

3. Maximal module runs give a unique word and ordered internal path covers. Every path edge belongs to exactly one internal run or one specific boundary port, so global converse followed by restoring those arcs is literally nonpath complement. Scalar cover counts cannot determine the port corrections. For block-contiguous paths there are precisely |B_i||B_j|-1 incorrect cross pairs at each adjacent block boundary in the naive uniform substitution, proving the stated named-decomposition boundary. The target need not remain a uniform substitution; the assertion is an exact marked construction, not closure of the substitution class.

4. The six signed triangle frames are a regular S3 carrier, not C6. Cross-block arcs impose one global sign on the quotient action; each internal triangle permutation must have that same sign. A single independent local reflection is invalid. Intrinsic recovery is sound: a cyclic triple meeting only two blocks cannot be cyclic, while a cyclic triple meeting three blocks is not a module because an unused vertex in one of its blocks distinguishes the other two triple vertices. The given bottom blocks are therefore exactly the cyclic three-vertex modules. No general module-recovery conjecture is needed.

5. The labelled conservation law is correct. For odd closure every pair occurs positively and oddly among path edges, so k(n-1)>=binom(n,2). The complement flips (n-1)(n-2)/2 pairs, odd exactly for n=0,3 modulo4, giving bipartiteness. This assumes n>=3 as the packet explicitly does; n=1,2 would have fixed complements and are excluded.

6. For n=2k, odd k>=3, the zigzag path for index i exhausts the matchings with endpoint sums 2i-1 and 2i. All residues occur once among these paths, so their edge sets partition K_n. Before its unique occurrence an edge is flipped exactly i times; alternating the initial direction by i therefore makes every requested path directed. Each edge is flipped k-1 times in total. A proper segment cannot close: an even-length segment flips exactly its nonempty union of disjoint path sets; an odd-length segment flips their nonempty complement. The constructed cycle is simple and has length k. Consequently odd girth exactly n/2 is proved for n=4r+2>=6.

7. Signed quotient holonomy must retain the relabeling disagreement vector and the converse bit. Ordinary relabeling closure gives sum p=(k mod2)1+Delta_g, and for odd k the support inequality follows because edges outside Delta_g must appear at least once. Converse closure changes the all-pairs coefficient to k+delta. None of these implies the same odd-girth bound after ordinary or converse-merged quotienting. The packet properly states this limitation.

8. At n=4r+1, k=2r+1, endpoint incidence gives 2k-e_v=(n-1)+(2-e_v). Each edge's excess multiplicity is even, so e_v is even and at most two. Since sum e_v=2k, exactly k vertices have e_v=2; the remaining k-1=2r vertices each have excess two. These excesses can only be one tripled edge at each such vertex, so tripled edges form a matching of size r, and multiplicity five is impossible. For an edge at step i with j previous occurrences, exactly i-j earlier moves flip it. The displayed orientation constraints are therefore necessary and sufficient. When consistent they construct a valid closed walk; any repeated vertex would split off a shorter odd closed walk, contradicting the already proved lower bound. The converse is valid.

## Reproduction and independent controls

- `check-rerun.out`: full producer run; class-pair bijections 24, 384, 17,280, 2,257,920 for n=3,4,5,6; all ordered pairs and colors match pinned enumeration. All 2,702 selected substituted paths pass port restoration. Signed two-block group has 18 of 72 block permutations. This bounded census does not certify an all-order bridge law.
- `holonomy-rerun.out`: producer directed witnesses n=6,10,14,18,22, n5 triangle, exhaustive parity checks n3,4.
- `independent.py` / `independent.out`: reviewer implementation uses directed edge sets rather than the producer bit-vector machinery and imports no producer or repository code. It checks the root's distinct n5 words (3,4,1,2,5), (4,5,1,2,3), (5,3,1,2,4), with tripled edge {1,2}; all are compatible and form a directed labelled triangle. It checks the Walecki family at every n=6,10,...,102. Finally it builds every labelled n5 tournament and every valid path move: L5 has 1,024 vertices, 3,840 edges and exactly 80 triangles. Every triangle satisfies the tripled-matching and endpoint-incidence anatomy.

## Scope cautions retained

No literature-priority finding was performed or claimed. The finite checks support, but do not replace, the proofs above. Exact adjacency simulation is not an all-order graph-shape classification. The n=1 modulo4 general sharpness question remains unresolved beyond the explicit n5 and n9 witnesses. The labelled and merged graphs must remain distinctly named in any downstream write-up.

## n9 amendment audited before final approval

The producer added a finite n9 witness during review. The final hashes above supersede the initial five-file snapshot. The new holonomy.py output was rerun and matched exactly. The independent directed-set verifier also reconstructed every n9 state, checked all five Hamilton paths and complement transitions, verified the repeated matching {01,23}, and matched the five distinct integer states. Together with the proved lower bound, this certifies odd girth exactly five for L9.

The optional discovery probe is not a trusted proof dependency. Its source was inspected; the MILP proposes an undirected decomposition, then searches path order and reversals. Its archived witness, not solver optimality, is certified by the integer and independent directed-set verifiers. The discovery MILP was not rerun. Frozen optional files:

- probe_n9.py: 966ec5da667df8f06c04a2eec6e4c555d3c2e92d6064633d404ece4de8cc2e2b
- probe_n9.out: fe7f2b9962bc5622b2defc934ad8e27544765ff4c9110e8c12cdc42a05d299f3
