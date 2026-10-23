# Independent P7 two-Fano path audit

Status: PASS. Checker `check_fano_paths.py` is independent of the producer's checker; Hamiltonian path counts use subset dynamic programming. All tests use explicit RuntimeError checks and give byte-identical ordinary/optimized outputs. Evidence: `fano_paths_receipt.json`, `fano_paths_check.out`. No Git operations.

## Definitions

P7 has vertices Z/7 and arc a→b iff b−a belongs to {1,2,4}. H counts directed Hamiltonian paths. For a directed Hamiltonian path P=(p0,...,p6), call it blue iff reversing the positions gives a tournament anti-automorphism, equivalently the tournament in this path's ordered coordinates equals its reversed converse. Black means not blue. The tests use this full adjacency condition, not the proposed Fano criterion, to determine color.

Let L+ be translates of {0,1,3}, and L− translates of its negative. Let u± count lines containing no path edge. Every such unused line can have all three arcs reversed without losing P.

## Findings

Independently enumerating all 5040 vertex orders finds exactly 189 directed Hamiltonian paths, partitioned as follows:

- (u+,u−)=(1,1), blue: 21 paths.
- (2,2), blue: 42 paths.
- (2,3), black: 63 paths.
- (3,2), black: 63 paths.

The 14 cyclic triples are precisely the disjoint union of the two seven-line Fano systems. Therefore the number of admissible single cyclic-triangle reversals is u++u−. Black paths have exactly five; blue paths have two or four. Each admissible reversal preserves the fixed path and regularity and gives H=171. The original P7 has H=189.

The full automorphism group independently enumerated from all vertex permutations has 21 members and consists exactly of affine maps x↦ax+b with a∈{1,2,4}. Every directed Hamiltonian path orbit is free of size 21. Normalizing p0=0 and p1=1 gives nine representatives, explicitly recorded with their increments and profiles in the JSON receipt.

## Structural explanation of the increment certificate

Write di=p(i+1)−pi∈{1,2,4}. Two path edges on the same Fano line must be adjacent: nonadjacent edges have four distinct vertices, whereas a line has three. A line cannot contain three path edges. Thus if r± counts consecutive edge-pairs lying on a line of L±, the six edges occupy 6−r± distinct lines, giving u±=1+r±.

For a consecutive pair with ratio d(i+1)/di equal to 2, its triple is a translate and QR dilation of {0,1,3}, so belongs to L+. Ratio 4 gives a translate and QR dilation of {0,1,5}, a line of L−. Ratio 1 gives neither. Hence exactly

 u+=1+#{i: d(i+1)/di=2},
 u−=1+#{i: d(i+1)/di=4}.

These ordered, not merely unordered, equalities are independently checked for every path.

Negation is an anti-automorphism. Composing it with every automorphism lists all anti-automorphisms as x↦−ax+b, a∈{1,2,4}. The position reversal is an involution; its slope must square to one, so its slope is −1. Therefore blue means p(6−i)=p0+p6−pi for every i. This is exactly palindromicity of the six increments. The nine normalized increment sequences are:

111111; 112114; 112242; 114411; 121144; 122422; 144112; 144224; 144441.

Only the first, fourth and ninth are palindromic. Their ratio counts give profiles (1,1),(2,2),(2,2); each remaining sequence gives (2,3) or (3,2). This nine-row certificate proves the balance/color equivalence, with all 189 paths covered by the affine action. The complete vertex paths appear in the receipt, making the certificate directly reproducible.

## Unused-line subcubes and the essential restriction

For either single Fano system, its seven triangle blocks partition all edges. Reversals on a set of its unused lines commute. Every resulting orientation keeps P and is regular, because each line is cyclic and contributes one outgoing and one incoming arc at each of its vertices, before and after reversal. The independent test checks all 1932 (path, plane, subset) instances; every one of the 1554 nonempty instances has H=171. Only the empty subset gives the original H=189 P7.

This is not true for the unrestricted full seven-line cube. Independently testing both full cubes finds 16 H=189 masks and 112 H=171 masks in each. The 16 exceptional masks are exactly: empty, full, the seven three-line pencils through a point, and the seven complements of pencils. Thus the full-cube finite classification gives a compact structural explanation for the path-subcube result: u±≤3 excludes all full and four-line masks, and a three-line pencil cannot be entirely unused, because the Hamiltonian path has at least one edge incident with its central point, lying on one of its pencil lines. Only the empty exceptional mask can remain.

The displayed exceptional-mask classification is verified by exact enumeration of 128 masks per plane, not asserted from an unproved all-orders principle. H=171 distinguishes every nonempty unused-line state from P7; the origin itself is literally P7. No isomorphism classification of other tournament orders is needed.

## Hostile regular order-seven tournament

For outgoing-neighbor bitmasks [98,84,41,19,37,74,28], the checker independently confirms regularity and H=171. Both of the following are directed Hamiltonian paths and each admits exactly five path-preserving cyclic-triangle reversals:

- (0,1,2,3,4,5,6): blue by full adjacency reversal test.
- (0,1,4,2,5,6,3): black by the same test.

Therefore “black iff five admissible cyclic-triangle reversals” is a P7 theorem, not a theorem for every regular tournament on seven vertices. No all-orders generalization is claimed. The earlier general Singer audit was about c3; this separate finite P7 result genuinely computes H and links the two-Fano incidence structure to the original path color.
