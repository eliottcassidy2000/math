# Independent audit: complete small rooted profiles and necessary width eight

**PASS: finite-exact classification and all-order exclusion of the eight stated profile requests.** This audit proves neither F(8,137) nor an eventual Hamilton-path spectrum threshold.

## Audited inputs and conclusion

The frozen sharp-width proof is `reviewed-sharp/PROOF.md`, SHA256 `979737422f9d388a985d03089c7cb826293842d1558fc03383f69fafaa8be834`. The finite-box proof is `reviewed-complete/PROOF.md`, SHA256 `187120e24014800629c27f9cfa21cbc581f7abb21ee7bf93134c99b0d36ddea7`. All files from each producer's manifest are preserved under these two directories. `INPUTS.json` identifies the original source and snapshot of every entry.

For a nonempty rooted tournament, let h count Hamilton paths, s count those starting at the root, and L count ordered spanning path covers whose first path is nonempty and ends at the root, while the second may be empty. Throughout this audit, H(empty)=1 and B=h+L-s.

The complete profile box 1<=B<=137, 1<=s<=26 has exactly 480 realizable pairs, at unrestricted original source order. Its strong-source part has 454 pairs. At X=137, the successful shifts among 0<=j<=8 are exactly {1,2,3,5} for residue 1 and {0,3,4,7,8} for residue 2.

At X=153 and residue 1, these eight pairs are absent at every source order:

    (153,1), (151,4), (149,7), (147,10),
    (145,13), (143,16), (141,19), (139,22).

The next shift is physically realized: the root-zero rows

    [222,128,2,6,14,31,62,124]

give (h,s,L)=(77,25,85), hence (B,s)=(137,25). The actual four-core readout has 461 Hamilton paths, equal to 3*153+2.

Consequently F(7,137) is false, and a true F(J,137) must have J>=8. This is a sharp necessary width at this particular base and residue. It does not show that width eight works at every base. It is not an exclusion of the target 461 from the Hamilton-path spectrum; the explicit repair realizes that target.

## Why the negative claims have complete scope

1. For s>0 the root is in the initial strongly connected component C. If q is the Hamilton-count product of the later components, then h, s, and L all factor by q. Therefore (B,s)=q(B(C),s(C)). A later singleton may be removed without changing either coordinate. If C is a singleton, B=s=q. A nontrivial strong component has at least three Hamilton paths.

2. Each of the eight negative requests has B>s and gcd(B,s)=1. Thus q=1, all later components are singletons, and removing them leaves a nontrivial strong source with the same pair. Arbitrarily long weak sources cannot evade the negative certificate.

3. Add a vertex a with root->a and a->every other old vertex. The physical edge extension E(U) has H(E(U))=h+L-s: paths starting at a account for h-s, and all other paths split at the forced root->a edge into exactly the L covers. The optional-empty second path is included. If U is nontrivial strong, E(U) is strong and has one more vertex.

4. The cited Busch minimum f(n) has f(10)=125, f(11)=225 and is increasing thereafter. Consequently H(E(U))<=153 implies |E(U)|<=10 and |U|<=9. This cutoff is a cited input, not a theorem proved by the finite computation. Its frozen source is `reviewed-sharp/dependencies/BUSCH_CORE_BOUND.md`, from [the published Session197 corollary](https://github.com/eliottcassidy2000/math/blob/17c3b0bb5cfd820737fc39449cda2fd3df14df53/session197/unrestricted-target/PROOF.md). That source identifies Busch's 2006 paper, Theorem 1 and Corollary 1, and the primary-source URL.

5. Every strong tournament of order n>=4 has a strong induced subtournament of order n-1. The proof using a maximal proper strong set is valid. If its complement has at least two vertices, all outside vertices are uniform toward it, and strong connectivity supplies a pair a,d with a->S->d->a. Maximality forces precisely those two outside vertices. Deleting any vertex of S leaves a strong tournament. The version preserving a specified root also holds: start with a triangle containing the root and delete a different vertex of S.

6. Adding a vertex cannot decrease H. In an old Hamilton order, insert it immediately before the first vertex it beats, or append it if it beats none. Deletion recovers the input order. Thus pruning H>153 cannot hide any later qualifying strong graph. Mixed one-vertex extensions are exactly the strong extensions of a strong parent.

These facts reduce the problem to a complete, finite strong bank, not a heuristic census or a search inside a special construction family.

## Independent bank and profile verification

`check_closure.cpp` does not import the producer's enumeration, canonicalizer, or subset-DP counting routine.

- It validates and literally counts each archived graph by depth-first enumeration of directed paths, with an early stop only after finding 154 Hamilton paths. Qualifying graphs retain their exact count.
- It proves pairwise nonisomorphism by exact adjacency-preserving backtracking. Outdegrees and the degree histogram of out-neighbors only restrict possible maps. A reported isomorphism always comes from a full vertex bijection. Rooted comparisons additionally fix vertex zero.
- It checks the triangle base and every mixed one-vertex extension of every retained graph of orders 3 through 9. Every extension with H<=153 must match an archived graph by this exact isomorphism test. Combined with the deletion lemma and monotonicity, this closure certificate proves the bank complete independently of the producer's canonical codes.
- It directly computes h, s and L by enumerating all directed paths on all nonempty subsets. It inspects every indegree-one vertex in the bank and inverts the physical edge extension. All strong sources with s<=26 must occur this way.
- It recounts all 535 stored strong-source witnesses, all 792 graphs in the secondary rooted bank, and all 480 explicit full-profile witnesses through B=137. Each extracted source also matches a secondary-bank graph by exact root-preserving isomorphism.
- It obtains every tail multiplier through 26 by multiplicative closure of the complete strong-count list. The resulting set is {1,3,5,9,11,13,15,17,19,23,25}. Component factorization then proves completeness at unrestricted original source order. It independently matches the 564 full profiles through B=153 and their exact 480-pair restriction through B=137.
- It checks all eight gcd-one absences, the physical shift-eight repair, and the literal 461-path readout graph.

Both optimized and unoptimized `_GLIBCXX_ASSERTIONS` builds give identical receipts:

- 1,113 strong graphs of orders 3 through 10, with class counts 1,1,6,35,344,471,228,27.
- 281,628 mixed-extension trials and 6,204 qualifying extensions.
- 10,141 exact isomorphism calls and 64,520 search nodes.
- 593,963,009 literal directed-prefix visits.
- 535 strong profiles from 792 indegree-one source occurrences.
- The exact B<=137 subbank has 884 graphs, class counts 1,1,6,35,321,362,147,11, and 454 strong profiles from 660 source occurrences.
- All 480 explicit B<=137 witnesses pass; their maximum order is 10.

The producer's small Python sharp-width check was also replayed in normal and optimized modes from the frozen snapshot. Both outputs match its original receipt exactly. The full producer enumeration was not redundantly rerun for this audit: the independently implemented complete closure check supplies the independent exhaustive verification.

## Source-code review and boundary

The producer's primary canonicalization is mathematically exact: order degree cells, exhaust their internal permutations, and minimize the complete orientation code. The secondary canonicalization similarly fixes the root and exhausts degree-and-root-incidence cells. The capped primary DP uses only nonnegative additions, so its overflow marker preserves all qualifying values. The secondary cut-matrix insertion formulas correctly separate the cases in which the new vertex is first, last, or internal. Its root-start formula excludes the new-first case. Root-preserving mixed insertion makes both B and s monotone, so its pruning is safe.

The primary encoder uses a 64-bit upper-triangle code. The audited runs have order at most 10 and at most 45 encoded edges, which is safe. Its original broad command-line bound range could permit order 12 or larger, where that representation is unsafe. Therefore this audit certifies the declared 137/26 and 153/26 runs only, not arbitrary command-line bounds. A separate release hardening guard may restrict ncap<=11 and positive endpoint caps, provided the original is preserved and the delta is checked. This non-load-bearing limitation does not affect either frozen proof or computation.

No absence outside the audited profile boxes is claimed. The eight selected negative pairs do have unrestricted source-order scope because their finite reduction is proved. F(8,137), arbitrary-base coverage, and the eventual threshold 23 remain open.

## Reproduction

Run from this audit directory:

    g++ -O3 -Wall -Wextra -std=c++17 check_closure.cpp -o check_closure
    python3 run_check.py ./check_closure NORMAL
    g++ -O0 -D_GLIBCXX_ASSERTIONS -Wall -Wextra -std=c++17 check_closure.cpp -o check_closure_debug
    python3 run_check.py ./check_closure_debug DEBUG
    cmp NORMAL.json DEBUG.json
    python3 reviewed-sharp/check_sharp.py > PRODUCER_NORMAL.json
    python3 -O reviewed-sharp/check_sharp.py > PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json PRODUCER_OPTIMIZED.json
    cmp PRODUCER_NORMAL.json reviewed-sharp/SHARP_NORMAL.json

Compiler: Debian g++ 14.2.0, x86_64 Linux. The optimized checker took approximately 1.35 seconds and the debug checker 7.10 seconds in this run, each with recorded peak child RSS 8,832 KiB. These timings are observations, not performance guarantees. Resource receipts retain the actual values; the host CPU model was not exposed by lscpu.

Independent checker SHA256: `db534ad2989d7e62046192f9ab975b39dd864bc003ecc8c3b593aaad2b1dd231`.

Identical normal/debug receipt SHA256: `09ca2f5de895a984fb17fec6a271ee17cedcff1677bad0961166942a485ac77c`.

The outer `SHA256SUMS` fixes the reviewed inputs, independent source, and receipts. It deliberately excludes local executable binaries and obsolete failed timing attempts. No shared checkout or producer input was edited during this audit.
