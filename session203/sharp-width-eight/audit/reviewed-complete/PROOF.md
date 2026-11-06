# Complete rooted-profile occupancy through B=137

**FINITE-EXACT CANDIDATE, pending independent audit.** This packet determines every physically realizable pair(B,s) with1<=B<=137 and1<=s<=26, with no restriction on the original source order. It does not prove F(8,137), which quantifies over every odd base X>=137. The infinite tail remains OPEN.

## 1. Exact scope and published cutoff

For a nonempty rooted tournament (U,v), let h=H(U), s=start_v(U), and let L count ordered spanning two-path covers whose first path is nonempty and ends at v and whose second path may be empty. H(empty)=1. Define B=h+L-s. All profiles in this packet have s>0.

The physical edge extension E(U) adds a vertex a with v->a and a->U-v. The inherited exact identity is H(E(U))=B. If U is strong and nontrivial, E(U) is strong of order |U|+1. The rooted transfer and edge-extension conventions are frozen in dependencies/ROOTED_COVER_CLOSURE.md.

The inherited Busch minimum f(n) for nontrivial strong tournaments is

    f(3k)=3*5^(k-1), f(3k+1)=5^k, f(3k+2)=9*5^(k-1).

The cited source and its earlier source audit are preserved in dependencies/BUSCH_CORE_BOUND.md. In particular f(10)=125 and f(11)=225. The sequence is increasing, so every strong tournament with H<=137 has order at most10. No conjectured spectrum classification is used.

## 2. Complete strong generation with safe pruning

**Lemma 1.** Every strong tournament of order n>=4 has an induced strong subtournament of order n-1.

**Proof.** A strong tournament contains a directed triangle. Choose a proper strong vertex set S of maximum cardinality. If its complement is a singleton, the assertion holds. Otherwise no outside vertex is mixed toward S, since adjoining a mixed vertex would give a larger proper strong set. Partition the outside vertices into A, which dominate S, and D, which are dominated by S. Strong connectivity forces both parts to be nonempty, and forces an arc d->a for some d in D,a in A. Otherwise A has no incoming arc from its complement. Then S union {a,d} is strong, using a->S->d->a. Maximality forces this set to be the whole tournament, so the complement of S has exactly two vertices. Delete any one vertex x of S. The remaining S-x is nonempty, and a->(S-x)->d->a makes the entire remainder strong, whatever its internal orientation. QED.

**Lemma 2.** Hamilton-path count does not decrease upon adding a tournament vertex.

**Proof.** In each old Hamilton order, insert the new vertex just before the first old vertex it beats, or at the end if it beats none. If this first beaten vertex is not first, its predecessor beats the new vertex by the chosen first position. Thus the insertion is legal. Removing the new vertex recovers the original order, so the map is injective. QED.

Starting from the directed triangle, the enumerator therefore generates every strong tournament with H<=137 as follows:

1. At each order, extend every retained isomorphism representative by all mixed orientations of one new vertex. Mixed is necessary and sufficient for the extension of a strong graph to remain strong.
2. Count H exactly up to the overflow marker138 and discard only graphs with H>137. Lemma2 makes this pruning safe for every later extension.
3. Canonicalize exactly and retain one representative per isomorphism class.
4. Stop at order10 by the cited cutoff.

Lemma1 and induction prove completeness. The capped endpoint DP uses only nonnegative additions. Replacing every intermediate value x by min(x,138) retains the exact answer whenever H<=137 and otherwise records overflow, so capping cannot discard a qualifying graph.

Canonicalization is exact, not heuristic. Partition the vertices by outdegree in ascending order and exhaust every permutation within each degree cell. Encode the upper-triangle orientation for each such ordering and take the minimum integer code. Every isomorphism preserves the degree cells, so two graphs are isomorphic exactly when these minimum codes agree. Degree cells are merely an exact reduction of the permutation universe; no unproved color-refinement uniqueness assumption is used.

The resulting complete strong bank has the following class counts:

    order3: 1; order4: 1; order5: 6; order6: 35;
    order7: 321; order8: 362; order9: 147; order10: 11.

There are884 graphs in total. These counts concern H<=137, not all strong tournaments at those orders. BANK.tsv stores order, exact H, canonical code, and outgoing adjacency rows for each graph.

## 3. Recovering every strong-root profile

In E(U), the new vertex a has indegree exactly1, with unique predecessor v. Conversely, given any vertex a of indegree1 in a tournament T, delete a and mark its unique predecessor v. The remaining rooted graph U satisfies T=E(U) literally.

For every graph in the complete bank, inspect every indegree-one vertex a. Retain its deletion U only when U is strong, and compute s=start_v(U). If1<=s<=26, record the pair(H(T),s). The edge-extension identity is additionally checked by direct L counting.

Every nontrivial strong-root source with B<=137 occurs in this procedure because E(U) is strong, has count B, and has order at most10. Conversely every accepted source is physical and has the recorded pair. Isomorphism reduction loses no marked source: all vertices of every representative are inspected.

There are454 distinct strong-root profiles, obtained from660 indegree-one source occurrences. STRONG_PROFILES.tsv stores one literal root-zero graph per distinct pair. No claim of uniqueness of these source graphs is made.

## 4. Recovering weak sources without an order cutoff

If s>0, the root belongs to the initial strong component C. Hamilton paths visit strong components in dominance order. A path ending at the root cannot leave C. Consequently, if q is the Hamilton-count product of all later components, then

    B(U)=q B(C),  s(U)=q s(C).                 (1)

Delete every later singleton component; this preserves both coordinates. Since s>=q, every tail multiplier in the declared universe is at most26.

The attainable tail counts through26 are obtained exactly by multiplicative closure of the strong counts through26, together with1 for an empty tail:

    1,3,5,9,11,13,15,17,19,23,25.             (2)

The strong bank is complete for these factors as well. Any component of a tail count q has H dividing q and thus H<=26, so no factor can be missing. Conversely ordered joins realize every product used by the closure.

If C is a singleton, its pair is(1,1), giving(B,s)=(q,q). Otherwise its profile is in the strong-root list. Multiply each strong-root pair(B0,s0) by every q in(2), retaining qB0<=137 and qs0<=26, and add the singleton-root cases. This proves completeness for every original source order. The result contains exactly480 distinct pairs, listed in ALL_PROFILES.tsv.

Every listed pair also has an explicit literal root-zero witness in ALL_WITNESSES.tsv. compile_witnesses.py constructs these by ordered joins from the finite strong bank and recounts B and s with an independent append-at-end subset recurrence. These positive checks do not prove completeness; Sections2-4 do.

## 5. Exact staggered occupancy at the first base

At base X=137, the staggered profile requests are

    B=137-2j, s=r+3j, r in{1,2}, 0<=j<=8.

The complete profile list gives exactly

    r=1: j in{1,2,3,5};
    r=2: j in{0,3,4,7,8}.                     (3)

Every listed shift has a literal graph in WITNESSES_NORMAL.json. Every other shift in the declared range is absent at every possible source order. In particular (B,s)=(137,1) is absent, agreeing with the separately inherited complete unit-start grammar check. This is a cross-check rather than an assumption of the enumeration.

The positive existence of both residues at X=137 was already supplied by the earlier finite certificate bank. The new conclusion is complete occupancy, including negative answers for all other shifts and all other pairs in the declared box. This does not establish F(8,137) at arbitrary X, and does not turn the threshold23 into a theorem.

## 6. Reproduction and boundaries

Run from this packet directory:

    g++ -O3 -std=c++17 enumerate.cpp -o enumerate
    ./enumerate 137 26 > NORMAL.json 2> PROGRESS.txt
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 enumerate.cpp -o enumerate-debug
    ./enumerate-debug 137 26 > DEBUG.json 2> DEBUG_PROGRESS.txt
    cmp NORMAL.json DEBUG.json
    python3 compile_witnesses.py > WITNESSES_NORMAL.json
    python3 -O compile_witnesses.py > WITNESSES_OPTIMIZED.json
    cmp WITNESSES_NORMAL.json WITNESSES_OPTIMIZED.json

The enumeration writes BANK.tsv, STRONG_PROFILES.tsv and ALL_PROFILES.tsv deterministically; their exact-byte manifests must agree across builds. The witness compiler writes ALL_WITNESSES.tsv. The program rejects mismatched canonical count collisions and validates graph orientation, strong connectivity, canonical idempotence, and every edge-source identity. No removably asserted gate is used.

The declared computation is only(B<=137,s<=26). No absence conclusion is authorized outside that box. The proof of the finite universe and generation is all-order, but the resulting occupancy data are FINITE-EXACT, not an infinite-tail theorem.
