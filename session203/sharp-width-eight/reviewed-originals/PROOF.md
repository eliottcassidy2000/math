# Width eight is necessary at the staggered base153

**FINITE-EXACT PROOF CANDIDATE, pending independent audit.** The complete bounded enumeration below proves that F(7,137) is false. It does not prove or disprove F(8,137), and it does not prove an eventual Hamilton-path spectrum threshold.

## 1. Exact statement

For a nonempty rooted tournament (U,v), let h=H(U), s=start_v(U), and let L count ordered spanning two-path covers with a nonempty first path ending at v and a second path allowed empty. H(empty)=1. Define B=h+L-s.

The previously published condition F(J,B0) asks: for every odd X>=B0 and each r in{1,2}, is there an integer0<=j<=J and an actual rooted source with

    B=X-2j,  s=r+3j?

**Theorem.** At X=153,r=1, no shift j in{0,...,7} is physically realizable, at any source order. Shift j=8 is realizable. Therefore F(7,137) is false, and any true F(J,137) must have J>=8.

The absent pairs are exactly the following eight requests:

    (B,s)=(153,1),(151,4),(149,7),(147,10),
          (145,13),(143,16),(141,19),(139,22).   (1)

The repaired pair is(B,s)=(137,25). An explicit root-zero witness has outgoing rows

    (222,128,2,6,14,31,62,124),                 (2)

with h=77, s=25, L=85, hence B=77+85-25=137=153-2*8. The inherited actual four-core construction has count3B+2s=461=3*153+2; this positive target count is an interpretation of the repair, not a new spectrum exclusion.

The width-eight conclusion is a necessary boundary for the stated sufficient hypothesis. It is not a disproof of the eventual spectrum conjecture, and it does not show every base accepts width eight.

## 2. Why weak components cannot hide a missing request

For s>0 the root lies in the initial strong component C. Let q be the product of the Hamilton-path counts of all later components. Component concatenation, and the fact that a path ending at the root stays inside C, give

    B(U)=q B(C),  s(U)=q s(C).                 (3)

Later singleton components can be deleted without changing either coordinate. If C is a singleton, then B=s=q. For every pair in(1), B>s and gcd(B,s)=1. Thus q=1, every later component is a singleton, and deleting those components leaves a nontrivial strong rooted source with exactly the same pair. A nontrivial strong component has at least three Hamilton paths and cannot contribute to q=1.

Consequently excluding strong rooted sources for(1) excludes all source orders and all weak decompositions. No guessed total-order cutoff is used for the weak case.

## 3. The actual edge extension gives a nine-vertex source cap

Add a vertex a with v->a and a->U-v. Call this physical tournament E(U). Its path split gives

    H(E(U))=h+L-s=B.                          (4)

Paths beginning at a correspond to the h-s old paths not beginning at v; every other path has v immediately before a and corresponds to an L-cover. This proves(4) with the declared empty-path convention.

If U is nontrivial strong, E(U) is strong. The new vertex reaches U-v directly and v through U, while every old vertex reaches v and then a. Its order is |U|+1.

Busch's published sharp minimum for Hamilton paths in a strong n-tournament is

    f(3k)=3*5^(k-1), f(3k+1)=5^k, f(3k+2)=9*5^(k-1).

The verified primary source and the earlier source-audited corollary are frozen in dependencies/BUSCH_CORE_BOUND.md. In particular f(10)=125 and f(11)=225, and f increases thereafter. Since every B in(1) is at most153, any corresponding E(U) has order at most10, and U has order at most9.

This is a cited-bound reduction, not a new proof of Busch's theorem. It makes the negative requests finite and exhaustive.

## 4. Complete enumeration of the edge extensions

The primary enumeration generates all strong tournaments T with H(T)<=153, through order10, without assuming a preexisting tournament census.

**Strong deletion lemma.** Every strong tournament on n>=4 vertices has an induced strong tournament on n-1 vertices.

To prove it, choose a proper strong set S of maximum cardinality; a directed triangle supplies such a set. If at least two vertices remain outside, no outside vertex is mixed toward S, since that would enlarge S while leaving it proper. Partition the outside vertices into A dominating S and D dominated by S. Strong connectivity forces A,D nonempty and an arc d->a from D to A. Then S union{a,d} is strong, so maximality forces it to be the whole graph. There are exactly two outside vertices. Delete any x in S; the remaining nonempty S-x, together with a,d, is still strong by a->(S-x)->d->a. This proves the lemma.

Starting from the directed triangle, extend every retained strong representative by all mixed orientations of one new vertex. Mixedness is sufficient as well as necessary for an extension of a strong graph to stay strong. The deletion lemma proves completeness inductively.

Prune only H>153. This is safe because every old Hamilton order extends injectively: insert the new vertex before the first old vertex it beats, or append it if there is no such vertex. The old order is recovered by deleting the new vertex. The capped subset DP uses only nonnegative additions and overflow154, which preserves every count at most153 exactly.

Canonicalization is exhaustive within degree classes. Put outdegree cells in ascending order, try every permutation inside each cell, encode the upper-triangle orientation, and take the least code. Every isomorphism preserves the degree cells; equality of these minima is therefore equivalent to isomorphism. No heuristic color refinement or assumed rigidity is used.

The complete class counts under H<=153 are

    order3:1, order4:1, order5:6, order6:35,
    order7:344, order8:471, order9:228, order10:27.

There are1,113 actual graphs. run153/BANK.tsv stores the order, exact H, canonical code, and every outgoing row.

For each graph T, inspect every vertex a of indegree1. Delete a, mark its unique predecessor v, and keep the source only if it is strong. Then T is literally E(U), so its H is B(U). Recount s and L directly. Every strong-root request in(1) must occur in this search by Section3, and every accepted profile is physical. Inspecting every vertex retains all marked possibilities after unrooted canonicalization.

This yields535 distinct strong-root pairs with B<=153 and s<=26, from792 source occurrences. None of(1) is in that complete set. The positive pair(137,25) is present and is separately verified from(2). This proves the theorem.

## 5. A second complete route and positive controls

The second program, rooted_scan.cpp, works directly with rooted strong graphs rather than their unrooted edge extensions. It reproduces the same535-pair set exactly.

The deletion lemma can preserve a designated root: start with a directed triangle containing that root, choose S containing it, and in the two-outside case delete a different vertex of S. Every vertex of a strong tournament belongs to a directed triangle, by shortcutting a shortest directed cycle through it. Hence all rooted strong sources arise by root-preserving mixed extensions from the rooted triangle.

Both B and s are monotone under these extensions. The graph E(U) remains an induced subtournament of E(U+x), so B is monotone by path insertion. For s, an old Hamilton order starts at the root; mixedness supplies an old vertex beating x. If its last vertex beats x, append x. Otherwise the sequence of signs along the old order has a transition from beating x to losing to x, and x can be inserted at that internal gap without changing the first root. This gives an injection preserving starts. Thus pruning B>153 or s>26 is safe.

The second program uses exact two-path-cut matrices to evaluate all new-vertex orientations. If I consists of old vertices beating x and O of vertices beaten by x, an ordinary Hamilton count after insertion is

    sum_(i in I) end_i + sum_(j in O) start_j
       + sum_(i in I,j in O) C_ij,

where C_ij counts ordered spanning pairs of nonempty paths ending at i and starting at j. Each Hamilton path splits uniquely around x. For the new B count this formula is applied to E(U), fixing its auxiliary vertex in I. For s, the prefix path must begin at the old root; the same cut formula uses root-started prefix counts and omits the case x first. These quantities are computed from exact subset tables.

Rooted canonicalization fixes the root and exhausts every permutation in the other vertices' degree-and-root-incidence cells. It is exact for the same reason as the primary canonicalization. The source-order cap is9. The independently generated pair sets agree, but neither this agreement nor the positive controls replaces the primary completeness proof.

check_sharp.py checks all eight gcd reductions and absences, compares the two full profile sets, directly counts(2), and literally counts its four-core child461. No unlisted negative conclusion is inferred outside the declared bounded computation.

## 6. Reproduction and scope

From this directory:

    g++ -O3 -std=c++17 enumerate.cpp -o enumerate
    cd run153
    ../enumerate 153 26 > NORMAL.json 2> PROGRESS.txt
    cd ..
    g++ -O3 -std=c++17 rooted_scan.cpp -o rooted_scan
    cd rooted153
    ../rooted_scan 153 26 > NORMAL.json 2> PROGRESS.txt
    cd ..
    python3 check_sharp.py > SHARP_NORMAL.json
    python3 -O check_sharp.py > SHARP_OPTIMIZED.json
    cmp SHARP_NORMAL.json SHARP_OPTIMIZED.json

Unoptimized C++ runs with _GLIBCXX_ASSERTIONS reproduce the same JSON and TSV bytes. The exact data files and dependency hashes are retained. All validation gates use explicit failures rather than removable assertions.

The primary program also prints its inherited base137 control and constructs bounded weak profiles; those fields are supplementary. The theorem here is the eight all-order absences at base153,residue1 and their actual shift-eight repair. F(8,137), arbitrary-base occupancy, and the threshold23 remain OPEN.
