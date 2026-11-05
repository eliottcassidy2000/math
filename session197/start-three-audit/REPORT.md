# Independent review of the start-three grammar and complete C3 enumeration

**PASS, no correction.** The all-order exceptional-root argument, rooted orbit corrections and final unrooted generating function are sound under the named, already reviewed C1/C3 structural dependencies. All copied source replays and independent finite controls agree. The original source files were not changed.

## Frozen source binding

The files under `session197-lowcycle/` are:

- `START_THREE.md`: 5b0010444d0339741c1ad3ece83e3a30511457d9c3e9f74f037af18510887753.
- `verify_start_three.py`: 629d9b42f2932ada166faf4c0057932904c34ef670467ce327754bbf68868a14.
- `start_three_certificate.json`: bc8e2338813a02000303cb430401839c5ab7771c1eade41b15cbcd61744c97f3.
- Pinned code helper `verify_start_two.py`: ca788fa6489cbca7bed8f9f9a84a553ff9f5de7502833544f9b51b9f059142ae.

The proof uses `START_TWO.md` for the C1 consecutive start-one block, its rigidity and refined Fibonacci counts; `BRIDGE_TWO.md` for the indegree-one vertex outside that block and the cycle-preserving two-degree-one bijection; and `C3_COMPLETE.md` for the all-order degree-one disjunction. These named sources were read. The previous independent C3 finite-census review is `session197-bridges/C3_COMPLETE_CENSUS_REVIEW.md`. This audit does not replace that structural dependency with an assumed finite stabilization.

## 1. Canonical entry grammar and orbit handling

Deleting the distinguished root gives the intrinsic ordered SCC list. The exact root-start formula separates the initial entry sum from the product of later component counts. A total of three therefore gives exactly entry sum three with singleton tail, or entry sum one with one triangle and singleton tail components. A nontrivial strong component of count three is a triangle: every larger strong tournament contains the strong four-vertex graph with five Hamilton paths, and induced path-count monotonicity applies.

For entry sum three, the partitions 3, 2+1 and 1+1+1 give exactly the three listed entry types. A start-one entry forces the initial strong graph to have one Hamilton cycle. The only whole-component entry set is the three vertices of a triangle. An empty tail in that case makes the new root a source; its subtraction by z^4 is necessary and sufficient. In every other branch, the prescribed incoming terminal incidence supplies strongness.

The triangle-tail root-incidence counts must be orbit counts. A terminal triangle permits three rotation orbits (subset sizes zero, one and two), while a nonterminal triangle permits four. Rooted start-one inputs are rigid, including the triangle with its root fixed. Ordered SCC positions cannot be exchanged. This yields exactly (2k-3)2^(k-3) tail choices, not a labelled-vertex multiple.

The mixed C1 start-one/start-two marks and the unordered start-one triple are also correct. Cloning the first vertex of a consecutive unique-start pair gives the complete C1 start-two-root class. The triangle's three possible pair positions are one rooted isomorphism orbit, whereas C1 graphs of order at least four are rigid. These boundary corrections give M=z^4(3-4z)/Delta and T=z^3(1-2z)/Delta, and the canonical recursion gives the stated B3 numerator.

## 2. Exceptional roots, including the small boundaries

Here exceptional means no indegree-one vertex other than possibly the marked root. It does not mean bridge-free.

- In the single-entry branch, all non-entry indegrees in the old strong graph are unchanged, while the entry indegree increases by one. The new root has at least four inneighbors in the old start-three input; tail vertices receive arcs from the entire input. Thus exceptionality is exactly inherited from the marked input, contributing K E3.
- In the mixed-pair branch, the C1 outside-indegree lemma forces an indegree-one vertex outside the entire start-one set. For the output to qualify it must be the selected start-two vertex. The clone characterization makes that vertex have outdegree two, hence indegree |U|-3. Therefore |U|=4. In the unique strong four-vertex graph there is exactly one qualifying mixed pair, contributing z^4 K. No larger C1 input can qualify.
- In the triple branch, the same outside-indegree lemma excludes every input of order at least four. The triangle input qualifies only with a nonempty tail; its old vertices gain indegree two and all tail vertices have indegree at least three. This contributes z^3(K-z).
- The triangle-tail branch never qualifies. A singleton initial component has indegree one from the new root; a nontrivial rooted start-one component retains an indegree-one vertex distinct from its entry, including at the triangle boundary.

Consequently E3=K E3+z^4 K+z^3(K-z), giving E3=z^5(2-z)/Delta. This yields 2,5,13,34,89,233,610,1597 at orders five through twelve, independently verified below.

## 3. Forgetting a degree-one mark

An outdegree-one vertex in a C3 graph starts exactly three Hamilton paths. Its deletion grammar has no free tail incidence bits. The initial component either has a marked start-three vertex followed by singleton components, or a marked start-one vertex followed by one triangle and singleton components. The second tail has k-2 triangle positions. This proves the formula for the marked series P.

A strong tournament of these orders has one or two outdegree-one vertices. Every automorphism fixes each individually: exchanging a two-vertex set would reverse its mutual arc. Therefore the marked multiplicity is one or two even when the graph is nonrigid. The proof correctly does not import C2 rigidity.

The cycle-preserving two-degree-one insertion/deletion bijection gives two-marked-degree graphs at order n in bijection with marked graphs at order n-1. Thus A=(1-z)P. For the outdegree-only class, the deletion input must be exceptional. The same bijection preserves absence of indegree-one vertices: the only old vertex whose indegree changes has outdegree one and remains far above indegree one at these orders. Hence O=z E3.

Finally the complete C3 degree-one disjunction says that the classes not counted in A are precisely the indegree-only classes, which are converses of the outdegree-only class. Thus C3=A+O, with no double counting of graphs having both orientations. This produces exactly the stated numerator

    z^5(1+5z-17z^2+11z^3-z^4)/Delta^2.

## 4. Independent finite controls

`check.cpp` reads the 27,945-graph low-cycle bank through order twelve, whose completeness and exact isomorphism comparison were independently established in the cited C3 census review. The named CSV `session197-bridges/c3-source-records.csv` has SHA256 1d0159168c03fcd011e90caab5f0fe1dd170362ba414d2f74730fe9d644b480a. Its rows are not treated as a newly proved complete bank in this review.

On every graph, the new checker independently:

- verifies strongness and counts Hamilton cycles by forward endpoint subset DP;
- counts start paths by a separate reverse subset recurrence, capped only at four;
- enumerates the full automorphism group by isomorphism-invariant degree refinement and exhaustive compatible bijections;
- computes actual start-three root orbits, rather than counting eligible vertices;
- checks the exceptional-root condition in each orbit;
- partitions every C3 graph by the two degree-one orientations;
- physically deletes the first vertex in every two-outdegree-one chain and checks strongness, cycle count and preservation of the indegree-one absence condition.

This checks 37,019 automorphisms, 6,156 mixed C1 mark orbits and 7,287 physical inverse deletions. The complete start-three orbit counts are 8,31,106,343,1070,3252,9694,28469. The exceptional counts are those in Section2. The C3 counts are 1,11,38,124,390,1193,3574,10537. All marked/unmarked and one-sided identities hold at every checked order. Optimized and unoptimized C++ with library assertions agree byte for byte.

`check_small.py` adds a separate literal control using every one of the 532 unlabelled tournaments through order seven from `session197-prime-atoms/source-bank-order7.json`, including graphs outside the low-cycle bank. It enumerates all vertex words and all vertex permutations directly, counts cycles and start paths, forms full automorphism root orbits, and reproduces the start-three, exceptional and final C3 statistics. Normal and Python -O outputs agree.

The source verifier was copied before execution because it writes a certificate beside itself. Its normal and -O outputs match each other and the frozen certificate exactly.

## Audit artifacts

- `check.cpp`: 3c8904da24061436eb82b8abb6815e59e97f8e8fad644ff6fb86af70213247b1.
- `independent-normal.txt`: 5bec5723b6a86761d01c7812e112678c611ec830d72f0b61a5c349b043b32632.
- `check_small.py`: 1d660928f5ba1fd923c96a659e48367b0a70114c7ec43d0b8a60bfc617388863.
- `small-normal.json`: f521c94b9f39280ee97e070b35cf964e276ab0cb55b0ab102dd239c6ce0c7cd3.
- `source-normal.json`: bc8e2338813a02000303cb430401839c5ab7771c1eade41b15cbcd61744c97f3.

No source correction, unsupported global spectrum closure, or outstanding blocker was found.
