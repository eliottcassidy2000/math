# Independent audit: moving-root one-seam grammar

Session 195, 5 November 2026. **PASS after two minor boundary clarifications.** The universal transfer, graph-class description, minimal linear dimension, arbitrary-depth bounded membership algorithm, gcd obstruction, primitive equal-response obstruction, and five-gadget coverage supplement are valid. The unrestricted `{7,21}` spectrum remains open.

The exact reviewed source snapshots are recorded in `reviewed_sources.json`. Suggested publication destination for this audit is `session195/rooted-interfaces-audit/`.

## 1. Scope and physical interface

The rooted state is (H,e,R), where e counts Hamilton words ending at the mark and R counts ordered two-path covers whose second path starts at the mark. The first path may be empty. These are not the nonempty-bin F2 conventions used elsewhere. In particular, the two-empty-slot count of the empty graph is one, so a singleton gadget has B=1.

All three universal transfer rows are correct for arbitrary nonempty rooted input and gadget, whether strong or weak:

    H'=gH-se+LR,
    e'=EH-delta*e,
    R'=QH+BR.

The Hamilton-word decomposition uniquely separates paths using the backward seam from those avoiding it. The endpoint correction occurs exactly for singleton gadgets. In a root-starting cover, using the seam forces it to be the first edge of the second path, giving the independent R and B factors. Empty first slots create no missing case. The deletion product d'=H(G-b)H(A) is also correct.

The stationary-root and universal-sink four-coordinate formulas were independently checked as supplementary claims. The census in the baseline packet uses only moving-root operations and does not silently include these extra moves.

The ordered-block description is exactly the recursively defined class C_k. Its decomposition need not be unique, and the proof does not claim that all tournaments with an everywhere-reaching marked root belong to it. Moving the root preserves the stated reachability invariant. Repeated triangle gadgets yield genuinely strong graphs with t vertex-disjoint directed triangles, hence feedback vertex number at least t.

## 2. Minimality and arithmetic restrictions

The three physical inputs and three contexts in the minimality proof yield the stated continuation matrix with determinant -4. A two-dimensional rational state with fixed linear context maps and linear Hamilton-count output would factor that matrix through dimension two, which is impossible. The initial-state encoding need not be linear. This proves exactly the claimed rational linear lower bound; it does not exclude arbitrary nonlinear encodings.

Integer transfer matrices preserve gcd(H,e,R) as a divisor of every descendant coordinate. The strong hostile (9,3,9), with deletion coordinate3, was checked directly. Its divisibility obstruction also survives the combined four-coordinate operations because all four coordinates are divisible by3. No claim is made that triple divisibility alone controls an arbitrary fourth coordinate.

### Primitive equal-response obstruction

`EQUAL_RESPONSE_BARRIER.md` is valid for every physical initial state satisfying R=h, h>=5 and h-e>=2, with arbitrary finite moving-root gadgets and target-dependent choices. The complete descendant support through3h+4 is exactly

    {h,2h-e,3h}.

The proof exhausts every nontransitive gadget and every marked transitive gadget. A key boundary checked independently is that a first transitive gadget with marked position j=1 gives3h for **every** size m>=2; m>=3 has sufficiently large R to prohibit continuation below the ceiling. The final proof includes this case correctly.

The explicit primitive strong states (5,1,5) and (49,13,49) satisfy the theorem. Their exact low supports are {5,9,15} and {49,85,147}. This obstruction genuinely goes beyond gcd content. It applies only to the moving-root operation and is not extended to stationary-root moves or other changes of graph presentation.

## 3. Completeness without an arbitrary depth cutoff

For every physical step, H'>=max(H,R) and R'>=H+R>R. The count-bounded search is therefore correct when it records a final count before discarding a state whose new R exceeds the positive integer bound X. Every state preceding a successful final step has R<=X; no longer word can return after crossing that barrier.

An explicit hostile control uses only the singleton gadget. At X=5 the retained states are (1,1,1), (1,0,2), (3,1,3), while the final (5,2,6) must still contribute H5. The independent implementation checks this terminal case.

The sharper Fibonacci depth bound is also valid. With F1=F2=1, after t steps H_t>=F_(t+1) and R_t>=F_(t+2). Thus H<=20,000 requires t<=21, and a bank of maximum gadget order5 has a successful presentation of order at most106. This is a safe bound, not an assertion of optimal order. The state algorithm remains the sharper practical test.

The initial wording was clarified in two places during review: the endpoint-floor inequality e'>=2q requires a **nontrivial** strong old graph, and the search theorem requires a **positive** integer X. The final producer implementation explicitly rejects a nonpositive bound. No unresolved mathematical correction remains.

## 4. Independent inventory, census and witnesses

The independent C++ verifier enumerates every labelled tournament of order1 through5 and every mark, with no isomorphism or count filtering:

- 1,099 graphs;
- 5,405 rooted gadget profiles;
- exactly65 distinct seven-constant matrices.

Every full profile is checked both by subset endpoint DP and literal permutation/cut enumeration, including deletion counts and both empty-slot conventions. The independently generated matrix set and deletion metadata match the producer inventory exactly.

The universal laws are checked on 3,376 actual joined graphs. This includes all841 rooted pairs with both orders at most3, whose full output profiles are additionally counted literally, and every matrix on39 rooted physical inputs. The stationary-root and sink formulas receive78 further physical checks.

The independently implemented arbitrary-depth closure gives:

- X=1000:12,711 retained states,493 scalar values;
- X=20000:1,085,407 retained states,9,993 scalar values.

Both have exactly the missing odd values

    7,21,67,163,173,247,511.

These state counts refer to retained R<=X states, not to all possible final marked triples. The numerical exclusions are complete for this specified grammar at arbitrary depth, not for all tournaments.

A separate independently written forward-append DP reconstructs all ten baseline witness graphs, checks every intermediate recorded triple, and verifies final rows, root, H and triple. Targets checked are11,13,41,43,47,127,137,139,299,983.

## 5. Five-gadget extension

All five proposed order-six gadgets have independently verified literal profiles. Their complete output graphs were reconstructed and independently counted. They are strong and have:

- H67, output order7, state (67,4,134);
- H163, output order9, state (163,9,339);
- H173, output order8, state (173,8,192);
- H247, output order9, state (247,18,408);
- H511, output order10, state (511,7,877).

The extended bank has exactly70 distinct matrices. Its independently computed complete closure through1000 has12,863 retained states and498 odd values, omitting only7 and21. The separate artifact checker also verifies all five full eight-coordinate gadget profiles, output rows, roots, construction steps and four-coordinate output states.

The claimed coverage through20,000 follows rigorously without requiring a new70-matrix search at20,000: the original bank already covers every other odd there, the five explicit words fill its remaining five positive holes, and the complete enlarged-bank search through1000 excludes7 and21 at arbitrary depth. This proves the stated **finite** coverage and does not prove eventual coverage.

The primitive equal-response low-support controls pass on all distinct qualifying small physical states plus the H49 hostile, under both the65- and70-matrix banks. Their all-gadget theorem rests on the proof, not on these finite controls.

## 6. Reproduction and artifact integrity

The standalone C++ source does not import producer code. Optimized and unoptimized bounds-checked runs agree byte for byte. The independent Python artifact checker also agrees in normal and optimized modes. Tests use explicit exception gates rather than removable assertions. Compiled binaries are excluded from the release manifest.

From this audit directory, with the producer packet at the sibling path `../rooted-interfaces`:

    g++ -std=c++17 -O2 -Wall -Wextra -pedantic verify.cpp -o /tmp/rooted-audit
    /tmp/rooted-audit > /tmp/rooted-receipt.txt
    cmp receipt.txt /tmp/rooted-receipt.txt
    g++ -std=c++17 -O0 -D_GLIBCXX_ASSERTIONS -Wall -Wextra -pedantic verify.cpp -o /tmp/rooted-audit-debug
    /tmp/rooted-audit-debug > /tmp/rooted-debug.txt
    cmp receipt.txt /tmp/rooted-debug.txt
    python3 verify_artifacts.py --packet ../rooted-interfaces > /tmp/rooted-artifacts.json
    python3 -O verify_artifacts.py --packet ../rooted-interfaces > /tmp/rooted-artifacts-optimized.json
    cmp artifact_receipt.json /tmp/rooted-artifacts.json
    cmp /tmp/rooted-artifacts.json /tmp/rooted-artifacts-optimized.json

The C++ verifier rewrites its own65-row `inventory.tsv`; the artifact checker compares that independent result to the producer inventory. No shared repository checkout was edited during this audit.
