# Independent audit of the hereditary positive-start floor

6 November 2026. **PASS for the rooted insertion lemma and the specified interface obstruction.** Every positive start count in an induced extension of the stated colored-body interface, with body order at least three, is at least 27. The four special wrapped words have the stronger floor 119. Neither result refutes F(8,137) or excludes any unrestricted scalar Hamilton count.

## Frozen input and construction scope

The reviewed proof has SHA-256 c0cd6133ad4b118249574f1c788f2eb61ded4b0b5f8bbbc179f794be90e58863. The four-word verifier has hash 4de0f422007892a2c8414351bf850535b0990bac5634abadf1dd639fc078e91f; the colored-kernel verifier has hash 64dc25350e4b911098fe0a1b2934b2dd7c59d17dc4a9311278202aeb03927e6e. All nine manifest entries pass. Exact proof, programs, and both complete kernel files are retained.

The underlying core is (102,0,10,83,135,30,182,47). Vertex 5 is replaced by a tournament body Q. Only old vertex 3 has arbitrary incidences to individual body vertices; every other outside incidence is inherited uniformly. Old vertex 1 is a global sink. Removing it leaves six fixed core vertices and the body. These six fixed vertices, not five, are needed in every nine-vertex kernel.

The main hereditary statement requires body order at least three. It is not a classification of arbitrary colored-body contexts, nor does it cover changed incidences at other core vertices or other internal core arcs. The four-word construction is fully specified in the proof, so no external count theorem is required for this obstruction.

## The root-preserving insertion map is injective

Let S be any nonempty induced subtournament of T, and assume r reaches every vertex of T. Choose a simple directed first-entry path from r to S, ending at v. If r is in S this path is just (r). All earlier vertices lie outside S. Concatenating this fixed prefix with a v-starting Hamilton path of S gives a rooted directed path on a fixed vertex set K.

For any proper K containing r, global reachability supplies an edge from K to its complement: take the first exit on a path from r to an outside vertex. Fix the next outside vertex using only K. This produces one vertex-insertion schedule for every initial kernel path, independent of the order of that path.

Given the chosen vertex w, at least one current path vertex beats w. Take the rightmost such vertex x. If x is last, append w. Otherwise the next path vertex loses to w, by rightmost choice and the tournament orientation. Insert w between x and that successor. This always inserts after an existing vertex and therefore never changes the initial root. Repeating gives a Hamilton path of T beginning at r.

Deleting all vertices outside S recovers the original kernel Hamilton word exactly: the fixed entry prefix and every subsequently inserted vertex were outside S, and insertion never reorders existing vertices. This is a left inverse, so the map is injective. Consequently

    s_r(T) >= s_v(S) >= min_(x in S) s_x(S).

Strongness of S is not needed for the injection itself. It makes the minimum positive in the intended application. If s_r(T)>0, its Hamilton path already proves that r reaches all of T. Thus the correct unconditional conclusion is s_r(T)=0 or s_r(T)>=min_x s_x(S), not a positive bound at every root without qualification. The proof correctly retains this zero alternative and allows r outside S.

## All finite kernels checked literally

The independent checker places the body first and the retained old core labels (0,2,3,4,6,7) last. This differs from both archived labeling conventions. It builds every arc from the stated core and color rule, then verifies the permutation to each archived graph and every root count.

Counts are obtained by nonmemoized directed-path enumeration, not either producer's subset dynamic program. Every completed Hamilton word increments its actual starting vertex. Independent reachability searches test strongness at all vertices.

All eight labeled three-vertex body tournaments and all eight colorings are covered, with explicit verification that the body-code/color index set is exactly {0,...,7} squared. All 64 resulting nine-vertex kernels are strong, and their minimum start counts are at least 27. Every archived adjacency row and all 576 endpoint counts agree.

The displayed representative minima also agree exactly:

- Transitive body (6,4,0), colors 0 through 7: 27,35,31,39,29,37,33,41.
- Directed triangle (2,4,1), colors 0,1,3,7: 45,51,57,63.

The representative classification is not needed for completeness because every one of the 64 labeled inputs is checked. The value 27 is attained. Adding back the universal sink preserves every surviving vertex's positive start count, so sharpness also has a literal full-interface witness.

The four eleven-vertex special kernels independently have:

    seed 10: H=7213, minimum start 119;
    seed 14: H=7221, minimum start 127;
    seed 17: H=7227, minimum start 133;
    seed 21: H=7235, minimum start 141.

All are strong, and every archived root count agrees. The original root need not realize the minimum; all-root counting is essential because the theorem permits arbitrary rerooting.

## Why the finite certificates imply the stated all-order floors

For an arbitrary body Q of order at least three, choose any three of its vertices. Together with the six fixed nonsink core vertices they induce exactly one of the 64 counted kernels: the selected internal body arcs form one of the eight tournaments, and their incidences to old vertex 3 give one of the eight colorings. No transitivity, module structure, or uniformity toward vertex 3 is assumed. Other core/body incidences must remain as specified.

Every tournament containing the full interface as an induced subtournament therefore contains a strong nine-vertex kernel with minimum start at least 27. The insertion lemma applies at any root with a positive start count. Arbitrarily many added vertices, arbitrary arcs to those vertices, and arbitrary rerooting do not change that conclusion.

For a wrapped word 1^k c 0^k, retaining only the five middle body vertices and six nonsink old vertices induces the corresponding eleven-vertex k=0 kernel. This is a graph restriction, not a count or profile equivalence. The independent checker verifies 104 such literal embeddings for k=0 through 25. The selection rule itself proves the embedding for every k>=0.

F(8,137) permits starts r+3j, with r in {1,2} and 0<=j<=8, hence positive starts at most 26. These interfaces cannot provide those source witnesses while the relevant strong induced kernel is retained. To escape this obstruction, a construction must leave this hereditary class, for example by altering the kernel's internal arcs or deleting part of it. Such alterations are not ruled out. The criterion F itself and unrestricted target coverage remain open.

## Boundary and actual-injection controls

The body-size restriction is substantive. All eight colored two-vertex bodies are also counted independently. Their minimum starts can be 15,17,19 or 21; thus the asserted floor 27 cannot be extended to body order two by the same statement. These smaller bodies are not excluded from F by this result.

The generic lemma is tested through actual maps, not only inequalities. For every labeled host tournament through order five, every nonempty induced subset, and every root reaching the host, the checker constructs the first-entry path, freezes the vertex schedule, inserts every kernel Hamilton path, checks the full output word and root, verifies the inverse deletion, and checks that the outputs are distinct.

This gives 1,099 host tournaments, 115,211 root/kernel maps, 55,702 cases with the root outside the kernel, and 88,886 actual injected paths. Some maps have an empty domain because the arbitrary test kernel is not strong; this is permitted by the general lemma. There are 106,400 schedules containing later insertions, and entry paths reach length five. Together with the finite kernels, the direct path enumeration visits 1,382,986 directed prefixes. Normal and optimized receipts are byte-identical, with explicit failure checks active in both modes.

Both producer scripts are independently rerun normally and with Python optimization. Their receipts and regenerated complete kernel files reproduce the frozen originals exactly.

## Reproduction and citation packaging

In this directory:

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The producer replays use PRODUCER_VERIFY.py and PRODUCER_COLORED_VERIFY.py. They write redundant KERNELS.json and COLORED_KERNELS.json, which must match the retained PRODUCER_KERNELS.json and PRODUCER_COLORED_KERNELS.json respectively.

Before publication, replace the proof's scratch-only four-choice citation with the immutable repository citation verified by the author:

https://github.com/eliottcassidy2000/math/blob/be393666d567723fb30675ed0e2952e5c784384c/session203/four-choice-endpoint-fiber/PROOF.md

That is a citation-only packaging change; preserve the reviewed original and its hash. No mathematical claim depends on unpublished files or on completeness of a source bank. Status promotion and publication require their separate final packaging check.
