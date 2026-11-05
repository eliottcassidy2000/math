# Independent audit: complete endpoint-one recursion

Session 195, 5 November 2026. **PASS.** The unique inverse decomposition, converse construction, exact coordinate transfers, count-bounded completeness theorem, and stated obstructions are valid. No circularity or unsupported completeness claim was found. The unrestricted `{7,21}` Hamilton-path spectrum remains open.

The reviewed producer proof is `session195-endpoint-recursion/PROOF.md`, SHA-256 `ebddac9ad6e5e2c9d99b4f9f1d13f912891ef9e8ffaa56c2275f7aff659523df`. This review includes the self-contained binary-apex recurrence, count-10,000 census and H511 witness additions. Suggested publication destination for this packet is `session195/endpoint-recursion-audit/`.

## 1. Structural classification

For a marked strong tournament (T,v) with exactly one Hamilton path ending at v, delete v and let C be the sink SCC of T-v, with preceding union P. Then

    e_v(T) = H(P) sum_(u in C, u->v) e_u(C).

Every e_u(C) is positive. Equality to one therefore forces H(P)=1, P to be a transitive chain, and a unique u in C beating v, with e_u(C)=1. This independently reproduces the proposed inverse theorem.

The sink SCC, preceding chain, distinguished u and chain order are unique for the specified labelled rooted graph. The inverse recursively terminates at a singleton. This is not uniqueness from the numerical triple, nor uniqueness after forgetting the distinguished endpoint; the proof correctly makes neither claim.

The converse is correct. With a nonempty preceding chain, v must beat its first vertex, and that condition is sufficient for strongness. With empty chain, the old core must be nontrivial; a singleton would give a weak order-two graph. No additional arc restrictions are missing. The same endpoint factorization proves the new terminal count remains exactly one.

For a nonstrong endpoint-one tournament, all SCCs preceding its sink SCC have product Hamilton count one. Thus they are transitive singleton components, and the sink SCC has the same H and a unit terminal endpoint. The reduction of **scalar support**, not the full marked state, to the strong class is valid.

## 2. Coordinate and word conventions

The state is (h,s,L), where s counts Hamilton paths **starting** at the marked vertex and L counts ordered two-path covers whose first path is nonempty and **ends** at the mark. The second path may be empty. These conventions are essential: the singleton has (1,1,1), and L includes the full-path contribution.

The proposed h, s and L transfers are correct for every allowed chain word. The h formula partitions paths according to whether they use the unique reversed core-to-new-root edge. The L formula partitions by the predecessor of the new root: no predecessor, the old marked root, or a zero-bit chain vertex. All bin assignments are recoverable, so the counts are injective and exhaustive.

Words must begin with 1 when nonempty. This exactly encodes strongness and gives a canonical positive binary integer A. Allowing leading zeros while retaining the formula s'=h would be false: the new root cannot start a Hamilton path when the first chain vertex beats it. The proof correctly excludes this case. Integer A=0 denotes only the empty word, and that operation is allowed only for h>1.

The finite-alphabet linear implementation maintains its stated eleven coordinate invariants exactly. Its final projection gives the same three-coordinate return. It is a fixed-dimensional linear implementation with unbounded integer states, not a finite-state automaton; the draft distinguishes these correctly.

The empty-word matrix satisfies M^3=M^2+M+I, so the displayed Tribonacci recurrence is an exact all-order consequence. The initial values 3,5,9,17,31,57,105,193 are correct.

## 3. Completeness of count-bounded closure

The bounded search is genuinely complete at all vertex orders:

- For a nonempty word, h'>=(A+1)L and A>=1. Hence A<=floor(X/L)-1 for any child of count at most X, and every allowed word is represented once by its positive binary expansion.
- Such a child has h'>=2L>=2h>h.
- For an empty-word child, h'-h=L-s>=h-s>=|C|-1>=2.
- Therefore inverse recursion strictly decreases H, and the finite branch bound excludes arbitrarily long hidden words.
- Merging equal triples is sound because every future allowed transition and its output depend only on the triple. Every reached triple has a physical witness, and every possible graph state reaches a smaller parent by the proved inverse theorem.

The theorem does not rely on a vertex-order cutoff, a strong-tournament minimum theorem, an assumed list of missing counts, or the finite census it certifies. No circularity was found.

## 4. Independently reproduced obstructions

The independent closure has exactly the stated thirteen triples for H<=13. Its scalar values are {1,3,5,9,11}; H13 is absent. The physical tournament `[2,12,9,17,7]` has H13 and terminal counts `[4,3,2,2,2]`. Thus no tournament of any order with H13 has a unit terminal endpoint; reversal gives the same conclusion for a unit starting endpoint.

Through X=1000, the independent closure has 9,293 states and 475 scalar values. Its absent odds are exactly

    7,13,21,41,43,55,67,77,83,163,173,231,233,
    247,251,277,299,421,449,503,533,569,583,637,821.

Through X=10000, it has 284,844 states and 4,971 scalar values. The missing list gains exactly 1067,1339,1861,2647.

These are all-order absences **from the endpoint-one class**. They are not new absences from the unrestricted spectrum. A recursive grammar requiring endpoint-one at every return is therefore insufficient by itself to cover all the desired atoms.

The added singleton macro `111010111` was independently built as an eleven-vertex graph. It is strong, has endpoint-one at its new mark, and has state (511,1,553), confirming the stated H511 witness.

## 5. Independent exact controls

The standalone C++ verifier does not import the producer's implementation. It uses terminal-removal subset DP, independent all-vertex reachability for SCC extraction, direct literal permutation/cut enumeration, and physical adjacency reconstruction.

Its exhaustive universe contains every labelled tournament of orders one through six: 33,867 graphs, with no isomorphism or H filtering. All 16,159 strong endpoint-one marked cases pass unique sink-SCC decomposition, strongness boundary conditions, old-root recovery, exact graph reconstruction and the three-coordinate transfer.

All 1,099 labelled tournaments through order five have H, all start/end counts, and every marked L count independently checked by literal permutations and cuts, including empty second paths.

The forward test uses every labelled strong endpoint-one rooted input through order four and every allowed positive word of length at most five, plus valid empty words: 2,527 physical graphs. All pass strongness, endpoint count, direct transfer and finite-alphabet transfer; 183 small outputs additionally receive literal word/cut checks. The three linear basis states are tested on every nonempty allowed word through length twelve, for 12,285 further checks.

The independent census computes the apex count g from a direct boundary formula rather than the producer's binary recurrence. For bits b_1,...,b_m,

    g = b_1 + (1-b_m)
        + number of adjacent 01 patterns
        + sum_(j<i, b_j=1, b_i=0) 2^(i-j-1).

This enumerates the possible last vertex before the apex and first vertex after it. The arithmetic state closure then independently reproduces all stated counts and missing lists for X=13,1000,10000.

The verifier uses explicit exception gates, not removable assertions. Optimized and unoptimized bounds-checked executions agree byte for byte. The source and receipts are included; compiled binaries are not publication artifacts.

## Reproduce

From this directory:

    g++ -std=c++17 -O2 -Wall -Wextra -pedantic verify.cpp -o /tmp/endpoint-audit
    /tmp/endpoint-audit > /tmp/endpoint-receipt.txt
    cmp receipt.txt /tmp/endpoint-receipt.txt
    g++ -std=c++17 -O0 -D_GLIBCXX_ASSERTIONS -Wall -Wextra -pedantic verify.cpp -o /tmp/endpoint-audit-debug
    /tmp/endpoint-audit-debug > /tmp/endpoint-debug.txt
    cmp receipt.txt /tmp/endpoint-debug.txt

No edits were made to the shared repository checkout.
