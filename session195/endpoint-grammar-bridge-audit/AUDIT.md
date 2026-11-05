# Independent audit: endpoint-one / A,K,D bridge

Session 195, 5 November 2026. **PASS.** The initial-K source stripping, start-count invariance, literal converse compiler, and equality of scalar supports are valid. No mathematical correction was required. This is a classification of the specified recursive class, not a proof of the unrestricted `{7,21}` spectrum.

Exact reviewed producer snapshots are recorded in `reviewed_sources.json`. Suggested publication destination for this packet is `session195/endpoint-grammar-bridge-audit/`.

## 1. The initial-K exception is real and correctly handled

A changes the root to the new vertex whose only old out-neighbor is the previous root. Every new root-starting Hamilton word therefore has that old root second, giving a bijection with old root-starting words. D appends a universal sink, also preserving the marked start count.

For K, the new vertex cannot occur internally in a Hamilton word beginning at the retained root, since its sole outgoing old edge returns to that already-used root. It must be last. When the old order is at least two, every old root-starting word ends somewhere else and appends legally. At the singleton this argument fails exactly as stated: K makes the retained root a sink, changing its start count from1 to0.

For a word beginning with K^t, the initial added vertices form x_1->...->x_t->r. None can ever become the mark later, because A marks only a fresh vertex and K,D retain the current mark. Each prefix vertex beats every later birth. Thus these vertices remain a transitive source prefix of the entire final graph, regardless of the subsequent word.

Deleting the complete initial K prefix leaves exactly the graph produced by the remaining word, under an explicit relabelling that preserves the root. Its contribution to the SCC product for H is one. This proves scalar invariance at all future depths, not just immediately after the prefix. It also covers an all-K word, whose reduced word is empty.

After erasure, the word is empty or begins with A or D, so its output has marked start count exactly one. It would be incorrect to omit the erasure and claim that every original word preserves marked start-one. The producer proof retains this distinction.

## 2. The graph compiler is exact

For an endpoint macro on (C,u), begin with C's converse, marked at u. The initial A creates v with v->u and all other old C vertices beating v. These are precisely the converse core incidences.

Read the chain bits in reverse order. A bit1 uses K and creates z_i with z_i->v; a bit0 uses D and creates v->z_i. Every C vertex beats each new chain vertex, and each previously created z_j with j>i beats z_i. This reconstructs the converse chain order and every cross edge exactly. The mark remains v. The empty word uses only A, including the harmless graph-level singleton boundary when no strong-output assertion is made.

The compiler is valid as a graph identity even for arbitrary rooted cores and words not meeting strongness conditions. The endpoint-recursion theorem supplies precisely the conditions required when compiling its strong class. Recursive compilation terminates because that inverse theorem strictly reduces the core.

### Exact scope of both directions

The stated proof handles weak endpoint-one witnesses by deleting their neutral transitive source components, which preserves H. The literal converse compiler is asserted for strong endpoint-one witnesses. For the other direction, initial-K erasure supplies a start-one rooted representative with the same H; taking its converse gives the endpoint-one representative.

Thus the audited conclusion is scalar-support equality plus the stated strong graph-level compiler. It is not equality of all original rooted output graphs at their original marks, nor equality of the three- and four-coordinate marked state sets. The finite audit also reconstructs weak endpoint-one converse graphs by D suffixes as an additional control, without expanding the producer theorem's stated scope.

## 3. Coordinates, dimensions and transported exclusions

The coordinate explanation is correct. A resets the deletion coordinate to the previous H. Subsequent K and D operations retain it. Under converse, the old moving-root coordinates (H,e,R) correspond to the endpoint coordinates (H,s,L), with reversal and interchange of the two path slots supplying the R/L bijection. Thus a completed endpoint macro closes in three coordinates even though unrestricted primitive continuations need four homogeneous linear coordinates.

The cited primitive continuation matrix was independently reconstructed from actual A,K,D words and has determinant-16. This confirms the relevant four-dimensional lower bound and the absence of a dimension contradiction.

The endpoint-one all-order exclusions consequently transfer to A,K,D. In particular H13 and the full29-element endpoint omission set through10,000 are absent from this grammar at every depth. These are class exclusions, not unrestricted tournament-spectrum exclusions. The proof does not assume the two finite scalar censuses match; their agreement is a consequence and a control.

## 4. Independent exact controls

The standalone C++ verifier imports no producer code. It uses its own graph births, reachability, endpoint DP, source stripping and recursive inverse compiler.

- All9,841 elementary words of lengths0 through8 are checked as literal graphs. The3,280 words with a nonempty initial K prefix pass exact root-preserving deletion/reconstruction, permanent source orientation, transitive prefix order and equality of H before and after stripping.
- Every reduced output has marked start count1. Every original output has start count0 exactly when its initial K prefix is nonempty. All364 words through length5 additionally receive literal Hamilton-word and path-cut enumeration.
- All17,955 single-macro graph compilations are checked: every labelled rooted input through order4 and every binary word of length0 through5. This includes leading-zero words and nonstrong boundaries, where only the graph identity is claimed.
- A,K,D start-count laws are checked on all285 labelled rooted inputs through order4, giving855 tests and including the singleton K exception.
- Every labelled tournament through order6 is enumerated. All26,175 marked endpoint-one cases are recursively compiled into their exact converse graph with the correct mark. These comprise16,159 strong cases and10,016 weak cases. The weak cases explicitly exercise the restoring D suffix.
- An independent four-coordinate A,K,D saturation through H<=1000 has11,722 retained states and475 scalar values. Its25 missing odd counts match the independently reproduced endpoint-one census exactly. Different marked-state counts are expected and are not treated as a discrepancy.

The C++ tests use explicit exception gates. Optimized and unoptimized bounds-checked receipts agree byte for byte. The all-order conclusions follow from the proofs; these tests specify their finite supporting universe exactly.

## Reproduce

From this directory:

    g++ -std=c++17 -O2 -Wall -Wextra -pedantic verify.cpp -o /tmp/bridge-audit
    /tmp/bridge-audit > /tmp/bridge-receipt.txt
    cmp receipt.txt /tmp/bridge-receipt.txt
    g++ -std=c++17 -O0 -D_GLIBCXX_ASSERTIONS -Wall -Wextra -pedantic verify.cpp -o /tmp/bridge-audit-debug
    /tmp/bridge-audit-debug > /tmp/bridge-debug.txt
    cmp receipt.txt /tmp/bridge-debug.txt

No shared repository checkout was edited. Compiled binaries are excluded from the release manifest.
