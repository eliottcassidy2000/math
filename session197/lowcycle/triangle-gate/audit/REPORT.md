# Independent audit of the physical C3 triangle gate

Result: **PASS**, with no source correction. This review binds the following frozen files in `session197-lowcycle/`:

- C3_TRIANGLE_GATE.md:886f89fd30ea6098b36e37c7580be033cc0394f7a0fad8d124787d85173c59ce.
- verify_c3_triangle_gate.py:06f17fa3ed65c817403a8a7f776afd5da192f3b042e8d3e6912d8e9e33415cdd.
- c3_triangle_gate_certificate.json:b3b6d4c644547ec0e71864042fc8952a23bf52602d872e1aa2c54d08e21384a8.

The result is an all-order physical transfer and an exact inverse for the named two-parameter family. It is not a classification of all strong C3 tournaments, all feedback-two tournaments, or their scalar counts. The absence of821 applies only to this family.

## 1. All-order proof audit

For a sole-out-neighbor vertex a->b, the intermediate identities in the source are exact. A Hamilton path ending at a is a path of T-a not ending at b. Every other occurrence of a cuts a path immediately before b and is encoded by the b-starting two-cover. A path ending at b must end a,b. In a two-cover whose reserved path starts at b, a ends one of the two paths. Putting it on the reserved path excludes precisely its old singleton-b case, accounting for the subtraction k0. These arguments also cover order2.

In the triangle replacement, every module exit enters b, so there is at most one exit. A two-chunk Hamilton path must end in the module and its intervening S-path cannot be singleton b. This proves the scalar and deletion formulas with the stated six/two internal partition factors. Paths ending outside the module have exactly one chunk, so the endpoint formula follows from the free cyclic symmetry on the triangle vertices.

The three two-cover cases are exhaustive and disjoint. With no exit, the contributions are d0+2(d0-f). When the reserved path exits, its module prefix is nonempty and starts at x; the triangle has seven partitions into the other path's tail, that prefix and the reserved tail. If the intervening S-path is singleton b, its module tail must be empty, leaving three. When the other path exits, its middle module chunk and the module-only reserved path are nonempty, giving four partitions, or two if its S-suffix is singleton b. An empty old unreserved prefix is valid in both cases, since b is on the other S-path and imposes no further endpoint restriction. The sum is exactly the source R formula.

A Hamilton cycle has one module chunk, giving three choices after contraction. A path starting at x has three available module prefix/tail partitions when |T-a|>=2. At |T|=2 the S-path is singleton b, so the tail must be empty and the start count is1, not3. The source explicitly excludes this boundary from its start-count law.

The first binary phase is a literal fixed-port graph. Independently, its one-apex two-cover count evolves K->2K for a1 bit and K->2K-1 for a0 bit, and its apex-ending count is the final-bit complement. Appending a universally dominated vertex preserves H; appending a vertex that beats only the apex changes H by K-e. These facts independently recover the stated g recurrence and seed profile. The two gate letters2,3 physically make the old a and the two births a directed triangle with identical incidences to every exterior vertex.

The old-root O and Z laws follow by placing the newborn at a path end or immediately before the old root. In a root-starting cover the newborn must end one of its paths, with exactly the singleton reserved-path exclusion in the Z case. Both preserve the root-start count for inputs of order>1. A Z output has exactly s_a cycles, regardless of whether its predecessor is weak, so every empty or final-Z suffix in this family is strong with exactly3 cycles. Removing the two ports leaves the unchanged transitive birth order.

The suffix formula is checked algebraically as a linear functional with coefficients

    (h,e,R,d) -> h+(t-1)e+(v-1)R+(g(v)-v+1-t)d.

Its first-Z recursion v=2l and first-O recursion v=2l-1 agree exactly with the physical matrices and the g recurrence, including v=1. Substituting the seed gives E(u,v).

Finally all first coefficients are positive, so replacing g(u),g(v) by1 and parity bits by their adverse bounds is valid. For u>=4 this gives E-7uv>=(u-3)(v+1); the separately handled u2/u3 cases give the claimed bound as well. Thus uv<=(X-1)/7 is necessary for a target X, and testing the exact formula is sufficient because each pair has an explicit word. No larger-order witness from this family can evade the inverse.

## 2. Independent finite controls

The standalone `check.cpp` imports no source research code and has distinct counting algorithms:

- Literal full permutations and all possible cuts count Hamilton paths, rooted endpoints and ordered two-covers. They check all1390 sole-out-neighbor incidences in every labelled tournament of orders2 through5, all1388 valid start-law cases, both order2 exceptions, and814 strong cycle-law cases.
- The same literal method checks10808 O/Z transformations over all roots of every labelled tournament of orders2 through5, including weak inputs.
- Assigning each triangle vertex to a labelled slot, then separately permuting within each slot, independently recovers the exact partition counts7,3,4,2.
- For fixed-port graphs, every Hamilton word is counted by assigning each nonport to one of three increasing slots around either port order. Hamilton cycles use two increasing slots around the anchored cyclic port order. This independently checks380 family graphs with2<=u<=20,1<=v<=20, plus the67/533/1861 examples. These counts do not use endpoint subset DP or the scalar formula.
- Direct enumeration of all directed path words supplies the seed endpoint/two-cover profiles for2<=u<=20.
- Every g(u),1<=u<=714, is recomputed from its literal one-apex graph by two increasing slots around the apex, rather than by the source recursion. These physical values independently reproduce the recurrence and are used in the exact family image calculation.

The complete inverse tests4094 pairs, gives1581 distinct values through5001 and omits821. The independently generated sorted image has SHA2566bc127e15eadfa9ed998a20c79d522dcd18d6a530e5d2caeb4881085bb7412c8, exactly matching the source certificate.

The three example row arrays agree literally with the source certificate. In particular the H533 array is the already reviewed fixed-two-port graph in birth order; the conclusion about minimum C3 at H533/FVS<=2 relies on that separate negative theorem, not on this family census.

## 3. Reproduction and hash binding

The source verifier was copied into the audit-only `source/` directory before execution because it writes a certificate. Its normal and optimized outputs both reproduce the frozen source certificate byte-for-byte; the original files remain unchanged.

Run the independent checker with:

    g++ -O3 -std=c++17 check.cpp -o check
    ./check > independent-normal.txt 2> image.txt
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check.cpp -o check-debug
    ./check-debug > independent-debug.txt 2> image-debug.txt

Optimized and debug outputs/images are identical. The independent hashes are:

- check.cpp:d48991c4a5b49a464d2b4a59ceb24965c617247d51e86316f20253dae5b97d30.
- independent-normal.txt and independent-debug.txt:f29ddf9e521bb2ed235c7702761eaa9c44712ad439594a9642b37bdbf9d0869a.
- image.txt and image-debug.txt:6bc127e15eadfa9ed998a20c79d522dcd18d6a530e5d2caeb4881085bb7412c8.

The proof establishes all-order validity. The finite controls independently test its physical interpretation, small boundaries, examples and exact bounded inverse output; they are not used as an extrapolation to all orders.
