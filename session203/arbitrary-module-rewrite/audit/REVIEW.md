# Independent audit: arbitrary shared-module correction and returning count ray

6 November 2026. **PASS in the stated scope.** The three physical constructions satisfy H(child(U))=2H(parent(U))+H(U) for every tournament U, including the empty module. The C rebase and its returning count ray are valid for nonempty U. This is not an occupancy or eventual-threshold theorem.

## Reviewed inputs

The reviewed PROOF.md has SHA-256 7edd5fb876daf5692581b0d35593c1998402c2b6cd5bea1275b08014086ab81a and is preserved as REVIEWED_PROOF.md. Every entry in the producer's nine-entry MANIFEST.sha256 was checked successfully; its exact text is preserved as REVIEWED_INPUT_MANIFEST.sha256.

The retained transitive proof and audit match the independently reviewed promotion hashes 55db518c68502746f63b3cdea4671d609425151a47963d7b55f262aeaec258d3 and 1b7af2329657e9aa26a0d6eae629aea875a7a6ed1f0e273e6f49e550455776d3. These establish the physical cores and all finite path/cycle skeleton coefficients. The new all-module inference is reviewed below, rather than inferred from finite module tests.

## Arbitrary-module bijection

Deleting U from a spanning path leaves one permutation of the nonempty core and a set of occupied gaps. Each occupied gap holds one nonempty directed path in U. Since all module vertices have identical incidences to every core vertex, any ordered path partition of U can be inserted in the same allowed gaps, regardless of U's internal tournament structure. Deletion and insertion are inverse. Thus the same skeleton coefficients a_r multiply the counts A_r(U) of ordered partitions into r nonempty directed paths.

For cycles, fixing one core vertex as the starting point removes rotation exactly once. The wraparound gap remains a single gap, and no module run is split by the chosen start. Internal cycles of U do not create extra factors, because each occupied gap contains a path, not a cycle. The cyclic skeleton expansion is therefore equally valid.

Selecting the occupied labels among q bins proves B_U(q)=sum_r binomial(q,r) A_r(U); finite binomial inversion gives the displayed reverse identity. These coordinates describe actual coupled covers of the same U. They are not independently assignable scalar parameters.

## Empty boundary and exact formulas

For empty U, A_0=1, A_1=0, and B_U(q)=1 for every q>=0. For nonempty U, A_0=B_U(0)=0 and A_1=B_U(1)=H(U). Therefore the child-minus-twice-parent skeleton vector (1,1,0,...) gives exactly A_0+A_1=H(U) in both cases. There is no missing empty-module correction.

With h=B_U(1) and b_j=B_U(j), independent expansion agrees with every formula:

    A parent H = h+17b_2+8b_3-b_0
    A child H  = 3h+34b_2+16b_3-2b_0
    B parent H = 5h+20b_2
    B child H  = 11h+40b_2
    C parent H = 27b_2-h-b_0
    C child H  = 54b_2-h-2b_0.

Cycle counts are A: h+b_2 and h+2b_2, and B/C: h and 3h. In particular A's child cycles equal twice the parent cycles minus h. Positivity of h, including the empty convention h=1, gives a spanning directed cycle in every parent and final child, proving strongness.

## Physical edits and deletion

The module remains an induced copy of the original arbitrary U. The core-to-module masks are literal and unchanged between the retained old vertices. A reverses exactly core arc 1->4 and adds one vertex, retaining all old vertices.

B and C reverse the stated arc incident to old core vertex 0, then add their new vertex. Vertex 0 is a universal sink even when U is nonempty, because the relevant module mask includes its bit. Deleting it preserves H bijectively and leaves the same order as the parent. All arcs between retained old vertices, including every edge within U, remain unchanged. There is no edit-minimality assertion when deletion is allowed.

## Coupling bounds and returning C ray

For nonempty m-vertex U, each Hamilton path and one of its m+1 cut positions gives a distinct ordered two-bin cover. Concatenation and the cut length recover the source, so b_2>=(m+1)h. Swapping the bin labels has no fixed point for nonempty U, so b_2 is even. Both stated parent lower bounds follow.

In U->singleton, each cover is extended by assigning the new sink to one bin and appending it there. Consequently B_(U->singleton)(q)=q B_U(q), with H unchanged. Repeating gives q^k B_U(q); the q=0 empty boundary agrees with the combinatorial convention for zero added vertices. The new independent controls check q=0,1,2,3 with one and two added sinks.

For nonempty U, b_0=0 in both U and its padding. Hence C's final child and the next parent P_C(U->singleton) both have 54b_2-h paths. Their cycle counts are 3h and h and their orders differ by one, so the proof correctly calls this a count-preserving rebase, not an isomorphism or induced extension. Empty U genuinely fails: the final child has H=51, while the padded parent has H=53.

For fixed nonempty U, the explicit parent representations based on U->Tr_k therefore have N_k=27*2^k B_U(2)-h and N_(k+1)=2N_k+h. The correction h stays attached to that original physical module. This is a returning count ray, not a method of choosing a fresh arbitrary remainder at each step.

An additional independent hostile confirms that the rebase preserves H only, not the whole cover profile. For U a singleton:

    final C child:       (H,B_2,B_3) = (107,2600,26235)
    padded C parent:     (H,B_2,B_3) = (107,4804,73041).

Thus the two graphs cannot be substituted interchangeably as modules in later operations that depend on B_2 or B_3. The valid returning ray explicitly reconstructs its next parent from the padded original module.

## Independent direct controls

The standalone check.py imports no producer code. It enumerates every labelled module of orders 0 through 4: 1,1,2,8,64 modules, 76 in total.

- The coordinates b_0,...,b_3 are counted by enumerating a full vertex permutation and all weak cut-length vectors into labelled bins, checking every within-bin arc. This differs from the producer's subset-convolution calculation.
- All 228 physical rewrite pairs and 684 constructed graphs are checked by nonmemoized literal directed-path DFS. It visits 24,225,202 path prefixes. Cycles use fixed starting vertex 0 and are counted modulo rotation.
- Every path/cycle formula, physical edit, retained-graph comparison, sink deletion, strongness check, padding law, nonempty rebase, and empty hostile passes.
- The higher-profile hostile above uses another direct calculation: enumerate every simple path to obtain all subset counts, then assign each vertex explicitly to one of q labelled bins and multiply the path choices. It does not identify equal-H graphs as equal profiles.

Normal and optimized Python receipts are byte-identical. Checks use explicit exceptions, not removable assertions. The finite controls are evidence for the implementation; the arbitrary-order conclusion comes from the proved path-partition bijection and inherited exact core coefficients.

## Reproduction and scope

    python3 check.py > NORMAL.json
    python3 -O check.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The result neither realizes arbitrary triples (h,b_2,b_3) nor proves that every sufficiently large odd count belongs to a displayed image. In particular an E(u) forward sector is not H(E(u)); inserting E(u) as the module retains its full correlated count. No threshold or free-sector extraction is inferred.
