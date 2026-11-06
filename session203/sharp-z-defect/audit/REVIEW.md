# Independent audit of sharp Z-defect equality

**PASS after the explicit singleton-boundary correction.** The equality values D=s>0 are exactly the ordinary nonempty tournament Hamilton-path counts. This is an all-order equality classification, not a proof that every positive odd value occurs, a classification of every equality graph, or an eventual spectrum threshold.

## Reviewed version and correction

The original proof, SHA256 `97273c71f0fbfc71fff9acf62a0689bfdb3ebb8524ecebbfb2eaa99ce32a4830`, and its original manifest inputs are preserved under `reviewed-inputs/`.

The approved corrected proof is `REVIEWED_PROOF.md`, SHA256 `582687a3d966a928efb7355cca6d9f1181d76d09d4065da17a0815d1343f457d`. `BOUNDARY_CORRECTION.diff` records the exact change. Only the interpretation of the physical Z difference was corrected; the main theorem, graph criterion, construction and producer code are unchanged.

The original Section 1 defined D for nonempty U and then identified it with B(ZU)-2B(U) without excluding the singleton. That identification is false at the singleton: B(U)=1, B(ZU)=3, the physical difference is 1, but the formula-defined D=R+b-e-B is 0. The corrected paragraph restricts the identification to nontrivial U and states both singleton values explicitly. The independent checker directly verifies this boundary.

## Exact theorem and necessity

For a nontrivial positive-start rooted tournament, put

    B=h+L-s,  D=R+b-e-B,
    p=d-e,   f=R-h+e-d.

The cut-after-root and singleton-root families of L covers have sizes h and d, and their intersection has size s. Therefore L>=h+d-s. Root deletion gives s<=d. In particular L=h implies s=d, not s=1.

Deleting the terminal root from an L cover gives precisely the labeled two-covers of U-v whose first path is empty or ends at an in-neighbor of the root. The complementary covers have a nonempty first path ending at an outneighbor a, and are counted by L_a(U-v). Thus

    b-L=sum_(a: v->a) L_a(U-v).

Each summand is at least d=H(U-v), by cutting Hamilton paths after a. Positive start on a nontrivial graph guarantees at least one such outneighbor. The root-ending deletion inequality gives p>=0. The cut-before-root and singleton-root families of R have sizes h,d and intersection e, giving f>=0. The exact rearrangement

    D-s=f+2p+(b-L-d)

therefore proves D>=s with no formal-state realizability assumption.

Equality forces exactly one outneighbor a and L_a(U-v)=d. With that unique outneighbor,

    s=start_a(U-v),  p=end_a(U-v),  f=0.

The last identity follows from the precise complement counted by f: its nonsingleton designated root-start path must contain a, while its disjoint ordinary path would have to end at a. This is impossible. Consequently the three stated conditions are an exact iff test for an individual equality graph.

Applying L_a=d=H(U-v) to the first union inequality gives

    s=H(U-{v,a}).

The residual is nonempty. If it were empty, U would be the rooted directed edge and end_a(U-v)=1, contradicting the equality criterion. This proves the necessary ordinary-count condition at arbitrary source order.

The positive-start and nontrivial hypotheses cannot be dropped from the displayed floor argument. The singleton exception is now handled explicitly.

## Strong physical sufficiency

For any nonempty tournament T with t=H(T), form the cyclic three-block source

    v->a, a->T, T->v,

rooted at v. It is strong even when T is weak: every old vertex reaches v, v reaches every old vertex through a, and these routes connect the whole graph.

Let c be the count of two labeled spanning paths of T, allowing empty bins. The direct disjoint decompositions in the proof are correct:

    (h,s,e,L,R,d,b)=(t+c,t,t,2c-t,t+c,t,2c),
    B=3c-t,  D=t.

For h, paths starting at v and a contribute t each; paths starting in T contribute c-t. For L, the root singleton contributes t, and the two disjoint placements of a contribute c-t each. For R, the root singleton contributes t and the remaining root-start covers contribute c. For b, a begins either labeled bin, giving 2c. Empty ordinary paths are included in exactly the stated cases.

This constructs a strong equality source for every attainable ordinary t and completes all three directions of the main equivalence. It does not freely prescribe B: this source has the particular value 3C(T)-t. The inherited four-core readout gives 9C(T)-t, not arbitrary target coverage.

The directed-triangle example has t=3,c=12,h=15,L=21,s=D=3. The separate source a->T has L_a=H=3 and start_a=3, correctly refuting the tempting start-one inference.

## Additional audit observation

Every equality source is itself strong. Equality forces p=0, so e=d>0, while s>0 is a hypothesis. A positive start places the root in the initial SCC; a positive end places it in the final SCC. Those SCCs must coincide, so the whole tournament is strong. This short observation is not needed for the frozen theorem's existence equivalence and does not classify equality graphs up to isomorphism.

## Independent controls and dependencies

The independent C++ checker enumerates literal directed path prefixes on every subset. It imports no producer subset-DP routine. Optimized and debug receipts agree on:

- All 33,867 labeled tournaments through order 6 and all 202,013 root choices.
- 164,030 nontrivial positive-start inputs.
- All 3,670 equality inputs, their exact individual-graph criteria, their nonempty residual counts, and their strongness.
- 3,806 physical Z children through order 5.
- 33,867 actual cyclic constructions, recounting every coordinate in the displayed formula.
- The singleton formula-defined D=0 versus physical difference 1, and the false start-one inference hostile.

The equality counts by orders 3,4,5,6 are 6,24,280,3360. Observed equality starts through this finite range are 1,3,5; no exhaustive all-order conclusion about other values is drawn from that finite list.

A separate Python checker enumerates actual Hamilton orders and ordered two-cover sets on all 285 marked labeled inputs through order 4. It checks both union injections and exact intersections, the exact f complement, the L-deletion complement bijection, and the unique-outneighbor f=0 statement. Its 8,787 literal cover visits agree under normal and optimized Python.

The producer's optimized/debug verifiers were also replayed and match the original frozen receipts. No validation gate uses removable assertions.

The all-order proof of the floor and equality criterion is self-contained in its path-cover arguments. The physical Z transfer and four-core readout are inherited from the published, independently reviewed sources frozen under `reviewed-inputs/dependencies/`, with their exact immutable provenance. In particular the transfer source explicitly contains the singleton Z start exception that motivated the correction above.

The later cover-equality normal-form and density packet is separate and is not audited or promoted by this report. Neither F(8,137) nor an eventual Hamilton-path spectrum threshold follows here.

## Reproduction

From this audit directory:

    g++ -O3 -std=c++17 check.cpp -o check
    ./check > NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check.cpp -o check_debug
    ./check_debug > DEBUG.json
    cmp NORMAL.json DEBUG.json
    python3 check_bijections.py > BIJECTIONS_NORMAL.json
    python3 -O check_bijections.py > BIJECTIONS_OPTIMIZED.json
    cmp BIJECTIONS_NORMAL.json BIJECTIONS_OPTIMIZED.json
    g++ -O3 -std=c++17 reviewed-inputs/verify.cpp -o producer
    ./producer > PRODUCER_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 reviewed-inputs/verify.cpp -o producer_debug
    ./producer_debug > PRODUCER_DEBUG.json
    cmp PRODUCER_NORMAL.json PRODUCER_DEBUG.json
    cmp PRODUCER_NORMAL.json reviewed-inputs/NORMAL.json

Independent C++ source SHA256: `41745759fdf21d55685af53c58e9d74a2c940690ade6e24ff4b7587142645efb`.

Identical normal/debug receipt SHA256: `09e58e19f89809156d2bd5c663e5e5360bb8a9ed4510b0887a62bb871a9032d5`.

Literal bijection source SHA256: `05b889259cdf828369d4daaee6b820b3398f3c3d0ea46624f97eb2e136324568`.

Identical bijection receipt SHA256: `2db8ec2d9a1aee0b175dad3beaab7f70a2c8a37e643044ed30dcb6bfa1b4de96`.

The manifest freezes the corrected proof, original lineage, independent source, and exact receipts, excluding binaries. No shared checkout or producer input was edited by this audit.
