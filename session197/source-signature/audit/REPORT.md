# Independent review: source signature and sharp product pruning

**PASS. No correction.** This report binds PROOF.md SHA256d5a9b57b62f8ce54c36fca51829efeaa493339512e47395160dc78e5474b35b6 and every file in the copied eleven-file source/RELEASE_FILES.sha256. Both source receipts are f4048fc4d54ba20b7281dc4f8038513274d15f91bdee2fa5d5be7560ed125f49. Original source files were preserved.

## All-order proof review

The recovery proof exhausts the physical source alternatives. Weak sources have E=0, odd d=B and A>=3B. Strong d=E sources have even B. For a canonical singleton-last strong source, 3E>=D makes A-3B=D-3E-d<0, separating every odd-B strong candidate from the weak branch. For even B the sign of 2A-3B-2C separates the two strong sheets; at equality their candidate deletion counts agree. The singleton, triangle and transitive two-vertex boundary are handled correctly. Thus the result is injectivity of numerical response quadruples, not rooted graph classes.

The exact membership procedure correctly treats this piecewise inverse as only a candidate and checks the genuine deletion count against the independently reviewed weak-source oracle. It does not accept the ambient integer cone. The old (D,E,U) deletion ambiguity is not contradicted: the new observables change when d changes.

The original N matrix has the stated kernel, and the descended matrix has determinant-1, proving rank3. Its signature gcd is odd because A is odd; therefore division by2 in the even-B recovery branch does not lose common-divisor information. The content identity follows in both directions. The cubic inventory bound follows from injectivity into the three-coordinate box, and the source-cap enumeration uses a genuinely monotone U response and B(L)>=L. Weak tails are still retained, so a Hamilton-count bound is not mistaken for an order or cover bound.

All physical inequalities in Section4 have valid injections. In particular the source condition gives an out-neighbor of the root, allowing insertion before its first out-neighbor in every deletion path and hence D-E>=d. The nonconcatenating root-singleton covers give U-D>=d-E without any cancellation of signed counts. This last inequality in fact holds for arbitrary rooted tournaments, but the stronger insertion inequality requires the stated nontrivial source condition.

Both positive-slack expansions in Section5 are exact. Every slack coefficient is nonnegative, and N0/M0 decrease with the endpoint counts. Substitution of the odd-D integer bound E<=(D-1)/2 gives precisely the displayed additive inequalities. The leading constants15/4 and17/2 are physically sharp by the explicit two-gadget source sequence. The packet claims sharpness of these leading constants, not separate optimality of every additive coefficient. The two-source C2 meaning is a named pinned dependency, under separate full review; the algebraic inequalities themselves use only the displayed physical source facts.

No fixed linear three-state return, graph-isomorphism injection, free reachability cone, or scalar-coverage conclusion follows or is asserted.

## Independent computational checks

The audit's bank generator interleaves O and M directly, rather than generating ordinary cores first and then appending weak tails. It computes g from a direct zero-digit contribution and Q from minimal-width reverse-complement bits, rather than the source's gadget recurrence. Its entire capped response set at U<=2000 equals the author's set exactly:32549 quadruples, with sorted-set SHA25682223401b636d4b555f21f287f15854d5845d77237fbe02ae1411c933fcb5a61.

The independently generated bank tests the signature implementation and reviewed membership oracle as systems under test. All184512 candidate triples in the U<=31 rectangle agree. All10609 ordered physical pairs through U41 satisfy the quotient form and both product bounds. An additional57600 pairs in a relaxed nonnegative slack parameterization verify both algebraic expansions independently of the physical inventory.

physical.cpp imports no author code. Its subset endpoint/cover counters enumerate all1099 labelled tournaments through order5 and independently produce all1695 start-one root responses. Every response agrees with the signature reconstruction. It builds the three sharpness sources literally and both graft types for each. A separate walk inclusion-exclusion counter verifies every Hamilton count; an anchored cycle DP verifies C=2 for all six graphs, through order14.

C++ optimized and assertion/undefined-behavior-sanitizer runs agree. Python normal/-O independent controls agree. Copied source normal/-O runs reproduce the frozen receipt exactly. The independent Python comparison deliberately imports the copied signature/oracle under test and executes the copied source verifier to compare full state sets; its bank and physical data are independently generated, not claimed to be source-free.

## Reproduction

From this audit directory:

    g++ -std=c++17 -O2 physical.cpp -o /tmp/source-signature-physical
    /tmp/source-signature-physical . > physical-normal.json
    python check.py > independent-normal.json

For the checked run use -O0 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined, and run python -O check.py. The physical checker writes literal-profiles.csv and sharp-graphs.txt as reproducible outputs. check.py needs the copied source subtree with its pinned dependencies; it also runs that copied verifier, which rewrites only the copied VALIDATION.json. Preserve the original packet by running a copy.

Hashes:

- physical.cpp:4c713b96f2beee7512d0612401a71c4ddd8ee5b3fb955d7a169b5e899208c5d3.
- check.py:410c3cdb78a93960b9ec5be3639502eae2e0c55aa5c0ef98f2590dd5f87d18d2.
- literal-profiles.csv:822f897d9e17c6a9e2eacebbaf747450d3cc893f12faec39a316d4ff2f3b5def.
- sharp-graphs.txt:ec200bf3aad9c56aa164d78e0997a98190829849c73cba180d582bfd1272d385.
- independent-normal.json and independent-optimized.json:37b07925411d7950dddb59f2fdfed38fc13a7abe502b1212051377054d3d0a63.
- physical-normal.json and physical-debug.json:e051117c08328d71674632a92d1cb318a3989c759a5252ce49adff098d6357fe.
- source-normal.json and source-optimized.json:f4048fc4d54ba20b7281dc4f8038513274d15f91bdee2fa5d5be7560ed125f49.

No mathematical correction or unresolved blocker remains.
