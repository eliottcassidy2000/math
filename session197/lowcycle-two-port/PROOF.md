# A faithful two-port compiler and exact low-cycle value language

Session197. **All-order graph compilation and an exact finite-control inverse-language criterion.** The scalar coverage statement is proved only through its stated bound. The full Hamilton-path spectrum remains open.

Define

    L = {1} union {H(T): T is strong and C(T) is1 or2}.

The strongness restriction is essential: adding a universal sink preserves H but destroys every Hamilton cycle, so allowing arbitrary nonstrong C=0 graphs would trivialize the cycle condition. No claim that L contains every odd integer except7,21,67 is assumed.

## 1. Physical language and dependencies

Start from a forward edge with ordered ports[0,1]. At each right birth, a new vertex beats exactly a selected subset of the two current ports and loses to every other old vertex. Temporarily append it to the port list and forget one of those three entries, keeping the other two in inherited order. Forgotten vertices remain in the graph and can never be exposed again. There are twelve moves, indexed by a two-bit mask and a dropped position.

The exact nonnegative27-coordinate path-cover interface is inherited from the independently reviewed session197-hybrid-constructions/packet. It retains covers with at most one nonport start. The unchanged implementation is included as twoport_core.py, SHA256733155121ac7a5936a38a0ec1f3a602759f868c1e4e547fd4b6f89bced25f9cb. Hamilton paths are the one-component coordinates.

The structural inputs are the independently reviewed C1_COMPLETE.md, START_TWO.md, BRIDGE_TWO.md and C2_COMPLETE.md in session197-lowcycle. They supply the canonical start-one/start-two decompositions and the fact that every strong C2 graph except the regular five-tournament has a start-two root. No additional scalar coverage theorem is imported.

## 2. Preserving prescribed start-one vertices

**Theorem1.** Every strong C1 tournament can be constructed in this two-port language with any prescribed set of one or two distinct start-one vertices exposed at the end. The compiler preserves the entire tournament, all original vertex identities, and that final exposed set. It returns the actual order of the final ports; it does not promise an arbitrary separately prescribed port order.

Allow the singleton as the induction base, with its sole vertex exposed. Choose one requested vertex b. In the canonical start-one decomposition, deleting b leaves an initial strong component U, possibly a singleton, with exactly one entry a beaten by b, followed by a transitive tail W. All other U vertices beat b, and every U vertex beats every W vertex. Write a tail bit1 for b→w and0 for w→b.

Two consequences of the exact old-root start-count law in START_TWO.md are needed, not merely the unmarked Fibonacci enumeration:

- If W is nonempty, a requested old start-one vertex lying in U can only be a. Every x in U−{a} has new start count at least3.
- If W is empty, any requested x in U other than a that has start count one in the new graph already had start count one in U. Indeed the law is s_x−p_(x,a)+D_(x,a), with D>=1. Equality1 forces s_x=p_(x,a). In particular the Hamilton path obtained by cutting the unique old cycle at x ends at a, so a→x is its closing edge. Every x-start path therefore closes to that unique cycle, and s_x=1.

Here p counts Hamilton paths starting at x and ending at a, and D counts the ordered two-covers starting at x and a. For a singleton U there is no other x. These are precisely the rooted inheritance facts ensuring an old requested vertex is not forgotten and later illicitly re-exposed.

If the other requested vertex x lies in U and W is empty, recursively compile U exposing{a,x}, or just{a} when x=a. Otherwise compile U exposing a. If W is nonempty and x lies in U, the first bullet says x=a, so no second old port is needed.

Now add b beating only a. Retain b and any requested old x. Append W in its transitive order. Each new tail vertex beats only b when its bit is0, and beats no old vertex when its bit is1. Retain b and, as soon as it is present, the requested x. At most two vertices are ever required simultaneously. If x is a tail vertex it can be kept while later tail vertices are born and forgotten; their only possible old outneighbor is b. This constructs every required edge and finishes with the prescribed set exposed.

Every recursive call is on a strictly smaller U. Thus the procedure is an all-order induction, not an appeal to the finite tests. A single requested start-one vertex is handled by the same procedure without a second required port.

### Standard forward-edge seed

The recursive construction may initially produce its first two vertices in reverse port orientation. Before executing the protocol, order these two initial vertices so their edge points forward. Thereafter re-index each mask and dropped position using the same desired retained SET as in the original protocol. Immediately before each later birth the exposed set is unchanged, so its required old outneighbors are still available. The inherited port order is then respected at every actual move. This is a compilation-time normalization of the initial edge, not a later swap or re-exposure. compiler.py performs and verifies this normalization explicitly.

## 3. Start-two graphs and the regular-five exception

**Theorem2.** Every strong rooted graph with two Hamilton paths starting at r has an exact two-port construction ending with r exposed.

Delete r and use the complete START_TWO decomposition. Its initial strong component U has either one entry a with start count two, or an unordered pair of start-one entries a,b. All later SCCs are singleton vertices forming a transitive tail.

In the one-entry case recursively compile(U,a) by this theorem. In the pair case use Theorem1 to expose the exact pair{a,b}. Add r beating precisely the one or two entries, retaining r. Stream the tail as in Theorem1, with each new vertex beating only r or no old vertex. Every incidence is exact, and r stays exposed. The pair case supplies the induction base; the one-entry case reduces the order.

The completed C2 classification says that every strong C2 graph except R5 has such a root, for example any outdegree-one vertex. R5 itself has rows[6,12,24,17,3]. From its forward edge on0,1, retain those two ports and append vertices2,3,4 with masks0,1,3, always forgetting the newborn. This gives R5 exactly. It has start count3 at every vertex, so it is a genuine exception to the start-two compiler, handled separately rather than silently included in it.

Consequently **every strong C1 or C2 tournament occurs in the physical two-port language**. Construction uses exactly n−2 births after the initial edge. Intermediate graphs may be weak; no intermediate-strongness claim is made. Every final graph is the exact input graph, not merely another graph with equal H.

## 4. Exact Hamilton-cycle readout

Let v be the27-coordinate state before a birth with selected-port mask S. Delete the newborn from a Hamilton cycle of the output. The remaining object is one Hamilton path of the old graph whose start is a selected port and whose end is not a selected port. Conversely each such path closes uniquely through the newborn. Fixing the newborn first removes rotation ambiguity.

Therefore the new cycle count is the nonnegative readout

    C_new = sum v_sigma over one-component sigma=(s,e)
            with s a selected port and e unselected or nonport.    (1)

It is independent of which port is then forgotten. In the supplied signature order, the one-component coordinates are indexed by(s,e) with s,e in{0,1,*}, in row-major order. Thus

    S=0:0;
    S=1:v_1+v_2;
    S=2:v_3+v_5;
    S=3:v_2+v_5.                                      (2)

These are actual directed-cycle counts, not values inferred from H. The formula also gives0 at the two-vertex boundary. A positive result is equivalent to strongness once the output has at least three vertices.

### Why a C>2 prefix can be rejected

Suppose a prefix graph P has C(P)>2. It is itself strong because it has a Hamilton cycle. If a continuation ended in a strong graph T with C(T)<=2, then P would be a strong induced subtournament of T.

Hamilton-cycle count is monotone between such strong induced graphs. One can grow P inside T by mixed singleton insertions, which inject cycles, or by a two-vertex ear P→q→p→P, whose cycle count is H(P)>=|P|C(P). Repeating these strong extensions proves C(T)>=C(P), a contradiction.

This compares the strong prefix directly with the final strong target. It does **not** assert monotonicity through the actual intermediate graphs, which may be weak and have cycle count0. Thus rejecting C>2 prefixes loses no eventual strong C1/C2 output.

## 5. Uniform finite control, independent of a target bound

For any nonnegative integer vector z, let cap_c(z) clip each coordinate at c. Nonnegative integral transfers satisfy

    cap_c(M z)=cap_c(M cap_c(z)).                      (3)

This follows coordinate by coordinate from nonnegative addition and multiplication: once a positive summand exceeds c, the sum is at least c; otherwise all relevant terms are exact. The same statement applies to the nonnegative cycle readout.

Use c=3. Every birth's value0,1,2 or at least3 is therefore determined exactly by the capped state. Starting from the forward-edge state, reject a transition when (1) exceeds2 and otherwise update by cap3. The complete reachable control closure has1,369 vectors. Add an absorbing rejecting state. Label an edge1 exactly when the newly born graph has C1 orC2, and0 otherwise. Partition refinement by all twelve output/next-class pairs gives a39-class edge-output machine, including one rejecting class. Every other class can reach an edge labelled1.

control.json contains every raw vector, every transition and cycle category, the quotient map and all quotient edges. Refinement stabilizes after three rounds; control.py regenerates the complete closed certificate. This dimension concerns finite cycle admissibility, not the minimal dimension of the unbounded H representation.

The raw automaton is justified at arbitrary construction depth by (3), not a word-length cutoff. The quotient preserves every future accepting edge by its stable transition/output identities. Combining it with the exact27-coordinate integer transfers gives the uniform inverse-language criterion

    h in L iff h=1, or there is a word accepted on its final edge
    whose exact Hamilton readout is h.                (4)

Soundness follows because an accepted physical graph is strong C1/C2. Completeness follows from Theorems1–2 and the explicit R5 construction, together with safe prefix pruning. Thus this criterion describes the FULL strong low-cycle scalar class, not just the outputs of a proper gadget subgrammar. It remains a matrix/finite-control criterion; it is not a solved scalar backward recurrence.

## 6. Exact bounded-output closure

For a target bound X, cap all27 coordinates at X+1 and discard states whose Hamilton readout exceeds X. Hamilton paths are monotone under induced extension, by deterministic vertex insertion, so this discards no later target<=X.

For every retained state, the one-component coordinates contributing to H are exact: if any such coordinate had clipped, H would exceed X. Equation(1) uses a subset of these same coordinates, so the cycle count readout is exact there as well. There is no signed readout, cancellation, or hidden headroom assumption. Reject C>2 transitions as justified above, and record H whenever C_new is1 or2, before state deduplication.

Recording before deduplication matters: the current graph's cycle count is an output of its incoming edge and need not be stored in the successor vector. Equal capped vectors have exactly the same capped future transfers and cycle tests; their different histories need not be identified as graphs.

The finite search is complete at arbitrary depth by (3) and the finitely many capped vectors. A resource-limit exit is explicitly marked incomplete and supplies no exclusion certificate. The frozen runs are complete:

    X=173:386,077 states,1,384,593 admissible edges,
          maximum shortest word length19;
    X=501:3,455,583 states,12,232,554 admissible edges,
          maximum shortest word length21.

At both bounds the missing odd values from3 through X are exactly7,21,67. Every other retained value has an explicit physical word and final cycle count in the receipt. These finite bounds do not imply eventual or all-target coverage.

## 7. Conditional route to full two-port coverage

The independent strong H67 witness has rows

    [122,116,105,82,68,24,32].

Its width-two word is0011132, where a mask's low bit refers to the older of the two previous vertices. Keep the two most recent vertices, forgetting the oldest at every later birth. Literal enumeration gives67 Hamilton paths and4 Hamilton cycles. Thus it is a valid two-port witness outside the low-cycle class.

Consequently the following implication is proved, while its hypothesis remains OPEN:

    If L contains every positive odd integer outside{7,21,67},
    then the two-port language contains every positive odd integer
    outside{7,21}.

The extra67 witness fills precisely that excluded value, and H1 is the seed's count. No finite support computation is substituted for the hypothesis. The known unrestricted omissions7 and21 are not newly claimed here.

The remaining arithmetic problem is to solve (4) uniformly, for example by a descending inverse on actual reachable profiles or a finite family of certified inverse branches. Neither gcd-one observations nor the finite control count alone establishes such a descent. Replacing a graph by an equal-count graph also requires an actual physical construction or a separate compiler; it is not a port operation by itself.

## 8. Reproduction and exact scope

compiler.py implements the two simultaneous promised-input inductions, retains original labels, records every mask/drop, normalizes the initial edge, and replays every adjacency entry. verify.py checks all labelled tournaments through order6:30,036 prescribed start-one sets,13,584 start-two roots, and all24 labelled R5 presentations. The H67 witness is counted literally. Normal and optimized receipts agree.

control.py regenerates the raw cycle-control closure and39-class quotient. support.cpp performs complete bounded-output closure using the serialized audited transfer. verify_support.py reconstructs that serialized table from the combinatorial generator, tests cycle readouts literally on all small marked graphs, independently repeats the X21 closure in Python, reruns the complete X173 and X501 closures, and replays every stored word with uncapped integer transfers and literal graph births. It compiles only the adjacent C++17 source in a temporary directory. Normal and Python-optimized receipts agree.

The all-order results are the exact graph compiler, safe cycle pruning, uniform finite control and bounded-target decision criterion. The numerical coverage remains bounded by501. No unconditional cofinite low-cycle support, arbitrary-target two-port realization, or unrestricted spectrum theorem is proved.

## Frozen structural dependency hashes

The exact independently reviewed proof snapshots used above are:

- C1_COMPLETE.md: d38a53b1c475df3578d58e83d9e9f97785d9455084c0a094aaed752577d8b289.
- START_TWO.md: a746d2e33b0fb78a98252afa602dc8985c25f718225c2f9ce37923fd20b0ada3.
- BRIDGE_TWO.md: 9542e1a23488d2b3ee308ad2aa79cd028a1e28ea7482eed107b43f4f7bdef652.
- C2_COMPLETE.md: 4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d.

All four are in session197-lowcycle. The strong-induced cycle monotonicity argument is restated in Section4, so safe pruning does not depend on a finite census. The raw transfer generator is included unchanged; no external data or package is required by either verifier.
