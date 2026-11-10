# Why scalar depth-two equality does not supply an inductive normal form

Status: **PROVED graph-preserving adjacent-swap normalization and explicit
same-graph counterexamples**; **FINITE-EXACT same-order scalar obstructions**.
The fixed-order depth-two scalar-image conjecture is **REFUTED** at order 12:
H=97309 has no depth-at-most-two order-twelve history. This does not imply
absence at larger orders: a depth-one order-twenty witness is explicit below.
Arbitrary allowed-odd scalar coverage remains **OPEN**.

## 1. Graph flattening is already false at order nine

Consider the nine-vertex chronological forest history with mixed-radix code
26485. Its choices at births 1,...,8 are

    0, 0, 1, 3, 4, 3, 6, 7.

Its outgoing-neighbor bitmasks, in vertex order 0,...,8, are

    502, 508, 424, 465, 484, 328, 132, 288, 64.

The supplied history has maximum forest depth three and H=217. The same
scalar H=217 also has an order-nine depth-one history, code 723. Thus a
failure to flatten this graph does not imply a missing shallow scalar.

**Proposition.** The displayed tournament has no chronological forest
presentation of maximum depth at most two, even after arbitrary relabeling.

**Proof.** Reading any chronological presentation backwards removes a vertex
of current outdegree zero or one. If its outdegree is one, its unique
remaining outneighbor is its forest parent. In the displayed tournament the
first four deletions are forced:

    delete 8 with parent 6;
    delete 7 with parent 5;
    delete 6 with parent 2;
    delete 5 with parent 3.

Thus every presentation contains parent chains 8->6->2 and 7->5->3. A
depth-at-most-two presentation would force both vertices 2 and 3 to have
depth zero. On the remaining vertex set {0,1,2,3,4}, the only deletable
vertices are 2 with parent 3 and 4 with parent 2. Deleting 2 contradicts its
required depth zero, so 4 must be deleted. On {0,1,2,3}, the only deletable
vertices are 2 with parent 3 and 3 with parent 0. Each contradicts its own
required depth zero. This exhausts all choices. The supplied depth-three
presentation proves that the graph's minimum presentation depth is exactly
three. QED.

**All-order separation.** Append any number of sinks to this graph. Its
Hamilton-path count remains 217, and the displayed presentation still has
depth three. The depth-at-most-two presentation property is hereditary:
restricting a chronological order to an induced subset deletes parent edges
and can only decrease forest depths. Because every sink-padded graph still
contains the displayed nine-vertex graph, its minimum presentation depth
remains exactly three. Thus graph flattening fails at every order n>=9,
even on an infinite family with the fixed scalar 217.

As a separate finite control, the verifier enumerates all 9!=362880 vertex
orders. Exactly 16 are valid chronological presentations: eight have depth
three, six have depth four, and two have depth five.

## 2. Even scalar normalization does not preserve fixed-order transitions

For a tournament T, let C(T) be the set of H-values of its ten order-ten
children obtained by adjoining a sink, or a vertex beating exactly one of
its nine old vertices. The displayed parent has the following counts in
choice order 0,...,9:

    217, 2041, 2241, 653, 769, 1197, 539, 437, 501, 373.

There are exactly 176850 order-nine histories of depth at most two. Of
these, exactly 458 have H=217. These 458 histories yield 30 distinct C(T)
sets, whose union has 212 values. That union contains none of

    373, 437, 501, 539, 2041, 2241.

Consequently the single transition 217->373 supplied by code 26485, followed
by choice 9, cannot be reproduced by *any* depth-at-most-two order-nine
parent with H=217. This is stronger than failure of one chosen shallow
representative: choosing a different shallow parent for that specific root
still does not reproduce the transition at the same order.

The child scalar 373 does have a shallow order-ten witness in the main
census. Therefore this obstruction neither refutes depth-two scalar-image
equality nor excludes that scalar from the forest class. It blocks a proof
that normalizes the parent using only its scalar and fixed graph order,
then expects the same one-birth continuation to remain available. A viable
proof must retain extra root/attachment information, change the construction,
or explain why a different parent scalar/order suffices.

## 3. Verification and exact scope

`verify_normal_form.cpp` independently constructs every depth-at-most-two
history at order nine. It recounts each literal tournament by endpoint subset
DP; for every H=217 parent, it independently recounts all ten literal child
graphs. It uses neither the singleton graft identity nor parent DP-table
reuse. A separate integer recurrence verifies the total shallow history
count: if r vertices have depth zero and s have depth one, the next birth
has one sink choice, r depth-one choices, and s depth-two choices.

`NORMAL_FORM_CHECK.json` contains the full 30-signature certificate, the
exact absent-child list, and the independent all-permutations depth count.
Compile and run with the same strict flags as the main census. There are no
random inputs, heuristic filters, or assert-dependent checks.

The exploratory discovery pass tested child-set normalization through order
eight before finding this order-nine obstruction. That discovery history is
not needed for either proposition, and no claim of globally minimal order
or minimal scalar is based on it here.

## 4. The scalar conjecture itself fails at order twelve

The exhaustive 12!=479001600-history census has 53549 distinct H-values.
Exactly one of them, H=97309, has minimum history depth three at that order.
A positive witness is code 25112126, with its full literal graph, root
responses, endpoint counts, and deletion counts in DEPTH_EXCEPTION_CHECK.json.
Every other order-twelve scalar has a depth-at-most-two witness. All orders
through eleven satisfy depth-two scalar equality, so order twelve is the
first failure of that fixed-order statement.

The literal graph's minimum presentation depth three also has a short proof.
Its first eight reverse deletions are forced:

    11->1, 10->3, 9->2, 8->6, 7->0, 6->4, 5->0, 4->0.

Here v->u means delete v with its unique remaining outneighbor u. Thus every
presentation contains the parent chain 8->6->4->0, requiring depth at least
three. The supplied history attains three. This graph-specific argument
alone would not exclude a different shallow graph with the same scalar.

The independent scalar exclusion instead comes from shallow_census.cpp. It
constructs only depth-at-most-two histories and never stores depth minima.
Using the graph-preserving normalization in Section 6, it exhausts 23670073
normalized histories at order twelve and obtains exactly 53548 values.
Their complete image equals the main census's recorded shallow image and
omits exactly 97309. An independent integer recurrence checks the normalized
history count. NORMALIZED_SHALLOW_AUDIT.json records the comparison.

## 5. A necessary variable-order positive control

The scalar 97309 is not a global shallow obstruction:

    97309 = 31 * 43 * 73.

Literal chronological histories (order,code)=(6,50),(7,296),(7,58) have
respective path counts 31,43,73 and each has maximum forest depth one.
Their ordinal join, in that order, is a twenty-vertex history of depth one.
Every directed Hamilton path in an ordinal join traverses its components
consecutively, since no edge returns from a later component to an earlier
one. Concatenation is therefore a bijection with one Hamilton path from
each component, proving the product count 97309.

verify_depth_exception.py independently recounts each factor, checks every
cross edge of the literal joined graph, checks its depth-one presentation,
and verifies the product. DEPTH_EXCEPTION_CHECK.json supplies its complete
adjacency and mixed-radix code. It does not claim a direct full twenty-vertex
DP, nor minimum possible order for a depth-one witness.

Thus the failed statement is equality of scalar images at every *fixed*
order. Variable-order shallow constructions, multiplicative closure, and
prime/primitive scalar targets remain separate research questions. The
single finite counterexample gives no global missing allowed odd value.

## 6. A genuine all-order normalization that survives

Suppose consecutive vertices x,y in a chronological forest presentation
satisfy y->x. Then x is the unique earlier outneighbor of y. Swap x and y
while keeping the tournament graph fixed. Only their mutual relative order
changes, so the reversal edge y->x disappears and no other reversal edge is
added or removed. The vertex y becomes a forest root; x keeps its old parent;
all other parent incidences are retained. Consequently every vertex's depth
stays the same or decreases.

Repeat until no adjacent reversal remains. The reversal forest has at most
n-1 edges, and each swap removes one, so termination takes at most n-1 swaps.
The terminal birth order is itself a directed Hamilton path. At birth j,
its allowed parents are earlier vertices other than the immediate
predecessor, giving exactly j possibilities including the sink choice.
Hence there are (n-1)! normalized histories, and every original graph has
one without increasing its presentation depth. This preserves the graph
and all root responses; it is much stronger than scalar preservation but
much weaker than a uniform depth-two bound.

For shallow normalized enumeration track r depth-zero vertices, s depth-one
vertices, and the last vertex's depth l. The next birth has one sink choice,
r-[l=0] depth-one choices, and s-[l=1] depth-two choices. This yields the
independent expected total 23670073 at order twelve. Without normalization
the coefficients are instead 1,r,s, giving 123489888 shallow histories.

**General bounded-outdegree corollary.** The same adjacent swap works for a
chronological order in which every newborn beats at most k earlier vertices.
It removes only the mutual reversal edge: y loses one earlier outneighbor,
x retains its old earlier outneighbors, and every other vertex retains its
set. Thus the earlier-outneighbor bound is preserved, and repeated swaps
produce a Hamilton-path chronological order with the same bound. Its
normalized history count is

    product_(j=1)^(n-1) sum_(i=0)^min(k,j-1) binomial(j-1,i).

The forest case k=1 gives (n-1)!. The depth-nonincreasing assertion is specific
to the one-parent forest case; no analogous unqualified depth claim is made
for multiple-parent histories. This is an elementary adjacent-swap application
of the Hamilton-path ordering principle; no external priority claim is made.
