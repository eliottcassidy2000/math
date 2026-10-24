# Conditional all-size occupied-sheet gluing, with a native deletion hostile

Status: **PROVED** combinatorial statement for any declared binary supply;
**FINITE-EXACT** deletion test in the authenticated28-source native menu.
This supplies an induction interface, not the missing all-order native supply.

## 1. The correct deletion-relative sheets

Let S be any finite inventory of independent capacity-one binary source
choices, with additive labels g_s in an abelian group, positive event weights
c_s, and fixed-margin fibre

    F_S(m,q) = {A subset S : |A|=m, sum_(s in A) g_s=q}.

Two states are adjacent when |A triangle B|<=r and
sum_(s in A triangle B)c_s<=B. Extra literal native constraints may be
included only if they are retained in each restricted domain below; no
unverified quotient is substituted for them.

For h in S, the occupied sheet is

    F_h = {A in F_S(m,q) : h in A}
        = {T union {h} : T in F_(S minus {h})(m-1,q-g_h)}.

Adding h preserves every edge's support and work because h is unchanged.
Thus connectivity of F_h is exactly connectivity of the DELETED-inventory
fibre on S minus {h}, not of the full-inventory fibre at margin m-1.

Treat a negative requested subset size as an empty fibre. For distinct h,j
the intersection is

    F_h intersect F_j
      = {T union {h,j} : T in F_(S minus {h,j})(m-2,q-g_h-g_j)}.

Hence overlap is a concrete binary subset-sum existence condition in a twice
deleted supply. Integer charge-span coverage does not establish it.

## 2. Gluing theorem and precise converse

Suppose m>=1. Discard empty sheets. If every remaining sheet is connected at
(r,B), and the graph on sheets with an edge for nonempty intersection is
connected, then F_S(m,q) is connected at the SAME(r,B).

Proof: each state lies in an occupied sheet. Follow an overlap-graph path
between its starting and ending sheets. Choose a state in each successive
intersection, and join successive chosen states inside the current sheet.
Their actual lifted edges keep h fixed and retain the same event costs.
No simultaneous move between disjoint sheets is assumed.

If r<2m, the converse also holds under the sheet-connectivity hypothesis.
Indeed any permitted edge A--B has

    |A intersect B| = m - |A triangle B|/2 > 0.

It therefore lies in at least one occupied sheet. A path in the full fibre
induces an overlap walk, with shared vertices passing between sheet labels.
Under joint support8, this exact-cover regime starts at margin5.

The restriction is essential. At margin4 a direct disjoint eight-slot edge
belongs to no occupied sheet. The complete two-state sharp22 witness is
connected by that edge although its occupied-sheet overlap graph is not.
Thus overlap connectivity is always sufficient under connected sheets, but
is not necessary when permitted edges may have disjoint endpoints.

## 3. Exact component-nerve version

Without assuming sheets connected, replace each sheet label h by all pairs
(h,C), where C is an actual connected component of its deleted-inventory
fibre after h is restored. Join (h,C) and (j,D) iff C and D share an actual
state. If r<2m, connected components of this component-overlap graph are in
bijection with connected components of the full fibre.

Proof: each local component is connected in the full graph; intersections
join those connected pieces. Conversely every full-graph edge is inside a
sheet component by the common-occupied-source argument. Thus the full graph
is precisely the union of these component subgraphs, and intersection gluing
accounts for every possible connection. This is an exact recursion on
smaller supplied inventories, not a claim that their components are trivial.

For a quantitative sufficient condition, if every sheet has diameter at most
D and the sheet-overlap graph has diameter at mostd, the construction gives
full-fibre diameter at most(d+1)D. Neither D nor d is asserted to be independent
of native supply size. Actual literal resources, pins and any charged shared
owner constraints still have to be authenticated before applying this model.

## 4. A native counterexample to deletion-hereditary22 connectivity

The complete original28-source inventory has uniform22-clock connectivity.
Yet this property does not pass to every source deletion.

Delete source2071 and consider margin4, old-basis charge(45,11). The ENTIRE
remaining27-source fibre consists of exactly three states:

    A={1394,1921,1923,2905},
    B={1532,1841,2835,3360},
    C={1691,1923,1974,2837}.

A and C share source1923; their event changes six cost-two slots, so costs12.
B is disjoint from both. B has clock weight13, while A and C each have weight10.
Both incident events at B therefore cost23. The fibre is disconnected at22
and connected at23. All states are actual inherited source/path mixtures.

Equivalently, in the complete28-source margin5 fibre of charge(56,2), the
occupied sheet with source2071 fixed at V has these three states with2071
adjoined. That sheet is disconnected at22, although the complete margin5
fibre is connected. A successful full-fibre route may have to release2071;
full-inventory connectivity cannot be used as if it proved sheet connectivity.

The complete single-deletion margin4 probe enumerates all17,550 states in
each of28 deleted banks,491,400 states in total. Exactly seven deleted sources
produce a failing fibre at22:

    1477,2071,2087,2837,2945,3000,3369.

There are eight failing fibres because deleting2071 produces two. Each has
two22-clock components, with a23-clock crossing. Thus23 is the sharp uniform
clock bound over all single-deletion margin4 fibres in this supplied menu.
The source list and every complete hostile are saved in `deletion_probe.json`.
A separate direct27-subset/Kruskal implementation verifies every fibre
again, including exact bottleneck histograms; its86-check normal/optimized
receipts are byte-identical.

The conditional theorem identifies the true all-size obligation: control
connectivity in deletion-relative supplies and the charge-supported overlap
relations, or retain their nontrivial components in the exact quotient.
The finite full28 theorem alone proves neither hypothesis at other sizes.
