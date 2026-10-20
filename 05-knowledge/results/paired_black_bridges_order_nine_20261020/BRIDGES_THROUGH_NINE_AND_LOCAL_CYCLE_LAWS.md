# Paired-black bridges through order nine, local cycle laws, and the remaining all-order gap

Research date: 2026-10-20.

Status: FINITE-EXACT full order-nine classification; PROVED local sufficient criterion; CITED all-order combined connectivity and its black-component consequence. The assertion that every paired-black bridge is pendant at every order remains OPEN.

## 1. Three graphs must not be conflated

The vertices of the merged metagraph are tournament isomorphism classes modulo converse. With a fixed directed Hamiltonian path, let t encode the free arcs, C complement every free bit, and R reflect positions. A literal line is {t,Ct}. It is blue if Rt=t and black otherwise.

- K is the literal black multigraph. Every black line has a distinct reflected partner with the same merged endpoints. Consequently K has **no ordinary bridges at any order**.
- Q identifies each pair of reflected black lines. Parallel multiplicity and loop multiplicity in K are divided by two. A loop contributes two to degree, and deg_Q(u)=(M_u-B_u)/2.
- The actual blue graph retains the reflection-fixed lines. Its degrees are B_u. It is not the auxiliary hypercube on anti-involution frames, and its bridges are not Q bridges.

The present complete bridge census concerns Q. A separate order-nine blue calculation is included as a control, not substituted for the Q question.

## 2. Complete order-nine result

Every one of the 2^28=268,435,456 fixed-path tilings was classified, without sampling, score merging, or pruning. The exact universe contains:

| Quantity | Value |
|---|---:|
| Ordinary isomorphism classes | 191,536 |
| Converse-merged nodes | 97,144 |
| Self-converse nodes | 2,752 |
| Q edges | 67,092,480 |
| Q loops | 4,499 |
| Active Q vertices | 97,138 |
| Q components, including isolates | 7 |
| Q cycle rank | 66,995,343 |
| Q bridges | 7 |

There is one active Q component and six pure-blue isolates. All seven bridges are leaf edges. Every leaf is non-self-converse and belongs to the already proved nonpalindromic singleton/C3 word family. There are no mixed Q leaves at this order.

The complete bridge representatives in explorer coordinates are:

| Orbit representative t | Leaf node | Other node | Other Q degree |
|---:|---:|---:|---:|
| 64 | 58 | 96545 | 97 |
| 4096 | 2964 | 96655 | 141 |
| 131072 | 42880 | 96554 | 137 |
| 2097216 | 96524 | 43209 | 131 |
| 16777280 | 97113 | 48851 | 123 |
| 16781312 | 97114 | 50272 | 187 |
| 67108928 | 97138 | 48911 | 131 |

Node IDs are ranks by minimum merged tiling mask in this census. They are not the canonical IDs in historical atlases.

The prior complete n<=8 census and this result together classify **all Q bridges through order nine**: each is incident to a leaf. They do not prove that a first nonleaf bridge cannot occur at order ten or above.

### Blue graph control

The actual order-nine blue graph has 32,768 lines, no loops, and one active component on all 2,752 self-converse nodes. It has exactly 24 bridges, each incident to a blue-degree-one vertex. These leaves can be mixed rather than pure-blue. The 94,392 NS vertices are blue isolates, giving 94,393 total blue components.

## 3. A stronger finite classification: every nonleaf edge has a short cycle

The full Q adjacency was reconstructed independently from the tiling map, with multiplicities kept. Every edge falls into exactly one of the following categories, in the indicated priority order:

| Category | Number of Q edges |
|---|---:|
| Loop | 4,499 |
| Member of a parallel family, hence a digon | 8,461,023 |
| Remaining edge in a triangle | 58,599,427 |
| Remaining edge in a quadrilateral | 27,524 |
| Leaf edge | 7 |
| Uncovered edge | 0 |

The total is 67,092,480. All 27,524 quadrilateral witnesses are recorded explicitly in census_n9_quadrilaterals.tsv. The support graph has 62,645,527 nonloop endpoint pairs. A quadrilateral witness (u,v,x,y) means the four distinct vertices carry edges uv, ux, xy, yv.

The same short-cycle property was independently checked through order eight. It suggests the sharper research target:

**Short-cycle conjecture.** Every Q edge that is not a leaf edge lies on a cycle of length at most four, allowing loops and digons.

This conjecture implies the pendant-bridge conjecture but is not proved here. The finite data do not turn it into an all-order classification.

## 4. Complete enumeration and independent bridge certificates

### 4.1. Why Hamiltonian-path enumeration is complete

For each as-yet unseen tiling t, reconstruct its tournament T. Enumerate all directed Hamiltonian vertex words P of T, and re-express every free arc in P coordinates. The resulting encodings are exactly the ordinary isomorphism fibre: an isomorphism to a fixed-path presentation pulls the fixed path back to such a P, and every such P supplies a presentation.

Adjoin each reflected encoding to obtain the converse-merged fibre. Distinct discovered fibres are asserted disjoint; all masks are assigned exactly once. Every tournament occurs because every tournament has a Hamiltonian path. Thus exhaustion of the tiling cube proves completeness, not just completeness relative to an externally supplied class list.

The DFS accumulates each tiling mask incrementally when a vertex is appended. It also counts all Hamiltonian words H. The ordinary fibre mass is M for SC nodes and M/2 for NS nodes, so the computed automorphism order is H divided by that mass. Oddness, integrality, blue parity, degree identities, and total mass are checked exactly.

### 4.2. Exact black orbit selection

C and R commute. On a black tiling the orbit {t,Ct,Rt,CRt} has size four. In particular Rt=Ct cannot occur, because reflection fixes the apex coordinate whereas complementation reverses it. Keeping the least mask in each such orbit gives exactly one Q edge, including parallel lines and loops. This is not support-graph simplification.

### 4.3. Three bridge checks

1. census_stream.cpp inserts every Q edge into an online two-edge-connectivity algorithm.
2. sparse_bridge_certificate.cpp scans all Q edges in reverse mask order and greedily builds two edge-disjoint forests. Ordinary low-link DFS on their union finds the same seven bridge masks. The certificate has 194,267 edges, versus 67,092,480 in Q.
3. short_cycle_audit.cpp constructs all actual multiplicities and verifies that every nonleaf edge has a loop, digon, triangle, or quadrilateral witness. No bridge algorithm is used in this third check.

The producer map and all node records were matched pointwise against the independent prior Python census at n=7 and n=8 before the n=9 run.

**Two-forest certificate lemma.** Let F1 receive a scanned edge when it joins two F1 components. Otherwise let F2 receive it when it joins two F2 components. Discard all remaining edges. Then H=F1 union F2 has exactly the bridges and connected components of the full multigraph G.

Proof. Each discarded edge has its endpoints connected by a path in F1 and by a path in F2. Every cut crossed by this edge therefore has at least one crossing edge from each forest, and the two edges are distinct. Consequently a cut of H with size zero or one cannot be crossed by any discarded edge. Conversely a G cut of size zero or one is unchanged in H. These facts preserve connected components and all one-edge cuts, hence exactly all bridges. Loops never cross a cut and cause no problem. QED.

This is an elementary graph certificate, not an unproved tournament-specific theorem.

## 5. A proved local sufficient criterion, with two essential hostile controls

For a directed Hamiltonian path P in T, write T^P for the tournament obtained by reversing precisely the n-1 arcs of P. The complement partner of the fixed-path encoding is (T^P)^op. Thus, after merging converse, a metagraph move is exactly Hamiltonian-path arc reversal.

Let P and Q be arc-disjoint directed Hamiltonian paths in T. Reversing P leaves Q directed, and reversing Q leaves P directed. Hence the four tournaments

    T, T^P, T^(P union Q), T^Q

form a closed four-move walk in the combined merged metagraph.

**Black commuting-square criterion.** If all four moves are black and the Q-edge represented by (T,P) occurs only once in that projected closed walk, then this edge is not a Q bridge.

Proof. The four moves give a closed walk entirely in Q. A bridge must be traversed an even number of times in every closed walk, because every traversal switches shores. The specified edge is traversed once, a contradiction. This also proves existence of a cycle containing it of length at most four, after deleting repeated subwalks. QED.

Both hypotheses after arc-disjointness matter:

- At n=5, t=2, take P=(0,1,2,3,4) and Q=(0,2,4,1,3). They are arc-disjoint. Encodings for the four moves are 2,10,52,20; the middle two are blue. This does not supply a black square.
- At n=5, t=14, take P=(0,1,2,3,4) and Q=(1,3,0,4,2). All four moves are black, with encodings 14,25,19,42. But their Q orbit representatives are 14,19,19,14. The square folds to a backtracking walk and gives no nonbridge certificate for edge14.

Thus a proof that merely produces a second arc-disjoint path is insufficient. The required sidecars are the four blackness tests and the quotient edge identities. The computer probe finds valid commuting-square certificates for 7,303 of the 8,004 nonleaf, nonloop Q edges at n=7. It does not cover all of them, despite every such edge having some cycle of length at most four.

### A score-only certificate that prevents both failures

The following strengthening makes the blackness and nonfolding conditions automatic.

**Score-separated square theorem.** Let P,Q be arc-disjoint Hamiltonian paths in T with four distinct endpoints. Write

    a = d+(start P) - d+(end P) - 1,
    b = d+(start Q) - d+(end Q) - 1.

Suppose a,b are nonzero, a is neither b nor -b, and a and b both have parity n-1. Then the four reversal states form an embedded black Q quadrilateral.

Proof. Reversing a directed Hamiltonian path changes only endpoint scores: its start loses one outneighbor and its end gains one. Since c3(T)=binom(n,3)-sum_v binom(d+(v),2), the c3 increment is a for P and b for Q. The endpoint sets are disjoint, so one reversal leaves the other increment unchanged. The four triangle counts are therefore c,c+a,c+a+b,c+b, and the numerical assumptions make them pairwise distinct. Triangle count is invariant under isomorphism and converse, so the four merged nodes are distinct.

If a Hamiltonian path were blue, its positional reversal would be anti-automorphic. Its endpoint scores would sum to n-1, so its triangle-count increment would equal 2d+(start)-n, of parity n. The stipulated opposite parity rules out blue for every move. Thus the four moves give a genuine black quadrilateral with distinct merged vertices and edges. QED.

### Explicit black quadrilaterals at every order n>=8

The theorem has concrete seeds in every residue class modulo four. In each row P=(0,1,...,n-1), and t is the explorer mask:

| n | t | Q | (a,b) |
|---:|---:|---|---|
| 8 | 644772 | (2,4,0,6,3,7,5,1) | (3,1) |
| 9 | 84938829 | (1,4,2,7,5,8,0,3,6) | (2,4) |
| 10 | 29152537641 | (4,0,7,9,3,1,8,6,2,5) | (1,3) |
| 11 | 28298342874890 | (5,0,3,6,4,10,8,2,7,9,1) | (2,4) |

Each seed contains an internal P arc u->v and an internal Q arc w->z whose four endpoints are distinct and avoid all four global path endpoints. Insert new vertices x0,x1,x2,x3 with the two internal path orders

    P block: x0,x1,x2,x3,
    Q block: x2,x0,x3,x1.

Their six required arcs are a copy of the unique strong order-four tournament and are pairwise compatible. Replace u->v by u->x0->x1->x2->x3->v and w->z by w->x2->x0->x3->x1->z. Let both global starts dominate every new vertex, and let every new vertex dominate both global ends. The chosen internal vertices avoid global endpoints, so these requirements do not conflict with the four insertion arcs. All other cross arcs can be oriented arbitrarily.

The new paths remain arc-disjoint, keep their four global endpoints, and have increments a+4,b+4. Positivity, inequality, and the parity requirement are preserved. For the next insertion use P arc x1->x2 and Q arc x0->x3, whose four internal vertices are distinct. This proves the construction can be repeated indefinitely. The four seeds thus yield an explicit embedded black Q quadrilateral at **every order n>=8**.

The constructor verifies literal tournament arcs, both Hamiltonian paths, disjointness, all four blackness tests, and all four triangle counts at every order8 through100 and eight larger orders through203. Those checks validate the implementation; the insertion argument is the all-order proof. This supplies an infinite constructive nonbridge family, not a claim that every nonleaf edge fits the criterion.

## 6. Cited all-order combined connectivity, and a new consequence for Q

Kenneth B. Reid's theorem says that, up to isomorphism, every n-vertex tournament is related to the transitive tournament by a sequence of reversals of directed Hamiltonian paths. His definition of k-reversal explicitly allows an isomorphic copy at the endpoint; this is not a labelled assertion.

Primary source: Reid's own announcement, abstract 691-05-29, *Notices of the American Mathematical Society*, January 1972, page A-38:
https://www.ams.org/journals/notices/197201/197201FullIssue.pdf

The announced Theorem 1 is the needed Hamiltonian-path case. The journal publication is K. B. Reid, *Equivalence of n-tournaments via k-path reversals*, *Discrete Mathematics* 6(3) (1973), 263-280. The title, volume, and pages were cross-checked against the bibliographic record. The exact primary-source announcement was available in the indexed AMS text; direct retrieval of the PDF returned HTTP403. No claim is made to have independently checked the full published proof.

Combining Reid with the identification in section5 yields:

**Cited theorem consequence.** The combined blue-and-black converse-merged metagraph is connected at every order n>=3.

This is not a novel tournament connectivity theorem; it is a translation of Reid into the current metagraph conventions. Neither removal of blue edges nor division of black multiplicities can be silently inserted into Reid's conclusion.

A useful genuine consequence is nevertheless available:

**Mixed-anchor theorem.** Every active Q component contains at least one mixed self-converse node. Equivalently, every NS node has a black-only Q path to a mixed SC node.

Proof. If a Q component consisted entirely of NS nodes, no blue line could be incident to it, since every blue tiling represents a self-converse tournament. It would therefore also be a component of the combined metagraph. But the transitive tournament is a pure-blue Q isolate outside this active component, contradicting combined connectivity. An active Q component cannot contain a pure-blue vertex because such a vertex has Q degree zero. Thus its required self-converse node is mixed. QED.

More generally, contracting every Q component leaves a connected multigraph whose edges are the blue lines. This is a necessary wiring constraint. It does not say that the active Q component is unique.

## 7. The pendant part is now exact; only nonleaf bridges remain unclassified

A separately proved and independently audited session-fifteen result identifies the only strong self-converse deficit-two center as S4. Thus D(z)=z^4, not merely a finite unknown polynomial. See the companion packet's [complete center proof](../deficit_two_centers_20261020/session15-deficit-centers/DEFICIT_TWO_COMPLETE_CENSUS.md) and [exact leaf series](../deficit_two_centers_20261020/session15-deficit-centers/EXACT_LEAF_SERIES.md). This classification is imported here, not reproved by the n9 census.

Together with the inherited leaf reduction, the complete leaf series is

    L(z) = (2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)].

Every leaf at order n>=5 is nonstrong: an NS leaf is a nonpalindromic rigid SCC word, while a mixed leaf has a nonempty rigid wrapper around S4. Every complement partner of a nonstrong tournament of order at least three is strong. Indeed, a Hamiltonian path traverses its ordered strong components, so its first vertex beats its last. Complementing the nonpath endpoint arc closes the preserved path to a directed Hamiltonian cycle. Consequently no two leaves share an edge at n>=5. At n4 the NS leaf and S4 are the two ends of a single bridge.

It follows that the **complete pendant-bridge series** is L(z)-z^4. If E(z) counts bridges incident to no leaf, then the exact total bridge series is

    Bridge_Q(z) = (2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)] - z^4 + E(z),
    E(z) = O(z^10),

with nonnegative integer coefficients in E. The order-nine census proves the stated vanishing range. It does not permit setting E=0.

The present packet therefore establishes:

- the entire pendant-bridge contribution, using the separately audited center classification;
- no nonleaf Q bridge through order9;
- explicit short-cycle coverage of every nonleaf Q edge through order9;
- an all-order local sufficient criterion, plus black quadrilateral constructions at every order n>=8;
- an all-order black-only path from every NS node to some mixed SC node, using cited combined connectivity.

The all-order active Q connectivity statement, the universal short-cycle conjecture, and the determination of E(z) remain OPEN. This is now a genuine nonleaf-bridge gap, not an unresolved leaf-center alphabet.
