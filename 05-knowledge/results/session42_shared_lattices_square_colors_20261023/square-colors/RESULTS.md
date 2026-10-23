# Square-sum island ancestry, color ensembles and incidence-rank decoding

Session42, 23 October2026. Proofs below concern the stated birth DAG or a supplied finite family of paths. They do not assert that a finite coloring scan solves unrestricted Hamiltonian wiring. No priority claim is made.

## 1. The three islands are exact ancestry objects, but they overlap downstream

Let Q_N have vertices1,...,N and edges a--b exactly when a+b is a square. Direct every edge from its larger endpoint n to its smaller endpoint. The predecessor set at birth is

P(n)={s²−n: sqrt(n)<s<=sqrt(2n−1)}.

Its size is floor(sqrt(2n−1))−floor(sqrt(n)). The only vertices without predecessors are1,2,4. For n>=5, one can either check5..8 and then use the square immediately above n, whose difference is less than n; hence no later founder arises.

Let A(n) be the set of founders reachable by decreasing paths. Then

A(1)={1}, A(2)={2}, A(4)={4}, and A(n)=union_(p in P(n)) A(p).

The user's three islands are genuinely monochromatic under this set-valued ancestry:

- {12,4,5,11} has sole ancestor4;
- {2,18,7,9,16} has sole ancestor2;
- {10,6,3,1,8,17} has sole ancestor1.

However13 has parents3 and12, hence ancestry{1,4};14 has ancestry{2,4};22 is the first vertex with all three ancestries. Vertex15, which is absent from these particular island lists, has pure ancestry1. Thus the lists are meaningful selected branches, not an exhaustive three-way decomposition of the later graph.

### Exact all-size saturation theorem

**Every n>=71 has all three ancestors. Vertex70, with ancestry{1,4}, is the last exception.**

Proof: the explicit recurrence certifies every n=71,...,99. For n>=100 put s=floor(sqrt(2n−1)) and p=s²−n. Then p<n and

p > n−2 sqrt(2n−1) >71.

The second inequality follows from (n−71)²>4(2n−1): at100 the difference is45 and its derivative is2n−150>0 thereafter. Thus p is an already certified predecessor, and induction gives A(n)={1,2,4}. The receipt retains the finite base and last-exception check. This is a finite-base structural proof, not extrapolation of a long scan.

The complete pure-ancestry islands are therefore finite and exactly {4,5,11,12}, {2,7,9,16,18}, and {1,3,6,8,10,15,17,19,30}. The user’s first two listed islands are already complete; the third gains exactly15,19,30. This follows from the exact recurrence through70 and the proved saturation tail, not a claim that these remain connected components of Q_N.

Consequently an intrinsic *single* founder label cannot equal reachability after70. The correct larger universe must retain choices, mixtures, or all ancestral owners.

## 2. Parent-selector forests give the exact larger universe

Choose one p(n) in P(n) independently at each nonfounder. Strict descent guarantees a spanning forest with exactly three rooted trees, at1,2,4. Its founder coloring is a coherent hard allegiance.

A color assignment c(n) in {1,2,4}, with fixed founder colors, occurs from such a forest **iff** every nonfounder has at least one earlier neighbor of its own color. If

m_n(c)=#{p in P(n): c(p)=c(n)},

then exactly product_n m_n(c) parent-selector forests induce that coloring. Therefore

sum_c product_n m_n(c) = product_n |P(n)|.

This multiplicity is important: uniform forests do not induce uniform colorings. At N=25 there are1024 forests but only192 distinct founder colorings, checked exhaustively with all shared choices retained.

For variables x_(n,j), the exact joint coloring partition function is

Z_N(x)=sum_(legal c) (product_n m_n(c)) (product_n x_(n,c(n))).

This is a shared-incidence sum, not a product of independently chosen vertex marginals. It can be evaluated with a frontier retaining precisely the earlier vertices queried by future births. Forgetting a shared parent prematurely changes Z_N.

Under uniform independent parent choices, the marginal allegiance vector h(n) is uniquely determined by h(4)=(1,0,0), h(2)=(0,1,0), h(1)=(0,0,1) and

h(n)=(1/|P(n)|) sum_(p in P(n)) h(p).

It is a rational probability vector, and its positive support is exactly A(n). Examples in order(founder4,founder2,founder1):

h(13)=(1/2,0,1/2), h(14)=(1/2,1/2,0), h(22)=(1/4,1/4,1/2), h(23)=(1/4,1/2,1/4), h(25)=(3/4,0,1/4).

**Shared-owner hostile:** color(13)=1 and color(23)=4 are individually possible but jointly impossible, because23 has parents2 and13. The product of their positive marginal probabilities would incorrectly assign positive mass. This is the same precise information-loss mechanism as independently compatible tournament ports, without identifying their underlying graphs.

This birth-directed averaging is a random-ancestry/harmonic recursion. It is not the undirected electrical Dirichlet problem unless a separate conductance construction proves that equivalence.

## 3. Proper coloring, allegiance and path observation are three different tasks

The first non-bipartite graph is Q15, with odd cycle6--3--1--15--10--6 (sums9,4,16,25,16). Bipartite BFS certifies all smaller graphs. The first triangle is6--19--30--6, with sums25,49,36, in Q30.

Thus two-color *proper* graph coloring already fails at15, although two-color path observations can be very informative. Island allegiance is not proper coloring, and it is not forced/optional/unused Hamiltonian edge status. Session40 already supplied explicit counterexamples to both identifications.

The inherited interval-Hamiltonicity route and cycle/chord insertion establish eventual optionality of every fixed square edge. That result is not reopened here. Exact prescribed endpoints and simultaneous edge requirements remain different questions.

## 4. A systematic finite possibility space, rather than cherry-picked examples

The executable study fixes the following feature grammars before scoring. Vertex features:

(n, binary digit sum(n), floor sqrt(n), Zeckendorf digit count(n), smallest Zeckendorf index(n), sum of Zeckendorf indices(n)).

Edge features for a<b, a+b=k²:

(a,b,k, binary digit sums(a)+binary digit sums(b), floor sqrt(a)+floor sqrt(b), Zeckendorf digit counts(a)+counts(b), smallest indices(a)+indices(b)).

For q=2,3, every nonzero linear combination of the indicated features modulo q is tested. Global renaming of colors is deduplicated. These are observation colors, not claimed proper colorings. The ordinary Zeckendorf convention uses indices>=2 with F2=1,F3=2.

The training carrier is **all392 inherited Hamiltonian paths of Q32 up to reversal**, oriented from their smaller endpoint. A label retains the actual endpoint pair and the ordered vertex- or edge-color word. Every imported path is independently checked for vertex coverage, adjacency, orientation and duplication; completeness is inherited from the pinned Session40 census, not independently reenumerated here.

Results:

| family | raw formulas | distinct colorings on Q80 vertices/edges up to names | decoders on Q32 |
|---|---:|---:|---:|
| vertex2 |63|63|56|
| vertex3 |728|364|359|
| edge2 |127|64|60|
| edge3 |2186|1093|1074|

Total:3104 raw formulas,1584 distinct candidates,1549 finite decoders. All1549 also separate the two supplied old birth-color collision paths at36 and the two supplied monochromatic birth-color hostile paths at80. Those are **pair-only holdouts**, not complete Q36/Q80 tests and not an all-size theorem.

The success rate is the important negative finding: many simple colorings pass these modest carriers, so passing them is weak evidence of a uniquely natural arithmetic coloring. The receipt gives the full score histogram, exact features and representative successful formulas. More demanding held-out families must be chosen before promoting any candidate.

## 5. A rigorous probabilistic explanation: shared-incidence rank

Let P be any finite family of distinct vertex paths with orientation fixed, and compare only pairs with the same retained endpoint label. Color either the graph's vertices or its edges independently and uniformly from q colors.

For a pair P,Q, form the **equality-constraint graph** on the colored objects: join the object at position i in P to the object at position i in Q. Loops impose no constraint. Put r(P,Q)=number of its vertices minus number of its connected components, ignoring unused objects. Then

Pr[ordered color words agree] = q^(−r(P,Q)).

Proof: within each equality component all colors must agree. Out of q^v assignments, exactly q^c survive. Repeated objects and shared incidences are automatically counted correctly; counting disagreeing positions independently would be false.

Therefore if sum_(same-endpoint pairs) q^(−r(P,Q)) <1, there exists a q-coloring separating every path in the supplied family. This is an all-size sufficient criterion for a **supplied path family**, not a method that avoids constructing it.

For Q32 the exact union bounds are:

- vertex2:74070791/268435456;
- edge2:14553243/134217728;
- vertex3:1665754372450/68630377364883;
- edge3:2943585822692/617673396283947.

All are below1. Random two-color observations therefore already have a substantial guaranteed success probability on this finite carrier. No special Fibonacci interpretation is needed to explain the high success rate above.

An exact counting refinement follows by inclusion-exclusion. For each collision event e, retain its equality constraints. For a subset S of collision events, let r(S) be the rank of the union of all their constraints. If M objects are colored, the number of separating q-colorings is

sum_(S subset events) (−1)^|S| q^(M−r(S)).

This is a partition-lattice/graphic-rank invariant of realized path pairs. It retains shared equalities; replacing r(S) by a sum of individual ranks generally overcounts independent restrictions. The finite direct enumeration used here evaluates specific feature families, not this exponentially large entire expression.

### Exact control and an all-size warning

The separate `check_collision_partition.py` exhausts all six paths from0 to4 in K5 and all15 pair-collision events. Inclusion-exclusion over32768 event subsets agrees with direct enumeration of every two- and three-color assignment. The vertex count polynomial is q^5−3q^4+2q^3, giving0 two-color decoders and54 three-color decoders. The edge count polynomial is q^10−6q^8+3q^7+10q^6−12q^5+4q^4, giving192 and30942 respectively.

More generally, for all endpoint-fixed Hamiltonian paths in K_n, a vertex-color word separates the paths iff all n−2 internal vertices have distinct colors. The exact number of successful q-colorings is q²(q)_(n−2), where the falling factorial vanishes for q<n−2. If two internal vertices share a color, exchanging their positions makes an invisible distinct path; distinct internal colors conversely recover every position.

An edge-color decoder must at least satisfy q^(n−1)>=(n−2)!, by counting words. Thus no fixed finite palette decodes every complete-graph path family at all sizes, even if the coloring is redesigned at each order. This is a general-graph boundary, not an assertion that square-sum graphs attain the same path counts.

## 6. Berggren connections: retain the proved map, reject the false growth identification

The numbers3,4,5 are the sizes of the three components at12, but after13 the component sizes are3 and10; after14 there is one component of size14. The later listed islands have induced-edge counts3,4,5, while their vertex counts are4,5,6. Changing the counted object is legitimate only when stated, not as a preserved Pythagorean invariant.

The birth DAG has three roots, merges, and vertices with arbitrarily many predecessors. A Berggren tree has one root and a unique parent per nonroot, with exactly three forward branches. Thus there is no bijection preserving all birth vertices and all directed predecessor incidences between these carriers. Choosing one parent produces a forest and erases the very merge data represented by the mixed allegiance vectors. Relabeling its three founders as Berggren branches is not an incidence-preserving map.

There *is* an inherited exact Fibonacci--Berggren map. THM-3339, with the later chamber correction retained, uses consecutive Fibonacci Euclid parameters and mandatory odd/odd primitive normalization to place triples on three interleaved rays (BA)^r, A(BA)^r, C(BC)^r. THM-3382 retains the ray object but distinguishes Fibonacci time from ternary heap address h(wA)=3h(w)+1, h(wB)=3h(w)+2, h(wC)=3h(w)+3; these change density and harmonic weighting. Neither theorem identifies those rays with the square-sum islands. They provide a genuine candidate target for a future map, which must retain predecessor merges or explicitly name the quotient it discards.

A primary literature reference for the Berggren-tree carrier is [Janičková and Csókási, Metric properties in Berggren tree of primitive Pythagorean triples](https://arxiv.org/abs/2304.05230). The established three-matrix generation theorem is used only for its stated tree structure here; no new theorem about all triples is claimed.

## 7. Reproducibility and provenance

Run `python3 explore.py` and `python3 -O explore.py`, and likewise `check_collision_partition.py`. Paired receipts are identical. `inherited_paths.json` contains the required path subset plus the SHA256 of its full source receipt at commit5258d52ac95924d524ee816a7dffa7329ddcf095. The two hostile files are extracted from the same pinned Session40 directory. The published manifest excludes the large local full receipt used for extraction.

Important next obligations: extend the finite coloring stress tests to independent complete carriers; study the all-size ancestor partition function with shared frontier states; and specify any island-to-Berggren map by actual preserved incidences rather than three matching labels.
