# Session57 synthesis: what the connections preserve

## 1. The precise target and three different readings

The target is nine **distinct positive integer squares** arranged in a 3-by-3 grid whose three rows, three columns, and two diagonals have one common total. Its rational counterpart is equivalent by multiplying every square root by one common denominator. Independent rescaling of cells, lines, or triangles is not allowed.

For every rational magic grid the center C determines the total T=3C. With C=e² nonzero, T cannot be an integer square. That elementary obstruction answers the stronger reading “square entries and a common square total,” not the usual magic-square-of-squares problem. Allowing different square line totals is weaker: the [quaternion-derived witness](square-total-variant/RESULT.md) has root grid

    2655  2976   320
    2624  2520  1665
    1440   895  3624

Its squared entries are distinct and positive. Six row/column totals are 4001², but the diagonals are 5151² and 2920². It is not magic. The search stopped at its first witness and proves neither minimality nor exhaustion of its parameter box.

## 2. Why the original graph fails, and the replacement works

A square-sum graph on the original nine values asks whether pairs sum to squares. Magic imposes eight **triple** sums with one shared total. These predicates are different. The [fixed-total hypergraph](carriers/THEOREM.md) retains the right predicate exactly; eight synchronized additive paths are an equivalent representation. The line-intersection graph alone loses the arithmetic.

The [finite-field controls](character-graphs/THEOREM.md) are concrete failures of an unqualified original-entry Hamiltonian criterion: F17 admits the reported path without the corresponding magic arrangement, while an F59 magic entry set has an isolated vertex. These are scoped finite-field examples, not integer nonexistence theorems.

The positive repair moves the square values from vertices to **edge sums**. In the [decorated K4,4](biclique-bridge/THEOREM.md), let S={0,x,y,x+y} and impose cross labels C+t−s for s,t in S. The sixteen edges carry precisely nine forms C+i x+j y, i,j∈{−1,0,1}, with prescribed multiplicities. Reflection, parallelogram, positivity, and distinctness conditions identify the intended magic grid. Two suitable Hamiltonian paths inspect all nine formal labels, not all sixteen edges and not the original nine vertices. The nine-edge union has two cycle equations; [four further arithmetic constraints](two-path-frame/THEOREM.md) remain necessary to recover the magic relations.

The [AP-K3,3 model](ap-biclique-bridge/REPORT.md) uses six carrier vertices and nine cross edges. Two paths sharing one edge cover them all, but the AP orders and diagonal conditions cannot be inferred from the bare skeleton. Positive integer gauge/scaling yields the stated induced normalization; this must not be silently transferred to arbitrary finite-field carriers.

## 3. Bare integer bicliques exist; the decoration is the bottleneck

The [common-neighbor factor bijection](integer-bicliques/REPORT.md) gives exact bounded enumeration: for a<b, write b−a=pq with p<q of the same parity, then x=((q−p)/2)²−a, retaining only positive in-range x different from a,b. Independent adjacency and factor routes establish minima 274 for K2,4 and 673 for K3,3, and absence of K4,4 through 10000.

A published family attributed to Jean Lagrange through Guy, reproduced by Choudhry, supplies bare K4,4 points. The [session's normalized bounded-family witness](integer-bicliques/PRIMARY_CONSTRUCTION.md) is

    L={8648,44943,242063,725903},
    R={1,62641,242353,1112833}.

Its sixteen cross sums are squares; both sorted opposite-pair defects are nonzero. The exact conclusion is 10001≤N_min(K4,4)≤1112833. The optimized upper witness is best only in the explicitly scanned primitive parameter box after the stated normalization/gauge procedure. Neither this bound nor a bare K4,4 supplies a magic square.

The [family obstruction](lagrange-obstruction/proof.md) is stronger than that finite search in a different direction: homogeneous factorization and exhaustive projective tests modulo 47 rule out every nonzero rational member of this particular family satisfying both parallelogram balances. Primitive-coordinate reduction makes the argument cover arbitrary rational denominators. The surviving A-balance boundary v=0 fails B-balance.

The [AP strengthening](lagrange-obstruction/ap-addendum.md) exhausts twelve indexed AP equations. All residual factors have finite-field projective certificates; the sole rational exceptional line gives three equal A values. Thus no original A triple consists of three distinct AP values. This excludes choosing original rows and columns for the AP-K3,3 bridge, including gauge, nonzero square scaling, permutations, and part exchange. It does **not** exclude selecting and regrouping other nine square values, finding fresh additive coordinates, using another family, or embedding into a larger graph.

## 4. Geometry separates the ambient room from the target locus

The [square-root variety](square-matrix-geometry/geometry.md) X_(r,s) is cut out by q_ij²+q_11²−q_i1²−q_1j²=0. Monic regular-sequence and multiquadratic arguments prove a geometrically integral normal complete intersection in characteristic different from 2. Its dimension is r+s−2, embedding degree 2^((r−1)(s−1)), and dualizing sheaf O(rs−2r−2s+2). Singular points correspond to zero rectangles, except the smooth X22 case.

In particular X44 has dimension 6, degree 512, and dualizing O(2). Its projective squaring-cover degree is instead 32768. Neither degree nor dualizing twist alone classifies a smooth resolution or proves the required rational point exists. A known numerical example independently proves bare-point existence on the relevant open subset.

The specified decoration restricts the seven-dimensional affine additive base to three parameters. Its pullback has 128 relative-sign components, each the six-quadric magic-square surface in P8, dimension 2, degree 64, dualizing O(3). These are not 128 different numerical magic-square constructions. The identification is an explicit coordinate derivation; the magic surface and its geometric study are already in the literature, including Auel–Singer. Its nonzero, nine-distinct-value open condition must still be imposed.

## 5. Repair a source reduction rather than discard its construction

For H(a,b)=a⁴−6a²b²+b⁴, the [quartic audit](norm-quartic-repair/repair.md) distinguishes cross matching, controlled by H(r,s)−m²H(p,q), from diagonal matching, controlled by (r²+s²)²−m²(p²+q²)². Choudhry's printed transition preceding (5.5) identifies the former with a formula for the latter. The PDF check confirms the discrepancy is not merely HTML extraction.

Explicit exceptional roots repair the rational construction on its defined domain. The separate [direct polynomial check](norm-quartic-repair/POLYNOMIAL_CHECK.md) verifies all sixteen identities of the printed polynomial formula over Z[u,v,x,y], including boundary specializations not covered by the rational substitution. Thus this is a correction of an intermediate reduction **and confirmation of the square-valued polynomial conclusion**, not a refutation of Theorem 3. Distinctness at every specialization is neither true nor claimed.

## 6. A classical elliptic bridge with the simultaneous guard retained

Fix rational d>0 and C_d={c>0: c−d,c,c+d are rational squares}. The [fixed-gap transpose](congruent-curve/PROOF.md) says a rational magic square with this gap d exists precisely when a center AP C−e,C,C+e lies in C_d with d,e>0, C>d+e, and e∉{d,2d,d/2}. These guards are the positivity and nine-value distinctness conditions.

Classically, c is the hypotenuse square of a positive rational right triangle of area N=d/4. On E_N: Y²=X³−N²X, the explicit triangle map and its inverse give c=4x(2P), with an eligible positive-branch point X>N, Y≠0. Hence the hard simultaneous condition is an AP among three eligible **doubled** x-coordinates on one fixed curve, with distinctness guards. Arbitrary x-coordinate APs, group progressions P,2P,3P, and torsion translations do not supply it. Torsion translates can repeat the same center.

For each fixed positive integer gap, the integer centers are exactly c=r²+s² square with d=2rs and r>s>0. This finite divisor description gives no guarded AP for d≤20000, and the 24 sampled E6 multiples likewise give none. Neither is a global exclusion over rational points. Common side scaling L sends areas N to L²N, centers and gaps to L² times themselves, and curves by (X,Y)↦(L²X,L³Y). It preserves *equality* of the areas, not their numerical value.

The triangle/curve correspondence, duplication formula, and magic-square connection are classical or recovered from earlier repository work. The session contribution here is a checked, guard-preserving synthesis and bounded certificates, not a priority claim for that connection.

## 7. From the AP subclass to every fixed-plane frame

Let A=C+U be a fixed nine-point affine F3-plane. A frame is an unordered pair {B,D} of three-element sets with the addition map B×D→A bijective. The [complete classification](affine-bare-frame-counts/THEOREM.md) now answers the earlier AP-only counting question.

After translation, differences of both factors lie in U. Unique addition means (B−B)∩(D−D)={0}. A noncollinear triple determines six of U's eight nonzero vectors; its complement is a punctured line. Therefore the other factor must be a full affine line. Both factors cannot be non-AP.

Counting a distinguished AP part gives 36q presentations: four directions, q/3 cosets, and 27 transversals of the three parallel target lines. The AP/AP presentations are counted twice on forgetting the distinguished part; mixed ones once. Thus there are 30q unordered frames, split as 6q AP/AP and 24q mixed. A gauge class has q distinct members and exactly nine with overlapping parts, yielding 30(q−9) disjoint frames. The AP/AP fraction is 1/5; among disjoint frames this fraction requires q>9.

Each K3,3 has 36 unoriented Hamilton paths. A six-vertex path recovers the unordered alternating bipartition and hence the frame, so these path counts are injective across frames. Consequently there are exactly 1080(q−9) distinct paths in this subclass, not merely incidences. Exactly 36 unordered path pairs cover a K3,3, each sharing one edge; the connected union recovers the frame. Hence the same total counts unordered two-path covers. These are not Hamilton paths of the full q-vertex graph or a census of all ambient paths.

The [independent audit](audit-affine-frames/AUDIT.md) supplies exhaustive triple and path checks, a 9! layout census, field irreducibility checks, and explicit two-path enumeration. F81 gives 2160 disjoint frames and 77760 paths; F729 and F2187 give 21600 and 65340 frames. Both the F81 and F729 planes contain zero; the F2187 witness avoids it. A graph interpretation requires the appropriate square convention. The [AP-only corrected edition](audit-affine-frames/THEOREM.md) and [preservation receipt](audit-affine-frames/PUBLICATION_CORRECTIONS.json) explicitly repair the former false suggestion of repeated path realizations. Numerical AP counts are unchanged.

## 8. Why the all-frame classification extends to odd primes

Let p be an odd prime, q=p^n with n≥2, and U a two-dimensional vector space over F_p inside the additive group of F_q. For one fixed A=C+U, take p-point parts B,D whose p² cross sums give A exactly once. This is an additive factorization theorem before imposing any square or graph condition.

The [prime-plane theorem](prime-plane-frames/prime_plane_frames.md) uses a necessary correction to the naive difference argument: for p>3, the presence of v in a difference set does not imply that every nonzero scalar multiple of v is present. Disjoint difference sets alone do not justify disjoint direction sets by that shortcut.

Instead apply a nonzero linear functional ℓ:U→F_p and a primitive pth root ζ. Unique addition makes the product of the two character sums zero. If Σ_(s∈S) ζ^ℓ(s)=0 for a p-element factor S, its projection-count polynomial has degree at most p−1 and vanishes at ζ. Divisibility by Φ_p=1+T+⋯+T^(p−1), together with total count p, forces every coefficient to equal 1. Thus one factor projects bijectively under ℓ. For every direction ker ℓ, at least one factor avoids that direction: their determined direction sets are disjoint.

The classical Rédei–Megyesi bound says a noncollinear p-point subset determines at least (p+3)/2 directions. Two noncollinear factors would therefore require at least p+3 disjoint directions in a plane with only p+1, impossible. One factor must be an affine line; the other is a transversal. The direction bound is an external classical input, with Somlai's modern proof cited in the packet. The Fourier/cyclotomic argument and counting synthesis are verified session derivations, not a claim of priority for the underlying factorization phenomenon or direction theorem.

Put C_p=(p+1)p^(p−1)−p(p+1)/2. Counting line-first presentations and correcting the double count of two-line frames gives C_p q unordered frames: p(p+1)q/2 have two line factors, and the rest exactly one. Each gauge class is free of size q. Its p² distinct cross differences give exactly p² forbidden gauges, since 2 is invertible. Hence there are C_p(q−p²) disjoint frames. Gauge equivalence of unordered pairs is intrinsic; a translation action requires a consistent orientation and must not silently depend on an arbitrary order.

The line/line class fraction is 1/(2p^(p−2)−1), giving 1/5 at p=3 and 1/249 at p=5. This is also the fraction among disjoint frames when q>p²; it is undefined there when q=p² because none are disjoint. In characteristic two the formula for disjointness fails: all 3q frames are disjoint when 0∉A, and none when 0∈A. No extension to nonprime subfield orders is claimed.

A K_(p,p) cross skeleton has (p!)² unoriented Hamilton paths, with alternating bipartition again ensuring injectivity across frames. Thus the exact prescribed-plane path subclass has (p!)²C_p(q−p²) paths. To interpret them inside a square-sum graph, separately require every value of A to meet its square convention. Neither inducedness nor vertex-square status follows from the additive theorem.

An elementary realization removes the square-plane hypothesis for one infinite family of zero-allowed examples: for every odd prime p, take q=p⁴ and U=F_(p²) as its subfield. Every nonzero u∈U satisfies u^((p⁴−1)/2)=(u^(p²−1))^((p²+1)/2)=1, so is a square in F_q. Thus A=U is a square-valued plane, yielding exactly C_p(p⁴−p²) disjoint square-sum K_(p,p) carriers with this cross-sum set. This is a classical subfield fact combined with the frame count, not a new finite-field existence breakthrough. Zero occurs in A and its edge label requires the zero-allowed convention. Carrier vertices need not be squares; nothing here lifts to integers.

A companion construction avoids zero: take q=p⁶, choose a two-dimensional F_p-subspace U of F_(p³), and C∈F_(p³)\U. Then A=C+U excludes zero and lies in the quadratic-extension square subgroup, since each nonzero u∈F_(p³) has u^((p⁶−1)/2)=(u^(p³−1))^((p³+1)/2)=1. The same theorem gives C_p(p⁶−p²) disjoint carriers with only nonzero square sums for every odd p. It still asserts neither inducedness nor square-valued vertices nor integer lifting.

Two Hamilton paths use at most 4p−2 edges, fewer than p² for p≥5. They therefore cannot cover all p² distinct cross-edge labels. The positive nine-label two-path result is genuinely small-size; it is not a general p-by-p magic-array construction.

Finally, affine parameterizations of a p-point line form AGL(1,p), of order p(p−1) and index (p−2)! in S_p. At p=3 this is all S3, explaining why all 72 ordered presentations of an unordered two-line frame are affine. At p=5 just 20 of the 120 single-line orders are affine. These are cyclic field parameterizations, not increasing real AP orders, and the observation claims no further representation-theoretic classification.

## 9. Productive follow-on questions, not completed results

1. Can a different published bare family meet the full decoration and distinctness locus? The Lagrange exclusions explain why repeating its original-row search cannot do so; they do not settle other families.
2. Can closure containers and label-preserving divisor filters improve the K4,4 minimum interval? Inspect four-subsets of closures rather than requiring the whole closure to have the decoration. Do not report the current upper endpoint as minimal.
3. Can the fixed-curve doubled-x AP constraint be attacked directly, deduplicating centers and retaining a common scaling? More isolated triangles or group progressions are insufficient.
4. What rational-point information is available on the nondegenerate decorated surface beyond ambient dimension and bare points? A resolution claim about X44 needs actual discrepancy analysis, not its dualizing twist alone.
5. The prescribed-plane frame classification and its distinct-path count are now proved above. What remains beyond this scope is a census of all ambient paths or frames with other cross-sum sets. Do not mistake the exact subclass count for that broader enumeration.

These are follow-on questions arising from the verified interfaces, not newly proved hypotheses. The session derives exact maps, counts, family-specific obstructions, repairs, and bounded controls. None resolves the nondegenerate integer magic-square-of-squares existence problem.
