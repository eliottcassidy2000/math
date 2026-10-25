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

## 7. Active frontier: affine frames and prime-plane counts

Further work studies how many disjoint additive frames realize one prescribed finite-field square plane, and whether this extends from characteristic three to prime-sized planes. Those producer packets and their evolving audits are deliberately reserved for a later checkpoint. This checkpoint makes no additional affine-frame enumeration or general-prime theorem claim. The distinction between AP parts, arbitrary parts, path incidences, and distinct ambient paths must be proved before identifying their counts.

## 8. Productive follow-on questions, not completed results

1. Can a different published bare family meet the full decoration and distinctness locus? The Lagrange exclusions explain why repeating its original-row search cannot do so; they do not settle other families.
2. Can closure containers and label-preserving divisor filters improve the K4,4 minimum interval? Inspect four-subsets of closures rather than requiring the whole closure to have the decoration. Do not report the current upper endpoint as minimal.
3. Can the fixed-curve doubled-x AP constraint be attacked directly, deduplicating centers and retaining a common scaling? More isolated triangles or group progressions are insufficient.
4. What rational-point information is available on the nondegenerate decorated surface beyond ambient dimension and bare points? A resolution claim about X44 needs actual discrepancy analysis, not its dualizing twist alone.
5. Can the finite-field frame count be refined beyond AP parts, or compared with actual ambient path counts? Specify what is counted and avoid overcounting representations.

These are follow-on questions arising from the verified interfaces, not newly proved hypotheses. The session derives exact maps, counts, family-specific obstructions, repairs, and bounded controls. None resolves the nondegenerate integer magic-square-of-squares existence problem.
