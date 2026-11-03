# The exact decorated square-sum biclique bridge

Status: PROVED elementary equivalence; FINITE-EXACT witness and negative-control checks. This does not construct an integer magic square of nine distinct squares and does not solve that open existence problem. No literature-priority claim is made.

## 1. Arithmetic setting and decoration

A magic square of squares (MSOS) here has nine distinct positive integer square entries, and its three rows, three columns and two diagonals have one common sum. Its entries are square values, not roots.

Over any field of characteristic different from 2, every 3-by-3 magic array has the form

    C+x       C-x-y     C+y
    C-x+y     C         C+x-y
    C-y       C+x+y     C-x.

Indeed, if its common total is T, the sum of all entries is 3T; adding the four lines through the center gives 4T=3T+3C, so T=3C, without dividing by 3. Opposite pairs then sum to 2C. These facts and two corner choices recover the displayed array. Conversely its eight lines all sum to 3C. This works in characteristic 3 too.

For integer arrays, signed swaps of x and y preserve the entry set and permit x>y>0. Nine distinct entries are then equivalent to x≠2y. Positivity is C>x+y. For a field, replace these inequalities by the explicit condition that the nine offsets 0, ±x, ±y, ±(x+y), ±(x−y) are distinct.

Let S be the indexed parallelogram (s0,s1,s2,s3)=(0,x,y,x+y). Given gauge g, put

    Li=C−si−g,     Ri=si+g.

The decorated bipartite graph has the sixteen edges LiRj, each labeled by its actual sum

    Eij=Li+Rj=C+sj−si.

The decoration means both the reflection identities Li+Ri=C for all i and the parallelogram identity R0+R3=R1+R2. Merely drawing a K4,4 or assigning arbitrary square labels is not this condition. Vertex values need not themselves be squares. Extra square-sum edges inside either side are allowed: this is a subgraph, not necessarily an induced K4,4.

## 2. Exact equivalence and reconstruction

The sixteen actual cross-sums have exactly the nine magic entries as their distinct values. Their multiplicities are: C four times; C±x and C±y twice each; C±(x+y) and C±(x−y) once each.

Conversely, suppose two indexed four-tuples of integer or field values satisfy the reflection and parallelogram identities, every actual cross-sum is a square, and there are exactly nine distinct cross-sum values. Set g=R0, x=R1−R0, y=R2−R0. The identities give Ri=g+si and Li=C−g−si, so the same formula holds. Reconstruct the magic array directly from edges as

    E01 E30 E02
    E12 E00 E21
    E20 E03 E10.

Every selected value is a square, the nine are distinct, and all eight line sums are 3C. In the integer setting, positive vertices force positive cross-sums. If integer vertices are not assumed positive, require positive cross-sums explicitly. This is the precise necessary decoration; a bare biclique is not an equivalent arithmetic problem.

## 3. Positive, disjoint integer embedding after one universal factor

Given any integer MSOS in the normal form C>x+y>0 with x>y>0, scale its entries by 4, equivalently multiply every square root by 2. Apply the construction with parameters (4C,4x,4y,g=1):

    Ri=4si+1,     Li=4(C−si)−1.

All Ri are positive. The least Li is 4(C−x−y)−1≥3, so all Li are positive. Each side has four distinct values. Every Ri is 1 modulo 4 and every Li is 3 modulo 4, so the sides are disjoint. All sixteen cross-sums are four times original MSOS entries, hence genuine integer squares.

Therefore an integer MSOS exists if and only if a positively embedded, reflected-parallelogram-decorated square-sum K4,4 with nine distinct labels exists. The forward construction changes the MSOS by square scale 4; it is not a pointwise bijection with unscaled arrays. Its image recovers the original array by dividing all nine labels by 4. Any integer square divisible by 4 has even root, so this division preserves integer-square status. Arbitrary decorated integer bicliques recover an MSOS without division. These are existence equivalences, not a method that guarantees a starting solution.

## 4. Exact gauge count in every odd finite field

Let F have odd cardinality q, with nine distinct offsets. Inside each side the four vertex values remain distinct for every g. A cross-side collision occurs exactly when

    Li=Rj  iff  2g=C−si−sj.

The set S+S has exactly nine elements. One especially short proof is the set identity

    S+S=(x+y)+(S−S),

because (x+y)−S=S. The nine-distinct-offset hypothesis says |S−S|=9. Since multiplication by 2 is a bijection, the collision equations forbid exactly nine gauges. Exactly q−9 gauges yield eight distinct vertices and hence an ordinary simple K4,4 subgraph of the square-sum graph on F.

This permits C=0, permits zero square edge labels, and is valid in characteristic 3. It does not apply in characteristic 2. If the intended graph excludes zero sums, the square array must additionally exclude zero; otherwise the zero-labeled edges are absent and no K4,4 claim follows in that convention. Coordinates in extension fields must use field arithmetic, not integer reduction modulo q.

## 5. Actual paths, and the information they must retain

A disjoint decorated biclique is a genuine square-sum graph carrier, not just an incidence metaphor. It has 4!·4!=576 unoriented Hamiltonian paths: each path alternates sides, and its unique orientation beginning on the left is given by independent permutations of the four left and four right vertices.

Three such paths cover all sixteen edges. Number left vertices 0,1,2,3 and right vertices 4,5,6,7. An explicit covering family is

    0,4,1,5,2,6,3,7
    0,6,1,7,2,4,3,5
    1,5,0,7,2,4,3,6.

Thus every labeled edge can be retained in a union of real square-sum paths. This is different from treating the nine square entries themselves as graph vertices and asking for a Hamiltonian path. A single seven-edge path cannot directly display nine different square labels. Any compressed path certificate must retain the global arithmetic decoration and justify that all nine labels are represented or otherwise proved square. The independent path-count lane supplies sharper two-path label covers; that is not required for this equivalence.

## 6. Exact checks and counterexamples

The verifier checks symbolic coefficient triples for every edge label, its multiplicity, the reconstruction, and all eight line identities. It exhaustively tests every gauge for pinned valid square arrays over F29, F59, F81 and F2187, obtaining respectively 20, 50, 72 and 2178 good gauges. F29 tests center zero; F81 and F2187 test characteristic 3. The extension-field polynomial and full witness are recorded in certificate.json.

An integer near-miss uses C=16,x=12,y=1. Its scaled positive disjoint decorated vertices are L=(63,15,59,11), R=(1,49,5,53). All reflection/parallelogram identities hold and there are nine distinct positive labels, but seven labels are nonsquares. Thus linear geometry and modular square evidence do not replace actual integer-square tests.

A bare square-sum K4,4 over F59 has L=(0,1,2,6), R=(3,15,19,20). Its eight vertices are distinct and all sixteen actual cross-sums are squares, but it has thirteen distinct labels, so no relabeling gives the required decoration. This controls the invalid inference that any existing square-sum biclique solves MSOS.

## 7. Inheritance and claim boundary

The repository already contains reflection matchings x↦t−x, their composition into translations, and square-sum path results in THM-4510, building on THM-4505. Wave-one Session57 material already distinguishes the wrong graph on square entries from incidence-preserving carriers. These are explicitly inherited background, not discoveries of this packet. Targeted repository searches for K4,4/biclique and arithmetic-progression/square-sum interfaces found unrelated K4,4 occurrences (Petersen/Heawood, cube and multipartite constructions), but not this exact reflected parallelogram equivalence. That is a bounded repository search, not evidence of global novelty. Source hashes and the repository commit are pinned in SOURCE_PINS.json.

Run `python verify.py --check` or `python -O verify.py --check`. All new verification conditions use explicit exceptions rather than optimization-removable assertions. No remote writes or modifications outside this packet were performed.
