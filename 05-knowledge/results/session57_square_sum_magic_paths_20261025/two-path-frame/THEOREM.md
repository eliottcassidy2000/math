# The exact decorated two-path frame

Status: PROVED over the rationals, with exact integer matrices and an exhaustive eight-vertex graph check. This is a representation/equivalence theorem conditional on the square labels; it neither constructs an integer magic square of distinct squares nor asserts global novelty.

## Frozen graph and indexing

Keep the indexed, oriented paths

    P = L0,R3,L1,R1,L2,R0,L3,R2
    Q = L0,R3,L2,R1,L1,R2,L3,R0.

Their union H has nine edges, in this order:

    e0=03, e1=11, e2=12, e3=13, e4=20,
    e5=21, e6=23, e7=30, e8=32,

where ij means LiRj. Their five shared edges form the three components {L0,R3}, {L1,R1,L2}, {R0,L3,R2}, of sizes 2,3,3. The symmetric difference has FOUR edges, not five: 13,20,12,23. Removing the leaf L0 and its edge to R3 leaves a theta graph with branch vertices L1,L2 and three internally disjoint routes

    L1,R1,L2; L1,R3,L2; L1,R2,L3,R0,L2.

Their lengths are 2,2,4. The whole union has exactly four unoriented Hamiltonian paths, recorded in the certificate; P and Q are two of them. The original all-biclique path-cover search, not this fixed-frame calculation, proves the global minimum of nine union edges for a two-path nine-label cover.

## The two automatic equations, and what they do not say

Let v=(L0,L1,L2,L3,R0,R1,R2,R3)^t, and let E be the 9-by-8 unsigned incidence matrix. The edge-value vector is z=Ev, so zj is an actual sum, not a freely assigned edge color. Connectedness and bipartiteness give rank E=7 and

    ker E = span(1,1,1,1,-1,-1,-1,-1).

The only linear identities automatically forced on z are the two alternating cycle identities

    z1-z3-z5+z6=0,
    z1-z2+z4-z5-z7+z8=0.                         (Zz=0)

Thus nine actual edge sums occupy a seven-dimensional space, not the three-dimensional magic-value plane. The common-edge forest alone does not add arithmetic equations.

Here is a constructive integral inverse modulo gauge. Set R0=0 and recover

    L2=z4, R1=z5-z4, L1=z1-z5+z4,
    R3=z3-z1+z5-z4, L0=z0-z3+z1-z5+z4,
    L3=z7, R2=z8-z7.

These reproduce the seven edges in P. They reproduce the remaining edges 12 and 23 exactly when the two cycle identities hold. This defines the integer 8-by-9 matrix T in the certificate. In particular TEv equals v plus R0 times the displayed gauge vector. Consequently im E=ker Z, over Q and also as integer lattices.

## Four independent extra constraints

Let A express the three reflection differences

    Li+Ri-L0-R0=0 (i=1,2,3)

and the parallelogram relation

    R0+R3-R1-R2=0.

The exact ranks are rank A=4 and A annihilates the gauge vector. Set B=AT. Explicitly, the FOUR additional edge equations are

    -z0+z3-z4+z5=0,
    -z0-z1+z3+z5-z7+z8=0,
    -z0-2z1+2z3-2z4+2z5+z7=0,
    -z1+z3+z7-z8=0.                            (Bz=0)

The identity A=BE follows immediately from TE's gauge formula. Therefore, for actual edge sums z=Ev, Bz=0 is equivalent to Av=0. The four equations are independent even modulo the cycle equations: rank [Z;B]=6, while rank Z=2. Their independence also follows because A has rank four and E identifies the potential space modulo its one-dimensional gauge with im E. They cannot be silently inferred from the two paths or their cycles.

For a simple negative control, taking vertex potentials (0,1,2,3,4,5,6,7) gives actual edge values satisfying both cycle identities but not the four extra equations. This control concerns linear implications, not square arithmetic.

## Exact identification with the magic plane

Under Av=0, define C=L0+R0, x=R1-R0, y=R2-R0, g=R0. Then

    R=(g,g+x,g+y,g+x+y),
    L=(C-g,C-x-g,C-y-g,C-x-y-g).

The nine edge labels in the specified order are

    C+x+y, C, C-x+y, C+y, C-y,
    C+x-y, C+x, C-x-y, C-x.

In particular the reconstructed array is

    z6 z7 z3
    z2 z1 z5
    z4 z0 z8,

namely

    C+x   C-x-y C+y
    C-x+y C     C+x-y
    C-y   C+x+y C-x.

Every line sums to 3C. Conversely every rational 3-by-3 magic array has this normal form, so its inverse permutation gives exactly this edge-label subspace. Equivalently: the seven equations obtained by subtracting the first row sum from each other row/column/diagonal sum have rank six, and their row space equals the span of the two rows of Z and four rows of B. The certificate records all seven equations, a greedy independent selection of six, and an exact rank-six check of the combined row space. No floating-point rank tests are used.

The geometry of the linear bridge is therefore

    ker A / (global bipartite gauge)  --E-->  ker Z intersect ker B

with dimensions 4-1=3 on the left and 9-6=3 on the right. This is an isomorphism, with inverse supplied by T, and remains an integral-lattice isomorphism when all values are integral. A fixed source choice C,x,y does not determine eight vertex values until a gauge is chosen.

## Square arithmetic, positivity, and clearing denominators

Imposing that these nine distinct labels are positive rational squares is an additional nonlinear arithmetic condition. Because H contains one representative of every label, it then follows that all sixteen sums in the uniquely completed indexed K4,4 are squares: every omitted edge repeats one of the nine values. This implication uses the decoration and the exact label coverage, not merely the presence of Hamiltonian paths.

For positive integer square labels, normalize by signed swaps of x,y so that x>y>0. Distinctness then requires x != 2y, and positivity is C>x+y. This normalization may permute the original labeled array; preserve the permutation if reconstructing that exact array matters. The frozen full-biclique bridge supplies a universal positive/disjoint realization after square scale 4:

    si=(0,x,y,x+y), Ri=4si+1, Li=4(C-si)-1.

All Ri are positive and 1 modulo 4; all Li are positive and 3 modulo 4. Values within each side are distinct, and the sides are disjoint. Restricting to P and Q gives two actual simple Hamiltonian paths on eight positive distinct integer vertices whose union carries the nine scaled square labels. The resulting magic array is four times the original, not the original pointwise; dividing the recovered labels by four recovers it and preserves integer-square status.

For positive rational square labels, choose d to clear the denominators of all nine rational roots. Scale labels by d^2, then use the same factor-four construction. The total square scale is 4d^2. This proves an existence equivalence, not existence itself. Conversely a positive disjoint integer decorated frame with nine distinct square labels reconstructs an integer magic square directly, with no division needed.

## Source, target, and loss ledger

- Source: indexed vertex potentials with reflection/parallelogram decoration, plus the frozen paths; for square-sum claims also require nine distinct positive square labels and a genuine distinct-vertex embedding.
- Target: a fixed permutation of nine edge values forming a magic array, retaining all square values exactly before optional scaling.
- Forgotten information: the one-dimensional bipartite gauge. Abstract/unlabeled paths also forget the Li/Rj indexing and therefore cannot recover the displayed array without extra data.
- Additional choices: signed coordinate normalization and optional denominator/positivity square scaling. These must be recorded if exact unscaled recovery is wanted.
- Arithmetic not supplied: rational or integer squarehood and nine-value distinctness. Linear rank calculations do not prove either.
- Not a bare Hamiltonian-path bijection: H has four paths, the full biclique has 576, while the mathematical correspondence is between decorated potential classes and magic value arrays.

The source path witness is copied to frozen_minimal_union.json; its SHA-256 is in certificate.json. The reflected full-biclique positivity construction is inherited from the sibling Session57 biclique bridge. Run python verify.py --check and python -O verify.py --check. Every check uses explicit exceptions, so optimized Python cannot disable verification. No remote writes or old-file modifications were made.
