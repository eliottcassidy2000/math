# Ordered root amalgamation and its positive binary inverse

Session197. This packet gives an exact physical product, its retained-root transfer, and a target-bounded inverse around one fixed rooted core. It makes no all-order Hamilton spectrum claim.

## 1. Product and inherited conventions

For rooted tournaments (A,a),(C,c), identify a with c, retain that root, and orient every A-a vertex toward every C-c vertex. Write A*C for this rooted tournament. This is an actual graph construction. The rooted singleton is its identity, and the product is associative: the nonroot blocks occur in the same ordered sequence under either parenthesization. Converse reverses the product order.

Use h=H(A), e=end_a(A), s=start_a(A), R and L for the labelled two-bin responses with the root at the start or end, respectively, of a designated nonempty path, and d=H(A-a). The other bin may be empty. All paths partition the vertex set; H(empty)=1. For C use h_C,e_C,s_C,R_C,L_C,H=H(C-c), and beta=B_(C-c)(2).

The O/Z singleton moves and rooted path-cover conventions are inherited from S195:

https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/combined-growth/PROOF.md

The one-edge stationary sewing specialization is the reviewed RETAINED_ROOT_SEWING.md, with its source attribution in RETAINED_SEWING_SOURCE_SCOPE.md. The scalar function F below is f(q)=g(q+1) from the published modular gadget theorem:

https://github.com/eliottcassidy2000/math/blob/48e10c6a4773de23f9ca19ef120b17012a0d7812/session197/modular-gadget-universality/PROOF.md

The proof below is a direct path-fragment derivation. No literature-priority claim is made for rooted amalgamation, four coordinates, or the inherited binary recurrence. The general retained-state product and its arbitrary-core inverse are distinguished from the one-edge seam and from C3-specific source-bundling applications.

## 2. Exact path-fragment laws

The Hamilton count is

    H(A*C)=(R-d)(L_C-H)+(h-e)H+(h_C-s_C)d+e s_C.  (1)

To prove this, put P=A-a and B=C-c, with P->B. A Hamilton word that has no B-to-P transition through the root has its nonroot P block before its nonroot B block. If the root lies before the end of P, its contribution is (h-e)H. If it lies after the beginning of B, its contribution is d(h_C-s_C). If it lies exactly between the two blocks, the contribution is e s_C. These cases remain valid when one nonroot block is empty.

Otherwise the word has the unique pattern

    P_0, B_0, root, P_1, B_1,

where B_0 and P_1 are nonempty, and the outer paths may be empty. The two A paths consist of P_0 and root,P_1; the two C paths consist of B_0,root and B_1. Their counts are R-d and L_C-H. All joining arcs away from the root point P->B. Cutting and concatenating are inverse maps, with no component-order or rotation factor. This proves (1).

The designated-root two-bin response satisfies

    R(A*C)=beta(R-d)+d R_C.                       (2)

Indeed each bin has its P vertices before its B vertices. If the designated root-starting path contains another P vertex, the A part has R-d choices and the two labelled B paths have beta choices. If it contains no P vertex, the entire P set is the other path, giving d choices, and the C cover gives R_C choices. This also covers empty ordinary bins. The dual L formula follows by converse.

If C is not the singleton, a Hamilton word ending at the root has all P vertices, then all B vertices, then the root, so

    e(A*C)=d e_C,  d(A*C)=d H.                    (3)

The non-singleton hypothesis is essential for the endpoint formula: with C the identity the answer is e, not generally d. Equations (1) and (2) themselves survive the identity boundary. The dual start-endpoint formula needs A non-singleton instead.

## 3. A positive retained-state representation

Set x=R-d and y=d-e. For every non-singleton C, equations (1)-(3) become

    h'=Hh+(L_C-H)x+(H-s_C)y+(h_C-H)d,
    x'=beta x+(R_C-H)d,
    y'=(H-e_C)d,
    d'=Hd.                                       (4)

All coefficients are nonnegative integers. Starting or ending Hamilton paths of C inject into paths of C-c, so H>=s_C,e_C. The root-singleton covers contribute H to each of L_C and R_C. Finally H(C)>=H(C-c): insert the root at the first legal gap of any Hamilton word of C-c; deleting the inserted root recovers that word. Every tournament word has a legal gap, by the usual endpoint-or-first-change argument.

More precisely,

    L_C-H>0 iff the root has an in-neighbor,
    R_C-H>0 iff the root has an out-neighbor.      (5)

For example an in-neighbor u supplies the path u,root and any Hamilton path of the remaining vertices as a non-singleton root-ending cover. If there is no in-neighbor, no such path exists. The dual proof gives the second statement.

Associate the identity matrix with the singleton core. Then M_(C*D)=M_D M_C, since the graph product is associative and the whole retained state is respected. This is a physical positive representation, not an invertible state action. Every nonidentity matrix has rank at most3 because the last two output coordinates are fixed multiples of d. Positivity also does not imply total nonnegativity; an explicit negative minor appears in Section6.

The one-edge seam is recovered by taking C to be a new root directed toward a disjoint B except for b->root. Its data give H=H(B), L_C-H=L_b(B), H-s_C=start_b(B), R_C-H=B_B(2)-R_b(B)-H(B), and h_C-H=L_b(B)-start_b(B). Thus (4) agrees with the reviewed seam transfer.

## 4. Arbitrary binary histories around a core

Let O append a universal sink and Z append a vertex beating only the retained root. They are products with the two possible rooted edges. Their exact laws are

    O(h,e,R,d)=(h,0,2R,d),
    Z(h,e,R,d)=(h+R-e,d,2R-d,d).

Define F(0)=1, F(2q+1)=F(q), and F(2q)=F(q)+q+(q mod2) for q>=1. A prefix from the singleton, reading q from most significant bit to least with O=1,Z=0, has state

    (h,e,R,d)=(F(q),(q+1) mod2,q+1,1).            (6)

For an arbitrary input state, a suffix reading t from least significant bit to most with Z=1,O=0 has count

    h+t(R-d)+d(F(t)-1)+(t mod2)(d-e).              (7)

Both formulas follow directly from the two state laws. One way to identify the coefficient of d in (7) is to apply the suffix to the singleton: its chronological bit word is the reverse complement of t, and actual converse of the one-apex graph preserves F. Leading prefix Z steps are full-state neutral only at the singleton; terminal suffix O steps are scalar-neutral only because the construction ends there. No other internal padding is removed.

For a fixed non-singleton core C, form the prefix, multiply by C at the root, and append the suffix. Put

    lambda=L_C-H, alpha=R_C-H,
    sigma=H-s_C, eta=H-e_C.

Equations (4),(6),(7) yield the symmetric positive formula

    N_C(q,t)=h_C+H(F(q)+F(t)-2)
             +lambda q+sigma(q mod2)+beta qt
             +alpha t+eta(t mod2).               (8)

Every q,t>=0 is physically realized. The canonical graph order is

    |C|+bit_length(q)+bit_length(t),              (9)

where bit_length(0)=0. Every binary prefix/core/suffix word has the same count as such a canonical graph.

There is an exact physical converse symmetry, not merely algebraic cancellation:

    N_C(q,t)=N_(C^op)(t,q).                       (10)

Reverse the outer chain orders and the root incidences. The old suffix becomes the new prefix and the old prefix becomes the new suffix with precisely these bit conventions. Algebraically, converse interchanges alpha with lambda and eta with sigma, and fixes h_C,H,beta.

## 5. Complete finite inverse, including zero slopes

Every term after h_C in (8) is nonnegative. Fix M>=h_C and write Z=M-h_C. For the interior q,t>=1,

    q <= floor(Z/(lambda+beta)),
    h_q=h_C+H(F(q)-1)+lambda q+sigma(q mod2),
    1<=t<=floor((M-h_q)/(alpha+beta q)).           (11)

Test (8) and retain only values <=M, or equality for an exact target. The denominator is strictly positive because beta>0 and q>=1. The bound follows from t>=1, and the inner bound discards only nonnegative remaining terms. Conversely each retained pair is an actual graph.

The axes require their own treatment. On q=0, if alpha>0, test 0<=t<=floor(Z/alpha). If alpha=0, (5) says the root is a universal sink; then h_C=H,e_C=H and the entire axis is exactly H F(t). Every value <=M has a canonical t=0 or positive even t<=2 floor(M/H)-2: delete terminal1 digits in the scalar recurrence until t is even, then F(2v)>=v+1. This chooses another actual graph of the same value; it is not a claim that arbitrary internal vertices can be deleted.

Dually, on t=0, if lambda>0 test 0<=q<=floor(Z/lambda); otherwise the root is a universal source and the axis is H F(q), with the same even-index normalization. Include the origin once. Thus every fixed core, including source and sink roots, has a complete finite target inverse. A finite bank is handled by union.

For an interior solution, beta qt<=M and beta>=1. Therefore bit_length(q)+bit_length(t)<=floor(log2 M)+2. A positive-slope axis has index <=M, and a zero-slope normalized axis has index <2M. Every retained count has a witness of order at most

    |C|+floor(log2 M)+2.                          (12)

This is a bounded construction certificate, not uniqueness of q,t or of the core. Its finite sums do not prove that a given bank realizes every odd target.

## 6. A physical repair and a positivity guardrail

Take

    C rows=[24,1,3,6,14], root=2.

Direct counting gives

    (h_C,e_C,s_C,R_C,L_C,H,beta)=(11,3,2,17,20,5,32).

Its retained-state matrix is

    [5 15 3 6; 0 32 0 12; 0 0 0 2; 0 0 0 5].

The minor in rows1,2 and columns2,3 is -96. Thus even this small physical positive core action is not totally nonnegative, and the rank defect rules out a global unimodular four-state inverse.

Nevertheless the scalar target inverse is constructive. At (q,t)=(24084,0), F(q)=20279, so

    N_C(q,0)=11+5(20279-1)+15*24084=462661.

The literal graph has20 vertices. Its core deletion has just four vertices and H=5. This uses arbitrary root incidences rather than the single reversed root-edge restriction of ONE_CORE_BINARY_INVERSE.md. That separate packet proves462661 is absent from its order-at-most5 right-block one-edge language, while also exhibiting a two-block repair. The present graph is a different repair using one general rooted core, not a contradiction to that narrower census.

The verifier checks the full law on8236 actual rooted pairs, literal associativity and converse controls, all singleton boundary formulas, small complete binary inverse sets including zero slopes, and the20-vertex target by subset Hamilton-path DP. It imports only the frozen retained-root verifier/helper; no exploratory bank is an input. Normal and -O receipts must agree.
