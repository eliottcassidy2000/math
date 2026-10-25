# Corrected publication edition of the AP-only theorem

This text corrects the original path-realization sentence, as required by AUDIT.md. The frame counts are unchanged. Both F81 and F729 pinned planes contain zero; F2187 does not. The complete all-frame extension is endorsed separately in bare-extension/THEOREM.md. Original and corrected hashes are recorded in PUBLICATION_CORRECTIONS.json; this is not a claim of byte-identical reproduction.

# Exact AP frame counts for one characteristic-three square plane

## The theorem and counting conventions

Let F=F_(3^n), q=|F|, and let A=C+U be a fixed nine-element affine plane, where U is a two-dimensional F3-subspace of F. Suppose A is contained in the square-value set Q={x²:x∈F}, with zero allowed. A disjoint AP frame means an unordered pair {B,D} of disjoint three-element affine F3-lines in F whose nine cross sums constitute A. This specifies both parts, not merely their six-vertex union, and does not require their square-sum graph to be induced K3,3.

There are exactly

    6(q−9) unordered disjoint AP frames,
    12(q−9) ordered pairs of unordered AP parts,
    432(q−9) ordered pairs of internally ordered AP parts.

The geometric count does not itself need the square-locus hypothesis: that hypothesis makes each frame's nine cross edges edges of the square-sum graph. If zero is forbidden, require A⊆Q\{0}; do not apply that stronger convention to a zero-containing witness. No existence claim for arbitrary q is made.

Each cross-edge skeleton is K3,3 and has exactly 36 Hamilton paths modulo reversal. Consequently there are 216(q−9) frame/path incidences, each a distinct six-vertex simple path in the ambient square-sum graph: the unordered alternating bipartition recovers its unique frame. Additional ambient edges do not change this injectivity. This counts exactly the AP-frame subclass, not all ambient paths and not Hamilton paths of the full q-vertex graph.

## Bijection and absence of hidden gauge multiplicities

U has four one-dimensional F3-subspaces. For each of their six unordered pairs, choose one orientation (L,M) once and for all. Distinct directions are complementary: U=L⊕M. Define

    B_g=g+L,       D_g=C−g+M,       g∈F.

Their cross sums are C+L+M=A, with every sum attained once.

For fixed oriented directions, equality of B_g and B_h forces h−g∈L; equality of D_g and D_h forces h−g∈M. Thus h=g. Equality as unordered part pairs cannot swap the two parts within this fixed orientation, because a three-point affine line has a unique direction, and L≠M. Different unordered direction pairs cannot yield the same frame for that same reason. In particular, despite each individual line having three choices of a base point, there is no additional factor of three or nine in this simultaneous parametrization.

Conversely, let B=b+L and D=d+M be AP parts with B+D=A. Their direction span is exactly U, since differences of cross sums span L+M and differences of A span U. Hence L,M are distinct directions in U. Choose the ordering of the parts matching the fixed orientation. Write uniquely b+d−C=l+m with l∈L,m∈M, and put g=b−l. Then B=g+L and C−g=d−m, so D=C−g+M. This proves surjectivity and the full bijection before imposing disjointness.

A common vertex is equivalent to

    g+l=C−g+m, or 2g=C+m−l∈C+U=A.

Because multiplication by 2 is bijective and A has nine elements, exactly nine values of g are forbidden. Intersections for these bad gauges have exactly one point, since the two directions are distinct. Each oriented direction slice therefore contributes exactly q−9 frames.

The side-swap operation sends (L,M,g) to (M,L,C−g). Its apparent fixed gauge g=C/2 is always bad, but it also swaps the ordered directions. It is not an involution inside the single chosen oriented slice, and one must not divide q−9 by two there. The factor of six has already selected unordered direction pairs.

## Changing the representative C

The affine set A, not a distinguished center, is fixed. Replacing C by C'=C+h, h∈U, does not create further frames. For a given oriented pair write h=l_h+m_h in L⊕M. The exact reindexing is

    g'=g+l_h.

Then g'+L=g+L and C'−g'+M=C−g+M. It is a bijection of all gauges and preserves disjointness. Thus the count is independent of the nine possible representatives C∈A.

## Why the 432 layouts and the 72-fold frame quotient agree

Every fixed nine-point affine plane has 9·8·6=432 labeled magic arrangements: choose its central entry c∈A, a nonzero u∈U, and v∈U\span(u). The array is

    c+u       c−u−v     c+v
    c−u+v     c         c+u−v
    c−v       c+u+v     c−u.

These choices are recovered from the center, top-left and top-right entries, so all 432 are distinct. The standard odd-characteristic parametrization of a 3-by-3 magic array gives necessity. For each such arrangement, its ordered additive edge matrix has exactly q potential solutions, one gauge parameter; exactly nine have coincident cross-part vertices. Hence there are 432(q−9) fully ordered disjoint frames.

Equivalently, every unordered frame has 2·3!·3!=72 fully ordered versions. Unlike over the real numbers, every permutation of a three-point affine F3-line is an AP order: for any ordering (x,y,z), x+y+z=0 implies x+z=2y. Therefore all 72 orderings satisfy the appropriate AP guards in the semimagic-to-magic permutation. Division by 72 gives the same 6(q−9).

Over the reals or integers, a distinct three-point AP has a unique midpoint, so just two of its six orders are AP orders in a fixed midpoint convention. Reversing each AP and swapping the two yields eight valid ordered choices, not 72. The independent real theorem that a distinct magic multiset has eight layouts uses ordering to recover its parameters; that theorem does not transfer to characteristic three. In particular the 432 finite-field layouts are not 432 real layouts or a claim about integer existence.

## AP restriction is essential

Bare square-sum K3,3 certificates need not have AP parts. Choose independent u,v∈U and set

    B0={0,u,−u},
    D0={C,C+v,C−v+u}.

D0 meets each of the three cosets of span(u) in A once, so B0+D0=A with unique representations. Its three elements sum to u≠0, so D0 is not an affine F3-line. Gauge B0→g+B0 and D0→D0−g makes the parts disjoint for q−9 choices, because D0−B0=A. These are legitimate bare bicliques but are excluded from the AP-frame theorem. Thus the theorem must not be described as counting all bicliques with cross-sum set A.

## Exact witness checks

The pinned characteristic-graph packet supplies square planes over F81, F729 and F2187. The local verifier independently implements base-three polynomial arithmetic, checks the square values and plane closure, enumerates every direction/gauge frame with duplicate detection, verifies all representative changes, constructs the non-AP counterexample, and checks the 432 explicit magic layouts. The exact unordered AP-frame counts are 432, 4,320 and 13,068 respectively. Fully ordered counts are 31,104, 311,040 and 940,896; frame/HP incidence counts are 15,552, 155,520 and 470,448.

For F81, an independent enumeration constructs all 1,080 affine lines of the ambient four-dimensional F3-space, examines all 582,660 unordered line pairs, and finds precisely the same 432 disjoint AP frames with the prescribed cross-sum set. This check does not assume in advance that the directions belong to U. The F81 pass additionally checks all 31,104 fully ordered frames under the explicit additive-to-magic permutation and obtains each of the 432 layouts exactly 72=q−9 times. The verifier also enumerates all six-vertex permutations and independently finds the 36 reversal classes of Hamilton paths in K3,3. All runtime guards remain enabled under python -O.

Source integers in extension fields encode polynomial coefficients in base three. They must not be treated as integers modulo q. The q=81 witness contains zero, so only the zero-allowed square convention is certified for it. There is no induced-biclique or integer-lifting inference.
