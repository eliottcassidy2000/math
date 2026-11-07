# A projective cell carrier for the reduced Hamilton count

Status: PROVED, conditional only on the inherited all-order odd-cycle collection formula where Hamilton counts enter; FINITE-EXACT standalone controls. No priority claim. This is a structured geometric realization of `(H(T)-1)/2`, not an independent family already known to have exactly the two numerical gaps 3 and 10, and not a new proof of tournament spectrum completeness.

## 1. Input, construction, and an actual graded correspondence

For a finite tournament T, let C(T) be the finite set of its directed simple odd cycles of length at least three, counted up to cyclic rotation. Different cycles on the same vertex set are different elements. Let K_T be the simplicial complex of pairwise vertex-disjoint collections of cycles. Write

    I_T(z) = sum_{r>=0} c_r(T) z^r,       c_0(T)=1.

The inherited OCF theorem, THM-002-ocf.md, says `H(T)=I_T(2)`. Here H counts directed Hamilton vertex words. The empty tournament convention H=1 is harmless; all conclusions about a nonempty carrier require T to have a cycle.

Give R^{C(T)} a coordinate for each cycle. Define

    X_T = {x : sum_C |x_C|=1, and support(x) belongs to K_T},
    Q_T = X_T / (x ~ -x).

The empty cycle set gives X_T=Q_T=empty. X_T is the usual octahedralization of K_T: a signed simplex for every packing, with at most one sign chosen for each cycle. Q_T is equivalently the union, in RP^{|C(T)|-1}, of the coordinate projective subspaces supported on packings. No metric or embedding of the tournament vertices is chosen.

**Theorem 1 (graded cells).** Q_T has a finite regular CW structure with

    f_{r-1}(Q_T) = 2^{r-1} c_r(T),       r>=1.             (1)

Its cells are in explicit bijection with pairs `(P,[epsilon])`, where P is a nonempty odd-cycle packing and `[epsilon]` is a sign assignment to its cycles modulo simultaneous reversal of every sign. A face is obtained by deleting cycles and restricting signs.

Proof. A packing of size r supplies 2^r signed open (r-1)-simplices in X_T. The antipodal map pairs them freely. A closed signed simplex contains no antipodal pair of vertices and has disjoint antipodal image, so the quotient is injective on that closed simplex. Thus every quotient characteristic map is an embedding, and its boundary is a union of lower-dimensional quotient cells. This proves regularity and (1). Different quotient cells can have the same set of vertices; Q_T is not generally an abstract simplicial complex. In particular two disjoint triangles give two vertices joined by two parallel edges. QED.

Define the positive-degree cell polynomial

    D_T(t) = sum_{r>=1} f_{r-1}(Q_T) t^r
           = (I_T(2t)-1)/2.                              (2)

Consequently

    total number of cells of Q_T = D_T(1) = (H(T)-1)/2,   (3)
    chi(Q_T) = -D_T(-1) = (1-I_T(-2))/2.                 (4)

Total cell count and Euler characteristic are different invariants. The correspondence in Theorem 1 is an actual cell/packing bijection. Formula (3) uses the cited OCF equality; this packet does not construct an additional bijection from individual Hamilton paths to cells and one distinguished object.

## 2. Packing number as a cohomological height

Let nu(T) be the maximum size of a disjoint odd-cycle packing. It is also the maximum size of a disjoint directed-triangle packing: each directed odd cycle lies in a nontransitive induced tournament, which contains a directed triangle; choose one within each disjoint support. The converse is immediate.

**Theorem 2 (height).** If nu>=1, then dim Q_T=nu-1. Let w be the pullback of the standard degree-one class from ambient real projective space (zero when that space is RP^0), equivalently the first Stiefel-Whitney class of the canonical double cover X_T -> Q_T. Then

    w^{nu-1} != 0,       w^nu = 0.                       (5)

In particular `dim_F2 H^j(Q_T;F2)>=1` for every `0<=j<nu`.

Proof. Formula (1) gives the dimension. A maximum packing P supplies an embedded coordinate RP^{nu-1} in Q_T. The composite inclusion into RP^{|C(T)|-1} is the standard linear inclusion. The standard calculation

    H^*(RP^d;F2)=F2[a]/(a^{d+1})

shows that w^j restricts to the nonzero a^j for 0<=j<nu. Dimensional vanishing gives w^nu=0. For nu=1 the assertion w^0!=0 says Q_T is nonempty; w=0. QED.

The standard cohomology input is recorded in [Hatcher, Algebraic Topology, Chapter 3](https://pi.math.cornell.edu/~hatcher/AT/ATch3.pdf), Section 3.2, Theorem 3.19; it is not a new topological theorem here. The new use is the explicit tournament packing carrier.

A maximum packing also contributes every signed subpacking, so

    k(T) = (H(T)-1)/2 >= (3^nu-1)/2.                     (6)

This recovers the familiar OCF packing bound geometrically; it is not a stronger Hamilton lower bound.

## 3. Dominance join and the missing relative sign

Let T -> U be a dominance join, with every vertex of T beating every vertex of U. Every directed cycle lies within one factor. Hence

    K_{T -> U}=K_T * K_U,
    I_{T -> U}(z)=I_T(z) I_U(z),
    D_{T -> U}(t)=D_T(t)+D_U(t)+2D_T(t)D_U(t).           (7)

When both cycle sets are nonempty, `X_{T -> U}=X_T * X_U` and Q is its quotient by the diagonal antipodal involution. A mixed packing has TWO possible relative signs after the signs in each factor have been taken modulo simultaneous reversal. This is the structural origin of the coefficient 2 in (7); replacing this operation by an ordinary join of Q_T and Q_U would lose that relative sign. If a factor is acyclic, its empty face-complex is a combinatorial join unit and the carrier is unchanged; no assertion about joining topological spaces with the empty space is needed.

On values, (7) is the shifted multiplicative law

    a diamond b = a+b+2ab,       2(a diamond b)+1=(2a+1)(2b+1).

It is not ordinary addition or multiplication of a and b.

## 4. A deletion attachment description

Fix a tournament vertex v. The old cycle complex is K_{T-v}. Every new cycle C contains v, so a packing contains at most one new cycle. All cycles compatible with C lie in T-V(C).

The carrier Q_T is obtained from Q_{T-v} by, for each odd cycle C through v, attaching a cone on X_{T-V(C)} along the map

    X_{T-V(C)} -> Q_{T-V(C)} -> Q_{T-v}.                 (8)

The first arrow is the antipodal quotient. The cone apex is the cell indexed by C. To see this, choose the sign of C to be positive; the remaining cycle signs are now actual signs, not signs modulo reversal. On the base, where the C-coordinate is zero, simultaneous reversal is again identified, exactly as in (8). Distinct new cones can meet only through their bases. When the complement has no cycles, the cone on the empty complex means an isolated point and has no attaching boundary.

The new cells in this cone consist of the apex and one cone cell for each signed nonempty complementary packing, so their number is I_{T-V(C)}(2)=H(T-V(C)). Therefore

    k(T)-k(T-v) = sum_{C containing v} H(T-V(C)).        (9)

This is a geometric refinement of the inherited OCF vertex-increment identity, not a new scalar identity or a realization induction. Attaching one cone can introduce many cells.

## 5. The inherited leaf grammar becomes projective space

A directed triangle has exactly one odd cycle. The unique strongly connected tournament S4 on four vertices has exactly two odd cycles, both triangles and intersecting. Singleton SCCs have no cycles. Thus SCC order and singleton padding are invisible to Q_T.

For a dominance word containing r triangle SCCs and no other cyclic SCC,

    I_T(z)=(1+z)^r,
    Q_T is homeomorphic to RP^{r-1},                    (10)
    f_{j-1}=binom(r,j)2^{j-1},
    k=(3^r-1)/2.

Here r>=1; r=0 gives the empty carrier. Indeed K_T is a full (r-1)-simplex, X_T is the boundary of the r-dimensional cross-polytope, and the quotient is projective space.

For one S4 SCC and r triangle SCCs,

    I_T(z)=(1+2z)(1+z)^r,
    k=(5*3^r-1)/2.                                     (11)

For r>=1, Q_T is homeomorphic to two copies of RP^r glued along their common coordinate RP^{r-1}. Each of the two S4 cycles extends the same simplex on the r triangle cycles, and those are the two maximal faces of K_T. Consequently

    Q_T is homotopy equivalent to RP^r wedge S^r.        (12)

For completeness, each RP^r is obtained from the common RP^{r-1} by attaching an r-cell along the usual sphere double cover. After the first cell is attached, the second attaching map is nullhomotopic through the first disk; attaching it therefore wedges on S^r. For r=0, Q_T is two isolated points. Over F2 its Betti numbers are one in degrees 0 through r-1 and two in degree r, with the r=0 interpretation b_0=2. Its Euler characteristic is `(1+3*(-1)^r)/2`.

Session209's audited paired-black leaf classification applies these statements as follows. Non-self-converse leaves are suitable nonpalindromic singleton/triangle words, giving (10). Mixed leaves are A -> S4 -> A^op with A a singleton/triangle word; here r is even, giving (11)-(12). Conversely, not every tournament with the same scalar H need have this topology. This packet does not assert a converse classification from H.

For the r=3 pure leaf count H=27, the carrier is RP^2 with cell vector `(3,6,4)` and total 13 cells. The topology does not depend on the placements of singleton SCCs; it therefore deliberately forgets information that the leaf classification retains.

## 6. Exact hostile: H and order do not determine the carrier

Let R5 be the tournament on 0,...,4 with every arc forward in the natural order except 4->0. All its directed cycles contain the backward arc. Its odd cycles are three triangles (choose one interior vertex) and the one 5-cycle, so I_R5(z)=1+4z and H=9.

Compare the two order-six tournaments

    A = R5 -> singleton,       B = C3 -> C3.

For A, Q_A is four isolated points, with cell vector `(4)` and Betti numbers `(4,0)`. For B, Q_B is RP^1, a circle with two vertices and two parallel edges, with cell vector `(2,2)` and Betti numbers `(1,1)`. Both have H=9 and k=4. Their Euler characteristics are 4 and 0.

This is the first possible same-order failure of recovering the graded carrier from H: below order six no two nontrivial odd cycles can be disjoint, so I_T(z)=1+c_1 z is determined by H. The explicit order-six pair proves sharpness. At order six, exact exhaustive replay finds 19 distinct H-values but 27 distinct packing polynomials. The H-values with more than one polynomial are 9,17,23,29,33,37,45.

## 7. What this says, and does not say, about the two holes

The inherited all-order exclusions H not in {7,21} imply that tournament-derived Q_T never have total cell count 3 or 10. Conversely, proving that these are their only absent total cell counts is exactly the open tournament Hamilton-count spectrum conjecture, because of (3). No simplification of that inverse realization obligation is claimed.

Removing tournament realizability destroys the holes: for an arbitrary simplicial complex K consisting of m isolated vertices, its octahedral antipodal quotient is m isolated points, for every m>=0. In particular m=3 and m=10 occur. The projection to an arbitrary independence complex is therefore not a proof of the two-gap spectrum.

The source is T; the target is the regular CW carrier with its canonical double cover. The map preserves the complete odd-cycle packing polynomial, its grading, disjointness dimension, and the shifted multiplicative law. It loses SCC order, singleton padding, and in general the full tournament orientation. Equality of total cell counts loses even the polynomial and topology, as Section 6 shows. The theorem offers a structured geometric correspondence, not a one-to-one correspondence of tournaments with integers or with unmarked carrier homeomorphism types.

## 8. Verification and reproducibility

Run `python3 verify.py` and `python3 -O verify.py` in this directory. The script is standalone, uses explicit checks unaffected by optimization, and compares Hamilton subset DP against independently generated directed-cycle packings for ALL labelled tournaments of orders 0 through 6. For a mask, bit b corresponds to the b-th pair (i,j), i<j, in lexicographic order; bit 0 means i->j and bit 1 means j->i. Cycles are normalized by starting at their least vertex.

The script separately constructs quotient-cell boundaries over F2 from signs modulo global reversal, checks projective-space Betti numbers for r=1,...,6, and the mixed-center Betti numbers for r=0,...,5. Independent Hamilton DP checks the SCC controls through r=3. These small computations test the all-order derivations; they do not establish their quantifiers.

RESULT.json contains exact outputs, the same-order hostile witnesses, and a SHA256 digest of every exhaustive row `(n,mask,H,packing_coefficients)`. The source dependencies and their status are recorded in DEPENDENCIES.json. An independent review is required before publication.
