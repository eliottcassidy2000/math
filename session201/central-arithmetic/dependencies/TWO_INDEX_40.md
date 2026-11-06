# Two integer parameters and forty physical source-pair templates

Status: the parameterization, finite inverse, graph-order bound, and conditional
reduction below are **PROVED** from the named reviewed source and graft laws.
The assertion that every odd integer greater than 691 occurs in this family is
**OPEN**. Exact finite verification through 1,000,001 supports that assertion; it
does not prove it. Session 198.

## 1. An explicit two-parameter family

Let g(1)=1 and, for h>=1 in the even case and h>=2 in the odd case, put

    g(2h)=g(h),             g(2h-1)=g(h)+h-(h mod 2).              (1)

For L>=2, the reviewed one-apex gadget, after taking its converse, is an actual
nontrivial source W(L) with rooted profile

    w(L)=(D,E,U,d)=(g(L), L mod 2, L, 1).                        (2)

Its root-deletion graph is transitive and its order is 1+ceil(log2 L).
An explicit adjacency recipe is: write the complemented binary digits of L-1
from low to high, put a transitive chain on those digit vertices, orient the arc
from root 0 to digit vertex i exactly when the corresponding complemented digit
is one, and then reverse every arc. This is the same physical compiler as
`dependencies/ONE_APEX_SOURCE.md`, Section 1.

The reviewed operation J adds a new vertex which beats only the old source root,
retains the new vertex as root, and has profile

    J(D,E,U,d)=(D-E+U,D-E,D+U,D).                               (3)

It is an actual start-one source operation, not an assignment of independent
profile coordinates. Define S_0(L)=W(L) and S_1(L)=J(W(L)). Every index in this
paper is at least two. In particular the singleton is not silently admitted.

Use all ten reviewed two-source grafts

    F in {N,M,A,B,C,T,J0,J1,J2,J3}.

N and M have exactly two Hamilton cycles; the other eight have exactly three.
Their sources may be any two actual nontrivial sources, so every choice

    F(S_i(L),S_j(K)),   i,j in {0,1},   L,K>=2                 (4)

is an actual strong tournament. This defines forty templates, each with only
two free integer parameters. The letter J for a source operation is different
from the graft labels J0,...,J3.

The apparent larger family obtained by optionally applying one of O,Z,J to
W(L) is exactly (4), because the physical retained-root moves satisfy

    O(W(L)) has profile w(2L),
    Z(W(L)) has profile w(2L-1).                                (5)

These are full four-coordinate identities. We need equality of source profiles,
not an unproved graph-isomorphism assertion, to replace a source in a graft.
No claim is made that all physical sources or all source pairs have form (4).

## 2. The complete finite list of bivariate formulas

For a nontrivial source use the reviewed nonnegative coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E).

Then U=a+f+2e+2p and d=e+p>=1. Put r(L)=(g(L),L,L mod 2,1)^t.
The two source coordinate maps are the following explicit integer matrices:

    C0 = [ 1  0 -1 -1 ]       C1 = [-1  1  0  0]
         [-1  1  1 -1 ]            [ 0  0  0  0]
         [ 0  0  1  0 ]            [ 1  0 -1  0]
         [ 0  0 -1  1 ]            [ 0  0  1  0].

Thus y(S_i(L))=C_i r(L). In particular Jw(L) has the especially simple positive
coordinates (L-g(L),0,g(L)-(L mod 2),L mod 2).

Let P_F be the ten explicit matrices in Appendix A. They are the reviewed
physical graft matrices in these coordinates. The path count in (4) is exactly

    Phi_(F,i,j)(L,K) = r(L)^t C_i^t P_F C_j r(K).                (6)

Equation (6), the two displayed C matrices, and Appendix A are a finite explicit
collection of bivariate forms, not an existential source-grammar specification.
`FORMS.json` supplies every expanded 4 by 4 matrix C_i^t P_F C_j, with row and
column order (g(index),index,index mod 2,1). A key concatenates the graft label
and the two source bits, e.g. `J201` means F=J2,i=0,j=1.
There are forty named templates and thirty-nine distinct expanded matrices:
J0 with bits 10 and J3 with bits 01 both give

    3 g(L)g(K)+12LK.

For another example the N template with bits 11 is

    g(L)g(K)+4g(L)K+4Lg(K)+6LK-(L mod 2)K-L(K mod 2).          (7)

Parity is prescribed by the two indices. Neither it nor either g-value is an
independently adjustable integer parameter.

## 3. Complete finite bounds and an arithmetic inverse

Let u=(1,1,2,2)^t and z=(0,0,1,1)^t. The matrices satisfy entrywise

    P_F >= P_N,       P_N >= 3uu^t+zz^t.                       (8)

For physical source vectors y,v this gives

    Phi >= 3 U_1 U_2+d_1 d_2 >= 3LK+1.                       (9)

Both source covers are at least two, so either actual cover is at most
floor((n-1)/6) in a representation of n. This justifies the complete finite
source bank used in the checker. No estimated or empirically chosen cutoff is
used.

The ten matrices are closed under transpose: B and T exchange, J0 and J3
exchange, J1 and J2 exchange, and N,M,A,C are symmetric. Swapping i,j as well
therefore permits L<=K without losing any value in the full family. A complete
index domain for target n is

    2<=L<=floor(sqrt((n-1)/3)),
    L<=K<=floor((n-1)/(3L)).                                  (10)

There is also a simple congruence inverse, useful for an independent exact
boundary check. Fix F,i,j,L and compute

    (alpha,beta,gamma,delta)=r(L)^t C_i^t P_F C_j.

Then the right-index equation is

    n=alpha g(K)+beta K+gamma(K mod 2)+delta.                  (11)

Here alpha,beta are strictly positive. Indeed every entry in the second column
of P_F C_j is positive; every entry in its first column is nonnegative and the
last two entries are positive. Dotting with an actual source vector y, whose
last two entries sum to d>=1, proves the claim.

For each q in {0,1}, write K=2t+q. Equation (11) requires

    2 beta t = n-delta-(beta+gamma)q  (mod alpha).             (12)

Divide by gcd(2 beta,alpha), if it divides the right side, to obtain one residue
class for t. Enumerate that progression only inside the interval (10), then
test (11) using the binary evaluation of g. This is a terminating complete
membership procedure for the forty-template family, and returns a literal graph
recipe for every positive answer. It asserts no polynomial-time worst-case
bound. In particular exact nonmembership is meaningful only for this family,
not for unrestricted C2/C3 sources or tournaments.

## 4. A uniform logarithmic graph-order bound

Let k_F be the number of vertices added to the two source orders by the graft.
The exact values and product lower constants are

    F:     N  M  A  B  C   T  J0 J1 J2 J3
    k_F:   1  2  2  2  3   2   2  2  2  2
    c_F:   3  6  7  7 14   7   6  6  6  6.

Appendix A directly verifies P_F>=c_F uu^t entrywise and
c_F>=3*2^(k_F-1). Hence a graph represented in (4) has order

    v=2+ceil(log2 L)+ceil(log2 K)+i+j+k_F
      < log2(LK)+k_F+6
      <= log2 n+7-log2 3
      < log2 n+6.

As v is integral, this implies the convenient bound

    v <= floor(log2 n)+6.                                   (13)

This is unconditional for every value represented in this family. It does not
depend on a cofinite-coverage conjecture, on a bound for a shortest arbitrary
source recipe, or on finite enumeration.

## 5. The precise remaining arithmetic lemma

Let Q be the union of the forty images in (6). The proposed sufficient assertion
is the following entirely arithmetic statement with no remaining realizability
qualification:

> **Open two-index lemma.** Every odd integer n>691 lies in Q.

If this lemma holds, every odd n>=3 other than 7 and 21 has a strong tournament
realization with at most three Hamilton cycles and at most floor(log2 n)+6
vertices. Proof: values n>691 use (4) and (13). The complete finite inventory
through 691 leaves only the twenty values in (14) below; all eighteen other
than 7 and 21 have the literal C1 repairs in `REPAIRS.json`. Those repairs also
satisfy (13), as checked directly. This is an implication, not an assertion that
the open lemma or the global realization statement has now been proved.

The number 691 is exactly absent from Q, as both complete enumeration and the
separate progression inverse show. Thus a proposed universal upper-threshold
statement for Q cannot lower its threshold below 691. Calling 691 the largest
exception over all integers would require the open lemma.

## 6. Frozen finite evidence and physical controls

The standalone verifier reconstructs the full Q image through 1,000,001. Its
missing odd values from 3 through that bound are exactly

    3,5,7,9,11,15,17,19,21,33,35,37,39,47,57,61,63,71,145,691. (14)

It generates 282,167 distinct source profiles and tests 4,540,802 unordered
source pairs using the proved cover cutoffs. Every one of its 499,980 selected
positive certificates is independently replayed with the bit evaluation of g
and the expanded bivariate formula, including the graph-order inequality.
The recursion and bit evaluation agree for all 166,666 needed gadget indices.

Additional literal graph checks comprise both source choices for L=2,...,64
(126 full rooted profiles), all forty graft templates on L,K in {2,3}
(160 graph path/cycle counts), and all eighteen C1 repairs. The largest selected
finite-range certificate has 23 vertices. Normal and `-O` receipts must agree;
all guards use explicit exceptions, not removable Python assertions.

The separate progression inverse certifies 691 absent. It also finds the
previously hostile two-gadget target

    913131 = J2(w(2),Jw(39102)) = 9g(39102)+21*39102,
    g(39102)=10221.

This is a genuine C3 source-pair witness. Its right profile is
(49323,10221,49323,10221). It illustrates why the present two-index family is
not the same restriction as two successive indexed source gadgets.

All finite evidence remains bounded evidence. The verified modular source
universality theorems do not prescribe the positive height in (11), and hence
cannot replace the open exact-target lemma. No universal binary-child identity
or complete arithmetic descent is asserted here.

## 7. Dependencies and scope

The reviewed `dependencies/ONE_APEX_SOURCE.md` supplies the actual gadget (2).
`dependencies/POSITIVE_SOURCE_GRAMMAR.md` supplies (3), nonnegative source
coordinates and the physical matrices in Appendix A. Its complete grammar is
not needed to assert that this smaller family is exhaustive among sources.
`dependencies/BILINEAR_C2.md` and `dependencies/C3_TWO_SOURCE.md` supply the ten
physical gluing laws, their source boundary conditions, cycle counts and vertex
costs. The C1 repair adjacency lists are independently counted and do not need
a separate unproved inverse claim.

This strengthens the preceding finite-arity route by isolating only two free
indices and improving the conditional/unconditional witness bound to +6. It
neither proves the all-odd spectrum theorem nor shows that the earlier larger
arithmetic family equals Q.

## Appendix A. The ten positive graft matrices

Rows and columns are ordered (a,f,e,p). These matrices and C0,C1 completely
specify every formula in (6).

    N =
        [ 6  4 10  9]
        [ 4  3  7  6]
        [10  7 15 14]
        [ 9  6 14 13]

    M =
        [11  8 20 18]
        [ 8  6 14 13]
        [20 14 34 31]
        [18 13 31 29]

    A =
        [14 10 21 20]
        [10  7 16 15]
        [21 16 33 31]
        [20 15 31 29]

    B =
        [14 10 26 23]
        [11  8 20 18]
        [20 15 35 32]
        [19 14 34 31]

    C =
        [23 18 44 40]
        [18 14 34 31]
        [44 34 80 73]
        [40 31 73 67]

    T =
        [14 11 20 19]
        [10  8 15 14]
        [26 20 35 34]
        [23 18 32 31]

    J0 =
        [12 12 24 24]
        [12 12 21 21]
        [15 12 30 27]
        [15 12 30 27]

    J1 =
        [12 12 21 21]
        [12 12 18 18]
        [18 15 32 30]
        [18 15 31 29]

    J2 =
        [12 12 18 18]
        [12 12 15 15]
        [21 18 32 31]
        [21 18 30 29]

    J3 =
        [12 12 15 15]
        [12 12 12 12]
        [24 21 30 30]
        [24 21 27 27]
