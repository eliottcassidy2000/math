# Exact seam transfer and the three coordinate quotient

**PROVED for words over the displayed or exactly certified atom alphabets.**
Finite atom checks are load-bearing certificates of explicit integer
identities, followed by induction for arbitrary word lengths. No universal
three-coordinate theorem for atoms of arbitrary size is asserted.

## T1 Literal sources and orientation gauge

A source is a dominance-ordered word of strong tournaments, each supplied
with the physical Hamilton path0,1,...,d-1. Concatenate these paths. Every
cross-component source arc goes from the earlier to the later component.
Complement every nonpath arc, retaining the path arcs, to obtain T(w).
Write h(w)=H(T(w)), counting directed Hamilton vertex words.

The vertex order used to process the target is the reverse of the supplied
source path. Within a new component, the relation is literally the reversal
of its path-complemented chart. Across a component boundary all old-to-new
arcs point forward, except the single arc from the first new vertex to the
last old vertex. This actual backward seam, not a free choice of attachment
orientations, is the operation being counted.

The golden alphabet is a singleton a, the directed triangle b, and the five
strong four-vertex charts q0,q1,q2,q4,q5. For qk, write k in three bits,
least significant first, indexed by(0,2),(0,3),(1,3). A bit1 means the arc
goes from the smaller to the larger endpoint; path arcs are forward. All
five diagrams are retained. Their sizes are1,3,4,4,4,4,4.

The golden word classification and its Fibonacci terminal flag are inherited
from [SCC arithmetic section2](../scc_basins_20261020d/arithmetic.md).
Its literal complement image and supplied width bound are inherited from
[golden complements C2 through C3](../fibonacci_parity_precision_20261103/COMPLEMENTS.md).
The general expiry mechanism comes from
[DEFECT_TRANSFER D2](../overnight_rank_return_20261029/DEFECT_TRANSFER.md).
Here endpoints are summed without separately querying n^2 endpoint pairs.

## T2 Six physical boundary states

During vertex introduction, retain each unfinished path component as its
ordered(head,tail) pair. Reject an attempted edge joining both ends of the
same component: that would close a directed cycle before a Hamilton path.
Each vertex has at most one chosen incoming and one chosen outgoing edge.

A processed vertex is active while a future vertex can send a backward
arc to it. When a head loses its last such future neighbor, it must be the
global path start; mark it S and reject a second such head. When a tail
loses its last exceptional future neighbor, forget its name and mark it *,
because it points to every future vertex. Its component association remains.
Do not choose an end in advance: the single remaining tail becomes the end.

At a component boundary only the final old vertex v can be active. There
are exactly the following six possible shapes:

    0: (S,*)                  3: (v,v)
    1: (S,v)                  4: (S,*),(v,v)
    2: (v,*)                  5: (S,*),(v,*)

State3 occurs only for a one-vertex processed prefix. A vertex cannot be
the head of one component and the tail of a different component. This
excludes the apparent extra pairings. The terminal observation is
`t=(1,1,1,1,0,0)^T`: states0 through3 finish one path; states4 and5 do not.

To introduce a vertex, independently choose an allowed existing tail as
its predecessor, or none, and an allowed head as its successor, or none.
Exclude choosing the same component twice; join the remaining components
through the vertex and expire their endpoint names as above. Add integer
multiplicities of equal resulting states. This is the recurrence implemented
by `advance`. The final vertex of a new block is retained as v for a future
seam; terminal observation handles the case of no further block.

**Correctness.** Each directed Hamilton path restricts on every processed
prefix to precisely such a path forest. Every edge is introduced exactly
once at its later endpoint in the processing order. Conversely, a history
accepted by t has indegree/outdegree at most one, no directed cycles and
one spanning component, hence one Hamilton path. Head expiry loses no
possible incoming edge; anonymous tails have identical future adjacency,
while their component pairing is preserved. The source decision is delayed
until forced, so it introduces no extra multiplicity. This proves both
directions, including n=1 and empty initialization.

This recurrence is valid for any supplied finite tournament block. Its
internal computation can be large; the next reduction is certified only
for the stated finite bank of blocks.

## T3 The seven positive matrices

Use row vectors. With columns and rows ordered(h,x,y), the matrices are

    P_a  = [[1,0,1], [2,0,1], [0,1,0]],
    P_b  = [[1,2,2], [4,3,4], [4,2,3]],
    P_q0 = [[1,4,4], [8,9,10], [8,8,9]],
    P_q1 = [[3,3,3], [12,6,9], [12,5,8]],
    P_q2 = [[5,3,4], [8,7,8], [6,6,7]],
    P_q4 = [[3,6,6], [6,8,9], [6,5,6]],
    P_q5 = [[5,4,5], [10,5,8], [8,4,5]].                (T1)

Initialize v(empty)=(1,0,0), and prepend a source symbol by

    v(cw)=v(w)P_c.                                    (T2)

Then the first coordinate is h(w). In particular all five Q4 phases have
different matrices. Their one-block vectors are respectively

    (1,4,4), (3,3,3), (5,3,4), (3,6,6), (5,4,5).

Equal one-block counts3 or5 do not identify phases: a left singleton gives
9 versus15 for q1/q4 and11 versus13 for q2/q5.

Here is a finite exact certificate proving(T2) for every word. Let M_c be
the6-by-6 transfer from T2 and let i_c be its empty-prefix initialization.
Set

    Z = [[2,0,0], [2,-1,1], [2,1,1],
         [2,0,2], [0,1,1], [0,1,1]].

The supplied matrices satisfy, entry by entry,

    M_c Z = Z P_c,             i_c Z = 2(1,0,0)P_c.    (T3)

The first column of Z is2t. Thus any six-state vector f projects to
v=fZ/2, preserving terminal H. Multiplying(T3) inductively proves the
recurrence for arbitrary word lengths. These identities are checked with
integers, without fitting a recurrence to initial terms.

From the matrix P_a,

    h(aw)=h+2x,       h(aaw)=h+2x+2y.

Hence the coordinates have the promised physical meaning

    x=(h(aw)-h(w))/2,     y=(h(aaw)-h(aw))/2.           (T4)

The apparent halves are integral and nonnegative: the positive matrices
and initial vector prove this directly. Redei parity supplies the same
integrality check. The recurrence uses O(number of components) arithmetic
operations for one H query; integer bit lengths grow and are not unit-cost.

## T4 Exact finite-bank extension

The same equations(T3), with the same Z, hold for every strong marked chart
of orders5 and6:50 and903 charts, respectively. The code exhausts all64
and1024 fixed-path charts, applies a literal strong-connectivity test, and
compiles each six-state transfer before solving and checking(T3).
All resulting P are nonnegative integer matrices with an odd top entry
and even lower entries in the first column. The full table, literal source
relations and matrices are in `atom_atlas.json`; initialization is rebuilt
and checked by the verifier.

This is a finite exact basis proof followed by all-length induction. It
proves a three-coordinate recurrence for words over **all strong atoms
through six vertices**, not just the original golden family. Equal matrices
may be merged for this scalar context observer: there are43 distinct
five-vertex matrices and544 six-vertex matrices. The full source chart label
must still accompany any request for an actual tournament relation.

There is no proof here that(T3) holds for every atom at every order.
Extending the finite table indefinitely is a separate problem. The general
six-state boundary mechanism remains valid, without this reduction.

## T5 Minimal left memory and exact modular precision

For any linear representation of dimension d, every matrix
`[h(u_i w_j)]` has rank at most d. Take u_i=a^i and w_j=a^j, i,j=0,1,2.
The resulting matrix is

    [[1,1,1], [1,1,3], [1,3,5]],       determinant -4.

Thus d>=3 over Q. The construction achieves3. This is a minimum weighted
linear carrier for the left-context language, not a minimum arbitrary set
encoding or a classifier of tournament isotypes.

For arbitrary weighted banks, let the state be an integer triple. All
future observation columns lie in

    L = Z e_1 + 2Z e_2 + 2Z e_3.

Indeed every P has its first column congruent to e_1 modulo2, and preserves
this column lattice. Empty, a and aa contexts already give columns
`(1,0,0),(1,2,0),(1,2,2)`, which generate L. Its Smith factors are1,2,2.
Therefore modulo m the exact future-observation quotient is

    Z/m x Z/(m/gcd(2,m)) x Z/(m/gcd(2,m)).             (T5)

It has m^3/gcd(2,m)^2 states for arbitrary weighted banks. The vectors of
a,aa,aaa form a unimodular row basis, so weighted combinations of actual
nonempty words span the whole integer triple lattice. Literal single words
are a narrower realizable set; their minimum modular observer can be smaller.
For example every such word has h=1 modulo2. This distinction retains
realization guards instead of treating any positive triple as a source.

## T6 Both ends require nine coordinates

Let

    A(w)=P_(last letter)...P_(first letter),
    h(w)=e_1^T A(w)e_1.

Prepending c right-multiplies A by P_c; appending c left-multiplies A by P_c.
Thus the full3-by-3 matrix is a nine-coordinate carrier for both-end contexts.
Three coordinates cannot suffice even on actual words:

    v(ba)=v(q5)=(5,4,5),
    h(baa)=9,       h(q5a)=13.                        (T6)

This hostile is already at source order4 and retains the same physical
complement definition. Both complemented tournaments have H=5; the boundary
placement is different. Neither arbitrary re-marking nor averaging is used.

The nine matrices indexed by

    empty, a, b, q0, q1, q2, q4, ab, aq1

have flattened determinant of absolute value128. They span M_3(Q).
The nine two-sided observations `(a^i, a^j)`, 0<=i,j<=2, separate that
space: the left and right Krylov matrices have determinants of absolute
values1 and4. The corresponding9-by-9 actual context Hankel determinant
has absolute value8192. Hence any linear representation for arbitrary
two-sided contexts has dimension at least9; the matrix carrier achieves9.
The certificate contains the full context matrix. No universal integral
Smith classification of the two-sided carrier is claimed.

## T7 Cubic clocks and a size recurrence for total mass

Repeating any fixed word w uses powers of the fixed3-by-3 matrix A(w), so
Cayley-Hamilton gives a scalar recurrence of order at most3. For individual
symbols the characteristic polynomials are

| Symbol | Characteristic polynomial |
|---|---|
|a|X^3-X^2-X-1|
|b|X^3-7X^2-9X-1|
|q0|X^3-19X^2-45X-1|
|q1,q4|X^3-17X^2-27X-9|
|q2|X^3-19X^2+23X-5|
|q5|X^3-15X^2-37X-21|

The a clock is the inherited March tribonacci interval family, not a new
scalar sequence. The two phases q1,q4 have identical repeated-word counts
but different mixed-context responses. A characteristic polynomial alone
therefore forgets relevant transport.

Let F_n=sum h(w), over **all marked golden source words of vertex size n**.
Set S(x)=xP_a+x^3P_b+x^4(P_q0+P_q1+P_q2+P_q4+P_q5). Summing every word gives

    sum F_n x^n = e_1^T (I-S(x))^(-1)e_1 = N(x)/D(x),

    N=1-x^2-6x^3-76x^4-72x^5+x^6+10x^7-7x^8,
    D=1-x-x^2-8x^3-95x^4-87x^5-64x^6-157x^7
        -629x^8-202x^9-27x^10-227x^11-585x^12.        (T7)

The exact polynomial determinant/adjugate calculation is independently
replayed. D retains degree12 modulo7, where gcd(N,D)=1. Hence N,D are
coprime over Q and the minimal scalar recurrence order is12. For n>=12,

    F_n = F_(n-1)+F_(n-2)+8F_(n-3)+95F_(n-4)+87F_(n-5)
          +64F_(n-6)+157F_(n-7)+629F_(n-8)+202F_(n-9)
          +27F_(n-10)+227F_(n-11)+585F_(n-12).

Initial values F_0,...,F_11 are

    1,1,1,4,32,154,465,1573,7501,31452,112644,411846.

This is Hamilton mass on a restricted marked family. It neither counts
ordinary classes nor supplies a new A000568 recurrence.

## T8 Fixed word wrappers cannot supply a nonzero constant offset

For any fixed source words u,v over the certified bank and any rational
alpha, there is no nonzero rational beta such that

    h(u w v)=alpha*h(w)+beta

for every sufficiently large singleton word w=a^n. In particular a fixed
pair of left/right word attachments cannot implement either3H+2 or3H+4
on the whole golden family, even if restricted to source counts H>=23.

**Proof.** Cayley-Hamilton applied to P_a gives

    P_a^3=P_a^2+P_a+I.

Both sequences h(u a^n v) and h(a^n) therefore satisfy the same tribonacci
recurrence for all n>=0. Their difference after multiplying the second
by alpha also does. If this difference is eventually the constant beta,
the recurrence forces beta=beta+beta+beta, hence beta=0 over Q.
Also h(a^n) is unbounded, so imposing H>=23 does not remove the tail.

The incoming
[Session151 ternary compiler target](../session151_spectrum_metagraph_20261103/strategy/TERNARY_COMPILER.md)
asks for a finite invariant bank and actual gadgets with outputs3H,3H+2,
3H+4. The argument here rules out fixed word concatenation on this entire
bank as the two offset gadgets. It does not rule out state-dependent
wrappers, a smaller invariant bank, changing the marked source relation,
or a different physical attachment. The conditional induction is inherited;
the obstruction is a direct consequence of the scoped seam recurrence.
