# Clone-stable readiness, binary carries and a literal infinite branch

**PROVED for all finite sources under the specified operations.** The small
checks are FINITE-EXACT support, not substitutes for the arguments below.
Notation: `H(empty)=1`; E_r(2) counts two ordered spanning path bins,
first nonempty ending at r, second possibly empty. A b-ready root has
`E_r(2)-H=b`. Dependencies are [READY](../two_entry_funnels_20261103/READY.md)
and [CL1/CL2](../return6325_completion_20261103/CLONE_LAW.md).

## S1. Complete classification of ready roots surviving cloning

Let T be strong and r be one- or two-ready. Replace r by a transitive
module r1->...->rm, m>=2. Among the **clone vertices**, a ready root exists
iff r originally had one in-neighbour p and `T-{r,p}` was transitive.
Then r was one-ready, r1 is one-ready for every m, and no later clone is
one- or two-ready. In particular no strong two-ready root survives this
operation at any clone mark. Re-marking outside the clone is not covered.

Proof. Deleting any one clone leaves a strong tournament, since at least
one clone remains. If a selected clone has indegree d, the inherited
defect decomposition is

    E-H = (q-start) + delta,
    q-start = sum of starting-path counts at its in-neighbours,
    delta >= d(d-1).

Every vertex of a strong tournament starts a Hamilton path, hence the
defect is at least d+d(d-1)=d^2. A later clone, or a first clone with two
old in-neighbours, has d>=2 and cannot be ready. For d=1, CL2 says the
first clone has defect `D=H(T-{r,p})`, an odd positive integer. It is1 or2
iff it is1, equivalently the deletion is transitive. The parent T-r is
strong by READY; p beats the first vertex of its transitive remainder,
so exactly one Hamilton path starts at p. Thus the original defect is1.
The converse and all m follow from CL2. QED.

## S2. An exact binary source grammar

For v>=1, write its binary digits least significant first as epsilon,
ending in1. Form a transitive chain with incidence word `w=00 epsilon`.
Add an apex p, with bit0 meaning p beats that chain vertex and bit1
meaning that vertex beats p. The first chain vertex r is the live root.
Call this rooted tournament R(v). It is strong, has `bit_length(v)+3`
vertices, and is one-ready. S1 says these are **all** clone-stable ready
roots, up to rooted isomorphism. The sole in-neighbour recovers p and the
transitive deletion recovers the word, so the parametrization is unique.

Cloning r into m vertices gives literally R(2^(m-1)v), after the explicit
chain-name transport: additional clones are the new initial zero bits.

Define

    c(v)=2-v
         + sum_(0<=i<t, bit_i(v)=0) floor(v/2^(i+1))
         + #{i<t-1: bit_i(v)=1, bit_(i+1)(v)=0}.

Then

    H(R(v))=4v+c(v),
    c(1)=1,
    c(2v)=c(v),
    c(2v+1)=c(v)-v-(v mod2),                  v>=1,
    H(R(2^(m-1)v))=4v*2^(m-1)+c(v).          (S2)

Proof of the counting formula. Deleting p from a Hamilton path leaves
two increasing chain runs, either possibly empty. Two paths place p
first or last. Otherwise let i be the last chain index before p, j the
first after p. Their bits are respectively1,0. If i<j, necessarily
j=i+1, giving one path for each adjacent10. If j<i, the i-j-1 middle
vertices choose either run freely, giving2^(i-j-1). The first two zero
bits contribute2v and v; the remaining zero-before-one pairs give the
displayed floor sum. Hence H=2+3v+sum+adjacencies. Appending a low zero
adds v to the floor sum and changes no adjacency count, proving the even
rule; appending a low one adds one adjacency precisely when v is even,
proving the odd rule. QED.

## S3. Two actual births giving3H+4 at every clone size

Take the word `0^(m+1)111`, the source R(7*2^(m-1)). Insert one new chain
vertex with bit0 between the first and second final1, then another bit0
between the second and third final1. Each new vertex has precisely the
incoming cut consisting of p and the chain prefix preceding it. These
are literal tournament births; all old relations and names remain.
Both intermediate and final tournaments are strong. The final word is
`0^(m+1)10101`, namely R(21*2^(m-1)). Since c(7)=-5 and c(21)=-11,

    28*2^(m-1)-5  -->  84*2^(m-1)-11 = 3H+4.        (S3)

The root and apex survive, and the output root is still one-ready.
`stable_words.json` contains literal cut receipts; the operation itself
is uniform for all m>=1. This extends a whole family, rather than only
a scalar identity. It does not give a second application to its output
family, or the +2 lane. Those remain open obligations.

## S4. Exact ternary carry and the obstruction to two fixed offsets

For a fixed v, suppose a target fixed core v' with a fixed shift of the
clone parameter satisfies H_target=3H_source+b for every sufficiently
large m. Coefficient comparison in S2 forces its dyadic coefficient to
equal3v. Since c is invariant under powers of2, the offset is uniquely

    b=c(3v)-3c(v).                                  (S4)

Thus this fixed-core dyadic grammar cannot supply both offsets2 and4
from one family. This does not exclude targets whose core changes with m
or a larger attachment grammar.

The actual stable18977 target has v4768=32*149, c=-95. At3v=14304,
c=-441, so its fixed-core offset is -156. Reaching the desired offsets
requires additional corrections158 and160. It is insufficient to replace
the old defect by a field phase or to solve the scalar equation only.

There is a finite integer carry recurrence available for that next search.
Let `(A,B,C)(v)=(c(3v),c(3v+1),c(3v+2))-3c(v)` and e=v mod2. Then

    (A,B,C)(1)=(-4,-2,-4),
    (A,B,C)(2v)=(A, A-3v-e, B),
    (A,B,C)(2v+1)=(B-2+4e, C+3v+3e, C-2+2e).         (S5)

Substitution of S2 proves these identities. Retaining v and e makes this
an exact binary digit transducer with unbounded integer registers. It is
not a finite-state realization of native returns. A finite search over
odd cores<=200000 finds no offset2 and only v7 for offset4; no all-size
exclusion is inferred from that search.

## S5. Converse is complemented reversal of the word

Fix the apex p. Global converse followed by reversing the chain replaces
w by complement(reverse(w)), with an explicit transport of every chain
name. It lies again in the displayed R grammar iff the last two original
bits are11. For t=bit_length(v), L=t+2, the new parameter is

    v* = [2^L-1-reverse_L(4v)]/4.                    (S6)

This is an involution on its domain and preserves H. It exchanges starting
and ending port responses and moves the selected root to the old last
chain vertex. It is not one native complement or a free physical re-marking.

An apex-preserving isomorphism to the converse must reverse the unique
chain, so it exists iff w=complement(reverse(w)). Consequently the number
of such words at order n is `2^((n-5)/2)` for odd n>=5 and zero for even n.
This counts apex-preserving self-converse sources; other maps moving the
apex are outside the claim. The involution, unlike an equality of numerical
counts, gives a concrete connection to mirrored-word constructions.
