# An exact inverse and a classification of ternary left attachments

**PROVED for the inherited 65-presentation bank. FINITE-EXACT certificates
are identified separately. Global Hamilton spectrum completeness is OPEN.**

## I1. Inheritance and the actual coordinates

Use the source words, five distinct Q4 phases and actual path complement
of [the seam theorem T1--T4](../golden_seam_recurrence_20261103/THEORY.md).
Write v(w)=(h,x,y), v(cw)=v(w)P_c, v(empty)=(1,0,0). The bank contains
a,b,q0,q1,q2,q4,q5, every strong five-vertex chart, and the eight six-vertex
charts s14,s85,s248,s261,s262,s268,s330,s453. Equal matrices do not erase
their distinct source presentations.

Incoming [universal block closure G1--G4](../golden_hamilton_transfer_20261103/RECURRENCE.md)
and [its coordinate transport](../golden_hamilton_transfer_20261103/INCOMING.md)
resolve the arbitrary-block three-coordinate closure question left open
in our previous packet. Credit belongs to that incoming proof. Its
coordinates are (A,B,C)=(h+x-y,y-x,x+y). The inverse here concerns a
specified finite alphabet; universal block compilation does not supply
the existence of an arbitrary requested target block.

For a nonempty target with the processed boundary port p, let E count
Hamilton paths ending at p. Let C count a path starting at p together
with a possibly empty second path, spanning all vertices; the first path
is distinguished. Then the incoming formula identifies A=H-E and B=E.
There is an injection from Hamilton words to these covers: cut immediately
before p, taking the suffix as the distinguished path. It is bijective
onto covers whose second path can be followed by the distinguished one.
The remaining covers have a backward joining edge p->last(second path).
Consequently the new coordinates

    d=x+y-h=C-H,     e=y-x=E,     f=h+x-y=H-E             (I1)

are all nonnegative. More concretely, d counts permutation words with
exactly one backward adjacent step, whose head is p. This is a counted
defect, not a signed relaxation. Summing d over all ports counts every
permutation with exactly one backward step once. It recovers the B_T
coordinate in the incoming [ordinal-pair ear formula](../session153_gadgets_20261103/cuts/README.md).

The inverse change is

    h=e+f,       x=(d+f)/2,       y=e+(d+f)/2.           (I2)

Thus d+f must be even, and h must be odd. These are necessary guards;
they do not assert that every such triple has a source word.

## I2. Positive defect transfer and descent

Let Q_c be the conjugate of P_c under (I1). For the seven golden symbols:

    Q_a  = [[0,0,1],[1,0,1],[0,1,1]],
    Q_b  = [[2,1,3],[4,1,4],[5,1,4]],
    Q_q0 = [[10,1,7],[16,1,8],[17,1,8]],
    Q_q1 = [[2,3,9],[4,3,12],[5,3,12]],
    Q_q2 = [[7,1,6],[9,2,9],[9,2,10]],
    Q_q4 = [[8,1,5],[14,1,8],[17,1,8]],
    Q_q5 = [[2,2,7],[5,2,11],[6,3,11]].                 (I3)

All matrices use row vectors. In particular a sends

    (d,e,f) -> (e,f,d+e+f).                            (I4)

Every non-a matrix in the bank has integer entries at least one. All
initial atom vectors have nonnegative defects; each non-a initial vector
has x+y>=4. These explicit integer assertions are checked for all 960
strong marked charts through six, a larger bank than needed. They prove
the claimed invariants for arbitrary word lengths, not arbitrary atom
orders. The singleton shift is an inherited tribonacci clock in new
physical coordinates.

The rank R=d+e+f=x+y increases strictly on every nonempty extension. Under
a it becomes R+h. Under a non-a symbol it becomes at least3R. The first
symbol increases the empty rank0. Thus every inverse step strictly lowers
a nonnegative integer. A sharper depth estimate uses the original
coordinates: every P_c is entrywise at least P_a, and every row sum of
P_a^2 is at least3. If L=h+x+y, every two letters multiply L by at least3.
Hence any word producing v has at most

    2 floor(log_3 L)+1                                (I5)

letters. Inherited C_16 gives L<=33h. Each atom has at most six vertices,
so all source words attaining a fixed h in this bank have vertex order
at most6[2 floor(log_3(33h))+1]. This is a bound within the family, not
an unrestricted minimum-order theorem.

## I3. A realization-preserving inverse, with finite fibre polynomials

All65 matrices have nonzero determinant. For a proposed v, each possible
first symbol c has at most one predecessor:

    u=v adj(P_c)/det(P_c).                            (I6)

Retain it only if every coordinate is integral, u is the formal empty
state or satisfies positivity, odd h, (I1), C_16, and its rank is smaller.
Then recursively require a source word for u. Stop successfully at the
formal empty state. Every successful branch returns the literal word
c followed by that predecessor's word. This proves soundness by actual
matrix multiplication and completeness by removing the first symbol of
any actual source word. Strict descent proves termination.

In particular divisibility, a cone, or a modular match alone never counts
as acceptance. For example (25,0,25) passes the nonnegative defect and
parity conditions but has no source word. Its scalar H=25 is attainable.

There is an exact uniform generating recurrence. Let R_v(t) count source
words of full state v, with t marking the number of source vertices:

    R_empty(t)=1,
    R_v(t)=sum over guarded c of t^|c| R_(v P_c^-1)(t). (I7)

Nonmembers have polynomial zero. This is a finite polynomial computation,
not a circular membership definition: rank decreases in each recursive
call. Repeated counting matrices contribute separately for their distinct
physical charts. The scalar fibre polynomial is the finite sum

    Phi_H(t)=sum_(0<=x,y<=16H) R_(H,x,y)(t),            (I8)

using the stated guards. For H=1 subtract the empty contribution if only
nonempty sources are wanted. There are O(H^3) possible memoized states
below H and O(log H) polynomial degrees, so a straightforward implementation
uses O(H^3 log H) arithmetic operations at fixed65 and fixed cone constant.
This is polynomial in H, not in the bit length of H. It decides whether
a target is realized by this bank; it does not prove that every allowed
odd target is accepted.

## I4. Complete classification of left ternary attachments

**Theorem.** For every nonempty bank word of state (h,x,y), h>=3, and
c in {1,2}, some left word attachment has count3h+2c iff exactly one of
the following applicable guards holds:

    x=h+c:       prepend a;
    x+y=h+c:     prepend aa.                          (I9)

No state supports both c=1 and c=2 by left word attachments, even with
unbounded source-dependent choice of the attachment word.

Proof: every non-a P dominates P_b entrywise. It sends H to at least
h+4x+4y>=5h, which exceeds3h+4 for h>=3. All later counts are nondecreasing.
Three singleton steps send H to3h+4x+2y>=5h, so an all-singleton attachment
has at most two letters. Their counts are h+2x and h+2x+2y, giving(I9).
For both offsets simultaneously, the only possible cross-case would
require y=1 with x=h+1, or y=-1. Both violate y>=x>=0 for h>=3. QED.

This is stronger than the prior fixed-wrapper obstruction: it allows any
selector, any memory and any unbounded menu, provided the operation is
left concatenation in the specified source language.

## I5. All sources for the two-singleton offset branch

The aa branch in(I9) requires d=1 or2. There are exactly **66 literal
source words, 41 distinct states**, with that property in the65 bank:
31 words have d=1 and35 have d=2. Their largest H is155. These are
all-size finite certificates backed by the following exhaustive reduction.

Three leading a symbols make the final d equal the preceding R. If d<=2,
that preceding word has R<=2 and must be empty, a or aa: every non-a
initial vector has R>=4, and all subsequent ranks increase. Otherwise
the first non-a symbol c has at most two leading a symbols. Its positive
Q_c makes every resulting coordinate at least the preceding R; hence
the preceding word again is empty, a or aa. It suffices to enumerate

    a^n, 1<=n<=5;      a^i c a^j, 0<=i,j<=2, c!=a.     (I10)

The receipt contains every accepted word and state. Independent inverse
fibre counts agree on their full multiplicities. Thus the aa branch is
completely exhausted, not merely tested up to an order or H bound.

The a branch x-h=1 or2 is different. Its integer functional can have
cancellation. A leading a gives2(x-h)=-d-3f and cannot work. For every
non-a leading symbol except r12 and r58, all three coefficients of
2(x-h) in the predecessor's defects have the same strict sign. Negative
ones cannot work; positive ones force R<=4. The formal empty predecessor
is checked separately. For r12 the coefficients are(9,0,-3), giving a
multiple of3 after the parity guard, so it too is impossible. Only r58
retains the mixed guard3(f-d)+2e=2 or4. The inverse decides every given
state, but no complete finite classification of this remaining mixed
slice is asserted here.

## I6. Two-sided and arbitrary-position obstructions

Reverse the source word, and reverse vertex labels while taking converse
inside every atom. This preserves H, interchanges the two ends, and sends
q1<->q4 while fixing q0,q2,q5. The full960-chart bank is closed under this
map, and satisfies the same positive defect bounds. Therefore right
non-a attachments or three right singletons also cost at least5h.
Any two-sided wrapper reaching3h+2 or3h+4 must consequently be

    a^l w a^r,       0<=l,r<=2.                       (I11)

Exact inversion gives Phi_25(t)=8t^6+3t^7: **all eleven** H25 source words
in the65 bank. None of their nine panels(I11) has count77 or79. Thus even
arbitrary reencoding to another H25 word in this bank cannot provide the
two offsets using end wrappers. Both77 and79 themselves are attained.
The failure concerns extending a retained source, not global existence.

There is also an infinite obstruction allowing insertion anywhere in a
source word. Write h_n=h(a^n). For n>=6, no insertion-only extension of
a^n in the bank has count3h_n+2 or3h_n+4. With at least two inserted
letters, matrix domination gives H>=h_(n+2). A single inserted non-a
symbol dominates a triangle, and

    P_b >= 2 P_a^2-I

entrywise, so its count is at least2h_(n+2)-h_n. The remaining cases are
no insertion and one singleton, giving h_n or h_(n+1). Finally

    h_(n+1) < 3h_n+2 < 3h_n+4 < h_(n+2),   n>=6.

The upper difference h_(n+2)-3h_n begins6,12,22 at n=6,7,8 and obeys the
positive tribonacci recurrence; the lower inequality follows from its
increasing positive terms. In particular the H31 source a^7 cannot
reach95 or97 by any such insertions. The positive terminal replacements
in TRANSPORT use a different H31 state and deliberately change a suffix.
