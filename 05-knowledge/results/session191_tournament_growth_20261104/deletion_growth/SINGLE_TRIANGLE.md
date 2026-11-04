# Reduce the single-triangle marked sector to strong endpoint responses

**PROVED reduction from the inherited one-arc-join law and Busch's cited
strong-tournament minimum; independently audited.** This does not prove
that a global minimum-increment root lies in exactly one triangle.

Write H for Hamiltonian vertex words, and Delta(T,v)=H(T)-H(T-v). For a
rooted tournament (A,a), let R_a(A) count ordered spanning two-bin path
covers with the second bin nonempty and starting at a. The first bin may
be empty. Write e_a(A) for Hamiltonian paths ending at a. The dual response
R'_b(B) has its first nonempty bin ending at b; write s_b(B) for paths
starting at b. These are the inherited Session29 one-arc-join conventions,
also used in INTERVAL_CHAINS.md.

## 1. Exact normal form and response

Suppose T is strong and v belongs to exactly one directed triangle. Put
A=N^-(v), B=N^+(v). Both sets are nonempty. The triangles through v are
in bijection with arcs B->A, so exactly one such arc exists, say b->a.
All other cross-arcs go A->B. Consequently a reaches every vertex of A,
and every vertex of B reaches b: otherwise the strongness of T would fail
on the source component of A or the sink component of B. Conversely these
two reachability conditions make the displayed construction strong.

Then

    H(T-v) = H(A)H(B)+R_a(A)R'_b(B)-e_a(A)s_b(B),
    H(T)   = H(A)H(B)+2R_a(A)R'_b(B),
    Delta  = R_a(A)R'_b(B)+e_a(A)s_b(B).                    (1)

The first line is the inherited one-arc-join identity. For the second,
apply that same identity to A'=A->v and B, reversing b->a. Here H(A')=H(A),
e_a(A')=0, and R_a(A')=2R_a(A): in a two-bin cover v may be appended to
either bin, and these are all possibilities. Subtraction proves (1).

If |T|>=4, T-v is automatically strong. Indeed b reaches a and then all A.
The vertex a also reaches b: either A has another vertex, reachable from a
and pointing to b, or B has another vertex, reached from a and itself
reaching b. Every B-vertex reaches b; every A-vertex other than a points
to b. Thus all vertices reach b, and b reaches all vertices. In particular:

**A weak deletion of a vertex from a strong tournament of order>=4 requires
at least two directed triangles through that vertex.**

## 2. Minimum source-root response and weak-side reduction

Always R_a(A)>=H(A): split each Hamiltonian path just before a to obtain a
unique counted two-bin cover. Suppose a reaches all A, and let C be its
source strong component, of size k. Write D for the ordered union of all
remaining components. Each of the two bins has its C part before its D part,
so

    R_a(A)=R_a(C) Z(D),
    Z(D)=sum_{S subset V(D)} H(D[S]) H(D-S) >= 2^|D|.       (2)

The first bin's C part may be empty, but its D part remains arbitrary.
Equality in the inequality holds for transitive D. If D is nonempty then
e_a(A)=0. If A is strong then e_a(A)>=1, by a Hamiltonian cycle.

Use the inherited, cited Busch minimum B(m) for strong order m>=3:

    B(3q)=3*5^(q-1), B(3q+1)=5^q, B(3q+2)=9*5^(q-1).

Set B(1)=1 and exclude m=2 from the strong domain. The source-root minimum
R over all order-m tournaments in which the root reaches all vertices is

    r(1)=1, r(2)=2, r(m)=B(m) for m>=3.                    (3)

For the weak ones it is

    w(2)=2, w(3)=4, w(m)=2B(m-1) for m>=4.                (4)

For strong A, (3) follows from R>=H>=B and is attained by an optimal interval
chain with its last spine endpoint marked; that endpoint has R=H. For weak
A, (2) gives B(k)2^(m-k), with k=1 or k>=3. The minimum is at the largest
permitted k<m: successive Busch ratios are strictly less than2, and moving
k=1 to k=3 reduces the value by the factor3/4. Thus k=m-1 for m>=4,
attained by an optimal rooted chain followed by one sink singleton. The
small orders2,3 give the other values directly. Formula (3) then follows
because the weak minimum is larger than B(m) for m>=3.

**Exact weak-side floor.** Among the single-triangle rooted strong
constructions of order n>=6 for which A or B is weak, the minimum Delta is

    2B(n-3).                                               (5)

Here a weak side has zero corresponding endpoint count, so (1) is exactly
the product of its R with the opposite response. Minimize w(a)r(b) for
positive a+b=n-1 and a>=2. If a>=4,b>=3, the product is
2B(a-1)B(b)>=2B(a+b-2)=2B(n-3). The product inequality follows by
concatenating optimal interval-chain compositions, yielding a strong chain
of order (a-1)+b-1. If a=2 or b=1, equality is attainable. The a=3 or b=2
cases follow from B(t+1)<2B(t), with the small boundary cases checked
directly. This also covers the possibility that both sides are weak.

For n=9+3j,10+3j,11+3j, the weak-side floors are respectively

    30*5^j, 50*5^j, 90*5^j,

strictly above the proved interval-chain candidates26*5^j,46*5^j,78*5^j.
At n=8, the same weak-side floor is18, above the chain minimum16.
Therefore any global minimizer at n>=8 that has exactly one root triangle
must have BOTH A and B strong. This is a genuine reduction within the
single-triangle sector, not a proof that arbitrary minimizers belong to it.

## 3. Remaining structural gap

There are two unproved steps toward the proposed all-root chain minimum:

1. Show that some global minimizer can be chosen with exactly one triangle
   through its root, or directly exclude every multi-triangle root below
   the chain candidate.
2. In the remaining two-strong-side construction, bound the response
   expression R_a(A)R'_b(B)+e_a(A)s_b(B) sharply enough to attain its minimum
   at the appropriate interval-chain endpoints.

Busch's scalar H bound alone does not control this response product and
endpoint correction. It is not a valid substitute for either step.

## 4. Greedy root-triangle removal is false

The missing first reduction cannot simply say that reversing some B->A
cross-arc always decreases Delta while retaining strongness. Take root0 in

    (48,5,9,3,14,30).

Here A={1,2,3} is a directed triangle, B={4,5} is transitive with5->4,
and all arcs go A->0->B->A. There are six root triangles. The tournament
is strong with H=33,H(T-0)=3,Delta=30. Reversing any one of the six
B->A arcs leaves the tournament strong and removes one root triangle,
but EVERY such reversal raises Delta:

* reversing4->a, a in A: H=37,Hminusroot=5,Delta=32;
* reversing5->a, a in A: H=43,Hminusroot=9,Delta=34.

The accompanying script checks this obstruction by literal enumeration of all 6! and 5! vertex permutations for all seven presentations, together with reachability checks. This is not a counterexample to the proposed
minimum law: its Delta30 is well above the order-six minimum6. It is a
counterexample only to an unrestricted greedy smoothing lemma.

The explicit witness can be checked with `python3 verify_smoothing.py`. Its finite computation is separate from the all-order reduction above. The strong minimum B is the cited Busch bound, as credited in the accompanying interval-chain note.
