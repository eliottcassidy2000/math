# Tensor parity compensates Pfaffian cancellation

Status: PROVED, independently audited. No new dimension cutoff is claimed
by this note.

## 1. An odd-coset tensor lemma

Let x_1,...,x_r be real with |x_i|<=1. Let b_S be integers indexed by
subsets of [r], with b_empty odd and every other b_S even. Define

    f_S(x)=sum_(J subset S) b_(S minus J) product_(j in J) x_j.

Then

    sum_(S subset [r]) f_S(x)^2 >= product_i(1+x_i^2).       (1)

Proof by induction. At r=0 an odd integer has square at least one.
Split the transformed coefficient vector at the final coordinate into
(a,b+x_r a). The coefficient vector producing a has odd empty coordinate
and even other coordinates, while the vector producing b is entirely even.
Thus a,a+b,a-b all satisfy the previous-dimensional odd-coset hypothesis.
Put L=product_(i<r)(1+x_i^2), A=||a||², C=||b||², D=<a,b>.
The induction hypothesis gives A>=L and A+C+/-2D>=L. Consequently

    |2D|<=A+C-L,
    A+||b+x_r a||²
      >=(1-|x_r|+x_r²)A+(1-|x_r|)C+|x_r|L
      >=(1+x_r²)L.

All coefficients used in the last comparison are nonnegative. This proves
(1) without bounding individual coordinates of the transformed vector.
In particular, the counterexample Pf=3x-2 does not refute (1): its odd
coefficient energy is retained in the companion coordinate.

## 2. Grouping all principal Pfaffians

Consider a nonempty weighted tournament A of order n. Suppose its nonbinary
unordered edges form a matching of size r; put x_i=K_(a_i,b_i),
K=A-A^T, q_i=(1-x_i²)/4. All other skew entries are +/-1. Let t=n-2r
be the number of vertices not in the uncertainty matching.

Arrange matching endpoints in consecutive pairs in the vertex order.
Group the even principal Pfaffians by fixing:
- which matching pairs contribute exactly one endpoint, and that endpoint;
- which of the t other vertices are included;
- all remaining matching pairs are then free to contribute neither or both.

Only groups whose fixed base has even size occur. The base is a binary
skew tournament, so its Pfaffian is an odd integer (empty Pfaffian=1).
The vector of Pfaffians in a group has exactly the form in (1): derivative
with respect to a full pair's adjacent skew entry deletes that pair with
positive sign. Its constant b_U is the Pfaffian on the principal subset consisting of
the fixed base together with all endpoints of pairs in U, with those full
pair entries set to zero. If U is nonempty, b_U is even. Indeed,
modulo2 its perfect matchings are those of an even complete graph avoiding
|U| disjoint specified edges. Inclusion-exclusion gives parity
sum_(J subset U)1=0, since every remaining complete even graph has an odd
number of perfect matchings, including the empty graph.

Apply (1) to every group. Exactly-one-endpoint choices contribute weight2
and odd parity, while a free pair contributes weight1+x_i² and even parity.
Thus the parity generating product is

    product_i[(1+x_i²)+2z] (1+z)^t,                     (2)

where only parity of the exponent of z is used.

The augmented-observer formula for u-v uses the same argument with a
fixed binary observer included in the base. It therefore takes the odd
rather than the even old-vertex parity in (2). The observer does not count
as an optional isolated vertex; it changes the parity selection only.
Here u=det(I+A), v=(-1)^n det A, and both Pfaffian sums have prefactor2^(1-n).

If t>=1, the parity sums are equal, and

    u+v >= product_i(1-q_i),
    u-v >= product_i(1-q_i).                           (3)

If t=0 (n=2r), the parity sums give

    u+v >= product_i(1-q_i)+(-1)^r product_i q_i,
    u-v >= product_i(1-q_i)-(-1)^r product_i q_i.       (4)

In particular, the sharp absolute matching-support margin floor is

    u-|v| >= product_i(1-q_i)-1_(n=2r) product_i q_i.   (5)

Binary tournaments r=0,n>=1 have floor1. The empty tournament is excluded.
The products in (4) are used only for r>=1. Ordinal weighted pairs and
singletons attain the appropriate bound, so the right side is sharp as a
function of the matching variances.

## 3. Exact nonlinear retained-baseline gap

Let Delta=H-u-|v| and B=E Delta(T), where T is independent binary rounding.
For a single currently unrounded edge e, the established interpolation
inequality has correction q_e (u_R-|v_R|), with R deleting its endpoints.
For a tournament whose uncertainty is a matching, apply (5) sequentially.
The positive parts of these correction lower bounds sum exactly to

    sum_i q_i product_(j later)(1-q_j)
       =1-product_i(1-q_i).

If initially t>=1, every remainder has an optional binary isolated vertex,
so there is no parity subtraction. If initially t=0 and r>=2, only the
first remainder lacks such a vertex: its subtraction is product_all q_i.
After that first rounding, its two endpoints are optional binary vertices
in every later pair-deletion remainder. For r=1,n=2, direct computation
Delta=B=0 yields the same formula without invoking an empty-remainder lemma.

Consequently, for every nonempty matching-uncertainty weighted tournament,

    Delta(A) >= B + 1-product_i(1-q_i)
                    -1_(n=2r) product_i q_i.           (6)

At r=0 use Delta=B. Ordinal pair/singleton constructions attain (6).
This is a nonlinear bound with the binary baseline retained; it is NOT
a false conditional inequality f(A)>=(1-q)E f(T). The full tensor grouping
and its integer parity are essential.

For an arbitrary weighted tournament, first round all edges outside a chosen
matching and discard their nonnegative corrections. Formula (6) then gives
an analogous lower bound using only the retained matching variances and
B for the complete original binary rounding. This does not by itself prove
the desired linear kappa_n bound at every remaining order.
