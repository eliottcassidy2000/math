# Cross-lane synthesis: exact filters from selector geometry

**PROVED:** return-join and triangle-wrap length-polynomial recursions below,
from the new selector decomposition and the inherited odd-cycle formula.
**FINITE-EXACT:** the compact 10010 realization of the eight-filter.
**OPEN:** an arbitrary-depth dyadic filter grammar and target coverage.

## 1. Do not compress the wrong sidecar

The selector lane retains the complete labelled odd-cycle list and the
induced-deletion records. The Fourier lane requires each cycle LENGTH, not
only total occupied vertices and number of cycles. A bivariate packing
polynomial Z(t,u) can forget the difference between lengths (3,7) and (5,5).
These have the same t²u^10 term but different factors
(1+z³)(1+z⁷) and (1+z⁵)². The multivariate cycle-address sidecar retained in
the return-join theorem is enough; a scalar or the coarse bivariate polynomial
is not silently promoted to enough information.

## 2. A positive recursive compiler for the Fourier response

Let R=R(A,alpha;B,beta) be the return join: all cross arcs A→B except
beta→alpha. Let p range over alpha-start paths in A and q over beta-end
paths in B; let ell=|p|+|q|. With K from `fourier/PROOF.md`, one has

    K_R(z)=K_A(z)K_B(z)
      +2 sum_(p,q with ell odd)
           (1+z^ell) K_(A−V(p))(z) K_(B−V(q))(z).       (1)

All coefficients are nonnegative integers. To prove (1), use the complete
cycle decomposition in `selectors/RETURN_JOIN_SELECTORS.md`: no cross cycle
contributes the first product, and one cross cycle of length ell contributes
its two colours and its two monochromatic-slot choices. The remaining
packing is on the two induced complements. This is exactly the length-
weighted OCF expansion; no enumeration extrapolation is involved.

Thus the same return-join grammar compiles both an exact address predicate
and an exact Fourier multiplier, provided the deletion sidecar survives.
Equation (1) is a genuine construction equation for a candidate filter:
vanishing modulo S_m can be tested on the ACTUAL recursively built block.
It is not permission to solve for a formal negative residual and call it a
new physical block.

There is also a min-plus companion. Define kappa(empty)=0 and let lambda(X)
be the positive leading quotient at z=-1 in the Fourier theorem. Then

    kappa(R)=min(kappa(A)+kappa(B),
                  min_(p,q with ell odd)[1+kappa(A−p)+kappa(B−q)]). (2)

At the minimizing terms, lambda(R) is the sum of lambda(A)lambda(B) for the
first branch, if minimal, and 2ell·lambda(A−p)lambda(B−q) for each minimal
cross branch. Positivity of every such leading term prevents cancellation.
This gives a recursively certified parity-jet resource in the same grammar.

## 3. The eight-filter is a triangle wrap of one five-bit apex word

For a nonempty tournament A, form W(A) by adding u,v with

    u→A→v→u.

A cross cycle is u,p,v,u, where p is a directed path on an odd number of
vertices of A. Summing by its vertex set gives

    K_(W(A))(z)=(1+z)² K_A(z)
      +2 sum_(nonempty S subset A, |S| odd)
            H(A[S])(1+z^(|S|+2))K_(A−S)(z).            (3)

The inherited ordinary count consequence is H(W(A))=H(A)+Q_2(A): setting
z=0 in (3) gives H(A)+2 sum_|S|odd H(A[S])H(A−S), and K_A(-1)=0 makes
the latter odd-subset sum Q_2(A)/2. This scalar identity is credited as an
inherited triangle-run law, not a new discovery.

Now let A be the one-apex tournament on a transitive five-vertex base with
binary in/out word10010, in the selector lane's convention. It has six
vertices, H(A)=15, and Q_2(A)=216. Its length polynomial is

    K_A(z)=15+36z+39z²+36z³+39z⁴+36z⁵+15z⁶.

The literal wrap W(A) is the eight-filter in the Fourier lane. Relative to
its displayed old vertex labels, the canonical order (apex, five increasing
base vertices, u,v) is

    (2,0,1,3,5,6,4,7).

Equation (3), evaluated by exact induced-subset arithmetic, yields

    K_(W(A))(z)=231(1+z)S_8(z).

The accompanying control checks that this is the SAME literal tournament
under that permutation, not just another block with the same polynomial.
This explains the finite filter's architecture: a very small one-apex
address gadget, followed by a physical cyclic wrapper. The primitive-
eighth-root cancellations are consequences of the complete residual
identity, not an analogy between an adjacency matrix and a division algebra.

An all-size question now has a concrete starting point: which one-apex
words w make K_(W(A_w)) divisible by S_(2^a)? Word10010 supplies a=3;
there is no claim that an untested doubling pattern works for higher a.
The lower bound on filter size and the wrapper's response-order cost both
remain constraints on such a programme.

## 4. Independent arithmetic addresses coexist with erased phase

The converse lane's [bilateral filtered code](converse/FILTERED_CODES.md)
combines the eight-filter with six valuation-independent atom species.
It proves that ALL nontrivial mod-eight length characters can be erased
while each original arithmetic address survives in the unnormalized
channel amplitude. Two sources can have identical normalized length laws,
identical H, and identical complete unmarked path-cover profiles, yet remain
exactly distinguishable and selectable by their rooted source amplitude.

This is the useful spin/flip lesson: specify which part of the response an
operation erases and which information remains to encode a target. It is
not a universal dimension threshold or a claim of arbitrary target control.
