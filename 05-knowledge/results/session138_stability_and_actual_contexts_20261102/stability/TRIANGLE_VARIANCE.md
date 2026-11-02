# Triangle variance controls

Status: PROVED, independently audited as a dependency of the sharp
order-six and order-seven graph certificates.

For a weighted triangle, write q1,q2,q3 for its edge variances and tau for
the sum of its two directed cycle weights. Then

    tau >= q1q2+q2q3+q3q1.                         (1)

## Exact three-variable Bernstein proof

Partition into orientation orthants K_e=s_e z_e, with z_e in [0,1] and
s_e in {+1,-1}. In the tensor degree-two Bernstein basis, q_e=(1-z_e^2)/4
has coefficient 1/4 at coordinate indices 0,1 and coefficient zero at2.
Call an edge active when its index is0 or1. Since tau is separately affine
in each K_e, its Bernstein coefficient is its value at
K_e=s_e k_e/2, k_e in {0,1,2}. Therefore the active magnitudes are at most
1/2, while the inactive magnitudes equal1.

With zero or one active edge, the right side of (1) has coefficient0,
while tau's coefficient is an actual nonnegative triangle probability.
With two active edges, its coefficient is1/16. The triangle identity
4tau=1+xy+xz+yz, minimized over |x|,|y|<=1/2, |z|<=1, is at least1/4,
so tau>=1/16. With three active edges, the right coefficient is3/16 and
min_{|x|,|y|,|z|<=1/2} tau=3/16. These elementary minima follow by separate
linearity and checking the endpoint signs. Thus every tensor Bernstein
coefficient of the difference is nonnegative, proving (1) for the entire
closed cube.

## Graph consequences

Summing (1) over triples gives, for any weighted tournament of order n,

    tau(A) >= W(q)
      := sum_i sum_{j<k, j,k != i} q_ij q_ik.

Put d_i=sum_{j!=i}q_ij. Then

    W(q)=(1/2)sum_i d_i^2-sum_e q_e^2
         >= 2V^2/n-V/4,

using sum_i d_i=2V and 0<=q_e<=1/4. Hence

    tau(A) >= max(0,2V^2/n-V/4).                    (2)

This explains the natural pair-block variance scale V=n/8. A maximal fair
matching at even order has V=n/8, tau=0, and saturates (2).

## Hostile to a too-coarse order-seven reduction

A possible route to sharp n7 stability was to retain the matching gap
P=V-M22+M222 and replace tau by W(q), then attempt

    P+(5/8)W(q) >= (37/48)V.

This polynomial is separately affine in q, so its endpoint check reduces
to undirected graphs with q_e in {0,1/4}. Exhaustive enumeration refutes it:
the graph C4 disjoint-union K3 is a witness. It has m=7 edges, W=7,
M22=14, M222=6 in unweighted matching counts. The scaled residual

    384 [P+(5/8)W(q)-(37/48)V]
       =22m-24M22+6M222+15W
       =-41.

This is NOT a counterexample to the weighted tournament stability
conjecture. It only shows that replacing directed triangle probability by
this coarse local variance bound discards essential orientation and
magnitude information. The same q-endpoint graph requires larger actual
triangle weights than the coarse bound can see.

`graph7.cpp` checks all 2^21 labeled undirected graphs on7 vertices exactly;
`graph7.json` records the minimum. Its role is a hostile control for this
failed proof route, not a positive stability certificate.
