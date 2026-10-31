# One-copy compatibility for every common square-class port

**Status: PROVED scoped all-size composition theorem; FINITE-EXACT consequence controls.** The universal tournament raw-H spectrum remains OPEN. The graphs here are selected square-sum carriers on finite sets of positive integers, not necessarily induced or initial intervals. Hamiltonian cycles are undirected edge sets, with rotation and reversal identified.

The [previous packet](../square_sum_spectrum_overnight_20261031h/PYTHAGOREAN-SEAMS.md) proved a two-seam product for matched ports and a four-copy synchronization for different shapes. We now remove that duplication cost. Its matched-port conic construction remains useful and is the degenerate boundary of the new argument.

## Theorem: exact product with one copy of each source

Let G,J be finite square-sum graphs with designated edges e,f. Suppose both endpoints of each edge lie in the same positive square class d, and the class d is common to the two edges. Then there are positive integers L,R such that:

1. L^2 V(G) and R^2 V(J) are disjoint;
2. deleting the scaled e,f and adding two explicitly paired cross-edges again gives a square-sum graph;
3. its Hamiltonian cycles correspond bijectively to pairs of input Hamiltonian cycles containing e,f;
4. maximum degree is preserved at the cut endpoints, so subcubic inputs give a subcubic output;
5. either new seam is an invariant edge whose endpoints again have square class d.

The output has exactly |V(G)|+|V(J)| vertices. If the input port-conditioned languages have respectively M,N cycles, the output has exactly M N. Source weights, existing port constraints and branch histories transport under the two declared square injections. The marked partition and scale factors must remain as sidecars for inverse recovery.

For different primitive port shapes, the arithmetic search needs at most |V(G)||V(J)|+1 frames. This is a bound on exact rational group operations, not a polynomial bit-complexity claim. A sharper source-dependent collision bound is given below.

## 1. Normalize the actual ports

Write the two edges, with orientations to be specified, as

    e={d(gu)^2,d(gv)^2},   f={d(hx)^2,d(hy)^2},
    gcd(u,v)=gcd(x,y)=1,
    u^2+v^2=d w^2,       x^2+y^2=d z^2.

Here d is squarefree and all quantities are positive integers. Existence of w,z follows from the original square-sum edges after dividing out the common root factors.

If d is odd, each primitive pair has one odd and one even member; orient it with u,x odd and v,y even. Then w,z are odd. If d is even, all u,v,x,y,w,z are odd; orient u<v and x<y. These assertions follow from reduction modulo4 and8 and coprimality. This arithmetic orientation is part of the port map.

If u/v=x/y, square-scaling G by h^2 and J by g^2 makes the ports equal. Apply the previous primitive norm-conic join. It uses one copy of each graph. The remaining case has different shapes.

## 2. The obvious frame is a collision point on an elliptic curve

Set

    A=(vx)^2, B=(uy)^2,
    E: Y^2=X(X+A)(X+B),
    P=(u^2 x^2, d u^2 x^2 w z).

A and B are distinct positive integers, so E is nonsingular. P is an actual rational point because

    X(P)+A=d x^2 w^2,   X(P)+B=d u^2 z^2.

Its three nonzero square classes are(1,d,d). A point Q in the same class has rational positive square roots

    S^2=X(Q), U^2=(X(Q)+A)/d, V^2=(X(Q)+B)/d.

Write t=S/(uv)=r/k in lowest positive terms. The square roots of the two primitive seam sums are

    sqrt(d[(ur)^2+(xk)^2])=kdU/v,
    sqrt(d[(vr)^2+(yk)^2])=kdV/u.

They are integers: their squares are integers and their displayed values are rational. The physical scale roots before a common cancellation are L=hr,R=gk. Multiplying both seam roots by gh gives the original port equations. Dividing L,R and both seam roots by gcd(L,R) preserves those equations and makes the physical roots primitive.

For Q=P, t=x/v. This makes the scaled second endpoint of e equal the scaled first endpoint of f. The obvious solution therefore fails source ownership. The task is to move the point while retaining its square classes.

## 3. The marked point is always nontorsion

The standard tangent formula, obtained by substituting a line into the cubic, gives

    X(2Q)=(X(Q)^2-AB)^2 / [4 X(Q)(X(Q)+A)(X(Q)+B)].

For any rational Q with v2(X(Q))<0 and integral A,B, the numerator has valuation4v2(X(Q)) and the denominator has valuation2+3v2(X(Q)). Thus

    v2(X(2Q))=v2(X(Q))-2.

No numerator or denominator vanishes in this regime. Repeated doubling gives infinitely many distinct points.

For odd d, X(P),Y(P) are odd, whereas AB is divisible by16. Hence v2(X(2P))=-2.

For even d, X(P),A,B are odd and v2(Y(P))=1. The nonzero difference

    X(P)^2-AB=u^2x^2[(ux)^2-(vy)^2]

has valuation at least3. It is nonzero because u<v and x<y. Consequently v2(X(2P))>=2. At that next point, X is divisible by4 while A,B are odd, giving

    v2(X(4P))=-2-v2(X(2P))<=-4.

In both cases the preceding denominator clock proves that P is nontorsion. No rank calculation, BSD assumption, numerical orbit or unproved rational-cycle assertion is used.

## 4. Odd multiples retain exactly the required square classes

For a line through rational points on E, evaluate its intersection polynomial at each root0,-A,-B. If its three intersections have x-coordinates x1,x2,x3, then

    (x1-e_i)(x2-e_i)(x3-e_i)

is a rational square. The tangent case uses the same identity with a repeated factor. Thus doubling gives trivial square classes, and adding2P preserves the three square classes of any odd multiple of P. None of these points is a two-torsion point, since P is nontorsion, so no zero-coordinate exceptional case occurs.

It follows that every Q_m=(2m+1)P has the required class(1,d,d), and therefore gives two integral square-sum seams by section2. These points have pairwise distinct x-coordinates: equality would mean Q_m=Q_n or Q_m=-Q_n, contradicting nontorsion unless m=n. The positive physical scale ratios L_m/R_m are consequently distinct.

This is the ordinary two-descent square-class mechanism, with the needed nonexceptional case proved here by the line identity. For background on the group law and square-class homomorphisms, see [Noam Elkies, Math223 notes](https://people.math.harvard.edu/~elkies/M223.24/index.html). The nontorsion point, valuation proof and source-collision application above are explicit parts of this packet; no novelty claim is made for the classical group law or descent map.

## 5. A finite forbidden-ratio set gives universal entry

For fixed source labels a in V(G), b in V(J), a collision can occur only when

    L_m^2 a=R_m^2 b.

This fixes one positive scale ratio. Since all the ratios are distinct, each pair(a,b) excludes at most one m. At most |V(G)||V(J)| frames are bad. Among the first |V(G)||V(J)|+1 odd multiples there is therefore a disjoint frame. This is an all-size finite-exclusion proof, not evidence extrapolated from successful small searches.

Once a disjoint frame is found, put in only the two declared seams. The reverse two-edge-cut argument from the previous packet recovers a Hamiltonian path on each source side; restoring e,f recovers the two input cycles. This proves the bijection, degree claim, source clocks and reusable output port.

For weighted cycles, use each input's generating polynomial after deleting its required port edge. The output polynomial is the product of those two polynomials and the two seam variables. This preserves actual individual source weights; it does not replace them by an average.

### Exact forbidden-ratio compression

A collision requires the two original vertex labels to have the same square class. More sharply, let K be the number of distinct positive rational numbers sqrt(b/a) arising from source pairs for which that square root is rational. The proof needs only K+1 frames, since these are precisely the forbidden physical scale ratios.

For intervals[1,n] and[1,m], the forbidden set is exactly

    {L/R : 1<=L<=floor(sqrt(m)), 1<=R<=floor(sqrt(n)), gcd(L,R)=1}.

Indeed a collision for primitive L,R forces a=R^2 t,b=L^2 t; conversely t=1 supplies one whenever both bounds hold. This is a primitive Farey rectangle. For n=m and q=floor(sqrt(n)), its size is2 sum_(j=1)^q phi(j)-1. At2812 and2861 both root bounds are53, so there are exactly1763 forbidden ratios and1764 frames suffice, improving the coarse8045133-frame pair bound.

After composition, vertex sets are disjoint unions of square-scaled original intervals. If the left source charts are(a_i^2[1,n_i]) and the right charts are(b_j^2[1,m_j]), the exact forbidden set is the union over i,j of(b_j/a_i) times the corresponding Farey rectangle. Thus the original source chart and scale sidecars retain the entire collision predicate without factoring the much larger output labels. Duplicate ratios are merged only after retaining their actual source owners; weights on those owners are not identified by this set representation.

`ratio_obstructions.py` checks320 complete interval pairs and a scaled-chart union against literal pairwise arithmetic. The second native composition has3436 distinct forbidden ratios in this representation. These are exact finite sets, not probabilistic coprimality estimates.

## 6. A faster explicit clock in square class1

The companion `one_copy_ports.py` gives an alternative arithmetic route when d=1. Put A=(vx)^2,B=(uy)^2 and start with

    S=|u^2x^2-v^2y^2|,
    U=u^2x^2+v^2y^2+2v^2x^2,
    V=u^2x^2+v^2y^2+2u^2y^2,
    T=2wz.

Then U^2=S^2+A T^2 and V^2=S^2+B T^2. S,U,V are odd and v2(T)=1. The homogeneous update is

    S'=|S^4-AB T^4|,
    U'=S^4+2A S^2 T^2+AB T^4,
    V'=S^4+2B S^2 T^2+AB T^4,
    T'=2SUVT.

Divide all four coordinates by their common gcd. The first three coordinates stay odd, so this gcd is odd, and v2(T) increases by exactly one. The two quadratic identities are preserved by expansion. Consequently the reduced denominator of S/(uvT) gains one binary digit per step and eventually exceeds the explicit separation threshold. This works also when A=B, so it does not require a nonsingular elliptic curve.

The integer and rational implementations agree on972 stages over324 frames from18 primitive Pythagorean triangles. The all-step proof is the parity identity, not the observed number of stages. This provides an example of source-preserving refinement whose clock is controlled algebraically at every depth; it is not a transferred Collatz convergence theorem.

## Consequence controls and scope

- The3:4 and7:24 native ports now join with one copy each: exactly2^50 cycles on5673 vertices, replacing the previous four-copy2^100/11346-vertex compatibility construction for this task. The earlier construction remains valid.
- Actual joins in classes1,2,5,10 use3P in the retained examples. Their infinite-order witnesses have the expected strictly decreasing2-adic valuations.
- A hostile pair of finite arithmetic source sets blocks the first eight frames; the ninth,17P, succeeds. Those hostile sets are not claimed to carry Hamiltonian cycles.
- Independent integer projective coordinates audit the rational affine group arithmetic across twenty square classes. A signed sixteen-cycle fibre independently verifies the weighted product identity.
- Other source ports and constraints survive by explicit maps; the original numerical labels are scaled. This does not prove arbitrary endpoint extension inside the full initial-interval square-sum graph.
- General odd-multiple search has a finite source-dependent operation bound. The size of its intermediate integers remains a separate quantitative question.

## Anchor / niche / wildcard board

* **Anchor:** exact same-class port compatibility is now closed. The remaining interval problem is ownership of all intervening labels, not the existence of arithmetic cross-edges.
* **Niche:** represent forbidden scale ratios by square-class root data or primitive Farey rectangles, retaining original source weights and labels.
* **Wildcard:** turn the positive source-cost automaton from packet h into a changing-word tournament coverage induction. The present elliptic argument escapes a finite collision set; a Collatz or tournament obstruction cannot be declared finite without its own proof.

The recurring move is precise: preserve a boundary predicate along an explicit infinite orbit, prove that orbit has distinct observations, and bound how many observations the fixed sources can forbid. The non-torsion denominator clock and the actual finite source sets are essential parts of that implication.
