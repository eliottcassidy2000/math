# Order-eleven stability: separated cycle payments

Status: **PROVED, computer-assisted, same-session proof audit**. The complete
extension certificate, independent pattern kernels, exact matrix controls
and source review pass. No independent-agent or formal audit is claimed.
The nauty generator is an explicit dependency.

Let A be a real weighted tournament of order11. Set H equal to its weighted
directed Hamilton-path sum, u=det(I+A), v=(-1)^11 det(A),
Delta=H-u-|v|, q_e=A_ij A_ji, and V=sum q_e. Let C_j count directed
simple cycles once up to rotation, retaining direction. Write tau=C3,
and M_ab... for unordered vertex-disjoint cycle collections. A2-cycle
has weight q_e. W is the weighted sum of unordered incident edge pairs;
P21 is a wedge together with a vertex-disjoint edge.

Theorem:

    Delta >= (781/1280)V + (153/640)tau
              + C4/64 + C5 + C6/16 + C7 + C8/4 + C9 + C10.   (1)

The variance coefficient is attained by five ordinal fair pairs and one
singleton. The triangle coefficient is a safe reserve, not claimed optimal.

## Inheritance and signed ledger

The [order-ten packet](../square_sum_archive_stability_20261102/STABILITY_10.md)
supplies the exact determinant and path-cycle conventions. The
[concurrent next-order target](../lower_frontier_complement_fibonacci_20261102/NEXT_11.md)
already identifies the necessary7/40 wedge charge on P3 disjoint union4K2.
Its graph polynomial and extension compiler are inherited, not newly invented
here. The new obligations are a complete certificate and a compatible payment.
Use the inherited uniform margins u_r+v_r>=2^(1-r),
u_r-v_r>=r*2^(1-r), with v_r=(-1)^r det A_R, and tau>=W.

For a cycle partition t, H has coefficient2^len(t) if every part is odd,
otherwise0. The coefficient in u is (-1)^(number of even parts).
At full support11, v has the negative of that coefficient. Thus the two
signed gaps have full-cover multiplier b=0 or2 in H-bu, and all proper
partitions have multiplier1. The complete finite identities are generated
in ledger.py; no sign of det A is assumed.

Put

    P=V-M22+M222-M2222+M22222,
    T=tau+3M33-M223-M2233+M2223-2M22223,
    D5=C5+M25+3M35-M225,
    D7=C7+M27-2M227,
    D8=C8-M28.

For b=0 take epsilon=+1, and for b=2 take epsilon=-1. Define

    D4=sum_(4-cycles C) w(C)(u_(C^c)+epsilon v_(C^c))+M44-M244,
    D6=sum_(6-cycles C) w(C)(u_(C^c)+epsilon v_(C^c))+M46.

For EACH signed gap the exact residual after subtracting

    P+T+D4+D5+D6+D7+D8+M23+C9+C10

has nonnegative coefficients. Selecting one cycle in a collection introduces
its actual multiplicity, checked by the compiler. In particular the
M44/M244 corrections are required. The remaining positive terms need not
be spent a second time.

## The nontriangle cycle groups

For D4, the complement of two4-cycles has three vertices and V<=3/4.
Thus M244<=3M44/4, and the seven-vertex margin gives D4>=C4/64.
For D6, the five-vertex margin and M46>=0 give D6>=C6/16.

For D5, a six-vertex complement has

    M22 <= V^2/2-W,  tau>=W,  W>=V^2/3-V/4.

Consequently M22-V-3tau<=-5V^2/6<=0, giving D5>=C5.
For D7, each edge in a four-vertex complement has one opposite edge,
so2M22<=V/4. Hence1+V-2M22>=1 and D7>=C7.
The complement of an8-cycle has V<=3/4, giving D8>=C8/4.

## The triangle reserve pays the missing wedge charge

Fix two triangles. Their complement has five vertices and M22<=15/16;
therefore M2233<=15M33/16. Fix three pairs and one triangle. Two vertices
remain, with q<=1/4; each four-pair/triangle collection is counted four
times. Therefore M22223<=M2223/16. These two estimates give

    T >= tau+(33/16)M33-M223+(7/8)M2223
      =sum_(triangles C) w(C)[1+(33/32)tau_R-M22_R+(7/8)M222_R],

where R has eight vertices. Using tau_R>=W_R, the bracket is at least

    B8(q)=1+(33/32)W-M22+(7/8)M222.

This polynomial is separately affine in all q_e in[0,1/4]. At a box vertex
represented by an eight-vertex simple graph G,

    512 B8=512+33w-32m2+7m3 >=212.                    (2)

The supplied exhaustive triangle8 certificate checks all12,346 unlabelled
graphs, with a second explicit edge-subset kernel on the whole universe.
The minimum212 is attained only on two disjoint K4s. Thus T>=53tau/128.
This is a universal real-weight conclusion from a finite endpoint lemma,
not extrapolation from weighted random examples.

The original M23 remains completely unspent. Applying tau>=W on the
complement of each selected pair gives M23>=P21. The final graph obligation is

    P+P21+(7/40)W >= (781/1280)V.                    (3)

Pay (7/40)W once from T. Since W<=tau, the remainder is at least
(53/128-7/40)tau=(153/640)tau. The certificate below proves(3), completing(1).

## The complete graph obligation

At q endpoints, 5120 times(3)'s left side minus its right side is

    J11=499m-320m2+80m3-20m4+5m5+80p21+56w.          (4)

Every eleven-vertex graph is a ten-vertex graph G plus a neighbourhood S.
The inherited compiler gives

    J11(G+S)=J11(G)+sum_(i in S)L_i
                       +sum_({i,j} subset S)(80Qij+56),
    Qij=m-d_i-d_j+Aij>=0,
    L_i=499-320m1(G-i)+80m2(G-i)-20m3(G-i)+5m4(G-i)
                       +80ell_i+56d_i,
    ell_i=W(G-i)+d_i(m-d_i+1)-sum_(j adjacent i)d_j.

All quadratic coefficients are positive. Removing a selected i with L_i>=0
cannot increase the expression, so a minimum has S contained in{L_i<0}.
Moreover J11(G)+sum_(L_i<0)L_i is a valid lower bound. The checker discards
only bases where this lower bound is strictly positive; all other bases
enumerate EVERY subset of their negative coordinates.

For speed, Z(S)=2^|S| sum_k m_k(G[S])(-1/4)^k is an integer with

    Z(S)=2Z(S-i)-sum_(j adjacent i in S)Z(S-i-j), Z(empty)=1.

Then J11(G)=5120-5Z(G)-781m+80p21+56w and
L_i=(5/2)Z(G-i)-781+80ell_i+56d_i. These are exact coefficient identities.
The complete run includes all12,005,168 unlabelled ten-vertex bases without
graph filters. Of these,12,005,149 have a strictly positive automatic lower
bound; the remaining19 require only1,214 subsets of negative coordinates.
The checker makes143,360 direct neighbourhood comparisons with a second
matching-count kernel. A separate Python edge-subset implementation checks
ALL1,024 neighbourhoods of EACH exceptional base,19,456 checks, without
either pruning rule. All values are nonnegative; the empty graph attains0.
The exceptional list is in audit_graphs.json. The producer's `zero_subsets`
counts only the pruned minimizer search, not all zero extensions.
Completeness rests on the declared nauty ten-vertex graph universe and the
proven augmentation reduction. Separate affinity extends(3) to the real q-box.

## Sharpness and equality

For five fair pairs plus a singleton, H=1, u=(3/4)^5, v=0, V=5/4.
Thus Delta/V=(4/5)(1-(3/4)^5)=781/1280.
If equality holds in the unreserved bound with V>0, the positive triangle
reserve forces tau=0 and hence W=0. Positive-q edges therefore form a matching.
Every outside vertex sees both ends of an uncertain pair in the same
direction: otherwise one positive-weight orientation of the pair creates a
cyclic triangle. The resulting blocks form a transitive ordinal quotient.
Thus H=1, v=0 and Delta=1-product_i(1-q_i). Pad the matching to five slots.
The separately affine polynomial
1-product_i(1-q_i)-(781/1280)sum_i q_i has endpoint values0 at zero and five
fair pairs, and strictly positive values at1..4 fair pairs. Independent
Bernoulli interpolation can be supported on those two zero endpoints only
when all q_i=0 or all q_i=1/4. This proves the equality locus without assuming
a classification of the graph lemma's zeros. The singleton may occupy any
ordinal position. For V=0, tau=0 forces the binary tournament transitive.

No assertion about completeness of the raw odd Hamilton-path spectrum,
Collatz convergence, or the plane chromatic number follows.
