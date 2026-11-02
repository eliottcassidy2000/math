# Recovering native meetings from odd ancestors and clocks

**PROVED for finite rooted inverse trees of odd Collatz states; FINITE-EXACT
checks on the selected SF5 family. Same-agent review, no independent audit.**
These are elementary reconstruction identities. No claim of literature
priority is made. They supply a positive refinement of the lattice obstruction:
an odd-state tree with its edge valuations reconstructs all joint native ports.

## J1. Shared halving tails are the missing identification

Use C(n)=n/2 for even n and3n+1 for odd n, and U=oddpart(3n+1).
Let V be a finite set of positive odd integers, with root r, carrying a
rooted tree whose edges u->q are actual U steps, stopped at r. Thus each
non-root has exactly one selected parent, all parent chains reach r, and
there is no selected cycle. Label the edge by

    a(u)=v_2(3u+1),    3u+1=2^a(u) q.

Its native expansion is

    u -> 2^a q -> 2^(a-1) q -> ... -> 2q -> q,        (J1)

of length a+1. Different children of the SAME odd q have different a.
Their even tails overlap, so subdividing every macro-edge independently
does not give an injective representation of the native integer graph.

For each q with at least one child, let A(q) be the set of incoming
exponents and M(q)=max A(q). Make one even chain

    2^M(q) q -> ... -> 2q -> q

and attach each child u to the vertex2^a(u) q on that chain. Chains for
different odd q have no common even vertex: every positive integer has
a unique odd part. They also contain no odd vertices. This construction
is exactly the union of the native paths of all selected sources to r.
It is a tree, with root stop at r, even if the un-stopped root belongs
to a cycle in the full arithmetic graph.

Consequently its vertex and edge counts are

    N = |V| + sum_(q with children) M(q),
    E = N-1.                                         (J2)

Independent subdivisions would give |V|+sum_edges a(u) vertices. The
exact duplicate count is

    sum_q [sum_(a in A(q)) a - max A(q)].             (J3)

For a binary odd tree this reduces to sum_q min A(q). These are actual
integer identifications, not a quotient chosen for visual convenience.

## J2. First meeting and the two clocks

Let u,v be distinct selected odd sources and q their lowest common
ancestor in the ODD tree. Let t_u,t_v be the sums of a+1 from u,v to q.
If one source is q, their first native common vertex is q, with clocks
t_u,t_v. Otherwise, let a,b be the exponents of the two distinct edges
entering q on these two paths, and put m=min(a,b). Then

    first native meeting = 2^m q,
    its two first-hit clocks = (t_u-m, t_v-m),
    native tree distance d(u,v)=t_u+t_v-2m.           (J4)

Proof: the two final even chains overlap precisely in2^m q,...,2q,q.
An earlier shared even vertex would have a common odd forward endpoint
strictly below q, contradicting that q is the odd LCA. An earlier shared
odd vertex gives the same contradiction directly. Therefore the top of
the shared tail is the first meeting, reached m steps before q by both
sources. The unique undirected path passes through this meeting. QED.

The smallest selected sibling example is

    13 -> 40 -> 20 -> 10 -> 5,
    53 -> 160 -> 80 -> 40 -> 20 -> 10 -> 5.

The odd LCA is5; the native meeting is40; the clocks are1 and3; the
distance is4. Treating the weighted odd tree as independent subdivisions
would give distance4+6=10. Its discrepancy is2 min(3,5)=6.

The finite checker tests every one of130305 distinct pairs among the511
odd nodes through SF5 depth8, comparing(J4) with literal native paths.
This preserves two asynchronous first-hit clocks. It does not silently
replace them with a synchronous meeting or declare them equivalent modulo
a clock. That distinction is inherited from
[the parity/tail quotient packet, section2](../collatz_tail_quotient_integration_20261009.md).

## J3. An exact shared-port generating recurrence

Let X_h(n) be the union of all native paths to n from the full SF5 odd
subtree of depth h rooted at n. Root n is a stopping boundary. Put

    Z_h(n;t)=sum_(x in X_h(n)) t^tau_n(x),

where tau_n is the native first-hit time, and let the two SF5 children
be (u,a),(v,b). Then

    Z_0(n;t)=1,
    Z_h(n;t)=sum_(i=0)^max(a,b) t^i
                + t^(a+1) Z_(h-1)(u;t)
                + t^(b+1) Z_(h-1)(v;t).             (J5)

The three terms are disjoint sets of INTEGER states: the root/even tail,
the child-u native subtree, and the child-v native subtree. Their odd
parts prove disjointness; the native clocks give the two shifts. Thus
this is an exact weighted recurrence without inclusion-exclusion debt.
It counts the shared tail once. In particular

    N_0(n)=1,
    N_h(n)=1+max(a,b)+N_(h-1)(u)+N_(h-1)(v).          (J6)

For n5,h1 the polynomial is

    1+t+t^2+t^3+2t^4+t^5+t^6,

and N_1(5)=8, whereas separate subdivisions count11 vertices. All
coefficients of(J5) are checked against independent path unions for
h=0..8. The recurrence is source dependent; replacing n by its mod9
class is not a closed all-depth update. The exact growing precision
requirement is proved in COLOUR_AND_CLOCKS.md, section3.

The same construction holds for ANY finite odd inverse tree of
U_p(n)=oddpart(pn+1), p positive and odd, with the same exponent labels
and native clock convention. It is local to the supplied tree. In
particular the5n+1 counterexample does not refute(J1)--(J6), because
none of them asserts that every source enters a chosen root.

## J4. The shared-tail Gram matrix and its exact rank

For m distinct sibling exponents0<a_1<...<a_m, represent a halving tail
by its edge indicator b_a=(1_(i<=a))_(i>=1). The shared-tail Gram matrix is

    K_ij=<b_(a_i),b_(a_j)>=min(a_i,a_j).

With L the lower triangular matrix of ones,

    K=L diag(a_1,a_2-a_1,...,a_m-a_(m-1)) L^T,
    det K=a_1 product_(i=2)^m (a_i-a_(i-1))>0.        (J7)

Thus it has rank m. The source-labelled sibling distance is
|a_i-a_j|+2 for i!=j: a line of halving stages with a private unit
edge for each odd source, rather than independent long arms at q.
These elementary nested-indicator identities are classical in form;
their use here identifies the exact port data hidden by the odd-tree picture.

For nonnegative edge weights c_i, replace a by S_a=sum_(i<=a)c_i.
The determinant becomes S_(a_1) product_(i=2)^m
(S_(a_i)-S_(a_(i-1))). Full rank holds exactly when each corresponding
stage block has positive total weight. Weight zero can erase a port
distinction; it does not certify that its arithmetic source disappeared.

Even a legal residue phase is insufficient for this response. At parent5,
exponents3 and9 agree modulo6 and give valid odd unit sources13 and853.
Against the fixed sibling53 of exponent5, their first meetings are40
and160, respectively. The two source clocks into those meetings are
(1,3) and(5,1). This is a source-specific, fixed-port separation of a
periodic exponent summary, within the already convergent basin.

More generally, take all leaves at odd depth h in the expanded SF5 tree
and their edge indicators along native paths to5. Their Gram matrix has
rank2^h: each leaf's first outgoing edge is used by no other leaf path,
so these columns contain an identity submatrix. The underlying finite
tree is contractible for every h. This supplies a literal arithmetic
instance of topology failing to bound the rank of a weighted response.
It complements, without identifying its operator with, the prior
[star-family return-Gram example](../lucas_mirror_gap_spaces_20261028/GAP_TOPOLOGY.md).

## What is and is not compressed

For a supplied odd tree, the root label and the exponents reconstruct all
odd labels by u=(2^a q-1)/3, subject to integral, odd, positive guards.
LCA plus the two final exponents supplies the earlier native meeting for
free; no independent table of all pairs is needed. The full tree and
the labels/clocks remain available. We do not claim the pair summary
alone is closed under arbitrary re-marking, new arithmetic operations,
or attachment to an unnamed source.

This gives a concrete instance of the same mechanism sought in tournament
composition: joint port data can sometimes be derived from a structured
carrier, but identifying independent branches without their shared native
tail changes the response. The81/27 tournament hostile is retained in
the verifier; there is no asserted map turning tournament automorphisms
into these arithmetic meeting clocks.
