# Exact all-size mean carriers for coloured clique attachments

**PROVED in the explicitly stated operation language; FINITE-EXACT controls.**
These are finite, named-clique observation modules. They are not new Ramsey
bounds, a classification of all colourings, or a claim to invent flag methods.
The broader counting framework is classical; compare Razborov's
[Flag Algebras](https://www.ias.edu/sites/default/files/math/csdm/05-06/arazborov_flag_algebras.pdf).
The exact carrier, integral observability proof and dimension below are the
specific deductions developed in this packet.

## C1. Universe, operations and terminal observation

Fix r>=3 colours and forbidden sizes t_1,...,t_r>=3. A source is an arbitrary
law on all r^binom(n,2) complete graphs with named vertices and exactly one
colour per edge. Allowed attachments add any finite number of named vertices,
fix every edge touching them independently of the source, and keep every old
edge. Observe the **separate vector** of monochromatic K_(t_i) counts, one
per colour. The theorem does not assume that only their sum is observed.

For each colour i and vertex set S write f_(i,S) for the0/1 indicator that S
is a colour-i clique. Retain:

* the constant coordinate1;
* for each edge, its indicators in r-1 colours; the last equals1 minus their sum;
* every named f_(i,S) with3<=|S|<t_i;
* each aggregate h_i=sum_(|S|=t_i) f_(i,S), when n>=t_i.

Take expectations for a probability law, or linear sums for signed weights.
The dimension is

    D_t(n)=1+(r-1) binom(n,2)
          +sum_i sum_(k=3)^min(n,t_i-1) binom(n,k)
          +#{i:n>=t_i}.                                  (1)

## C2. Exact closure and actual isolating probes

For a proper enlarged clique S, let A be its old vertices. If any fixed new
edge within S has the wrong colour, its indicator is0. Otherwise it is exactly
f_(i,A); for |A|<=1 this is1, and for |A|=2 use the edge coordinates. Because
|S|<t_i, the needed old flag is retained. For an enlarged terminal count,
the all-old contributions give the single old aggregate h_i; each remaining
term uses a proper old flag. This proves closure by an integer matrix under
every allowed attachment. Relabelling also merely permutes coordinates.

Every proper flag is observable by an **actual attachment**, not just by a
formal coefficient extraction. For |S|=k<t_i, append t_i-k new vertices
forming a colour-i clique. Give their edges to S colour i and every edge
to old vertices outside S one other colour. Then

    new h_i - old h_i = f_(i,S).                          (2)

Any new colour-i K_(t_i) must use all the added vertices and all of S: fewer
new vertices would require more than k eligible old vertices. A disjoint
new colour-i K_(t_i), with all cross edges another colour, similarly yields
the constant difference1. The original aggregates are already terminal
observations. Thus the entire carrier is in the integer span of actual
attachment queries. Other colours may create extra cliques, but the
observation vector keeps their components separate.

Unlike a spanning Hamilton path, a forbidden clique need not consume every
old vertex. This is why the isolating neighbourhood in(2) works; it is not
an all-size proof of minimality for the tournament endpoint carrier.

## C3. Independence and integral precision: no hidden modulus exception

For each old edge e use variables x_(e,i) for the first r-1 colours. The
last colour is1-sum_i x_(e,i). Functions on all colourings have an integer
basis of monomials choosing either1 or one x_(e,i) per edge. The one-edge
evaluation matrix has determinant+/-1; the full basis change is its tensor
power and is again unimodular.

For |S|>=3, the highest edge-degree part of f_(i,S) has support exactly
E(K_S). Different vertex sets give different supports. For each non-last
colour its leading monomial is pure in that colour. The last-colour flag
also has mixed leading monomials: choose different colours on two edges.
Since r>=3, such a mixed monomial exists and has coefficient+/-1, whereas
no other colour flag with this support contains it. Hence these r leading
polynomials are independent, with an integral unit triangular pivot block.

Sort by decreasing clique size and apply this argument to each support.
A terminal aggregate of colour i has the same leading coefficient on every
t_i-subset, but no individual flag of that same colour and size is retained.
Choosing any one such support gives its unit pivot. A mixed pivot separates
the last-colour aggregate even when other colours retain individual flags
at that size. At size2 the r-1 edge indicators and1 already form a basis.
This constructs a maximal minor+/-1 in the monomial coefficient matrix.

Consequences:

1. Equation(1) is the **minimum rational linear dimension** for unnormalized
   laws and this query language, at every n. Normalized laws have D_t(n)-1
   varying affine coordinates.
2. The coordinate-evaluation map on arbitrary integer source weights has
   Smith factors all1. The proof uses a unit coefficient minor and the
   unimodular evaluation transform; it does not require a chosen actual
   query submatrix itself to have determinant1.
3. For every integer m>=2 the smallest exact closed observation quotient
   on arbitrary source weights modulo m has m^D_t(n) states; with total
   mass fixed to1 it has m^(D_t(n)-1). Equality of carrier vectors is
   equivalent to equality of all future query responses by C2, and the
   coordinate map is surjective modulo m by the unit-minor argument.

These modular weights are not asserted nonnegative probability laws.
For positive laws the realizable carrier set is a specific finite convex
polytope; arbitrary coordinate assignments need not lie in it.

## C4. The requested Ramsey targets

For (3,3,s), s=5 or6, the formula specializes to

    D_s(n)=1+2 binom(n,2)+sum_(k=3)^min(n,s-1) binom(n,k)
                     +2[n>=3]+[n>=s].                    (3)

| n | D_(3,3,5) | D_(3,3,6) |
|---:|---:|---:|
|3|10|10|
|4|20|20|
|5|39|39|
|6|69|75|
|57|427466|4614572|
|91|2802349|49306807|

These are dimensions for **three-colour clique means**, not the tournament
dimensions7,28,126. In each case the retained species is dictated by the
terminal predicate and allowed attachment operations. The inherited
Ramsey cascade stalls where its neighbourhood counts cease to be saturated;
C2 retains the individual unsaturated clique coordinates instead of
replacing them by a forced degree.

This does not settle the Ramsey problems. A genuine probability law with
sum_i E[h_i]<1 contains a colouring with no forbidden clique, and a law
with all E[h_i]=0 is supported on such colourings. Proving an upper Ramsey
bound requires excluding zero violation for every actual colouring, or
equivalently for every realizable law. Linear independence, modular
surjectivity, and positive semidefinite moment matrices do not establish
that exclusion. No new numerical Ramsey bound is claimed.

The use of finite named clique indicators is related to, but not identified
with, asymptotic flag algebras. It also differs from a positive constructor:
equation(2) is a difference of observations on two genuine graphs.

## C5. Incoming tournament refinement

While this packet was developed, public commit2a836157af established the
pooled-parity tournament carrier. It improves the order-six upper bound661
to the exact minimum571; the earlier exact minima7,28,126 remain correct.
Its retained data are isolated labels, endpoint pairing, pooled interior
vertices and component parity. Subsequent Session146 commit14823946db
proves the exact order-seven minimum2668 with independent modular
certificates, and gives all-order support-lifting lower bounds. Its next
commit84154edbc1 proves minimum13028 at eight by certifying the2233
remaining directions. The next bounds are56296<=d(9)<=65908.
We inherit these advances; we do
not present661 as the tournament minimum.

Our all-size minimum theorem C3 is for the different, explicitly coloured
clique language. It neither supplies the missing all-order tournament lower
bound nor authorizes conditioning or native path complements on either
mean carrier. The earlier triangle-parity hostile remains a direct test of
those extra operations.
