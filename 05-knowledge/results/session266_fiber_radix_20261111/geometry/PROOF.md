# The electrical defect of a padded strong core

**PROVED from the inherited SCC and weighted-path laws.** The exact rational
examples are new-constructor controls. No external priority claim.

Session24 already proves finite core words and eventual fixed-H polynomials.
Session23 already gives a uniform-simplex gap limit for fixed H/|Aut|, with a
richer neutral alphabet. The fixed-H gap limit below is a specialization of
that mechanism, not a newly discovered probability principle. The new object
here is an additive electrical statistic which survives arbitrarily long
singleton padding and refines the multiplicity compiler.

## 1. Turn one core into a two-terminal network

For a tournament C of order m let H=H(C), let a_v and b_v count Hamiltonian
paths starting and ending at v, and let W_C(u,v) count paths using u,v
consecutively. Add an ordinary source s before C and sink t after C. Every
Hamiltonian path is s, then a path of C, then t.

Normalize the weighted network by H. Its conductances are

    c(s,v)=a_v/H, c(v,t)=b_v/H, c(u,v)=W_C(u,v)/H.

Define rho(C) to be its effective s-to-t resistance, and define

    delta(C)=m+1-rho(C).

This is a rational, converse-invariant number, determined by the retained
weighted path network and its two endpoint distributions. It is not determined
by H alone.

**Bounds and equality.** For every nonempty C,

    2 <= rho(C) <= m+1;
    rho(C)=2 iff a_v=b_v for every v;
    rho(C)=m+1 iff C is transitive.

*Proof.* Hold s at1 and t at0. The trial potential1/2 at every core vertex has
energy1/2, since both boundary conductance sums are1. Effective conductance is
the minimum energy, so rho>=2. The trial is the minimizer exactly when it is
harmonic at every core vertex, which is a_v=b_v.

The normalized network is the average of the H unit-conductance path networks
s,P,t. Each such path has minimum energy1/(m+1). The minimum of their average
energy is at least the average of their minima, giving rho<=m+1. Equality
requires one potential to minimize every path energy. Along a path the unique
minimizing potential records each vertex's position, in steps1/(m+1). Two
different Hamiltonian paths cannot have that same potential. Hence equality
is equivalent to H=1, and the inherited tournament theorem says this is exactly
transitivity. QED.

Thus a nontrivial strong core has 0<delta(C)<=m-1. Endpoint-balanced cores
have maximal delta=m-1. This does not require vertex transitivity, although
vertex-transitive examples are endpoint-balanced.

Exact examples:

| core | H | rho | delta |
|---|---:|---:|---:|
| directed triangle |3|2|2|
| unique strong four-vertex type |5|445/219|650/219|
| regular cyclic five-vertex type R5 |15|2|4|
| Paley seven-vertex type |189|2|6|
| transitive three-vertex type (hostile) |1|4|0|

The four-vertex example shows why endpoint imbalance is a genuine coordinate:
its resistance is slightly greater than2 even though its converse is isomorphic.

## 2. An exact series law for separated SCC words

Let T have nontrivial strong cores C_1,...,C_k and singleton gap lengths
(g_0,...,g_k), each at least1. Let n be its order, H=product H(C_i), and let
s,t be its first and last singleton vertices. Normalize W_T by H.

The inherited rank-one boundary law makes each C_i and its neighboring
singleton vertices precisely the two-terminal network above. Different such
blocks meet at articulation vertices. The remaining singleton chain edges
have conductance1. Therefore ordinary series addition gives

    H R_(W_T)(s,t) = n-1-sum_i delta(C_i).                 (G1)

Here resistance is measured in the unnormalized W_T network before multiplying
by H. Every g_i merely lengthens a unit-resistance wire; it cannot alter the
finite defect sum. This is a physical network identity, not subtraction of
abstract tournament counts.

For example a separated triangle/strong-four pair has H=15 and defect

    2+650/219=1088/219,

whereas a separated single R5 core also has H=15 but defect4. The invariant
therefore distinguishes two same-H families for arbitrarily large orders.
Reversing the core order or replacing every core by its converse preserves
this scalar; the full ordered interface contains more information.

The condition g_i>=1 matters. With adjacent nonsingleton cores the connector
is rank one but not a single articulation, so (G1) is not asserted without
an additional boundary network calculation. Terminal singleton vertices are
also required to specify the same two-terminal observation.

## 3. Resistance-resolved multiplicity and its limiting distribution

Fix an attained odd h>1 and retain the finite catalog of ordered strong core
words W with product h. Write k(W),s(W),D(W)=sum delta(C) for its length,
vertex total and defect. Exactly

    binom(n-s(W)-1,k(W))

padding choices have every singleton gap positive, with zero if
n<s(W)+k(W)+1. Consequently their resistance-resolved generating function is

    sum_W z^(s(W)+k(W)+1) u^[D(W)] /(1-z)^(k(W)+1).       (G2)

The bracket exponent is a formal marker in the additive monoid of rational
defects, not necessarily an integer exponent. Equation(G2) is an exact refined
classification of the separated families. For the labeled channel each word
is weighted by1/product|Aut(C)|, then the order-n coefficient is multiplied
by n!, just as in the unrefined compiler.

Let r be the maximal core length and let W_r be the finite set of individual
maximal words, not merely their lists of H values. At fixed h and n→infinity,
all gaps are positive and the word has length r with probability1-O_h(1/n)
under the uniform ordinary-isomorphism-class law. Moreover:

- each W in W_r has limiting probability1/|W_r|;
- conditional on W, normalized gaps converge to uniform Dirichlet(1,...,1);
- conditional on the all-positive-gap event, the observable
  n-1-H R_(W_T)(s,t) converges to the finite distribution of
  D(W) for a uniformly chosen W in W_r.

For uniformly sampled labeled tournaments conditioned on H=h, replace the
word probabilities by normalized weights1/product|Aut(C)|. The gap law stays
the same. The two sampling measures must not be conflated.

*Proof.* One word contributes binom(n-s+k,k)=n^k/k!+O(n^(k-1)). Requiring
positive gaps changes it to binom(n-s-1,k) with the same leading coefficient;
the excluded boundary has degree at most k-1. Shorter words are lower order.
Uniform weak compositions, divided by n, converge to uniform simplex volume
by elementary Riemann sums. Weighted counts differ only by a constant attached
to W. Formula(G1) then gives the last statement. This proof also explains why
word sizes disappear from the leading mixing weights at fixed H, whereas
Session23's exponential fixed-H/Aut ensemble has Boltzmann size weights. QED.

In particular, for h=3^a the unique maximal word is a triangles. Thus the
limiting defect is the deterministic number2a, and the strong-core positions
converge to the order statistics of a uniform points in[0,1]. For h=15 the
two maximal words are (C3,C4) and(C4,C3), so the limiting defect is1088/219,
not the defect4 of the strong R5 realization. Growth frequency chooses one
network geometry even though both have the same total Hamiltonian count.

## 4. Connection contract and boundary

Source: literal tournament with separated strong cores and singleton gaps.
Target: series electrical network plus finite core-defect word.
Preserved: H, total order, two-terminal resistance, core order and normalized
endpoints before the final scalar projection.
Lost by D alone: order of cores, gap locations, chirality and internal arcs.
Repair: the ordered core and gap word, or Session265's recovered SCC interfaces.
Transfer: equations(G1),(G2) convert exact fiber multiplicities into a finite
limiting resistance distribution; the program constructs both sides literally.

No Collatz step acts on these circuits here. This is a new sidecar for the
multiplicity problem, not a proof of the odd-count spectrum conjecture.
