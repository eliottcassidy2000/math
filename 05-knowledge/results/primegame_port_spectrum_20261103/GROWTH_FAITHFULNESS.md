# Actual growth recovers every Boolean derivative and every labelled source

**PROVED, all finite orders.** Elementary proof below, supplemented by exact
finite controls. Same-agent audit; no literature-priority claim. This is an
observability theorem, not a proof of Hamiltonian-count spectrum completeness.

## 1. Specify the operations before taking a quotient

Let H(T) count directed Hamiltonian vertex words. Retain the named vertices
of T as ports. Two allowed operations add one vertex without changing any
old arc:

* S appends a universal sink.
* B_x appends z with z->x and v->z for every other currently present v.

Here x is a specified original port. Compositions act rightmost first.
The operation B_x is a singleton-cut specialization of the inherited
[insertion calculus](../zero_heavy_ports_20261020g/profile_growth.md), not a
new general attachment operation. Names and creation order are retained.

For distinct original x,y define the four-query response

    L_xy H(T) = H(B_y S T) - H(B_x S T)
                 - H(B_y B_x S T) + H(B_x B_y S T).            (1)

Every term counts an actual Boolean tournament. If T comes with a specified
Hamiltonian path, each term keeps it and appends the newly created vertices
in their creation order. The first sink ensures that no later chosen port
is the current path endpoint. Thus these are legal extensions of the
supplied path, retaining the original labels, creation order and attachment
length as explicit data.

## 2. Growth-only derivative theorem

Write T^(x->y) and T^(y->x) for setting just that old edge in the two ways.
Then at every order n>=2,

    L_xy H(T) = H(T^(x->y)) - H(T^(y->x)).                     (2)

The right side is the signed Boolean derivative, independent of the current
orientation of xy. It need not be positive or nonzero. The theorem recovers
a numerical edge-edit response from growth observations; it does not execute
that old-edge edit by appending vertices.

**Counting proof.** For a parent A, let e_x count its full paths ending x,
and let D_x count ordered covers (P,Q) of all its vertices, with Q nonempty
and starting x; P may be empty. Let C_xy count such covers with P nonempty
and ending x, Q nonempty and starting y. Deleting the new vertex gives

    b_x := H(B_x A) = H(A)-e_x+D_x.                          (3)

Paths ending in the new vertex contribute H(A)-e_x; all others continue
immediately to x and contribute D_x. No joining arc between P and Q is assumed.

Put a new z at x and w at y. Swapping their creation order changes only
the edge zw, after aligning their names. Paths not using that edge cancel.
Those using z->w have either an old prefix ending outside x and an old
suffix beginning y, or an empty old suffix. Reversing z,w gives

    H(B_y B_x A)-H(B_x B_y A)
       = b_y-b_x+2(e_y-e_x)-C_xy+C_yx.                       (4)

The endpoint correction is essential. For the transitive two-tournament
0->1, b=(3,1) but the two double births both have H=5.

Now take A=S T. Every full path ends at the new sink, so e_x=e_y=0 for old
ports. Also C_xy-C_yx equals the edge derivative of H(A): concatenate the
two paths across xy, or delete that adjacent edge from a full path. This
derivative equals the derivative of H(T), since a sink leaves H unchanged.
Equation (4) rearranges to (1)-(2). QED.

A useful independently checked marginal law, for arbitrary A of order n, is

    sum_x H(B_x A) = (n-2)H(A)+E(A),                         (5)

where E counts ordered spanning two-path covers allowing either slot empty.
Indeed sum e_x=H and sum D_x=E-H. This scalar law loses the port allocation.

## 3. All mixed derivatives are observable by actual growth

Let f_n(a) be the specific multilinear polynomial

    sum over vertex permutations v of product_i p(v_i->v_(i+1)),

where a_ij=p(i->j) for i<j and p(j->i)=1-a_ij. Its degree is at most n-1;
at every Boolean a it equals H of that labelled tournament. This formula,
not an arbitrary interpolating polynomial, defines the continuous extension.

For an edge set E of size r, applying the operators L_e for its edges
recovers partial_E f_n(a). To justify iteration, every growth context in
(1) copies every old edge unchanged and its construction does not depend
on those edge values. Therefore differentiation in another old edge
commutes with that context. Induction on r proves the assertion. The
result is a signed sum of at most 4^r actual H queries, each adding at most
3r vertices. Repeated differentiation in one edge gives zero as required
by multilinearity. The original ports are never renamed or erased.

## 4. Why the full derivative response identifies the Boolean source

Let m=binom(n,2), n>=3. Fix any orientation of the complete graph to define
its n-by-m signed incidence matrix D, with -1 at the tail and +1 at the head.
At the fair point a=(1/2,...,1/2), the Hessian of f_n is exactly

    Hess f_n = c_n (2 I_m - D^t D),
    c_n = 2^(4-n) (n-2)!.                                  (6)

**Direct proof.** Diagonal entries vanish by multilinearity. For disjoint
edges, reversing either ordered two-vertex block cancels their signed
contributions among the path words. For two incident edges, only words
having their three vertices consecutive contribute. A directed chain and
its reversal give +2(n-2)! words; a fork gives the negative. The remaining
n-3 path factors each equal1/2. These are precisely the entries of (6).

Since DD^t=nI-J, the eigenvalues of (6) are

    2c_n          with multiplicity (n-1)(n-2)/2,
    (2-n)c_n      with multiplicity n-1.                      (7)

Neither is zero. This is a cut/cycle decomposition of a particular response
operator. It is not a claim that ordinary graph cycles encode every
arithmetic dependence, or that a nonsingular Hessian makes f_n surjective.

**Faithfulness theorem.** Suppose two labelled n-tournaments have equal H
after every growth context built as above, adding at most 3(n-2) vertices.
They are the same labelled tournament.

**Proof.** Their Boolean vectors a,b have equal mixed derivatives through
order n-2 by Section3. Derivatives of order n-1 are constant polynomials,
and higher ones vanish. Thus their entire Taylor polynomials agree:

    f_n(a+z)=f_n(b+z) for all real z.

Let d=b-a. A polynomial invariant under translation by d is invariant under
translation by td for every real t: it is invariant for all integer t, and
a one-variable polynomial with infinitely many equal values is constant.
Differentiating this identity implies Hess f_n(u) d=0 for every u. Evaluate
at the fair point and use (7) to get d=0. QED.

**Converse improves the even-order bound.** Reversing each path proves
f_n(1-a)=f_n(a). In centered variables a-1/2, all odd total degrees vanish.
Thus a degree bound is d=2 floor((n-1)/2), and the proof only needs
derivatives through order d-1. The attachment bound sharpens to 3(n-3)
when n is even and n>=4; it remains 3(n-2) when n is odd. In particular,
H together with its first derivatives is already proved faithful at
orders three and four. The finite checks at orders five through seven
go beyond that degree bound.

There is an explicit small-order decoder. For n=3 or 4 the centered
polynomial has degree two, so for the signed derivative vector g=grad f_n(a),

    a = (1/2)1 + (1/(2c_n))(I - D^t D/(n-2)) g.

Indeed DD^t=nI-J gives this inverse of the Hessian, and the centered
linear term vanishes by converse symmetry. This reconstructs a source
from its whole derivative vector; it does not invert the single scalar H.
All 72 Boolean sources at these two orders are independently decoded.

This supplies a finite decoder in principle: compare the measured derivative
responses with the finitely many Boolean candidates. No efficient decoding
bound is claimed. At n=2 the derivative polynomial is constant, but direct
one-port births distinguish the two orientations (counts1 and3 at a named
port); at n=1 there is only one orientation.

## 5. Consequence for the smallest closed carrier

Any deterministic observer retaining H under all these named attachments
must distinguish every labelled orientation at fixed n. Actual adjacency
is itself sufficient. Consequently its minimum fixed-order cardinality is

    2^binom(n,2)                      for all labelled tournaments,
    2^binom(n-1,2)                    for the fixed-Hamilton-path chart.

The fixed-chart statement holds for n>=3 by the path-preserving queries in
Section1; at n=1,2 the chart has one point. Adding native complement or legal
re-markings cannot shrink this lower bound. Symbolic storage, shared factors
and equivariant transport of command names can still reduce implementation
cost; they do not merge sources distinguishable by fixed named commands.

The previous [finite relocation result](../two_edge_transport_20261102/ESCAPES_AND_CARRIERS.md)
used a different language. This theorem does not prove its open all-order
legality-only conjecture. The restricted C,R carrier can remain small
because it omits these growth probes.

The incoming [depth-three separation theorem](../primegame_path_orbits_20261102/THEORY.md)
already disproves the extrapolation from depth-two legality tests: two
eight-vertex sources with H=171 and H=173 agree on all legality words of
length at most two, but a word of length three separates them. Ordered
strong-atom compositions amplify that obstruction. Eventual legality-only
faithfulness remains a separate open question.

## 6. Individual faithfulness does not reconstruct arbitrary source laws

There is a complementary all-order obstruction. Let the old n-vertex
orientations be uniform on even, respectively odd, parity of their
m=binom(n,2) edge bits. For n>=3 these are different positive probability
laws. For **every fixed source-oblivious attachment context**, their mean
Hamiltonian counts are equal, even if arbitrarily many new vertices are added.

Here source-oblivious means all new arcs are specified independently of the
old orientation, using named old ports; old arcs are preserved. All S,B_x
words qualify. A context chosen adaptively after querying the old source,
or conditioning on a legal re-marking, is outside this assertion.

**Proof.** In a Hamiltonian word of the enlarged tournament, the old-old
adjacencies form a linear forest on the old n vertices. At most n-1 such
arcs occur. Therefore its probability product has degree at most n-1 in
the old edge bits, regardless of the number of new vertices. Summing words
preserves the degree bound. Uniform parity laws agree on every monomial
in fewer than m bits, and n-1<m for n>=3. QED.

More generally, every shared-source product of r fixed context-count
queries has degree at most r(n-1) after Boolean reduction. The parity laws
therefore agree on all such joint moments whenever r(n-1)<m. A possible
separating moment must have order at least ceil(n/2); sufficiency at that
order is not asserted in general. Independent copies instead lose the
shared-source product and retain only products of the equal means.

The threshold is attained at n=3. Use the context B_1 B_0, with no initial
sink needed for this example. For old masks0,...,7 its actual H values are

    5,5,13,9,9,11,15,13.

The two parity laws both have mean10, but second moments113 and111,
variances13 and11, and probabilities1/4 and0 of attaining H15. All fixed
growth-context means nevertheless agree by the theorem. This is an
all-context obstruction, not just agreement on one tested list of means.

Thus a family of scalar functions can separate every individual Boolean
source while its linear span fails to separate probability laws. Recovering
an unknown orientation from joint exact responses is different from
recovering an ensemble from separately averaged responses. Nonlinear joint
tests or source-conditioned information are required for the latter.

### Shared copies repair the rational ensemble inverse

Allow a further explicit operation: take several copies of the SAME old
source T, apply specified growth contexts independently to those copies,
and concatenate the resulting tournaments as an ordinal sum. The count
of the result is the product of the context counts. This realizes every
shared-source moment by an actual tournament; it is not the same experiment
as drawing fresh independent old sources. Copying old arcs into new copies
also falls outside the source-oblivious hypothesis above.

For any finite universe X of labelled tournaments of a fixed order, with
the same named ports, the means of these shared-copy counts, together with
the constant1 query giving total source mass, determine every rational
source weight. The constant is the zero-copy convention, or equivalently
the count of a fixed one-vertex tournament. For probability laws that
mass is already1. To prove this, for each
distinct x,y in X choose a growth response phi_(x,y) separating them, as
supplied by faithfulness. The function

    delta_x(z)=product_(y!=x)
      (phi_(x,y)(z)-phi_(x,y)(y)) /
      (phi_(x,y)(x)-phi_(x,y)(y))                            (10)

is1 at x and0 at every other source. Expanding it expresses each source
weight as a rational linear combination of shared-copy mean counts, using
at most |X|-1 copies in any term. Consequently their rational span is the
entire |X|-dimensional function space. This is finite interpolation applied
to actual contexts, not an efficient inverse or a new general interpolation
theorem. Normalized laws have affine dimension |X|-1.

The coefficient field matters. By
[THM-001 / redei](../../../01-canon/theorems/THM-001-redei.md),
every actual tournament H is odd, and every
ordinal product of such counts is odd. Thus every query row reduces to the
all-one row modulo2; an H-only response matrix has rank at most1 over F2.
For |X|>1 it cannot give a unimodular inverse over all moduli. Division or
extra integer-normalized observations are required. This differs from the
incoming legality-conjunction zeta inverse, whose 0/1 query rows are
unimodular. No arbitrary integer denominator is silently inverted modulo2.

In the n3 example, two shared copies of B_1 B_0 T have H equal to the square
of the displayed count. Their mean counts113 and111 separate the parity
laws using actual ten-vertex tournaments. The construction still supplies
products, not arbitrary target primes: a nontrivial ordinal product cannot
close the prime realization lane.

## 7. Precise consequence for the spectrum problem

Every proposed inverse can now be tested on exact old-edge derivatives
using legal growing witnesses with preserved path ports. All-order memory
compression by forgetting labelled orientation is unavailable for this
attachment language. Neither fact supplies an orientation with prescribed H.

The remaining global obligation is still actual target entry for the prime
and7p lanes. The four-query identity is a signed numerical identity, not
four positive branches whose counts can be subtracted inside a tournament.
In particular it cannot manufacture the forbidden values7 and21.
