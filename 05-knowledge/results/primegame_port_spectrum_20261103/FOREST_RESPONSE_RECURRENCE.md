# A closed carrier for mean responses: retain path blocks and their endpoints

**PROVED all-order upper bound and closure / FINITE-EXACT sharpness at
orders three through five. Same-agent proof audit and independent rank algorithms.**

**Subsequent refinement:** [PARITY_CARRIER](../positive_source_jacobian_prime_clocks_20261103/PARITY_CARRIER.md)
shows that internal allocation across path blocks can also be forgotten
when endpoint pairing and individual path parities are retained. Its closed
dimensions571/2668 improve661/3928 at orders6/7. The theorem below remains
valid as an upper bound; its all-order sharpness candidate is refuted.

For a law on labelled tournaments, every fixed source-oblivious attachment
response is determined by a vector indexed by set partitions with two
designated endpoints in each nonsingleton block. The dimension bound has
the exact exponential generating function

    B(x) = exp(x + x^2 exp(x)/2).

At orders three through five the bound is attained: the minimum rational
linear dimensions are 7, 28 and 126. These are dimensions of averaged observations,
not cardinalities of deterministic states. Individual named sources still
require 8, 64 and 1024 states by GROWTH_FAITHFULNESS.md. The universe here is all
labelled orientations; a fixed-Hamilton-path chart is a different restriction.

## 1. The operations and coefficient field

An attachment retains every old edge and fixes all new arcs independently
of the old orientation. It retains the old vertex names and records the
birth order. This includes S and B_x, and also arbitrary fixed cuts for a
new vertex. There is no adaptive choice, conditioning on legality, or
copying the unknown old edge bits into new arcs in this operation class.

The state of a source law will be its rational expectations of specified
polynomials. All statements about minimum dimension are over Q, with linear
updates and readouts. They do not assert a finite-state cardinality or a
minimum over arbitrary nonlinear encodings.

Let a_ij be the old Boolean orientation bit for i<j and put

    y_ij = a_ij - 1/2,       y_ji = -y_ij.

Thus each directed old-arc probability in the weighted Hamiltonian
polynomial is 1/2+y_ij. No relaxed orientation is used as a constructor.

## 2. The first bound: only linear forests occur

In any Hamiltonian word, its old-old adjacencies form a linear forest:
an undirected acyclic graph of maximum degree two. Expansion in the old
y variables only deletes edges, so every nonzero monomial also has linear-
forest support. This holds regardless of how many external vertices exist.

Let l_n count labelled linear forests on n vertices, including isolates.
A nonsingleton component on k labels has k!/2 path orders. Decomposing at
the component of one distinguished vertex proves

    l_0=1,
    l_n=l_(n-1)+sum_(k=2)^n binom(n-1,k-1)(k!/2) l_(n-k),
    sum_n l_n x^n/n! = exp(x + x^2/(2(1-x))).

Hence all separate mean attachment queries span at most l_n dimensions.
This is an upper bound on observations, not tournament enumeration.

It strengthens the parity obstruction. If an old edge set F contains a
cycle or a vertex of degree at least three, fix all other old bits and
compare uniform even/odd parity on the bits of F. Every fixed attachment
has the same mean under these two laws. In particular both a triangle
and a three-edge star are invisible parity directions. The star explains
why acyclicity alone is insufficient.

For r shared-source products, every monomial is supported on a union of
r linear forests. Therefore parity on F remains invisible whenever F
cannot be partitioned into r linear forests. This is the linear-arboricity
constraint, used here only as a necessary condition. For all edges of K_n,
the bound of n-1 edges per forest gives r>=ceil(n/2), as in the earlier
degree proof. Sufficiency of that moment order is not claimed.

## 3. A stronger bound: internal orders combine with their signs

For a vertex set I of size at least two and endpoints u<v in I, define

    P_(I;u,v)(y) = sum over all orders w of I\{u,v}
                      y_(u,w1) ... y_(w_last,v).

For |I|=2 this is simply y_uv. A partition pi of the old vertex set,
with an unordered endpoint pair in every nonsingleton block, defines

    Q_pi(y) = product over nonsingleton blocks (I;u,v) of P_(I;u,v)(y).

Singletons contribute 1. Every linear-forest monomial belongs to exactly
one such Q_pi. Its coefficient there is +1 or -1, according to traversal
from the smaller designated endpoint to the larger one in each component.
Different Q_pi have disjoint monomial supports, so they are independent.
Uniqueness of multilinear interpolation makes them independent as functions
on the full Boolean orientation cube as well.

**Endpoint-block theorem.** Every fixed attachment count polynomial lies
in the span of the Q_pi, at every order and every attachment size.

**Proof.** The nonforest coefficients vanish by Section 2. Fix a forest F.
A word contributing to its coefficient must use every edge of F, so each
component of F occurs consecutively, traversed in one of its two directions.
Replace its internal order by any other path on the same vertex set with
the same endpoints. This is a bijection on contributing words. Arcs between
that block and the rest of the word are unchanged, as are all fixed new
arcs and all factors 1/2 contributed by unselected old edges.

The only change in the coefficient is the product of directed signs along
the replaced block. Reversing a block changes both candidate signs by the
same factor (-1)^(|I|-1), so their ratio does not depend on which direction
the word traverses. Thus, after multiplying by the canonical traversal
sign, the coefficient is constant over each Q_pi support. This is precisely
membership in their span. QED.

The signs are essential. This is not permission to replace an actual
tournament path by an arbitrary internal order: the proof compares formal
coefficients at the fair old-edge point. The external arcs stay fixed.

## 4. Closed updates and the exact carrier recurrence

The same replacement proof applies to the pullback of Q_pi on an enlarged
vertex set. Its terms are products along disjoint spanning paths with fixed
endpoints. After substituting the fixed new arcs, their old supports are
linear forests. Reordering an old component keeps its vertices in the same
enlarged block and preserves every attachment to its endpoints. Therefore
the pulled-back polynomial lies in the old Q span as well.

Consequently the expectation vector (E[Q_pi]) is closed under every fixed
source-oblivious attachment. It is also closed under a fixed relabelling
and ordinary converse. A source-dependent legal re-marking, native fixed-
path complement, or shared-source copying needs its own pullback check;
none is silently included. These are rational states and updates, not
nonnegative probabilities or a positive target-construction rule.

Let b_n be their number. A singleton component has one decoration; a
k-element component has binom(k,2) endpoint choices. Hence

    b_0=1,
    b_n=b_(n-1)+sum_(k=2)^n binom(n-1,k-1)binom(k,2)b_(n-k),
    sum_n b_n x^n/n! = exp(x + x^2 exp(x)/2).

The initial sequence is

    1,1,2,7,28,126,661,3928,25628,181828,1396261,... .

This is an exact count of the defined carrier coordinates, not yet an
all-order theorem that every coordinate is observable from H alone.
The upper bound is uniform and the update closure is proved. Minimality
must additionally establish the span of actual response queries.

## 5. Exact minimum dimensions at small orders

For each source, evaluate all words of S and the original named B ports
through the stated length. Fraction-free integer row elimination gives
the rational rank; reductions modulo 101 and 103 independently agree.

| Old vertices | Allowed word lengths | Number of queries | Rational rank | All-context upper bound |
| ---: | ---: | ---: | ---: | ---: |
| 3 | 0 through 2 | 21 | 7 | 7 |
| 4 | 0 through 3 | 156 | 28 | 28 |
| 5 | 0 through 3, plus 70 specified length-four words | 329 | 126 | 126 |

Since these actual queries attain the proved upper bounds, no larger fixed
attachment can enlarge the mean-response span at these orders. The minimum
linear dimensions are therefore exactly 7, 28 and 126 for unnormalized laws.
For normalized laws the constant coordinate is known, leaving 6, 27 and 125
independent measured coordinates. There are respectively 1, 36 and 898
dimensions of source-law ambiguity on 8, 64 and 1024 orientations.

For four old vertices there are 34 linear-forest monomials. The six extra
relations arise because each pair of endpoints has two spanning path
orders but only one signed coefficient. For endpoints 0,3, the paths
0-1-2-3 and 0-2-1-3 have opposite signs. Their edge masks are 41 and 26
in lexicographic pair-bit order. Writing chi_F(a)=(-1)^(sum_(e in F) a_e),
the two strictly positive laws

    mu_+(a) = (1+(chi_41(a)+chi_26(a))/4)/64,
    mu_-(a) = (1-(chi_41(a)+chi_26(a))/4)/64

have identical means after every fixed attachment. All weights are
1/128, 1/64 or 3/128. This ambiguity survives even after accounting for
the missing cyclic and branching monomials. Joint source moments can
detect it: a single B_1 birth has mean 15/2 under both laws, but second
moments 285/4 and 279/4. The receipt reconstructs every weight and response.

## 6. What this changes for the spectrum problem

We now have different exact objects for different tasks: a labelled source
for individual future behavior; a closed endpoint-block moment vector for
separate mean responses; and shared-copy products for reconstructing an
arbitrary source law. Their dimensions and allowed operations differ.

The recurrence supplies no source law concentrating at a requested H, and
its signed path polynomials are not tournament constructors. In particular
neither b_n nor a nonsingular response matrix proves the required prime
and 7p realization lanes. The numerical inverse still needs actual source
membership and compatible attachment ports.

## 7. Joint moments and the next positive certificate

Write q_n(T) for the column of endpoint-block coordinates. For a fixed
attachment C, the closure theorem gives a rational matrix M_C with
q_(n+k)(C(T))=M_C q_n(T). If an actual source probability law has
m=E[q_n] and R=E[q_n q_n^t], its exact updates are

    m' = M_C m,       R' = M_C R M_C^t.

If ell_n is the H readout, then for any requested odd h its error is

    E[(H-h)^2] = ell_n R ell_n^t - 2h ell_n m + h^2.

Thus the new closed carrier also supplies exact recurrences for weighted
joint-port control and variance, by retaining the second-moment matrix.
Combining it with the inherited
[Boolean energy criterion](../metagraph_primegame_20261102/INCOMING.md),
an actual positive source law with this error below 4 certifies a target
tournament, with target mass at least 1-error/4. This follows because two
different odd H values have squared separation at least 4.

The missing proof obligation is now precise: construct such actual laws
uniformly on the prime and 7p lanes and prove their error below 4.
Positive semidefiniteness of R alone does not ensure realizability by
Boolean sources. Nor does this fixed-operation recurrence authorize
conditioning a source law using only its first two moments. Those extra
operations require their own retained information. This corollary combines
the new closure with the inherited criterion; it does not claim that the
criterion itself is new or that the uniform bound has been achieved.

## Reproduction and status

The inherited mechanism is the full named-path transfer in
[Session32 / ports/PORT_TRANSFER](../session32_square_monster_20261023/ports/PORT_TRANSFER.md).
Its later [cyclic-product correction](../square_sum_spectrum_overnight_20261031l/README.md)
requires all linear-forest boundary states: restoring cut ports can produce
several local cycles. We retain all forest types here. The new compression
is of formal signed coefficients across internal path orders with the same
vertex set and endpoints; it is not restriction to cuts of one local cycle,
nor contraction of a branched graph to one path vertex. General boundary
dynamic programming is inherited, not claimed as a new method.

    python -B 05-knowledge/results/primegame_port_spectrum_20261103/forest_observations.py
    python -B -O 05-knowledge/results/primegame_port_spectrum_20261103/forest_observations.py --output forest_observations_optimized.json
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/forest_observations.py --five --output forest_observations_five.json
    python -B 05-knowledge/results/primegame_port_spectrum_20261103/audit_forest_rows.py

All edge sets through six vertices independently check both component
recurrences and the partition of forests into endpoint blocks. All tested
response rows check every forbidden Fourier coefficient and every endpoint
sign relation. The rank audit uses original Boolean evaluation rows at
orders three and four, avoiding a circular rank reduction. The optional
five-vertex probe uses the fully checked Fourier relations to compress
rows before elimination; its depth-four stopping rule is recorded explicitly.
An independent C++ starting-endpoint counter reconstructs all 347,048 actual
graphs in the three sharp query banks and obtains ranks 7, 28 and 126
directly from full Boolean evaluation rows over both prime fields. This
independent lower bound does not use the endpoint compression. All 2240
pullbacks of basis polynomials under arbitrary cuts from order 3 to 4 and
4 to 5 separately verify update closure.
No literature-priority claim or independent-agent review is made.
