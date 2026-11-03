# Rooted inventories between the tiling cube and A000568

**PROVED elementary identities, with finite controls.** Ordinary tournament
classes, rooted classes and Hamilton-path presentations are different
objects. This note gives the exact maps and the extra inventory required
by a recurrence; it does not claim full metagraph classification.

## 1 The correct fibre weights

Fix the directed Hamilton spine `0->1->...->n-1`. Its n-1 arcs are fixed;
all remaining arcs are independent bits. Hence there are

    2^((n-1)(n-2)/2)

fixed-spine charts. For an ordinary tournament class T, the automorphism
group acts freely on its ordered Hamilton paths: an automorphism fixing a
path fixes its vertices one by one. Thus its fibre over the fixed-spine
chart model has exactly `H(T)/|Aut(T)|` elements, and

    sum_[T] H(T)/|Aut(T)| = 2^((n-1)(n-2)/2).           (18)

For a rooted class (T,r), with r at the first position, the fibre instead
has `s_T(r)/|Aut(T,r)|` elements. Summing these over vertex-root orbits
recovers(18). This is the precise connection to the user's fixed-path
tiling model. Counting charts once does not count A000568 once per class.

At start budgets1 and2 the rooted stabilizer is trivial, so the chart
counts are exactly one and two times the class counts. With D=1-3x+x^2,
the possibly weak rooted strata have series

    A_1(x)=x(1-2x)/D,
    A_2(x)=x^4(1-x^2)/D^2.

The corresponding first-position tiling strata have series A_1 and2A_2.
For strong charts they have series U and2V. At budget3 the strong wrapper
`C3->a->b->C3`, rooted at b, already has an order-three root stabilizer.
Therefore the no-weight shortcut is
special to budgets at most two. This is a reason to enumerate the next
layer with cycle indices rather than extrapolate the Fibonacci scalar law.

## 2 A universal entry-subset inventory

For a strong core C, including the singleton, set h_C=H(C), G=Aut(C), and
s_v=s_C(v). Define

    B_C(z) = (1/|G|) sum_(g in G)
              product_(O a vertex cycle of g) (1+z^(sum_(v in O) s_v)). (19)

Burnside's lemma says this counts **orbits of entry subsets J** by their
total starting response `sum_(v in J) s_v`. The empty set contributes1,
the whole set contributes `z^h_C`, and `b_C=B_C(1)` is the number of
vertex-subset orbits. An orbit with two entries is joint data; two
independently averaged vertex orbits do not determine it.

The definition is exact for all cores and all budgets. To preserve later
operations, an enumerator must also retain the actual subset representative,
its stabilizer and ordered physical port transport. Formula(19) by itself
is only a class-count inventory, not the full attachment carrier.

## 3 Exact root-deletion identity with multiplicative count clocks

Delete the distinguished root of a strong tournament. Its remaining
components are an ordered sequence of strong cores. Let J_i be the subset
of component i beaten by the root. The exact strongness guards are

* one component: J_1 is nonempty and proper;
* at least two components: J_1 is nonempty, J_last is proper;
  intermediate J_i are arbitrary.

Necessity follows because no vertex outside the first component can lead
into it, and no vertex in the last component can leave it except through
the root. Sufficiency follows by strongness within each component and
the directed chain of components. Independent component automorphisms
act on their own entry subsets. This is a bijection on rooted classes.

For a series f(z), put `D_h f(z)=f(z^h)`, and define formal x-adic operators

    I(x,z) = sum_[C] x^|C| (B_C(z)-1),
    Z(x,z) = sum_[C] x^|C| (B_C(z)-1-z^h_C),
    L = sum_[C] x^|C| b_C D_(h_C),
    R = sum_[C] x^|C| (b_C-1) D_(h_C).

The ordered SCC path factorization in TWO_ENTRY(1) gives the universal
strong rooted start-budget inventory

    V_all(x,z) = x [ Z(x,z) + R (1-L)^(-1) I(x,z) ].   (20)

Here `(1-L)^(-1)=sum_(j>=0)L^j` is well-defined by positive x-degree.
Each later component multiplies the number of starting paths by h_C,
which is exactly the dilation D_h. Intermediate components have b_C
subset choices; the terminal component has b_C-1. This proves(20).

Projecting to z^1 or z^2 removes all nontrivial later components because
their h_C is at least three. The remaining singleton operator is just
the binary tail clock. Moreover `[z]I=x+U` and `[z^2]I=V+P`.
Equation(20) therefore gives exactly the U and V equations proved in
TWO_ENTRY. Higher coefficients bring in strong cores, subset stabilizers
and multiplicative factors at their actual first possible budgets.
This is a systematic extension of the funnel strata, rather than an
ad hoc list of small examples.

## 4 Ordinary and rooted class recurrences

Let `C(x)=sum_[strong C] x^|C|`, including the singleton, and
`B(x)=sum_[strong C] b_C x^|C|`. The unique SCC decomposition gives the
classical ordinary tournament identity

    A(x)=1/(1-C(x)),                                  (21)

whose coefficients are A000568, including the empty tournament. If an
arbitrary vertex is rooted, delete it; every component now has an arbitrary
entry subset and no strongness guard. Thus the ordinary **rooted** class
series is

    R_all(x)=x/(1-B(x)).                              (22)

This yields the convolution recurrence
`r_(n+1)=sum_(j=1)^n b_j r_(n-j+1)`, with r_1=1, where b_j is the sum of
subset-orbit counts over strong cores of order j. The elementary structure
behind(21)-(22) is SCC decomposition and Burnside pointing. The weighted
start-budget refinement(20) identifies exactly how the two-entry grammar
sits inside it.

Finite independent class canonicalization gives ordinary counts
`1,1,1,2,4,12` through order5, rooted counts `1,2,4,12,48`, and(22) gives
296 at order6 from the complete strong inventory through5. The strong
rooted totals from(20) are `1,4,24,194` at orders3 through6. The first
three are checked against direct root orbits; the order6 result is a
recurrence evaluation, not a new A000568 term.

Equations(20)-(22) are specifications with unresolved input inventories.
They do not close A000568 on finitely many scalar registers. In particular
the full vertex cycle index and the core library do not follow from H,
from a single rooted count, or from the first two strata. The inherited
[THM-4845, Gaussian/cycle-index A000568
formula](../../../01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md)
and its non-P-recursive boundary remain the broader counting context.

## 5 What a full metagraph recurrence must additionally enumerate

A fixed-spine chart retains a path, its orientation gauge and every
off-spine bit. The literal native complement reverses those off-spine
arcs. After quotienting charts to ordinary classes, the chosen path orbit
and its endpoint transport determine the actual target. Therefore a full
recursive metagraph inventory needs:

* source class plus vertex and entry-subset stabilizers;
* joint ordered path components, multiplicities and selected path orbit;
* the full complemented target relation and its normalization transport;
* integral realization guards and the phase/size clocks used by the next move.

The new two-entry grammar enumerates a complete low-budget part of the
source inventory; (19)-(20) specify its next layers. TRANSPORT gives actual
complement programs on a useful infinite subfamily. Neither result removes
the target-reconstruction obligation. This is compatible with the inherited
[THM-4807, legal-path reconstruction and endpoint
decoder](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md):
the full legal path language is much richer than a class-count coefficient.

For the twelve-vertex hostile pair, |Aut|=81 or27 implies every vertex
stabilizer has size at least |Aut|/12>2. Its free action on starting paths
excludes all one/two-start roots. Equality of labelled cyclic triples
therefore cannot turn this restricted grammar into a general classifier.
Their inherited two-run ranks24/22 and three-run ranks430/435 modulo3
remain meaningful joint controls, with their original universes and scope.
