# Ten arithmetic vertices: prime-incidence cycles and exact rank compression

**Status: PROVED, independently audited by the root coordinator and geometry
lane; FINITE-EXACT controls in the explicitly stated universes.** The two arithmetic inputs are actual integer values and their
multiplicative relations. The finite graph retains unbounded integer labels;
it is not a reduction to a fixed residue modulus. Collatz and G2 remain open.

## 1. Inherited mechanisms and the question that survives

The closest proved mechanism is
[THM-4493, gcd and height rank](../../01-canon/theorems/THM-4493-collatz-gcd-height-rank-certificate.md).
It certifies independence from private primes or strict gcd dominance and
gives a pointwise no-dip source-height theorem. Its canonical hostile
`(6,10,15)` is independent but has neither private primes nor strict gcd
dominance. The corrected near miss is replacing valuation rank by positivity
of the log-gcd kernel: `(2,4)` defeats that replacement. The underused sidecar
is the *incidence and endpoint exponent* of a shared prime.

An earlier graph operation in
[THM-2521, signless-potential module](../../01-canon/theorems/THM-2521-k13-drift-k14-potential-module-bridge.md)
uses edge measurements `c_i+c_j`; here the exact prime equation supplies
`u c_i+v c_j`, with u,v positive integers. The equal-weight case is literally
the same signless incidence operation. This is a map between linear
operators, not a transfer of the older LRC conclusion. Another useful guardrail
is [THM-3387, q=2 gcd graph](../../01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md):
an intrinsic gcd relation can be symmetric and should not be made into a
tournament. Its sheet-blocking predicate is different from the rank predicate
here. The Sun reflection graphs in
[THM-4246](../../01-canon/theorems/THM-4246-sun-reflection-orbit-two-adic-tower-and-odd-cycle-firewall.md)
are bipartite; unlike those carriers, the arithmetic graph below genuinely
contains odd cycles on actual Collatz paths.

The concept board is: coprime atoms; prime/node incidence; weighted cycle
consistency; actual Collatz carries; and separate coordinates of a pair.
The anchor is a finite certificate for the exact rank target. The niche is an
odd-cycle certificate that improves THM-4493's simple tests. The wildcard is
whether the prime pair 223,233 imposes a ten-vertex structure; its exact local
arithmetic is recorded below, without identifying an unsupported map.

## 2. A finite exact carrier, constructed without prime factorization

Let `a_0,...,a_(N-1)>1` be positive integers. Begin with the set of their
distinct values. Whenever two distinct current factors b,c have
`g=gcd(b,c)>1`, replace b,c by `g,b/g,c/g`, discarding 1 and duplicates.

This process terminates. The product of the current distinct factors strictly
decreases at every split: before deduplication the replacement product is
`bc/g`, and coincidences with other factors can only decrease it further.
All original inputs remain products of nonnegative integral powers of the
current factors, by substitution. At termination the atoms `b_1,...,b_s>1`
are pairwise coprime, and

    a_i=product_h b_h^E_(h,i), E_(h,i)>=0.                 (1)

The exponents can be recovered by repeated exact division. The atoms need
not be prime or prime powers. Any prime p dividing b_h divides no other atom,
so its valuation row on the inputs is

    (v_p(a_i))_i = v_p(b_h) (E_(h,i))_i.                 (2)

Every atom has such a prime, hence the atom matrix E and the complete prime
valuation matrix have exactly the same rational and integral kernels. In
particular multiplicative rank equals `rank_Q(E)`. This is an exact
relation-preserving carrier, not a probabilistic factoring substitute.

Coprime bases and efficient factor refinement are established algorithms:
see Bernstein, *Factoring into coprimes in essentially linear time*,
[author-hosted primary paper](https://cr.yp.to/lineartime/dcba-20040404.pdf),
*Journal of Algorithms* 54 (2005), 1--30,
[publisher record](https://doi.org/10.1016/j.jalgor.2004.04.009).
The simple implementation here is justified by the decreasing-product proof;
it does not claim Bernstein's running-time bound or any factoring novelty.

## 3. Weighted graph reduction with an exact residual matrix

Each one-support row of E forces its coefficient c_i to vanish. Mark that
vertex as grounded. Each two-support row gives

    u c_i+v c_j=0, hence c_j=-(u/v)c_i, u,v>0.           (3)

Retain a labelled edge for each such row, including parallel edges when
different atoms share the same two endpoints. Traverse each connected
component, assigning a root weight 1 and transporting rational weights by
`-u/v`. A component is forced zero if it contains a grounded vertex or if a
cycle returns a weight different from the root weight. Otherwise its entire
solution space is one dimensional, `c_i=t_C w_i`, with nonzero rational
weights w_i fixed by the traversal. Isolated ungrounded vertices count as
one-dimensional components.

Let b be the number of these surviving components. For every original row
with at least three nonzero entries, substitute their component expressions:

    R_(h,C)=sum_(i in C) E_(h,i) w_i.                   (4)

Terms in forced-zero components contribute zero. Then the full kernel of E
is exactly the image of `ker R` under `t -> c`. Therefore

    rank(E)=N-b+rank(R).                                (5)

This is an equality, not merely a sufficient certificate. In particular no
common-prime row may be discarded. Its support may be large before graph
reduction but become decisive in R. All operations are rational and all
denominators may be cleared when an integral relation is required.

An odd cycle necessarily forces zero: its gain product has negative sign
and cannot be 1. An even cycle also forces zero unless its product of positive
endpoint exponent ratios is 1. Parallel edges give the same check at length
two. Thus graph parity alone is a sufficient special case; endpoint exponents
are necessary for the exact result.

The source is a labelled integer family, and the target is a graph with
positive endpoint weights plus a residual matrix. The map preserves its
entire rational relation kernel and therefore multiplicative independence.
The unweighted graph alone forgets exponent ratios; the graph without R
forgets higher-support primes. Neither quotient is claimed exact.

## 4. A factor-free graph-only certificate and iterative removal

One can recognize the useful unweighted graph using gcds alone. For a vertex
i, start from a_i and repeatedly divide out gcds with all other a_j until
each such gcd is 1. A nonunit remainder means some prime occurs only at i.
For a pair i,j, start from `gcd(a_i,a_j)` and strip all primes shared with
any remaining a_k, k distinct from i,j. A nonunit remainder means some prime
has support exactly `{i,j}`. Draw precisely those edges.

For any multiplicative relation, a private vertex has zero coefficient.
On a pair-exclusive edge its endpoint coefficients are both zero or have
opposite signs. Consequently every connected component with a private
vertex or an odd cycle is forced zero. If all components have one of these
features, the family is independent. Once coordinates are proved zero,
remove those vertices and recompute: formerly higher-support primes may
become private or pair-exclusive. Each removal remains valid because the
removed coefficients already vanish. This iterative test is still only
sufficient; equation (5) gives the exact completion.

This strictly improves a requirement that every node have a private prime.
It does not imply the numerical gcd-dominance condition from THM-4493.
Both are sound certificates with different failure boundaries.

## 5. Genuine ten-vertex arithmetic and Collatz examples

Assign a different prime to every edge of a graph, and let a_i be the product
of primes on incident edges. There are no private primes when the graph has
no isolated vertices. Its prime matrix is the signless edge-incidence matrix.

- For the ten-cycle C10, rank is 9. The product of the even-positioned values
  equals the product of the odd-positioned values.
- Square one endpoint's prime factor on one C10 edge. The unweighted graph
  is unchanged, but the even-cycle gain becomes unequal to 1, giving rank10.
- For the ten-vertex Petersen graph, the connected graph has an odd cycle,
  so rank is10, with no private prime. No Hamiltonian property is used.

These are constructive integer models, not claims that arbitrary graph
realizations form chronological Collatz paths. Actual paths also contain
the required configuration. On the odd orbit of27,

    91=7*13, 175=5^2*7, 325=5^2*13.

These three nodes have no private prime within their family. Each node is
exactly the product of its pairwise gcds with the other two, so strict gcd
dominance fails everywhere. Their exclusive-prime triangle nevertheless
forces all coefficients to vanish. Seven additional nodes from the same
orbit give the following actual ten-node selected family, written as
`(odd index,value)`:

    (1,41),(2,31),(3,47),(4,71),(5,107),(7,121),
    (8,91),(9,137),(13,175),(35,325).

The seven other vertices have private primes. The three ungrounded triangle
vertices are certified by the odd cycle. Applying the same certificate to
the first slots `3m_i` retains rank10: the common prime3 is irrelevant to
these particular non-3 certificates, but remains in the exact matrix.

There are also consecutive ten-node examples. Source199 gives

    199,299,449,337,253,95,143,215,323,485.

Here `299=13*23`, `253=11*23`, `143=11*13` form an ungrounded exclusive
triangle inside the actual ten-vertex family. In the finite census below,
source4347 is the first no-dip prefix with an odd component:

    4347,6521,4891,7337,5503,8255,12383,18575,27863,41795.

Every displayed node is at least4347. Its odd cycle has node indices0,3,6;
the edge primes are7,23,29. This cycle need not be essential to certification
because its component also contains private vertices.

**A minimal useful boundary.** Iterative unweighted removal leaves the pair
`(117,39)` unresolved. Its atom rows are `(2,1)` at3 and `(1,1)` at13, so two
parallel weighted edges force zero and rank is2. These are actual first slots
of odd nodes39 and13 in the ten-node prefix from39. Retaining only the single
unweighted edge loses this decisive information.

## 6. Keep both coordinates for G2

For actual pairs `(-3m_i,3m_i+1)`, construct one coprime basis across all
`2N` positive integers `3m_i` and `3m_i+1`. Keep the atom exponent rows of
the two coordinates **separate**, then stack them into one matrix. Its rank
is N if and only if these pairs are multiplicatively independent: full rank
rules out a relation; a rational kernel vector can be cleared to integers
and doubled to remove the first coordinate's possible minus sign. Thus the
finite torsion does not change this independence criterion.

One must not add corresponding coordinate rows or merely check the first
slot. For the injective odd nodes `(3,5,1)`, first slots `(9,15,3)` have rank2,
while the actual paired matrix has rank3. The ten-node experiment has96
first-slot rank failures but no paired-rank failure. This is FINITE-EXACT
evidence, not a universal G2 theorem.

The factor3 itself cannot be dropped without a sidecar: the raw positive odd
tuple `(5,75)` has rank2, while its scaled first slots `(15,225)` have rank1
because `225=15^2`. This is a typing hostile, not a claim that5 and75 are
consecutive Collatz nodes. The implementation retains3 in every first slot.

A stronger hostile is realized on a single actual orbit:

    507 -> 761 -> 571 -> 857 -> 643 -> 965 -> 181 -> 17 -> 13.

The selected raw odd nodes `507=3*13^2` and13 have rank2, but their first
slots are `1521=39^2` and39, of rank1. The script directly replays all eight
odd transitions, checks the arithmetic identity and both ranks, and confirms
that the two actual *pairs* still have rank2. Thus this is an actual failure
of raw-node-to-first-slot rank transfer, not a G2 counterexample.

For Collatz, every prime common to m_i,m_j divides the exact subword carry
`C_(i,j)` from THM-4493. Thus the graph edges also have an arithmetic
realizability constraint; arbitrary edge-prime assignments do not inherit a
Collatz realization. Conversely, this carrier uses the actual unbounded
integer labels and exact gcds, so it does not fall under the complete finite
residue-lift obstruction of the earlier223 session.

## 7. What 223 and 233 do, and do not, share

Exact computations give:

| p | ord_p(2) | ord_p(3) | log base3 of2 | index of <2> |
|---|---:|---:|---:|---:|
|223|37|222|180|6|
|233|29|232|72|8|

The identities `2^37-1=223*616318177` and
`2^29-1=233*1103*2089` give genuine residue clocks. At both primes the
quadratic characters satisfy `chi(2)=+1, chi(3)=-1`. But `-1` is nonsquare
mod223 and square mod233. Therefore quadratic-residue differences give the
Paley tournament at223 and an undirected Paley graph at233. These are
different native objects. The difference233-223=10 provides no map from
either residue object into the ten-node valuation carrier above. The latter
is useful because its edges encode actual equations, independently of the
numerical gap.

## 8. Exact universe, outcomes, and stopping point

Run the [standalone script](../../04-computation/experiments/crossroads10_20260926_arithmetic.py)
normally and with `-O`. Its [output](crossroads10_20260926_arithmetic.out)
is exact stdout. All checks use integers or Fractions, not floating point.

The complete small universe consists of all8835 multisets of lengths2,3,4
with entries2..20. All agree between direct prime-valuation elimination,
coprime-atom construction, and exact graph compression. Refinement uses15488
gcd splits. The unweighted graph certificate succeeds for2578 families.

The chronological universe is every odd source3..20000 whose first ten odd
nodes are distinct, giving9446 prefixes. Of these:

| Test or property | Count |
|---|---:|
| Every first-slot node has a private prime |6487|
| One pass of the graph-only certificate succeeds |8481|
| Iterative graph-only removal certifies full rank |9061|
| Exact compressed first-slot rank is10 |9350|
| Exact compressed paired rank is10 |9446|
| Graph has an odd component |45|
| All ten odd nodes stay at least their source |793|
| Both preceding properties |2|

The first-slot rank histogram is `{7:1,8:4,9:91,10:9350}`. Independent
dense prime-row elimination audits both coordinates separately on all retained
sources<=2001 and every deficient first-slot case,885 cases total. The synthetic ten-node models,
the selected orbit27 family, and the counterexamples are checked separately.

The decisive new restriction is that a surviving relation must live in the
ungrounded, cycle-consistent graph components and satisfy every compressed
higher-support row. The exact remaining target is injectivity of that
compressed matrix for all relevant chronological paired families. A uniform
proof is not supplied by the finite examples. The graph changes when another
node is appended: a pair-exclusive prime may become a higher-support prime.
Consequently a ten-node certificate cannot simply be repeated to certify an
entire infinite orbit, and full rank alone still does not force descent.
