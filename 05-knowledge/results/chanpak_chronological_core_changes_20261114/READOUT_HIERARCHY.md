# Chronological paths: infinite hierarchy and a 7-adic two-path strip

Let H be the number of directed Hamilton paths of a tournament. Its
chronological budget κ is the least number of directed Hamilton paths
whose union meets every directed cycle. Equivalently there is one rank
of the vertices such that the paths cover every backward arc; order
the acyclic graph remaining after deletion. The
[proof for 63](../chanpak_collatz_exact_addresses_20261114/READOUT_63.md)
introduced this budget and proved that it composes by maximum over
strong components. The results below are about the κ of particular
tournaments, not generally the minimum κ over all realizations of the
same H.

## The budget is finite, and can grow linearly

For a strong tournament, Busch–Jacobson–Reid's arc-traceability
structure theorem implies that if x→y lies on no Hamilton path,
some vertex z separates them, with x in the initial SCC and y in
the terminal SCC of T−z, and at least four SCCs in all. See the
[authors' paper](https://doi.org/10.1002/jgt.20173) and the explicit
statement of Theorem 4.12 in [Busch's thesis](https://digital.auraria.edu/files/pdf?fileid=de2a8eac-a8fb-44fa-b93e-e60575695962).
In T−z every out-neighbor of y is in its terminal SCC and is
therefore beaten by x; x also beats y and at least two vertices in
middle SCCs. Adding z can change the out-degree difference by at
most one. Hence every untraceable arc goes strictly down out-degree.
Order by decreasing out-degree, breaking ties arbitrarily, and
select one Hamilton path per backward arc. All such arcs are traceable.
In particular κ≤binomial(n,2) in a strong n-tournament. The SCC law
extends finiteness to all tournaments. This does not assert that the
rank given here attains κ.

The first vertex in any residual rank originally has at least
minimum in-degree incoming arcs, and all must be deleted. A path can
supply at most one incoming arc at this vertex. The same argument
at the last vertex uses minimum out-degree. This yields a lower
bound by the greater of the two degree minima (useful within each SCC).

* **Explicit variable-size hierarchy.** For transitive blocks in a
  directed triangle C3[TT_r,TT_b,TT_r] of sizes r,b,r with b≥r≥1,
  κ=r. Put the blocks A,B,C in that rank; the only backward arcs
  are C→A. Partition these into r perfect matchings by cyclic
  shifts. Fix r nonempty consecutive increasing runs of B. The
  path using successive groups (c_i,a_(i+s),B_i), with indices
  modulo r, covers the entire shift matching s. The lower bound
  is min indegree r, attained by the first vertex of B or C.
* **Regular tournaments with a Hamilton decomposition.** If a
  regular n=2r+1 vertex tournament decomposes into r directed
  Hamilton cycles, then κ=r. For any rank, delete one forward
  edge from each cycle, which always exists, to get paths
  covering all backward edges. The degree bound gives equality.
  Thus [Kühn–Osthus's large-order Kelly theorem](https://arxiv.org/abs/1202.6219)
  gives κ=(n−1)/2 for every sufficiently large regular tournament.
  For any circulant prime-order tournament no deep theorem is
  needed: the chosen steps are already disjoint Hamilton cycles.
  Their translation orbits on Hamilton paths also show q|H at
  prime order q. These explicit κ results concern the family,
  not an assertion that a count cannot occur at a lower budget.

## An exact two-path construction parametrized by a smaller core

Write P_D(j) for the number of ordered j-tuples of disjoint directed
paths whose vertices partition D, allowing empty paths. In particular
P_D(1)=H(D). For an ordinal sum R→TT_m, restriction to each summand
and gluing in each corresponding slot are inverse, so

    P_(R→TT_m)(j)=P_R(j) j^m.

Indeed every one of the m transitive vertices can be assigned to any
slot: once its slot is chosen its order inside that slot is forced.

For a directed three-block tournament C3[TT_2,B,TT_2], color a Hamilton
path by block. After contracting each maximal monochromatic run, the
color word advances deterministically A→B→C→A. Each A or C has either
one or two runs; one run gives one internal assignment, while two
runs give two assignments. The table records the total weights of
the admissible start-color and run-length pairs, grouped by the
number of B runs:

| Total runs | one B run | two B runs | three B runs |
|---|---:|---:|---:|
| 3 | 3 | 0 | 0 |
| 4 | 4 | 1 | 0 |
| 5 | 4 | 4 | 0 |
| 6 | 0 | 12 | 0 |
| 7 | 0 | 0 | 4 |

The B runs themselves must be a cover of B with nonempty ordered
slots. Their counts for one, two, and three slots are respectively
P_B(1), P_B(2)−2P_B(1), and P_B(3)−3P_B(2)+3P_B(1), by
inclusion–exclusion for empty slots (B is nonempty). Hence exactly

    H(C3[TT_2,B,TT_2]) = 4P_B(3)+5P_B(2)−11H(B).

Choose any nonempty core R with κ(R)≤2, set B=R→TT_m with m≥1,
and use the two-shift construction above. In each matching path,
make the first B run a selected Hamilton path of R and the second
the increasing TT_m; if only one R path is needed, it may be reused.
The two paths cover the C→A arcs and all backward arcs of R for
a common rank. The minimum in-degree is at least two (each block
has size at least two). Thus these are *strong*, exact κ=2
tournaments with path count

    N_R(m)=4P_R(3)·3^m+5P_R(2)·2^m−11H(R).

This is a physically closed move for κ≤2, even though count recursion
does not close on H alone. For example a further iteration needs
higher path-cover statistics of the new core. If R is a singleton,
one recovers 4·3^b+5·2^b−11 for B=TT_b, b≥2, including primes
137, 9377, and 81281 at b=3,7,9. This transitive-middle strip never
has 7 as a divisor, by the periods of 2 and 3 modulo 7.

## A concrete core that unlocks the 7-adic clock

Use the order-five core with outgoing adjacency bitmasks
(18,8,3,5,14). A single path (0,4,1,3,2) covers all arcs backward
in the rank (4,2,3,0,1). It has (H,P_R(2),P_R(3))=(13,104,471)
by the direct cover replay below. Its κ=2 extensions thus have

    N(m)=1884·3^m+520·2^m−143.

For m=2 this is exactly 7·2699, with 2699 prime by trial division
to its square root. Modulo 7 the six values for m=0,...,5 are
(0,4,0,5,5,3), and repeat. The two lifts of the zero modulo 7
have genuinely different behavior:

* **m=6u+2:** N(m)/7≡4 (mod 7) for every u≥0, so the
  valuation is exactly one.
* **m=6u:** N(6u)/7 is 7-adically bijective, in the following
  finite sense. For every s≥0 its values induce a permutation
  modulo 7^s. Thus there are infinitely many valid physical m
  with each prescribed positive exact 7-adic valuation.

For the proof, 3^6=1+7·104 and 2^6=1+7·9. For m=6u+r,
writing A=1884·3^r, B=520·2^r, the binomial theorem gives

    [N(6(u+7^s t)+r)−N(6u+r)] / 7^(s+1)
        ≡ (6A+2B)t (mod 7)

for nonnegative u,t and s≥0. This uses the elementary fact that
(1+7c)^(7^s t)≡1+7^(s+1)ct (mod 7^(s+2)). For r=2
the coefficient is zero mod 7, and N(2)/7=2699≡4.
For r=0 the coefficient is 3 mod 7, so each choice of next
u digit yields each next output digit exactly once. The argument
also gives explicit nested addresses for exact valuation rather
than merely an existence lift. The exponent m=0 is used for
congruence arithmetic, not claimed as a κ=2 witness; any needed
class of u has arbitrarily large positive representatives.

## Even a bank of all one-apex sources is too small

Allow R to range over **every** one-apex tournament, meaning there is
a vertex whose deletion leaves a transitive tournament (irrespective
of whether this R itself satisfies κ≤2). This is a larger strip
family than the defined κ=2 one, so bounding its size suffices.
Every assignment of core vertices to three slots admits a path in
each occupied slot, by Rédei, so P_R(3)≥3^d, where d=|R|.
Also H(R)≤P_R(3) by using only a single fixed slot. Therefore
for m≥1,

    N_R(m) ≥ (4−11/3^m)3^(d+m) ≥ 3^(d+m)/3.

A marked one-apex core of order d is specified by the binary win
mask of the apex on the uniquely ordered transitive deletion, so
there are at most 2^(d−1) possible unmarked cores. Write
L=floor(log_3(3X)). Even counting duplicates, all outputs at most
X have at most

    sum_(1≤d<L) (L−d) 2^(d−1) = O(2^L) = O(X^(log_3 2))

presentations. Thus arbitrary changes to a one-apex source in
this two-path construction still give density zero in the odds; a
more general core geometry is necessary for a density conclusion.

Even arbitrary ordinal sums of these source blocks are too small.
The three P_R values multiply component by component under ordinal
sum, so order of the blocks has no effect on the output. With
a_d≤2^(d−1) block types of order d, count possible multisets by
F(z)=prod_(d≥1)(1−z^d)^(-a_d). For 0<z<1/2,

    log F(z) ≤ sum_d a_d z^d/(1−z^d) = O(1/(1−2z)).

The choice of m≥1 contributes the bounded factor z/(1−z).
For total order d+m≤L the number of presentations is at most
z^(-L) F(z) z/(1−z), because its series has nonnegative
coefficients. Take z=(1−1/sqrt(L))/2 for L>1. Thus allowing any
number of one-apex blocks, arbitrary orders, widths, and masks
still gives at most

    X^(log_3 2) exp(O(sqrt(log X)))

distinct positive counts up to X. This conclusion counts sorted
signatures, not ordered physical tournaments: counting all ordered
lists instead would bring in the misleading growth base 3. More
general source libraries with exponential type growth 1<λ<3 have
the same power-sparsity mechanism, replacing 2 by λ (or by a
slightly larger λ still below 3 if using only an exponential
upper bound). The lower bound on the arithmetic size still
comes from three path-cover slots rather than from a particular
tournament's Hamilton count.

**Limit:** this strip gives an infinite family of two-path
tournaments with independently controllable 7-adic precision.
It supplies specific 7p counts, not all such counts or infinitely
many prime quotients. It does not prove Chan–Pak. A fixed R still
gives only logarithmically many exact values up to a height; the
problem is to choose R so the higher covers in this exact
equation select the absolute target, not only its 7-adic address.

Run `python3 replay_readout_hierarchy.py` for direct count controls
and verification of the concrete rank, certificates, and 7-adic lifts.
