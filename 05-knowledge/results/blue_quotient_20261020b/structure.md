# Black closed paths and cycle-cut mobility

Status: PROVED for Sections 2–4; FINITE-EXACT and OPEN as marked below.
Scope: the converse-merged **complement-line** metagraph, n>=3.
The operation below reverses all arcs of a directed Hamiltonian path.
It gives the same merged line as complementing every nonpath arc,
by [THM-4729, Proposition K](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
It is not the single-arc-flip metagraph.

## 1. Inheritance and concept board

The closest proved mechanisms are
[THM-4729, closed/open line incidence](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md);
[THM-4731, black sigma-twins](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md);
the explicit odd-circulant prefix obstruction
in `../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`;
and the black-component bridge quotient in
`../boundary_support_20261020_metagraph.md`.
The prefix proof repairs the older affine-only anti-automorphism gap.
The companion `../boundary_support_20261020_incoming.md` supplies
the rigid-strong deduction and exact pure-blue leaf/branch inventory.

The new distance-moment packet
`../merged_metagraph_distance_moment_constraints_20261020/PROOF.md`
retains compatible Walsh Gram matrices and exact compression losses.
Those are necessary wiring constraints, not a proof of connectivity.
Its hostile even leaves the entire complement graph unchanged.

Concept board:

1. Black core / path reversal / irreducibility; lost path order is restored
   by a literal Hamiltonian word.
2. Closed paths / cycle cuts / the exact coset of blue cuts.
3. Score current / triangle descent / maximum versus average.
4. Ordered rigid SCC words / black leaf basins / local boundary transport.
5. Distance moments / cut indicators / missing compatible Boolean fibre data.

Anchor: black connectivity. Niche: classify symmetry along a fixed cycle.
Wildcard: use endpoint currents to build a terminating descent.
The results below close the niche and provide an explicit operation for
the anchor. The wildcard's tempting averaged sign is refuted.

## 2. Every strong black-active tournament has a black closed path

For a Hamiltonian path P, let r_P be its positional reversal.
Call P blue when r_P is an anti-automorphism, black otherwise.
Call P closed when its last vertex beats its first.

**Theorem.** A strong tournament has every closed Hamiltonian path blue
if and only if it is C3 or the regular tournament R5. Consequently every
strong black-active tournament has a black closed Hamiltonian path.

The singleton is omitted from this statement: the line model starts at n=3.

Choose a Hamiltonian cycle C=(v_0,...,v_(n-1)), using the classical
strong-tournament Hamiltonicity theorem. A primary source already routed
by [CORE-PAPERS](../../reference/CORE-PAPERS.md)
is J. W. Moon, [On Subtournaments of a Tournament](https://doi.org/10.4153/CMB-1966-038-7),
Canadian Mathematical Bulletin 9(3) (1966), 297–301, Theorem 1:
its vertex-pancyclicity theorem includes the required Hamiltonian cycle.
See the [source intake](../../reference/CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md).
The concurrently inherited
[endpoint-transport note, Section 3](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md)
also supplies the elementary proof, checked in this session: take a
longest directed cycle. An outside vertex with both arc directions to
the cycle inserts at a consecutive transition. Otherwise all outside
vertices dominate or are dominated by the cycle. Strong connectivity
forces an outside edge from the latter set to the former, inserting
two vertices and contradicting maximality. Thus the cycle spans.
Its cut starting at v_j has reversal

    r_j(v_i)=v_(2j-1-i),                    indices modulo n.

If one cut is black we already have the required path. Suppose all cuts
are blue. Then r_(j+1) r_j is the rotation rho^2, hence an automorphism.

### Odd order

For odd n, rho^2 generates rho. Thus T is circulant in the C coordinates.
For n>=9, the inherited prefix obstruction produces a black path
whose last vertex remains 2m and whose first vertex is respectively

    F1: 2;       F2: k-1;       F3: 1.

The required closing increments modulo 2m+1 are respectively 3, k, 2.
Each is explicitly positive in the hypotheses of its family.
Unit reindexing preserves this closure. Thus the inherited construction
actually produces a **black closed path**, a stronger consequence than
was needed for its original purity theorem.

The four literal order-seven paths in that proof also close: their
closing increments are 4,4,3,2, all in their respective connection sets.
At orders three and five, normalized circulants are respectively C3
and R5, whose every Hamiltonian path is blue by the inherited direct
positive-control enumeration. This handles all odd orders.

### Even order: an elementary replacement for the crossing dependency

Write n=2k. The automorphism rho^2 has order k. Tournament automorphisms
have odd order, so k must be odd. The case k=1 has no strong tournament;
therefore k>=3.

Set u_i=v_(2i), w_i=v_(2i+1), with i modulo k.
Rotation by two makes each induced U and W tournament circulant.
The anti-automorphism r_0 maps u_i to w_(-i-1). Applying its arc-reversal
rule to u_i->u_j shows that the two induced circulants have the
**same** connection set.

Choose epsilon in {1,-1} such that u_i->u_(i+epsilon); the corresponding
arcs also hold in W. Define the cross word

    a_d = 1[u_i -> w_(i+d)].

It is independent of i. The original Hamiltonian cycle gives a_0=1
and a_(-1)=0, so this cyclic binary word is nonconstant.
Since 2epsilon is a unit modulo odd k, there exists t with

    a_(t+epsilon)=1,        a_(t-epsilon)=0.             (1)

Indeed, otherwise the absence of a 0-to-1 transition around the
single cyclic orbit of addition by 2epsilon would make the word constant.

Now concatenate the two directed paths

    u_0,u_epsilon,...,u_(-epsilon),
    w_t,w_(t+epsilon),...,w_(t-epsilon).

The first equation in (1) joins the U endpoint to the W start.
The second joins the W endpoint back to u_0. This is a Hamiltonian cycle.

All U vertices have one common score alpha and all W vertices one
common score beta. Moreover alpha+beta=2k-1, so alpha!=beta.
Cut an internal U edge in the new cycle. The resulting closed
Hamiltonian path starts and ends in U. An anti-automorphism sends
score s to score n-1-s, so its two endpoint scores would have to sum
to 2k-1. Instead their sum is 2alpha, an even integer. The cut is black.

This contradicts the assertion that all closed paths are blue.
The even argument is elementary once a Hamiltonian cycle is available;
it does not use Lu's crossing theorem. It also yields a black closed
path with **equal endpoint scores** whenever the starting cycle has
all cuts blue.

## 3. Exact blue-cut cosets on a Hamiltonian cycle

Fix any Hamiltonian cycle C, without assuming all its cuts blue.
Let B_C be its set of blue cut indices in Z/nZ. If i belongs to B_C,
then

    B_C = i + H_C,
    H_C = {j in Z/nZ : rho^(2j) is an automorphism of T}.          (2)

Proof: r_(i+j) r_i=rho^(2j). With r_i anti, r_(i+j) is anti
exactly when that product is an automorphism. H_C is a subgroup.

Consequences:

- At odd n, a nonempty blue-cut set has odd cardinality dividing n.
  If it is proper, it has at most n/3 elements.
- At n=2k, its nonzero cardinality is 2d, where d is an odd divisor
  of k: the image of j->rho^(2j) has order k and its automorphism
  subgroup has odd order d.
- A cycle with any black cut has at least two black cuts. In fact
  at odd n at least 2n/3 cuts are black; at even n at least n/2 are
  black. At n=2 mod4 a proper nonempty blue set again has size at
  most n/3.

The subgroup concerns automorphisms that are rotations of this particular
cycle; it is not asserted to describe the whole automorphism group.

## 4. A literal two-black-step operation

Let P=(v_0,...,v_(n-1)) be an **open black** Hamiltonian path of T,
so v_0->v_(n-1). Reverse its path arcs, obtaining U. Then

    C=(v_(n-1),v_(n-2),...,v_0)

closes in U using the untouched shortcut v_0->v_(n-1).
The cut Q_0=(v_(n-1),...,v_0) is black: path reversal preserves
the blue/black condition, since the reversed path-arc set itself
is invariant under its positional reversal.

Choose any other black cut Q of C; Section 3 guarantees one.
Reverse Q in U. The resulting tournament T' differs from T in
exactly two arcs:

    the shortcut v_0->v_(n-1), and
    the original P arc at the cut omitted from Q.                (3)

Every other path arc was reversed twice. Thus (3) is realized by
two actual black lines, with the intervening strong tournament U.
No appeal to a quotient degree or abstract V4 realization replaces
the literal paths.

**Important failure boundary.** This is a labelled two-arc move, not
necessarily progress between distinct merged nodes. At n=4 the
non-self-converse merged node {C3->1,1->C3} has only one black
neighbor. Its two-step cut move returns to the same merged node.
Thus coset mobility alone does not prove black connectivity.

## 5. Endpoint descent: a positive exact frontier and a failed average

Write s(v) for outdegree, and reverse a Hamiltonian path from a to b.
Only the endpoint scores change:

    s'(a)=s(a)-1,        s'(b)=s(b)+1.

Hence

    sum_v s'(v)^2 - sum_v s(v)^2 = 2(s(b)-s(a)+1),
    c3(T')-c3(T) = s(a)-s(b)-1.                              (4)

**OPEN descent lemma.** Every strong black-active tournament has a
black closed Hamiltonian path with s(b)>=s(a).

This holds on every merged strong black-active class through n=7:
the class counts are 0,1,5,22,216 at n=3,...,7. It follows already
from Section 2 for regular tournaments, and for any tournament that
has a Hamiltonian cycle with all cuts blue. The remaining case
requires choosing among mixed-color cycle cuts, possibly changing
the cycle.

**PROVED limited descent.** The descent lemma also holds for every
strong non-self-converse tournament: every Hamiltonian path is black,
and the endpoint score differences over the cuts of any Hamiltonian
cycle sum to zero. At least one difference is nonnegative.
Thus the unresolved family is nonregular self-converse tournaments
whose cycle cuts need not all be blue.

The weaker theorem in Section 2 does not imply the score inequality.
Nor would the descent lemma alone prove black connectivity: iteration
can land at distinct nonstrong black-active basins.

An attractive average shortcut is false. Summing s(b)-s(a) over all
closed path-orbits gives zero, by grouping rooted cuts of Hamiltonian
cycles. But the blue summand need not be nonpositive. At n=5, canonical
merged code 10 has

    blue closed current +2,     black closed current -2,

while still having two black closed paths with s(b)>=s(a).
The ascending fixed-path tiling x=0 is a representative, with scores
(1,1,2,3,3). Thus a negative black mean does not remove its successful
maximum. Positive blue currents occur in 2,4,38 merged classes at
n=5,6,7.

Even equal-score endpoint selection is too strong: at n=5 code 12
has black closed paths but its only equal-score closed path is blue.
At n=7 code 310 has no equal-score closed path of either color.
Both nevertheless admit the required strict black descent.

### A literal order-nine obstruction to staying on one cycle

The stronger claim that **every** Hamiltonian cycle of a black-active
tournament has a descending black cut is false. On vertices 0,...,8,
take the out-neighbor bitsets

    (98,92,409,273,225,206,140,267,115).

Equivalently this is ascending fixed-path tiling 48440728 with free
pairs in lexicographic order. It has the Hamiltonian cycle (0,1,...,8),
score sequence

    (3,4,5,3,4,5,3,4,5),

and blue cut indices exactly {0,3,6}. All six other cuts have
`score(end)-score(start)=-1`. Their path reversals leave c3 unchanged;
none gives the strict triangle descent requested in (4).
The blue cuts instead have endpoint difference +2.
This realizes the periodic strictly-increasing runs allowed by the
coset theorem; the coset theorem cannot exclude that configuration.

Nevertheless the literal closed path

    (0,1,2,3,8,4,5,6,7)

is black and has endpoint scores 3 and 4. Its reversal lowers c3 by
two. Therefore this is a hostile to fixing the cycle, **not** a
counterexample to the still-open global descent lemma.

## 6. The SCC normalization obligation

After descent, nonstrong tournaments remain. Their Hamiltonian paths
traverse the ordered SCCs. In particular, a nonpalindromic ordered word
in singleton and C3 components is black-active but has one ordinary
path orbit. Its converse-merged fundamental black degree is one.
A proof must connect these basins through strong intermediates.

The formal local rewrite `C3,1 -> 1,C3` is deceptively simple:
on vertices a,b,c,d, the first block has a->b->c->a and all three
beat d. Flipping only c->a and b->d produces a singleton a followed
by the cycle b->c->d->b. But a local path construction on this block
changes both block endpoints. Inserting exterior SCCs destroys the
global Hamiltonian boundary arcs needed by that construction.

This is exactly the boundary-sidecar obstruction: a local two-arc
identity is not yet a realization by full Hamiltonian path reversals
in every exterior. A useful sufficient next lemma would realize these
adjacent rigid-word moves by uniformly bounded black paths while
retaining the exterior attachment data.

## 7. Cut formulation and distance-moment limitation

Let h(x)=1_A(pi(x)) indicate a union A of merged nodes in the fixed-path
cube, and let tau complement all bits and sigma reflect coordinates.
Then A is a union of black components exactly when

    h(x)=h(tau x) whenever sigma x != x.                         (5)

Thus the complement-antisymmetric function h-tau h is supported on
the blue subspace Fix(sigma). This exact Boolean formulation is a
potential interface with the new distance-moment packet.

However its Walsh Gram positivity and compression-loss identities
hold for many partitions; they do not classify the tournament-fibre
indicators satisfying (5). Connectivity would require proving that
the only such indicators on the black-active locus are constants.
That is not a consequence of the moment identities currently proved.

The missing invariant is literal path normalization and its compatible
boundary data. Sections 2–4 retain it; the false mean shortcut discards it.

The concurrent endpoint-transport note, equations (2.1)–(2.3), further
prescribes each node's black-neighbour triangle-current polynomial and
folded-score capacities from its black Hamiltonian-path endpoint
inventory. A proposed black-component cut must contain the required
destination shells. Regular black-active nodes have black neighbours
in the almost-regular shell, which is strong for n>=5. These are useful
necessary constraints, but the inventory omits whether the path is
closed and the destination's identity within its shell. It therefore
does not prove the closed-path descent claim or connectivity. Its
regular-end-SCC characterization of current rigidity applies to the
total path inventory; no such iff is inferred for the black restriction.

## Reproduction

Run `python3 structure.py` from this directory. The script imports only
the inherited independent fixed-path canonicalizer, then checks every
tiling at n=3,...,7, closed-path currents, the descent frontier, the
displayed hostiles, blue-cut cosets through n=6, and the constructive
even two-orbit proof on all parameters at n=6,10 and selected larger
orders. It also verifies the literal order-nine same-cycle hostile and
its successful changed-cycle descent. It writes only its own `structure.out` when shell redirection
is explicitly used.

No all-order connectivity or non-leaf-bridge classification is claimed.
