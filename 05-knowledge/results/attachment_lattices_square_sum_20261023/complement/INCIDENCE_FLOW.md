# Incidence flow: exact conditional variance, and a mixture that loses its source

23 October 2026. **PROVED elementary:** the reversible finite Markov chain,
stationary law, electrical/variance identity and constant-energy kernel.
**FINITE-EXACT:** the order-four mixture collision. The later incoming
coupling theorem supplies the exact spectral gap 8/9, as cited below.
There is no lossless decoder for mixtures or native metagraph-path lifting claim.

## 1. The incoming incidence is exactly the one in this packet

The incoming commit `cf7c4fc327`,
[correspondence/CORRESPONDENCE_AND_OBSTRUCTIONS.md, sections 4--5](../../session42_shared_lattices_square_colors_20261023/correspondence/CORRESPONDENCE_AND_OBSTRUCTIONS.md),
proves all-size connectedness of the bipartite graph whose vertices are ordinary
graph/complement types and tournament/converse types and whose weighted edges
count backward words. It also proves that the empty/complete graph pair is its
only left vertex of degree one. Its four-order obstruction independently gives
the classification in our README; our two-order construction additionally
computes the full unequal-neighbor triangle-count increment.

For unmerged graph and tournament types, its word counts are precisely our
b(G,T): normalize an ordered graph to the reference observer and obtain one
backward word. Thus there is an exact common finite carrier, not merely an
analogy between two tables. Summing rows and columns over complement/converse
orbits gives the incoming merged weights. Normalizing a merged row gives our
distribution mu_C. Its connectedness proof also works without merging, since
the ordinary graph and affine tournament relabelling actions already generate
all single-edge translations.

## 2. An all-size reversible chain on graph types

Fix n>=2, let M=binom(n,2), and let B=(b(G,T)). Write

    r_G = sum_T b(G,T) = n!/|Aut(G)|,
    c_T = sum_G b(G,T) = n!/|Aut(T)|,
    Z = sum_G r_G = sum_T c_T = 2^M.

Starting at G, choose T with probability b(G,T)/r_G; then choose G' with
probability b(G',T)/c_T. The resulting graph-type transition matrix is

    P(G,G') = sum_T b(G,T)b(G',T)/(r_G c_T).             (1)

It is stochastic and reversible for

    pi(G) = r_G/Z,
    pi(G)P(G,G') = (1/Z) sum_T b(G,T)b(G',T)/c_T.        (2)

Thus pi is exactly the graph-type law of a uniformly random labelled graph.
Every diagonal entry is positive, since every row has positive mass. The
incoming connected incidence graph makes the two-step graph chain irreducible.
It therefore has a unique stationary law and converges to it from every
initial distribution. The elementary incidence argument alone supplies no
uniform spectral gap. The stronger incoming theorem identified below does;
an efficient native sampler or enumeration algorithm does not follow from it.

Equivalently, D_r^(1/2) P D_r^(-1/2) is the Gram matrix

    (D_r^(-1/2) B D_c^(-1/2))
      (D_r^(-1/2) B D_c^(-1/2))^transpose.

Its eigenvalues lie in [0,1], with eigenvalue one simple. Zero eigenvalues
are permitted and measure lost linear information; they do not obstruct the
connectedness of the chain. The identical construction works after merging
graph complements and tournament converses, using the summed weights.

## 3. Electrical energy equals the conditional variance that was discarded

For any real graph statistic f, put C_(G,G')=pi(G)P(G,G'). Its Dirichlet form is

    E(f) = (1/2) sum_(G,G') C_(G,G') (f(G)-f(G'))^2
         = (1/Z) sum_T [sum_G b(G,T) f(G)^2
                  - (sum_G b(G,T) f(G))^2/c_T].        (3)

The bracket is c_T times the variance of f under the conditional graph law
b(G,T)/c_T. Expanding the square proves the identity. The conductances in (2)
therefore give a literal electrical network, and (3) measures exactly how much
source variation remains after the tournament observation is fixed.

The energy vanishes iff f is constant on every graph fiber over T. Incoming
connectedness then makes f globally constant. This is the weighted form of
the common-invariant obstruction, rather than a metaphor about information
flow. It concerns the observer-incidence carrier. It does not identify these
edges with Hamiltonian-path complement moves or transport an arbitrary walk
to an actual port-preserving tournament wiring.

## 4. Distinct source distributions do not make arbitrary mixtures reversible

Our finite census distinguishes every individual graph/complement pair through
order seven by its triangle histogram, hence also by mu_C. A stronger proposed
decoder for arbitrary mixtures fails already at order four.

Let E be the empty/complete graph pair, S the pair containing the three-leaf
star, and D the pair containing two disjoint edges. Their three coordinates
below are the transitive tournament type, the source/sink-triangle pair, and
the strong four-tournament:

    mu_E = (1,   0, 0),
    mu_S = (1/2, 0, 1/2),
    mu_D = (1/3, 0, 2/3).

Therefore

    mu_S = (1/4) mu_E + (3/4) mu_D.                    (4)

A known star source and a randomized mixture of the empty and matching sources
have exactly the same tournament-type observations. In triangle polynomials,
Z_E=24, Z_S=12+12z^2 and Z_D=8+16z^2 satisfy Z_E-4Z_S+3Z_D=0.
All six individual rows are distinct; the collision is genuinely in mixtures.

The dimension obstruction is general: six merged source rows at n=4 feed only
three target coordinates; at n=7, 522 source rows feed 272 target coordinates.
Thus the order-seven observation map has at least 250 independent linear
dependencies, even though its 522 individual rows are pairwise distinct. Any
nonzero left null vector of a row-stochastic matrix sums to zero, so its positive
and negative parts normalize to different probability mixtures with identical
observations. Exact ownership or a source selector is essential if mixtures
are to be reversed.

## Checks and scope

`check_flow.py` counts backward words directly and compares the unmerged counts
to the independently constructed ordered-object incidence in `receipt.json`.
Through n=4 it checks orbit row/column totals, detailed balance, stochasticity,
self-loops and connectedness. It checks the quadratic-form identity on basis
vectors and their pairwise sums, which determine the full symmetric quadratic
form. Exact rational elimination produces (4); no two distinct rows collide,
and the displayed three-source linear dependence has minimum support in this
census. All arithmetic uses integers and `Fraction`.

```text
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/check_flow.py
python -B -O 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/check_flow.py
```

The all-size flow follows from the proved incidence balance and incoming
connectedness. The finite matrix supplies an explicit information-loss witness,
not a conjectural explanation of a large computation.

## Later incoming spectral completion

Commit `313ad831f1`, [the coupling theorem](../../session42_shared_lattices_square_colors_20261023/coupling/COUPLING_OPERATOR.md),
diagonalizes this same P by graph-shaped Walsh sums and proves exact gap8/9
for every n>=3, including the merged carrier. This is an inherited stronger
theorem, not a new claim established by our order-four matrix tests.
[INCOMING.md, section 5](../INCOMING.md) proves the explicit bridge to our
triangle statistic: its centered order-average is the P3 Walsh mode and
contracts by exactly1/9 under P. Shared-port and fixed-path restrictions
change the carrier; the incoming gap cannot be imported to them without proof.
