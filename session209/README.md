# Session 209: merge Sessions 14, 15 and 27

The principal interpretation of “14, 15 and 27” is the prior sessions:
Session14's paired-black leaf reduction, Session15's complete S4-center
classification, and Session27's automorphism-weighted SCC enumeration and
marked response experiments. A separate packet checks the numeric reading
at order14 with counts15 and27. No literature-priority claim is made.

## Exact leaf spectrum and two different typical count laws

[Proof, standalone checks and independent review](leaf-spectrum/README.md).

The fundamental paired-black converse-merged multigraph has leaf count
spectrum exactly

    {3^k : k>=1} union {5*9^k : k>=0}.

Minimum orders **among leaf realizations** are3k+1 and6k+4, respectively.
The joint order/count generating function gives exact multiplicities.
No15^k is a leaf count. This does not prohibit tournaments with those
counts; in particular Session208's order14 H=15^6 construction remains valid.

Mixed S4-centered leaves are exponentially rare, yet attain the largest
leaf H at every order congruent to0 or4 modulo6. A uniform quotient leaf
has typical log(H)/n about0.213410; a uniform labelled tournament conditioned
on its quotient image being a leaf has rate about0.146747. Exact algebraic
constants, exponential-moment asymptotics and normal limits are proved.
These are different measures, not conflicting estimates.

## Optimal binary reconstruction inside the known atom grammar

[Proof, exact experiments and independent review](cyclic-port-decoder/README.md).

For a known dominance join with s independently reversed C3/C5 factors,
every one-port response table is blind to the s labelled orientation bits.
Exactly s two-port queries of H modulo4 recover them, with a matching
information lower bound for that binary oracle. The queries are fresh
experiments, not cumulative insertions. General module hiding and visibility
are inherited from later repository work and explicitly attributed.

For any labelled tournament representing a quotient leaf with r triangle
SCCs, its COMPLETE one-port response fibre, among all labelled tournaments
on that vertex set, consists exactly of their2^r independent reversals.
After full one-port preprocessing, exactly r additional two-port modulo4
queries suffice and are necessary. The fibre presentations are isomorphic
as unmarked tournaments; this does not create2^r distinct quotient nodes.

[The exact bridge](QUERY_SPECTRUM_BRIDGE.md) is that r=v3(H) equals both
the logarithm base2 of the one-port fibre size and the optimal number of
additional binary queries. Thus the leaf-spectrum normal laws above also
describe residual query complexity. The full one-port table can require2^n
entries; its preprocessing cost is not included, so this is not a linear-
time reconstruction claim.

## A numeric order14 comparison

[Proof, pinned inventory and standalone replay](numeric-14-15-27/README.md).

At order14 the H15 stratum has55 quotient nodes, with exact paired-black
degree histogram {3:5,5:41,15:9}. The H27 stratum has105, with histogram
{1:28,9:65,27:12}. Thus H27 has28 leaves, only4/15 of that stratum, while
H15 has none. The exact first order at which H27 leaves become a majority
is23:336 of660 nodes, versus280 of561 at order22; they remain a majority
at every later order. This is a population threshold inside the fixed H27
stratum, not a threshold for realizing all odd counts.

Degrees are taken in the full quotient graph, not the induced fixed-count
subgraph. The checker also independently verifies the inherited
strong order14 witness with H=15^6 and its7/3 construction-rate upper bound;
neither minimum order nor global optimality is claimed.

## Boundaries and verification

The leaf classification is inherited from audited all-order mathematics;
small computations test the new interpretations and formulas. The numeric
packet inherits completeness of its order7 atlas and cites Busch's strong
count bound to exclude larger low-count atoms. Each packet includes its
precise dependency pins, source code, deterministic output and independent
review. Normal and optimized Python replays agree.

Paired-black multiplicity degree, distinct-neighbor degree, and labelled
native degree are different observables. The two quotient leaves at order4
share one edge, so leaf counts are not automatically bridge counts. No
complete all-tournament Hamilton-count spectrum, eventual odd threshold,
global minimum-order formula, or LRC(14) proof is established.
