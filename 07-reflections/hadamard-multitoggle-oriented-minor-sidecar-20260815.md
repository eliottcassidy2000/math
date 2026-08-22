# Hadamard multi-toggle response: from pair curvature to oriented minors

**Historical research reflection, 2026-08-15.**  Current truth is
[THM-3407](../01-canon/theorems/THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance.md),
not this note.

## Inheritance and concept board

The closest proved mechanism was THM-3403's polar inverse and uniform one-bit
tariff.  Its canonical hostile was `H_4`, where positive two-toggle
plaquettes become singular and closed quadruples sit exactly at the small
boundary.  The corrected near miss was the temptation to let a binary rank or
pair-plaquette packet classify every higher response.  The least-used relevant
sidecar was not another rank: it was the oriented event matrix, equivalently
the packet of signed core minors.

The live board was:

| object / representation | predicate | operation | lost coordinate | cheap test |
|:---|:---|:---|:---|:---|
| binary Hadamard core | maximal determinant | sparse bit toggling | orientation | compare direct and low-rank determinant updates |
| event matrix `Q` | exact determinant response | principal-minor expansion | toggled matrix itself | all `3x3` masks |
| plaquette signs | two-toggle shell | pair quotient | directed cycles | find equal-pair-data triples |
| rank-one sign rectangles | equality/singularity | full rectangle flip | higher-rank trades | row swaps and closed quadruples |
| odd full binary circuit | dependency rigidity | pair-product map | real sampled signs | compare circuit injectivity with triple response |

The wildcard was to regard the determinant on the toggle cube as a Boolean
interaction model rather than a scalar statistic.  This paid immediately:
its Mobius coefficients are exactly signed minors, and repeated row/column
events disappear because they are not matchings.

## The precise quotient failure

The pair quotient is exact at rank two.  Each event pair contributes one
plaquette sign `chi`, and `chi` determines the whole two-toggle response.  At
rank three there are two directed 3-cycle gains.  Their product is the product
of the three plaquette signs, so the loss is especially small:

- odd plaquette parity forces the signed `3x3` minor;
- even plaquette parity leaves exactly one directed-cycle bit;
- that bit splits the two possible triple shells.

This is stronger than saying merely that “pair data fails.”  It identifies
the first failed implication, the minimal missing coordinate, and the exact
condition under which the quotient remains sufficient.  The order-8/12/20
hostiles show the missing bit is realized inside genuine Hadamard cores, not
only in an abstract sign matrix.

The connection contract is therefore:

~~~text
source: labelled sparse toggle set in a normalized Hadamard core
target: pair-plaquette packet
map: retain every two-event Mobius coefficient
preserved: all responses through toggle rank two
destroyed: orientation of cycles of length at least three
sidecar: directed cycle gains, or response-completely all signed core minors
decisive test: equal labelled plaquettes with unequal triple determinant
~~~

## What the equality wall teaches

The rank-one rectangle formula `rho=1-2ab/N` turns three familiar phenomena
into one response line:

- `ab<N` is strict loss;
- `ab=N/2` is the singular wall;
- `ab=N` is a Hadamard trade.

Row swaps and closed-quadruple switching are the same area-`N` mechanism with
different rectangle shapes.  The energy proof of the trade floor then says
that no higher-rank geometry can cross back to equality with fewer than `N`
changed entries.  Thus the response formula supplies constructions while the
energy estimate supplies the global exclusion.

The odd full circuit removes closed quadruples by making pair-product
signatures Sidon-like, but it cannot remove the universal row-swap trade.  This
is a useful warning: eliminating one equality mechanism is not classifying the
equality locus.

## Incoming-work comparison: analogy, not bridge

The concurrent MISTAKE-389 correction has the same proof-engineering shape but
no object-level map to Hadamard matrices.  In that LRC lane, half-grid
integrality was pairwise-compatible data, but a mode-divisibility sidecar was
needed before it represented a common centre.  Here plaquette signs are
pairwise-complete data, but an oriented-cycle sidecar is needed before they
represent a higher determinant response.

This is a **method analogy only**:

~~~text
local compatibility + omitted higher coherence != global object.
~~~

It does not transfer an LRC predicate to a determinant predicate.  Its value is
to suggest the same hostile protocol: after a pairwise or scalar quotient,
write the smallest coherence cell and search for two lifts with the same
quotient but different target.

THM-3404's principal-part lane gives a second method-level comparison.  Separate
valuations detect poles but lose coefficient cancellation; principal parts
restore it.  Pair plaquettes detect quadratic interactions but lose oriented
cycle cancellation; signed minors restore it.  Again there is no theorem map,
only evidence for the existing controlled-forgetting method card.

## New questions exposed

1. The full determinant response polynomial is a signed determinantal rook
   polynomial.  Which principal-minor/Pluecker relations characterize the
   response polynomials realizable by Hadamard cores?
2. At toggle rank four, can cycle type separate the response into a bounded
   packet of two-, three-, and four-cycle gains, and what is the cheapest
   genuine Hadamard hostile to a triple-only quotient?
3. For odd quarter-order, the full binary circuit forbids closed quadruples.
   What are the next equality-trade shapes after universal row swaps, and can
   their support hypergraph be bounded using the signed-minor packet?
4. The order-668 core now has an exact sparse determinant oracle whose cost is
   the number of touched rows or columns, not `667`.  This makes targeted
   searches for low-support near-trades feasible without dense determinants.
5. The equality masks form a minimum-distance code in the labelled toggle
   cube.  Determine its local distance enumerator from the response polynomial,
   while keeping Hadamard equivalence distinct from labelled Hamming geometry.

The stopping boundary is honest: THM-3407 classifies two toggles and the exact
information needed at arbitrary rank, but it does not classify higher shells
or equality trades.  The next useful move is operation-first: enumerate cycle
types or trade supports, not another scalar determinant histogram.
