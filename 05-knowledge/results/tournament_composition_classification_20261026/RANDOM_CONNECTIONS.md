# Reproducible random recovery and connection contracts

**Status:** recovered results are cited by path; transfers below are proved
only where a map is supplied. Historical reflections are idea provenance.

`random_recovery.py` freezes commit `fc3768ac5fcd`, samples 24 Markdown files
from 12,883 candidates with seed `812720261026`, and records source hashes in
`random_sample.json`. It reproduces the actual Windows `rg --files` sampling
order. Later additions cannot alter this sample. The excluded immediately
preceding norm-seven packet was already known, not a fresh random discovery.
Five sampled objects received the deeper probes below; the remaining 19 are
recorded rather than represented as having been fully audited.

## Coprime encoding versus lost decoration

Source: [Sylvester, Mersenne, and Lossless Encoding](../../../07-reflections/sylvester-mersenne-lossless.md),
a historical reflection. The exact survivor is that the Sylvester terms
`a_1=2`, `a_(k+1)=product_(j<=k) a_j+1` are pairwise coprime, so exponent
vectors map injectively to products of these terms. A prime dividing `a_i`
cannot divide any later `a_j`, proving both claims by valuations.

Target/map: replace a decomposition tree by the valuation vector of its
automorphism order. Preserved predicate: divisibility and multiplication of
group orders in a known extension. Destroyed information: child types,
quotient stabilizers, and permutation actions. Needed sidecar: decorated
tree and cycle index. Cheapest hostile: rigid tournaments of distinct types
all map to order one; even the 12-vertex example has identical internal
`3×3×3` factors but a different external factor. The Sylvester injectivity
argument cannot transfer because independent coprime factor coordinates
have been discarded. This stopping reason is structural.

The sample also exposed an arithmetic error: iterating `x²−2` from 4 gives
`4,14,194,37634`, not `4,7,47,2207`; the latter tail starts from 3. A correction
is attached to the historical source and the mistakes ledger. The reflection's
global Mersenne-intersection claim is not used as a proved input here.

## Equal group orders, different actions

Source: [THM-4764 — modular cancellation depth](../../../01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md).
Its regular `C25` and `C5×C5` actions have the same degree and group order but
different subset-orbit polynomials:

    C25:    (1/25) sum_(d|25) phi(d) (1+t^d)^(25/d),
    C5²:    ((1+t)^25 + 24(1+t^5)^5)/25.

At `t=1` these give 1,342,184 and 1,342,208 orbits. Target: automorphism
actions on decorations/ports. Map: apply the permutation action to subsets
or to colourings of its vertex set, then use Burnside. Preserved predicate:
fixed colouring. Lost under group order: the cycle partition of each group
element. Sidecar: the cycle index, or the full action when needed. Cheapest
test: the displayed integer sums already separate the two cases. The repo's
50-vertex realizations establish that this distinction occurs in tournament
consumers; this packet does not claim a new realization.

## Partial transport versus independent modules

Source: [odd cell gluing](../odd_cell_gluing_20261023/README.md).
For a connected graph of chart domains and partial-bijection constraints,
global sections inject into each chart and hence number at most the smallest
chart. This requires actual transport constraints. Uniform tournament arcs
between two modules impose no such bijection on their internal choices.

Target: modular substitution. Map: restrict an automorphism to its induced
quotient permutation and internal child maps. Preserved predicate: arcs.
Destroyed by a presumed diagonal transport: independent child symmetries.
Sidecar: the kernel `product Aut(T_i)` and the colouring stabilizer. Cheapest
hostile: two uniformly joined cyclic triples have `3×3=9` independent
automorphisms, not at most three. The repaired exact sequence is (2.1) of
`CLASSIFICATION.md`. For path presentations one must additionally transport
the marked path: `f(T,P)=(f(T),f(P))`.

## Finite observers and their first invisible direction

Source: [THM-4044 — sixty-clock Hasse alias](../../../01-canon/theorems/THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall.md).
The source identifies an observer kernel involving `(P^60−1)^k`; its
marked characteristic and degree cap are essential. Target: the cyclic-triple
observer and metagraph composition quotient. The transferable procedure is
to compute an observer fibre and then its first operation-visible ambiguity.
This is a method transfer, with no map from polynomial residuals to tournament
arcs. The decisive tests here are the sharp `4p` and order-five hostiles.
No Jacobian-conjecture consequence or shared algebraic kernel is claimed.

## Why binary lift choices do not imply a Kummer theory

Source: [THM-3496 — marked graph Kummer degree](../../../01-canon/theorems/THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction.md).
Its graph-to-Kummer cohomology comparison retains coefficient, field and
degree data. Candidate target: the binary prime-node converse choices in
a cyclic-triple fibre. There is no supplied field or cohomology map here;
linear nodes also contribute factorial order choices. Moreover unlabelling
has nonconstant stabilizers, witnessed by orbit sizes `16,48,48,16`.
Therefore a Galois-torsor identification fails the map/sidecar gate. The valid
replacement is the explicit finite action groupoid and mass identity (3.1).

## Portfolio and concept-board updates

The anchor is the user-assigned classification; LRC(14) remains open. The
niche is arithmetic of permutation actions; the wildcard is marked transport
and information loss. Five live concepts were compared after each pull:
modular tree, observer fibre, local group action, boundary port, and size scale.

| New evidence | Effect on the five concepts |
|---|---|
| 81/27 fibre | tree needs child colours; fibre needs stabilizers; action loses a 3-cycle; metagraph multiplicity changes; first scale is 12 |
| base-`p` blocks | tree carries repeated towers; fibre changes its `p`-part only at `4p`; action attains Sylow bounds; does not recover ports; size is controlled by digit sum |
| order-five composition hostile | same factor classes do not determine boundary composition; retain rooted presentation; scalar degree also fails; earliest size is five |
| random C25/C5² | equal group order does not fix an action; cycle index is a necessary next coordinate; no automatic metagraph reconstruction follows |
| transport/Kummer probes | independent kernels must not be diagonalized; observer groupoid is available; a geometric equivalence remains unsupported |

No new universal meta-pattern is promoted: the observer/sidecar discipline is
already in the protocol, and these results refine its concrete triggers.

## Incoming connection, separate from the frozen random sample

Commit `36346d6f19` arrived during close-out. Its
[orientation/monodromy note, Section 3](../orientation_monodromy_rule30_20261026/README.md)
constructs a signed subset algebra from tournament edge signs and explicitly
rephases its subset basis to an exterior algebra. Consequently its unmarked
algebra forgets the individual edge signs. The relevant transfer is precise:
source = the algebra with its specified subset basis; target = its abstract
isomorphism class; map = forget the basis; lost predicate = the original arc
orientation. Retaining that basis restores the edge coefficients. This is a
second concrete reason to retain markings in a classification interface;
it supplies no equivalence between the algebra and metagraph composition.
The same incoming note's four-object converse dictionary transports only
the involution, not a native size-changing product. The hinge product here
is defined on marked path presentations and does not silently fill that gap.
