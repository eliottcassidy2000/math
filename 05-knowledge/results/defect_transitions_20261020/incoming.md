# Concurrent results integrated after the local proofs

The session began from `d9674b6a00`, read incoming results while working,
and rebased onto `7eea99e6eb` before its first push. Both sides of the
concurrent corrections in `01-canon/MISTAKES.md` were retained.
This addendum distinguishes inherited computations from our own experiments.

## Endpoint rows constrain the proposed quotient

Commit `8dec6684b0` supplied the
[endpoint-transport theorem](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md).
For a tournament representative, the colour-resolved row of neighbour
triangle counts is determined by its Hamiltonian-path endpoint inventory:
the current is \(d(\mathrm{first})-d(\mathrm{last})-1\), with the appropriate
automorphism and converse-fibre factors.

The root read its proof. Our [blue checker](blue_check.py) then applied the
row constraint to our abstract rewiring and rejected it. At node zero its
blue current row changes from \(z^{-3}+2z\) to \(z^{-3}+z+1\).
Thus the fake demonstrates loss in degree/periphery/common-hub data only;
it does not survive the stronger row constraint. The actual local black
detours are literal path certificates and are unaffected.

## The independently certified order-nine black core

Commits `126830b14b` and `7eea99e6eb` supplied
[the blue-quotient continuation](../blue_quotient_20261020b/README.md).
Its [audit](../blue_quotient_20261020b/AUDIT.md) records a complete
191,536-class ordinary tournament census, 97,144 converse-merged nodes,
and a certified black spanning tree on all 97,138 black-active nodes.
All 65,536 blue tilings are separately checked.

**Inherited FINITE-EXACT result; not rerun here:** at order nine the full
metagraph has five bridges, all leaves. After black contraction there is
one black-active hub, five single-edge leaves, and one triple-edge
pure-blue attachment. Blue-only bridges are not full-metagraph bridges.

Our [trivalent hub theorem](blue.md) proves the last attachment pattern
uniformly in the order, without a finite-census extrapolation. It does
not prove that the black-active hub is unique at every order.

The incoming [closed-path proof](../blue_quotient_20261020b/structure.md)
also gives every strong black-active tournament a black closed Hamiltonian
path and describes blue cycle cuts as a coset. The root checked the even
two-orbit cycle construction and the two-black-step arc identity; the odd
case uses the independently repaired prefix theorem.

Two different obstructions now point to the same next operation:

- The incoming order-nine cycle can have black cuts but no descending
  black cut. A different cycle repairs the example.
- Our \(R_4\to Z\to R_4^{op}\) family can have black paths but no black
  replacement with the chosen two endpoints. A method must change them.

The combined target is therefore a black-path normalization that can change
both cycle and endpoints while keeping an exact progress measure.
Even a universal strong triangle-descent lemma would leave distinct
nonstrong SCC-word basins to connect. Neither packet closes that obligation.

## Overlapping and complementary Petersen results

The concurrent [snark proof](../blue_quotient_20261020b/snarks.md) independently
obtains the same isochromatic two-wire linkage and iterated Petersen
insertion mechanism as [our pole proof](petersen.md). This is overlapping
research, not two separate new constructions or a literature priority claim.
Both retain boundary support and a rooted minor certificate.

The incoming proof adds the complementary heterochromatic reduction:
two disjoint Kempe paths plus a shortest connector yield the rooted
two-adjacent-junction minor. The root checked the support-forced pairing,
the connectedness argument and the connector construction.
Its four-cut proposition therefore handles both colour-open sides.

The graph class remains essential. Reducing Petersen's pole in its original
closure gives a looped dumbbell; that is a support-preserving minor but not
a smaller simple bridgeless snark. Our insertion theorem avoids this error
by restoring a specified predecessor exactly.

The incoming packet also constructs a **connected** 17-vertex seven-pole
with profile \(4J\), with two authenticated reductions to the repaired
junction. This complements our explicit triangle-times-Petersen pole
with profile \(2J\). Both are constructive examples. Neither proves that
every arbitrary pole of support \(J\) contains the simultaneous rooted
claw-plus-two-wires forest; our separately forced Kempe linkages leave
that simultaneous-disjointness question open.

The incoming [112-closure measurement theorem](../blue_quotient_20261020b/junctions.md)
recovers all 91 seven-port extension coordinates by adding seven physical
tree probes to the 105 junction closures. This is an inherited theorem
about extension tensors, not graph reconstruction or a minor certificate.
No computation from it is used to establish our defect coefficient.

## Connection ledger after the pull

The late [marked switching-class correction](../../../01-canon/theorems/THM-4726-level-branch-as-a-weighted-count-over-two-element-decorations-and-the-decorated-flip.md)
was also read. It changes the order-eight decorated fixed count from nine
to seven: an isomorphism of the underlying object need not commute with
the retained involution. Its exact weights use the relevant centralizers.
This is a different carrier from ordinary weighted tournament cores and
changes none of our counts. It does give a concrete parallel to our rule:
retain marks through the operation, and quotient only with the equivalence
appropriate to those marks. A bijection between these carriers is not asserted.

| Source | Target and exact map | Preserved | Still lost |
|---|---|---|---|
| Path endpoint inventory | Coloured neighbour-current row | Every current coefficient | Destination identity within a score shell |
| Local black path certificates | Black-component quotient | Explicit connectivity and nonbridging | A terminating all-order normalization |
| Cubic pole colourings | Kempe paths and rooted replacement | Labelled support and the certified minor | Graph-class validity unless separately checked |
| Deficiency-two descriptions | Canonical heavy-gap owner | Exact orbit count and third coefficient | General successor-overlap quotient |

The useful common lesson is operational: exact local response data must be
paired with a legal move and its boundary. Equality of aggregate counts,
extension tensors, or degree inventories alone does not supply that move.
