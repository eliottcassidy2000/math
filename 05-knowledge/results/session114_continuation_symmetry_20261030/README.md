# Symmetry under continued operations

**Session 114, 30 October 2026. PROVED scoped results and independent exact audits.** The major conjectures remain open. No literature-priority or proof-assistant certification claim.

## 1. Lonely Runner: an exact all-size structured family

For q>=3, set C_q={1,...,q−1,q+1} and let s(q) be the least prime not dividing q. The [Farey-tent proof](lrc/GENERAL_CORE.md) determines the maximum clearance after adjoining every positive speed v:

- q does not divide v: maximum 1/q.
- v=qm: maximum max{1/(q+2), m/(qm+s(q))}.

For m>=2, the family is tight at the conjectured (q+2)-runner bound precisely when every prime at most 2m divides q. The strict branch has exactly two maximizers and boundary slope 2(1/s(q)+1/(qm)). The tightness criterion is inherited from Goddyn–Wong. The uniform non-tight spectrum closes the earlier THM-1286/1271 gate-completeness obligation and yields explicit positive excesses whose runner-count-squared normalization tends to zero. Previously computed rungs retain their credit. This is a theorem for a structured family, not arbitrary speeds.

The [q=12 certificate](lrc/RESONANT_TAIL.md) separately closes the all-m ladder left open in THM-1230. It retains the essential m=1,2 exceptions, both with maximum 1/14. The proof uses exact directional binding owners; a modular orbit count alone would lose them.

## 2. Collatz: exact ordinary all-layer periodic means

[Periodic content rates](collatz/PERIODIC_CONTENT_RATES.md) give a closed rational formula for every fixed positive valuation controller and entry, with an O(log N) error for N repetitions. A fixed exit does not change the mean. The proof controls all layers simultaneously using cyclic affine subgroups, commutator depth, primitive word roots and LTE; it does not assume that a Haar mean equals a trajectory mean.

If the primitive controller's affine order grows fivefold from modulus 5 to 25, the asymptotic content per valuation clock is at most 5/8, sharply. This is an all-layer theorem for a broad class of fixed controllers. Stalled controllers, arbitrary changing words, the universal finite-word bound and Collatz convergence remain open.

## 3. Keller maps: one genuinely symmetric coordinate suffices

[One-coordinate symmetry](jacobian/README.md) proves that one coordinate with genuine mixed-sign torus weight forces every polynomial Keller companion into an explicit triangular automorphism. The companion needs no degree or symmetry bound. Finite-order certificates and the first aliasing shell are classified. Actual symmetry remains a hypothesis; averaging cannot manufacture it for free. Planar JC remains open.

## 4. Tournament charts: singleton output, exponentially many owners

[Parent multiplicity](charts/PARENT_MULTIPLICITY.md) counts and phase-classifies all regular cyclic parents inside C_m wr C_r. Their number is m^(r−1)φ(gcd(r,m))/gcd(r,m). The inherited singleton selector retains every one of these owners simultaneously; for m=3,r=3^k this is 2·3^(n/3−2).

Distinct named parent owners are not necessarily distinct support sets. The source-cover tradeoff is not a universal lower bound on compressed representations or unrestricted tournament enumeration.

## 5. Aliquot: a complete quotient can fail after one new operation

[Aliquot continuation](aliquot/ALIQUOT_CONTINUATION.md) classifies the complete binary response under prime dilations by source parity and the odd squarefree kernel. Adding one actual aliquot step splits a single such class into infinitely many distinguishable classes. No bounded finite-valued sidecar repairs this quotient when every prime dilation is allowed.

For a fixed finite dilation alphabet there is a different, exact finite response count. No infinite lower bound for one fixed finite alphabet, nor aliquot termination result, is asserted.

## What connects these results

A useful symmetry quotient must preserve the target under the actual next operation. Here that requires directional endpoint owners, return phase, polynomial bracket channels, named parent multiplicity, or successor factor support. The carriers differ; no isomorphism between the conjectures is asserted.

[Checks and reproduction](CHECKS.md) specify independent audits and finite universes. The proofs state their all-parameter scope separately from those tests. Prepublication arithmetic/indexing slips were repaired before they entered the published results.

Next mathematical boundaries: arbitrary LRC cores lack this complete Farey envelope; arbitrary Collatz words lack a fixed controller; arbitrary Keller maps lack certified coordinate symmetry; chart ownership does not supply an efficient full census; full mixed aliquot-language equivalence is unclassified.
