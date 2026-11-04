# Independent audit of the exposed-anchor noncancellation theorem

**PASS in the stated exposed-anchor scope.** The exact singleton-port averaging lemma, general-order anchor battery, proper nested-pair theorem, and exposed-anchor laminar theorem are correct. Fully covered laminar families and arbitrary crossing families are not promoted by this result. No exact full-response collision or universal reconstruction claim follows.

## Proof review

The weighted run-reversal involution in the averaging lemma is well defined. Every external old name and the uniformly attached second port see the chosen module uniformly. Ignoring only the first port's module contacts therefore leaves a genuine bijection after reversing each module run. The singleton-attachment multiplicities are m,1,m-1,1 for no contact, incoming only, outgoing only, and both contacts. In the last case the two module names are distinct because a Hamilton word cannot repeat a name. Reversal preserves these contact types and distinctness. This proves exact integer equality, not just a residue or an individual-query invariance assertion.

The pinned incoming E1 identity is stated and proved for every n>=3. Its actual mixed quadratic/shared cubic finite differences yield the proposed battery: the coefficients of f(e_v) and f(0) are d and n-1+d. Thus the odd-order extension, optional queries, and n-1/n+1 query counts are correct. The n=2 boundary is excluded from E1 and handled by its endpoint parity. Division by two is legitimate because every tournament Hamilton-path count is odd; consequently a sum2 modulo4 forces an individual contrast2 modulo4.

For a proper nested pair, an anchor in B minus A exists. The whole-B certificate lifts with the odd H(T/B) factor by the pinned incoming modulo-four transfer. On the proper module A, the battery queries indexed by its vertices form exactly one complete singleton average with a common exterior attachment. All other queries have uniform port incidences on A. The restoration of A consequently preserves the battery sum exactly, and preserves the separate endpoint query individually. It is unnecessary, and generally false, to assume B minus A is itself a module.

For a laminar family, the exposed anchor avoids every selected proper module. Their reversals do not change arcs incident to that anchor, so the original battery's partition and attachment choices stay valid at every step. Laminarity preserves each selected set's module property throughout the sequence. Selected sets outside B lie in disjoint maximal modules and are uniformly seen by both lifted ports. Their individual count invariance is exact. This gives the claimed noncancellation and no more.

The disjoint theorem in the companion PROOF.md is also consistent with this reasoning: isolate one selected module's incoming separator, and the remaining disjoint module reversals preserve that query exactly by ordered-path-cover reversal. The six-vertex fully covered example remains a boundary control rather than an application of the exposed-anchor theorem.

## Independent finite audit

`check.py` imports no producer code and reads no polynomial response bank. It enumerates original modules and laminar families independently, reconstructs every exposed anchor, computes endpoint parity separately, and counts actual two-port Hamilton paths by literal subset recursion. It verifies the source's exact physical masks, simultaneous edge-XOR targets, integer query sums and endpoint counts, exact aggregate invisibility, and modulo-four separator.

Normal and optimized executions agree byte for byte:

-15,818 exposed laminar families and all33,808 exposed-anchor records through five old names;
-486,648 literal query occurrences, comprising82,332 distinct actual attachment queries;
-all3,828 singleton-port averaging configurations through four old names;
-1,064 fully covered families counted but explicitly not promoted by the theorem.

The independent family/anchor set equals the complete source record set, so these are not sampled records. All checks remain active under `-O`. The receipt pins ANCHORED.md, both source checkers and the complete11,356,449-byte anchor receipt.

The incoming dependency was read at commit76b3b1a77d124100f03026420f6b821533cdec15, path `05-knowledge/results/two_port_parity_separation_20261104/PROOF.md`, blob3ef22092e86c0028e162943fc6cddb3e561aea68. Its E1 and E4 scopes match their uses here.

Run `python3 check.py` and `python3 -O check.py` after restoring the sibling simultaneous_response receipt package. The audit deliberately uses direct counts for every certificate query, preserving the effect of all cubic and quartic terms without truncation. The unbounded conclusions rest on the reviewed bijections and parity proof, not finite exhaustion.

## Complete larger-response controls

`check_controls.py` separately evaluates every attachment assignment for each retained larger control: all4,096 queries on each side of the three six-vertex comparisons and all1,024 on each side of the five-vertex crossing example. It then applies the full integer Boolean Mobius transform to those actual counts, without using the producer's polynomial routine or cover formulas.

All26,624 literal counts and complete coefficient vectors match. Every off-basis coefficient is zero, and all nonzero degree counts, including every cubic and quartic, agree. The fully covered six-vertex example has exactly36 quadratic,132 cubic and132 quartic differences; the disjoint and nested six-vertex controls retain108 and84 quartic differences respectively. The uniform-block query records and their modular blindness agree as well. Normal and optimized receipts are identical.

One uncited historical census comparison was removed from the companion proof and search-receipt scope during review. The final wording makes only the current mechanism-specific finite claim. No mathematical data or theorem changed.
