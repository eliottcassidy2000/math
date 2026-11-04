# Independent review of the exact tripler and symmetry boundary

2026-11-04. Accepted after proof review, inherited audit replay, and independent literal controls. No mathematical correction required.

The multiport formula immediately multiplies H and root defect by3 under TT2-to-C3 prefix replacement, since their (H,P(2)) pairs are(1,4) and(3,12). Empty-allowed two-slot profiles multiply under ordinal sums: each slot is uniquely split into its block pieces, including empty pieces. Thus the stated L->TT2->R generalization is correct. The count identity needs no strongness assumption. Strong suffix/nonempty W suffice for strongness, as stated.

The homogeneous triangle rotation fixes all outside vertices and is an actual order3 automorphism. It acts freely on ordered Hamilton lists, because fixing a list would fix every physical vertex. Divisibility by3 therefore persists under every alteration that retains this homogeneous cyclic module, regardless of exterior size or strongness. Both desired offsets are nonzero modulo3. The scope allows edits that break uniform incidences, internal orientation or the module itself; it does not exclude arbitrary multi-edge repairs.

The edge-reversal identity K(v,u)-K(u,v) is correct: changed Hamilton paths are exactly the ones using the changed adjacency, and cutting that adjacency gives the two nonempty spanning slots. The backward endpoint recurrences in triple.cpp and forward extension recurrence in audit.cpp implement different evaluations. Root index12 is correctly physical39 for the frozen input. At n18 the bound19*18! is below signed64 capacity, so all counts and two-slot sums are safe.

Recompiled and replayed the author's audit.cpp on all154 literal graphs, checking the source and153 reversals. This is a replay of the author's independent implementation, not a new implementation of that large computation. The receipt is replay.json.

The separate new verify.py checks actual name/induced-edge preservation and homogeneous triangle orientation; confirms the delta file contains exactly all153 distinct directed edges, not duplicates or omitted pairs; and checks every delta difference, no+2/+4, unique zero6->8, and minimum nonzero absolute change28062 at22->2 with sign negative. It also directly enumerates permutations for every suffix through order3 and every nonempty W, totaling63 TT2/C3 comparisons, including nonstrong suffixes. Both H and E_root(2) triple. Normal and optimized runs byte-match.

The supplemental51-hit census is explicitly scoped as the same bounded menu as the closure packet, not independent evidence. Inspection confirms safe pruning: any uniform target has h=c+J<=c+F<1537279, and F<=1537279-h, so the target checks are not discarded. Pruned port subsets cannot have unpruned supersets, since every marginal F is positive. Clock stopping follows the stated bound. This review did not rerun that entire supporting census.

Later cross-packet review found one supporting-census physical-label typo: in the fixed suffix order, W6144 denotes ports33,40, not40,54. The author was notified to correct that prose. The actual count, candidate matrices, tripler theorem and153 reversal certificates are unaffected. This corrects the original review's omission of that label check.
