# Independent review: pair-entry and exposed-root bridge

**PASS, no correction.** The pair transfer, adjacent clone-root switch, nonadjacent two-birth graph isomorphism, retained-root operations and finite target criterion are valid at their stated scopes. None is an all-target coverage theorem.

## Exact source snapshot

The reviewed packet is session197-hybrid-constructions/pair-entry-bridge. Its source was not modified.

- PROOF.md: fcff7ca192c4055243aeab2344582686214aa7aa25504d1233f3b6b33465a856.
- core.py: dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042.
- verify.py: ada513bcfc1bfb9adfa089a81438eeed696a90d5403527de41517e2224c9a292.
- control_graphs.json: d4603d925bf5d647d98fd6e828ea4d64189e09b762ba8a7aceaaecda60f227d9.
- VALIDATION.json and VALIDATION_OPTIMIZED.json: 69fb3aff826d548fcac47e5b8ee6baf292bce80aed918949beda8c15e9dfcbc0.

The frozen structural dependencies listed in the source have already been independently reviewed. This review checked the new uses and their scope rather than substituting finite inventories for those all-order inputs.

## Proof review

The scalar pair correction K is exactly the forbidden old-prefix/gadget-prefix join when the gadget prefix is the singleton new root. It carries no missing tail multiplicity because the remaining gadget vertices form one transitive path. A path cannot use both reversed arcs, since both leave the same vertex. The ending row and reserved-start two-cover row are independent bijections and give the stated formulas.

For an adjacent pair a->b with a's sole out-neighbor b, the K identity follows by deleting a from the reserved path. The physical first/last-clone isomorphism and both rooted response rows are correct. Both marked vertices really remain exposed. The positive difference for a subsequent gadget is therefore an authorized choice of root on the same graph, not an implicit replacement by another equal-count graph.

In the nonadjacent case, deleting the first vertex of the three-start-one block leaves a strong C1 graph: the cycle predecessor shortcuts to the middle vertex because the latter's sole out-neighbor is the block's last vertex. The displayed permutation checks every local and external incidence. In particular it sends the desired new root to the LAST clone. The other clone generally has a different profile; the explicit H13 hostile example correctly preserves this distinction.

The J, Z and O rows respect their different retained roots and deletion responses. The formal inverse hostile is genuine: an H3 tournament is an ordered join of a directed triangle with singleton components, so endpoint count2 cannot occur. This does not invalidate the transfer, only the claim that its formal inverse always produces a physical state.

For the bounded inverse criterion, every strong C2 graph other than regular5 reaches a genuine two-entry event by following its canonical start-two decomposition backwards. The complete C1 source grammar supplies all actual marked pairs. Source order is bounded using induced-path monotonicity and the cited strong-order inequality. The bound K<=R is termwise, and the lower bound h'>=h+(L-1)R for L>=2 is valid. Recording terminal values before discarding excessive R is essential and is stated correctly. The remaining ordinary closure has finitely many bounded profiles and loses no future H<=N. These conclusions are exact for each target; the source never infers cofinite scalar coverage.

## Reproduction and new independent checks

The copied source verifier ran normally and with Python optimization; both outputs match the frozen receipt byte-for-byte.

The separate independent.py imports no source counter or subset DP. It inventories Hamilton paths and ordered path covers by literal permutations and cuts, constructs every tested adjacency matrix independently, and checks:

- 208 general pair-gadget profiles, including arbitrary weak source graphs of orders2 and3;
- 54 adjacent first/last-clone profiles on all labelled C1 graphs through order4;
- 54 complete nonadjacent graph isomorphisms with the exact distinguished root;
- 306 J/Z/O rooted response rows.

Normal and optimized independent receipts agree.

- independent.py: 319586f30c693922fb616a15d31d6d9b1e7f8380f68cc8c4e2e41a6c55a21abf.
- INDEPENDENT.json and INDEPENDENT_OPTIMIZED.json: 853720fede90beba1cb174df0a546899c7f5180ccfd88841a93a2fc60efe61de.

The review does not cover any later proposed cone corollary or two-parameter family added outside this frozen snapshot.
