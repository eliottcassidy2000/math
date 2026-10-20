# Regular score-fibre transport and icosahedral face incidence

**Status: corrected V3, rechecked after a literal path-table typo was caught before publication.** All-order statements are proved in their specified labelled score-fibre setting. The regular order-five return graph and incidence geometry are proved by structural arguments and explicit finite certificates.

## Main results

- `regular_return_theorems.md`: the same-endpoint return criterion; sharp six-edge support and block-interchange classification; the two-icosahedron graph at order five; exact endpoint formulas for local Hamiltonian-path count changes.
- `gram_return_addendum.md`: the literal incidence factorization `W=BB^T`, with raw Hamiltonian-path counts on the diagonal; actual intermediate tournaments representing all forty icosahedral faces; private-column/face factorization, two sampling rules, and exact finite spectra.
- `AUDIT.md`: amended independent proof review, finite-check universes, reproduction instructions, additional precursor attribution, and precise scope boundaries.
- `QUOTIENT_STATIONARY_MASS.md`: the reversible regular-return walk and its exact ordinary/merged stationary masses h and kappa h.
- `CORRECTION_HISTORY_V3.md`: the displayed-word correction, why the initial audit missed it, and the source-aware regression checks.

## Computation

`check_regular_transport.py` and `check_gram_return.py` are the producer checks. The second imports the first. Their JSON and normal/optimized output receipts are included.

`independent_audit.py` and `independent_gram_audit.py` use a separate arc-set/literal-word implementation. The second imports the first. Their JSON receipts and identical optimized-mode receipts are included.

`FROZEN_SHA256SUMS` preserves historical V1/V2 proof/code hashes and staging paths. `PUBLICATION_SHA256SUMS_V3_1` verifies the current corrected package. V1 and V2 remain separately preserved; `V2_TO_V3.patch` records the changed prose exactly.

`check_literal_examples.py` reads the current proof examples directly and verifies them against literal paths and tournament codes, including the original erroneous word as a hostile control.

## Figure

[The two labelled icosahedral components](regular_return_icosahedra.svg) is a projected drawing of the derived labelled return graph, not the original merged metagraph or a unit-distance embedding. All 24 displayed states belong to one ordinary tournament isomorphism class.

## Boundaries

This is a derived labelled two-complement return graph, not the original tournament metagraph. The work proves neither all-order return-component classification, black-metagraph connectivity, nor the ordinary raw Hamiltonian-path spectrum conjecture. Raw H, automorphism-normalized path mass, return-presentation mass, and simple graph degree are distinct quantities.

The additional precursor recorded in the audit is `merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md`, Section 4: the almost-regular shell's high-to-low endpoint kernel. No priority claim is made for general Gram positivity.
