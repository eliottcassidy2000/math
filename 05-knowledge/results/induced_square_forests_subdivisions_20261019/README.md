# Induced square-sum forests and graph subdivisions

Read [induced forests and finite-stage attachments](INDUCED_SQUARE_FORESTS.md), [the mod-seven quartic construction](SUBDIVISION_MOD7_CONSTRUCTION.md), and [independent review](review/INDEPENDENT_REVIEW.md). All finite forests have induced square-sum embeddings at near-quadratic height; every finite graph has an induced exact two-subdivision at O(K^4) height, where K is the subdivided vertex count. Countable extensions retain only finite-stage existence and the stated scope.

Reproduce from this package directory:

    python3 -B review/independent_exact_review.py
    python3 -B -O review/independent_exact_review.py
    python3 -B review/verify_mod7_subdivision.py
    python3 -B -O review/verify_mod7_subdivision.py

Compare outputs with the matching normal and optimized JSON receipts under review/. Together the scripts run 8,278,445 explicit checks. These finite controls corroborate the all-size proofs; they do not establish an unstated global conjecture. review/FROZEN_SHA256.txt records pre-publication reviewed-input hashes under their original source-root names; SOURCE_PROVENANCE.json gives the published hashes and presentation changes.
