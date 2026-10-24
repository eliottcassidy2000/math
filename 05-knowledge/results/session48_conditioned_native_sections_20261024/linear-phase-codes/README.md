# Native affine C3 phase codes

Read AFFINE_PHASE_CODES.md for proofs and scope.

- Canonical phases are k=owner(8r)-8r modulo3. Raw owner coordinates require the right-side shift (8r)A1, separately from the converse-head phase shift2A1; other owner embeddings require an explicitly transported phase chart.
- Arbitrary explicitly declared affine linear predicates on actual suffix-owner phases give exactly affine ternary codes, with no discarded helper phases or invented resources.
- Budget components are cosets of the span of low-weight kernel vectors. An original-generator basis gives a constructive route in at most its rank events, without enumerating the full state space.
- Hamming distance distributions are exact code weight enumerators; one changed coordinate costs one slot and three actual owner-clock changes.
- A precise DAG characterizes per-slot full-record non-return. A native five-resource control is budget3-connected but needs budget4 for a specified no-return endpoint route.
- The three-mode pair/helper family has3^(L+1) states, three budget2 components for L>=2, and its first mode bridge at L+1 slots.
- The27-state/five-resource case has full budget3 Laplacian spectrum0^1,9^8,12^6,15^12. A selected all-size phase graph has an exact Fourier spectrum and sharp normalized gap27/40, with growing move budget and exponential degree explicitly retained.

This is not arbitrary pairwise-owner-equality realization, one global observer, or a spectrum of an older metagraph. Full guards must truly factor through fixed incidence or be retained separately. Converse and native transports are authenticated, not inferred from coordinate notation.

Reproduce with Python standard library only:

    python -B check.py
    python -O -B check.py
    python -B verify_manifest.py

The deterministic receipt is frozen; --write deliberately regenerates it. Sourcearchives.txt identifies vendored provenance. MANIFEST.json covers the entire standalone package with LF-normalized SHA256 hashes, excluding only itself.
