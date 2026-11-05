# Independent review: retained-root unit-clock reduction

**PASS. No correction required.** The theorem is an exact graph-level reduction for an affine return beginning at the original state. It does not classify all scalar counts, arbitrary port moves, or clocks after an ignored first transient.

## Frozen binding

Source directory: `session197-hybrid-constructions/retained-unit-clock-reduction/`.

- PROOF.md: cd68ee0b90616fd058b8e1d861c54da37f5dd23afb615babbc182ace9b40847b
- core.py: dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042
- verify.py: 12949d64dd6a66b44dbeeb910b09123373274241001891bdd3143de06c6b5261
- controls.json: 79c39b39f1a63a52b26acbf68497a875924371f130879102f1d5f8837179aa4a
- VALIDATION.json and VALIDATION_OPTIMIZED.json: eafc4eb4038900d8620f9c74595d50facb24aafabe8721736e5398444798f552

The original files were not edited. A copied source replay under normal and optimized Python reproduces the frozen receipts exactly.

## Proof review

The Z recurrence gives different initial and tail affine offsets unless e=d. Their difference is exactly 2(e-d), so the iff criterion really includes the original first step. The singleton convention H(empty)=1 is consistent. Other constant unit examples need not be singletons; they do not invalidate the stated nonconstant theorem or the structural result for order at least two.

For an arbitrary appended sink module S, a Hamilton path has at most two module fragments because its only exterior exit is r. A two-cover whose reserved path starts at r cannot exit a module fragment, so it has at most one fragment at each path end. In both settings the internal lifts are exactly H(S) and B_S(2)-2H(S). The two transitive controls determine their actual skeleton coefficients, rather than extrapolating arbitrary module structure from sample counts. Substitution gives all three displayed rows and the exact factorization of the offset by H(S). The e=d row follows independently from the ordered deletion (U-r)->S, including U a singleton.

The terminal-SCC step is sound. Every vertex of a strong terminal component can end a Hamilton path, by the fixed-end insertion argument. Hence e_r=d_r forces every terminal vertex to point to r, while every earlier component points to that terminal component. This supplies exactly the physical module hypothesis, with U=T-S still containing r. Unit offset forces H(S)=1; a nontrivial strong tournament has at least three endpoint-distinguished Hamilton paths. Thus the terminal SCC is a unique singleton v whose only out-neighbor in T is r.

Every Z newborn has the same exterior incidences as v and follows all previous clones, so the literal clone identity and the converse implication preserve the actual old root. The amplitude is positive precisely in the nonconstant case under discussion. No scalar rebasing or replacement graph is hidden in this argument.

I checked the named prime-amplitude dependency and its audited hash. Its Theorem3 excludes both amplitude-seven signs for arbitrary tournament clone sources, not only strong sources. The signs transfer correctly: H13 with affine +1 would give 7*2^m-1, while H15 with affine -1 would give 7*2^m+1. The H13 opposite-sign example and the H5 first-transient example are consistent with the scope.

## Independent physical controls

The new checker uses literal vertex permutations and ordered cuts, with no imports from the source graph DP. It verifies319 arbitrary module constructions, covering every labelled rooted source through order three and every labelled nonempty module through order three. Empty deletion counts and both empty-cover slots are included.

Across all5,405 rooted labelled tournaments through order five it checks the exact first-step versus eventual-offset boundary. Among them1,598 states have e=d. A separate subset-defined terminal-SCC finder verifies its incidence and offset divisibility in every case. All1,070 unit states have the required singleton, and3,210 explicit clone-chain adjacency maps check three successive Z births. For the smaller inputs the resulting scalar counts are also counted literally. Normal and optimized independent runs agree.

Independent artifacts in this review directory:

- independent.py: 991ee40db927e311f6f25d68986eb847d6fd439933503ec7cb272d7a5070e138
- independent.json and independent-optimized.json: ffb0ce7b5316ba25e7082642e7e6ca8ecf3ffef43134926698b8154387b8da84

The finite checks corroborate the fragment and SCC proof; no census completeness is used as an all-order premise. No unresolved issue remains.
