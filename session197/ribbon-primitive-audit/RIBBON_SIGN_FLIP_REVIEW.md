# Independent audit of the span-two sign-changing interface

Result: PASS, with frozen source files unchanged.

Reviewed inputs and SHA-256:

- `session197/ribbon-sign-flip/PROOF.md`: `a2ffb030a46d021064bc3d68a94c3e739e1115712d25e5bd02c5c7ed7c7d18d5`
- `session197/ribbon-sign-flip/verify.py`: `a31362b079f04cadd4f553c72b4c73695990378cf90a25c41465e7d0dbd6673b`
- `session197/ribbon-sign-flip/VALIDATION.json`: `b5cb62ad8c11e49c5a85f0a35afd1c7a7764d6caf1538b6040d08024ed06a199`

I checked the reordered A graph, transitive D graph, exact width-two word conventions and seam distances. Strongness does not assume D is strong: its displayed entry-ending Hamilton path splices with the old Hamilton cycle, and the forward join cannot be the exceptional seam. The next exit may then be used as the root of the constructed cycle.

The two profile matrices, charge reversal, D-squared identity, reduced affine maps and inverse formulas all agree algebraically. The A cone is strict under either charge. The D lower boundary is genuinely inclusive, equals R-h, and correctly accepts the seed's first D step. The growth inequalities prove strict inverse descent. Oddness, residue preservation, observable-gcd identities and zero-density count all hold with the stated limited language scope.

I reran the full normal and Python-optimized checker; both match the frozen receipt byte-for-byte. These reproduce the 15 literal graph words and 2,047 state words. I added a separate inverse test that does not use the cone-selection code or the displayed inverse formulas: it inverts both 2-by-2 affine matrices by their determinants, recursively tests all possible predecessors, and compares with exhaustive forward words. It agrees on every one of 125,500 admissible charge-recorded states with h<=501, including the D boundary and all accepted states. This auxiliary result is recorded in `ribbon-sign-flip-independent.json` beside this report.

No correction or unresolved blocker was found. This is an exact state-language decoder and sparse physical construction, not a scalar completeness theorem.
