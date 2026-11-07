# Session221 positive-profile cap audit

**Verdict: PASS, scoped to the corrected `session221-provenance/THEOREM.md` and positive-profile compiler.** The mathematical all-length conclusion is a proved application of inherited transfer models. The reported enumerations are exact completed finite computations, not extrapolations. No publication was performed by this auditor.

## Independent argument and boundaries

For `B>=1`, coordinate saturation at `B+1` is a homomorphism of the nonnegative integer semiring. Addition and multiplication both commute with saturation, including the zero-factor case. Consequently each nonnegative profile vector, every future transfer, and the nonnegative final readout remain exact up to the cap. This does not hold for the signed rank-14/62 representations.

Starting from the singleton tournament, at most `(B+2)^d` vectors exist. A BFS first predecessor gives a shortest nonempty word for that vector; all words reaching the same capped vector have identical bounded future outputs. A shortest accepted target therefore has a shortest representative in this finite graph. Prefix pruning at `H>B` is valid only because vertex insertion cannot decrease `H`. It must not be imported into arbitrary edit reachability.

The optional product with a regular-language automaton is mathematically valid, but is not implemented or included in these receipts. An arbitrary language can distinguish the first letter even when all first letters encode the same singleton, so all first-letter DFA successors must be initialized. This caveat is present in the corrected theorem. The strong-connectivity sidecar has no first-letter ambiguity.

The two-sided local-edit observer `cap(R M_b L)` is sound. A bit changes a physical arc only when its offset exists. Matrix singularity prevents assuming head movement updates the right context by inversion. The stated block-product-tree complexity is a standard unimplemented observation, not a tested feature.

## Independent executed reproduction

`audit_cap.py` implements a separate BFS, separate sparse matrix multiplication, and independently written backwards endpoint DP for every selected output witness. It imports only the inherited profile generator and physical word decoder from `ribbon_inherited.py`, pinned to main commit `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`, original path `04-computation/experiments/ribbon_spectrum_20261104/ribbon.py`. Thus this is an independent reachability and physical-output verification, not an independent proof or regeneration of the inherited profile transition model.

Completed runs:

- Span 1, `B=99`: 65 nonempty capped vectors, 130 nonempty transitions; longest representative length 13.
- Span 2, `B=99`: 5,102 nonempty capped vectors, 20,408 nonempty transitions; longest representative length 12; all positive odd values at most 99 except 7, 21, 47; longest shortest output witness length 8.
- Span 3, `B=47`: 22,399 nonempty capped vectors, 179,192 nonempty transitions; longest representative length 13; all positive odd values at most 47 except 7, 21; shortest witness for 47 has length 7.

The auditor starts one level earlier at the empty profile, so `cap_receipt.json` includes one additional state and respectively 2, 4, or 8 additional transitions. The empty profile has readout zero and is not recorded as an accepted tournament. Removing it gives the counts above and exactly matches the author's span-2 `B=99` run. Alphabet numbering in the independent receipt is explicitly listed and differs from bitmask order at span 3.

The independent span-3 witness for 47 has letters `[], [], [], [3], [2,3], [2], [3]`. Its reconstructed seven-vertex graph has 47 Hamiltonian paths by the independent backwards endpoint DP.

## Corrections and originality

The author applied the audit's requested positive CLI parameter validation, arbitrary-DFA first-letter clarification, and existing-offset qualifier before this PASS. There is no unresolved correctness blocker in this packet.

Session135 `acceleration/THEORY.md`, section 4, already supplies the abstract nonnegative-polynomial cap identity. THM-4880 supplies the physical ribbon transfer. THM-4881(C) already proves all-length span-2 `H<=99` saturation at 2,865 uncapped reduced states. Session216 `inverse/PROOF.md` already gives the general positive-automaton cap BFS, shortest witnesses, all-target termination, and insertion-based safe pruning, citing Session203 saturation. Thus the present contribution is the concrete ribbon specialization and marked-local-write observer, not a new termination theorem, spectrum exclusion, finite-state technique, or universal odd-count construction. Numeric-budget polynomiality at fixed width is not polynomiality in the input bit length. No external originality claim was audited.

## Reproduction and release

Run from this directory:

    python3 -B audit_cap.py

Compact audit release files:

- `CAP_AUDIT.md`
- `audit_cap.py`
- `ribbon_inherited.py`
- `cap_receipt.json`
- `CAP_SHA256SUMS`

The remaining audit scripts and receipts concern separate physical constructions and are not required for this cap checkpoint.
