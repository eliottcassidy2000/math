# Final synthesis scope audit

Status: APPROVED after correction. Scope: Session37 synthesis sections 1–6 and the model/operation caveats in section 10. Sections 7–9 are not recertified by this scoped review.

Approved synthesis SHA256: `170061540ab3e2ba4a0719624df19fde47a9c4ef8b673f9be960c98d450e0c96`.

## Corrections verified

- The two connected-space dimensions do not first differ at order eight. Exterior graph connected dimensions are 1,0,1,2,8,44,398,6414; the free Lie-super dimensions are 1,1,1,2,8,44,398,6413. They differ at order two, agree at orders one and three through seven, and differ again at eight. Both the synthesis and corrected producer note now say so. Their parity assignments already differ at order one.
- The exterior-generator minima 2n and 2n-1 are explicitly restricted to n>=1.
- Both loop multiplicities are explicitly set to one: odd cycles receive weight two and even cycles cancel. Cycle lengths themselves do not double.
- The graph determinant condition is explicitly a condition on invariant/coinvariant orbit lines, not on whether a monomial exists in the full exterior algebra.
- The exact graph Fock algebra is distinguished from the PBW associated-graded algebra of the SCC word model.

## Accepted scope

The edge-parity and vertex-parity models remain separate. Character equality is ungraded; tournament converse is not exterior parity. The graph superalgebra has a genuine symmetric/exterior component factorization, whereas PBW does not identify the original word algebra with a supercommutative algebra.

The old real-root threshold is correctly attributed to THM-020/025. The new one-back-arc family formula and H specialization agree with the prior proof audit and independent full-polynomial checks through n=15. The obstruction applies to an entire cycle-fugacity polynomial being det(I+zK) for finite Hermitian K>=0; it does not exclude interacting, signed, nilpotent or non-Hermitian models or a determinant representation of a single scalar.

The section-2 matrix distinction is sound. An additional direct 5! check of the C3-plus-two-singletons transitive and cyclic quotients gives H=3 and H=15 respectively and the same skew characteristic polynomial x^5+10x^3+21x.

The section-10 caveats properly avoid claims of physical quantum dynamics, operation-preserving equivalence, literature priority, metagraph classification or unrelated conjecture resolution.

## Correction lineage and frozen artifacts

The initially reviewed synthesis had SHA256 `56a186ecd57a1f34184a1e43e169cf5a82af64cd2d95f822c96e270d62499663`; its false agreement-through-seven sentence was reported before approval. The accepted hash is above.

The previously recorded producer-note hash was `9e6d2b16c1b5e3cabe8226e3a8e02b3106083444764b1f45285881f66154bb32`. The current corrected `SUPER_PBW_AND_SIGNED_CURL.md` has SHA256 `380f5f62eddcbe3479f3240b6afed1b3bcf51c8edc79c5ac5c574ee9e18ab5a0`. Its comparison paragraph is independently consistent with this audit's existing exact receipts. This lineage records successive audited versions, not an assertion that the intervening byte difference consists only of that paragraph.

No mathematical computation changed. `AUDIT.md`, all six audit sources and their existing computation receipts remain unchanged. The input-hash manifest is refreshed for the corrected producer note and the newly reviewed synthesis; the publication manifest and checksum list now include this scope report.
