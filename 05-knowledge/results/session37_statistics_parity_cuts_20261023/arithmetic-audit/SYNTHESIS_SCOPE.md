# Session37 synthesis arithmetic scope review

Status: **APPROVED** for the arithmetic/category scope below. The requested precision fixes and disambiguations were verified in the final bytes.

Reviewed file: `session37-root/SYNTHESIS.md`.
Reviewed initial SHA256: `56a186ecd57a1f34184a1e43e169cf5a82af64cd2d95f822c96e270d62499663`.
Scope: sections 2, 8 and 9, and arithmetic/category caveats in sections 1 and 10. No new mathematics, source changes, Git actions or publication authorization are part of this review.

## Required precision fixes

1. Section 8 must not imply that spectral compression is the only information-loss stage. Full Cayley U recovers S, but the separate squaring map S→K already identifies S and -S. Taking spectra loses additional information.
2. Section 2 should describe adjacency determinant/permanent data as different, not proved richer than skew/energy spectra. No refinement relation between these two observers was established.

## Small disambiguations requested

- Say skew/energy-spectrum-only estimator in section 9, avoiding an unintended claim about adjacency-spectrum estimators.
- Define pc_T(2) as unordered spanning covers by two nonempty directed paths, with singleton paths allowed.
- Attribute the odd/even sector decomposition to the half-sum and half-difference of the logarithms, not ambiguously to the rational cut bounds. The positive even sector is (log B-log F)/2.

## Substantive conclusions already accepted

- Rational/irrational and algebraic/transcendental are distinct partitions of the real numbers; Dedekind cuts represent all real types.
- Nonzero energy μ is a totally positive algebraic integer, Q(μ) is totally real, and Q(μ,sqrt(-μ)) is CM. Zero modes are excluded.
- The Cayley/charge maps preserve their specified fields; the full Cayley matrix is invertible data, unlike spectral quotients.
- Irrational energies can give torsion phases; order five is the first possible transcendental normalized angle θ/π, whereas raw nonzero θ is already transcendental at smaller orders.
- Rational fugacity values in the cited positive/convergent domain give rational determinant functions and positive nonzero transcendental real logarithms. Exact rational enclosures, rather than decimal evidence, define their cuts.
- P163 arithmetic does not use or imply additional Heegner, modular, VOA or RH results.
- All eight regular-base block extensions share skew/energy/phase/charge spectra and their determinant/log observables, not full Cayley matrices or asserted adjacency characteristic/permanental data.
- The unbounded H ratio and √N symmetric multiplicative-estimation obstruction follow from the exact path-cover inequality. The cyclic extension is strong and its transitive counterpart reducible. No fixed-strong, fixed-score or regular-target-class impossibility is asserted.
- The summary keeps classical input and repository-relative synthesis separate from literature-wide novelty and does not overclaim major open-problem conclusions.

## Final approval

The revised sections 1, 2, 8, 9 and 10 pass the arithmetic scope review. Every requested precision fix above is applied. The separate statistics-comparison revision in section 5 is present; its substantive approval remains with that lane's auditor.

Approved final synthesis SHA256: `170061540ab3e2ba4a0719624df19fde47a9c4ef8b673f9be960c98d450e0c96`.

The principal arithmetic source remains frozen at SHA256 `5318cc2a27c43587287b783982e1d81858c28ebdc3225eb5cca3762c522ec196`. No mathematical source or executable changed during this synthesis review. The existing 13,865-check independent arithmetic receipt is unaffected.

This approval covers the specified synthesis bytes and arithmetic scope, not publication authorization or subsequent edits. No Git mutation was performed.
