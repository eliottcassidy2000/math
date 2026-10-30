# Odd eighths: uniform formulas with their carriers retained

**Session 115, 30 October 2026. In progress.** Scoped proofs and exact checks; no major conjecture closure or literature-priority claim.

## Verified phase classifications

[The phase classification](quadratic/ODD_EIGHTH_PHASES.md) gives all sixteen relative Fano gauges, with linear-coordinate orbit sizes1,7,7,1. Their eight cubic classes yield exactly the odd-eighth Walsh measurement profiles. The profile is relative to a chosen base table and is not an invariant of abstract octonion isomorphism.

The uniform extension characterizes, in every binary dimension d>=3, precisely the Boolean phases with a Walsh coefficient of magnitude2^d−2: they are single-point defects of affine phases. Their character-twisted probabilities lie in {1/2^d,1−1/2^d,1/2−1/2^d,1/2+1/2^d}. Exhaustion of all odd-weight phases is special to dimension3; dimension4 supplies a hostile.

The producer passes11,796 explicit gates in normal and optimized modes. Independent proof review and a separate fast Walsh transform checker exhaust all256 dimension3 and all65,536 dimension4 Boolean functions, with66,103 gates and matching normal/optimized output. Scripts and receipts are included.

[Word-address transport](unitary/WORD_UNITARY_CLASSIFICATION.md) gives an exact all-depth classification. The multiplier change q→p preserves products, real values, positivity and pullback structure on depth-A character algebras iff p≡q modulo2^(A−1). Distinct multipliers first fail at depth v₂(p−q)+2. The signed half-character defect explains why3→5 first fails at the odd eighths.

The unitary producer passes89,001 gates; the separate literal-decoder audit passes44,313, normally and optimized. [Independent proof review](audit_phases/AUDIT.md) covers both packets. Unitary equivalence is not a positivity-preserving probability transport outside the precise congruence gate.

## Remaining session work

Other lanes classify affine-group even-order probabilities, harmonic-tail LRC families and Collatz one-reset entry phases. They will be added after independent audits. The same fractions are not automatically the same invariant: each connection must retain its source, map, predicate and lost information.
