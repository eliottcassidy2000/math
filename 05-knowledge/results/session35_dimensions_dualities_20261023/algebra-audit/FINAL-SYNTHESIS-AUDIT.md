# Final synthesis scope audit

2026-10-23. **PASS after corrections.** Scope: tournament/recursion, finite-count, generating-function, contextual-Walsh, and spectrum/zero terminology. Lie sections have a separate reviewer. This is a synthesis scope check against existing audited results, not a new replay of every session computation.

Reviewed:
- `session35-root/SYNTHESIS.md`.
- `session35-root/OBSERVABLE_COMPLETION_SYNTHESIS.md`, including section6.

## Corrections verified

1. The synthesis now separates completed Epstein reflection from the negative-even zeros of the **uncompleted** Epstein zeta. Completed Lambda is finite nonzero at those points; the earlier wording incorrectly assigned it the zeros.
2. Fixed-path restoration is now conditional on a whole-line flip actually changing path edges. Restricting to unused lines already preserves the path, while still defining an auxiliary move relation rather than automatically reproducing the original off-path-complement adjacency.
3. Deterministic contextual quotient size is distinguished from linear translate-span dimension; no unsupported positive-realization minimality claim remains.

## Checks against the available audits

- P7 nine-orbit counts21/42/63/63, unused-line degrees2/4/5, and1932 path/plane/state incidences agree with the independent Fano-path audit. Full cubes versus restricted unused-line cubes are expressly distinguished. Same-order H171 hostile prevents universalizing the criterion.
- The affine-reflection parity proof is valid in its stated scope. An extra independent direct check confirms the displayed Paley11 path has legal arcs, degree32, and fails the full mirror predicate. No arbitrary-order anti-automorphism classification is silently assumed.
- Spectral E/D closure, equivalence with adjacency cospectrality, the minimal order6 H29/25 collision, and odd-cycle difference4 agree with the independent algebra audits. No universal H discrimination after D is claimed.
- The positivity horizon now explicitly repeats t>=0. Augmentation census conclusions explicitly say “Through odd order7,” retaining the actual enumerated scope. Both final qualifiers were reread before the hash pin.
- THM341 correction, s6=903, formal-series nature, and conditioned fixed-k morphology match the independent good-cut audit. Recovery and inherited formulas are attributed rather than presented as new discoveries.
- Projective triangle-observer counts2^56 and1597 are not transferred to Hamiltonian counting. The new fixed-path caveat is compatible with the P7 restricted-move theorem.
- Conference theta collapse and adjacency/skew modulus correction agree with the separate zeta audit. No RH or metagraph-realizability conclusion is asserted.

The checked wording makes no claim that publication succeeded. Actual Git state remains a separate receipt matter; this audit does not certify a push or approval.

## SHA256 pin

SYNTHESIS.md:
385c90bf5cd29f7c5cc6b8fca253fb58b969ec1ac8d6913817fb12ba3cc84ff4

OBSERVABLE_COMPLETION_SYNTHESIS.md:
3a28b04509e58bf38b2bc7ab00af8ee25b543611cdd9f0603346535bc58fa139

## Stronger-result scope repin, 05:37 UTC

The preceding pin records the earlier checked revision and is retained as correction lineage. The new synthesis was compared with the completed independent reports `UNGRADED_AUDIT.md`, `SPLIT_REAL_AUDIT.md`, `CHARGE_DYNAMICS_AUDIT.md`, and `FOURTH_ORDER_AUDIT.md` in `session35-lie-arithmetic-audit`. This is a scope cross-check against those mathematical audits, not a duplicate proof audit or script replay.

- The ungraded minimum248 is explicitly over C, for finite-dimensional Lie algebras with nondegenerate invariant symmetric form and an injective embedding retaining E7⋉H57. Equality forces simple E8 and recovers the contact grading. No grading assumption or positivity of the ambient metric is silently restored.
- The general bound retains simple complex s and nontrivial irreducible symplectic V. The constructive standard-symplectic family has m>=1, so its Heisenberg dimensions are odd and at least3; the original phrase “every odd Heisenberg dimension” was flagged for that endpoint correction.
- The split-real conclusion is restricted to the specified split E7 input and its absolutely irreducible56-dimensional real module. It does not classify arbitrary real forms or marked embeddings.
- The stronger quartic exclusion is exactly the GL orbit of the standard Freudenthal quartic, without a preselected alternating form. The audited full stabilizer quotient C2 and odd order57 justify it. Other invariant quartics remain explicitly allowed.
- The38/5 elliptic paragraph specifies positive integer powers, distinct rational charges, arcsine equidistribution on(1,25), and reduced denominators10^m/2. It excludes a dyadic-density inference and VOA construction. These match the charge audit; negative powers would repeat charges, so the positive indexing is essential.

All earlier Epstein, fixed-path, deterministic-state, and finite-scope corrections remain in place. The observable-completion note remains a valid weaker graded special-case summary; it does not contradict the new ungraded theorem. Publication remains separately receipted. The parent's reported42 normal/optimized executions are not recertified by this scope-only review, and no published-hash update or successful push is inferred here.

The endpoint qualifier “at least3” was applied and reread. Final stronger-revision verdict: **PASS**.

### Current SHA256 pin (supersedes the earlier revision pin)

SYNTHESIS.md:
bed2336bcff1ab10e5d39fdc058a97c277cdc5645109234264d926f05ce3c545

OBSERVABLE_COMPLETION_SYNTHESIS.md:
3a28b04509e58bf38b2bc7ab00af8ee25b543611cdd9f0603346535bc58fa139

