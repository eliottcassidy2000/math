# Symmetry under continued operations

**Session114, 30 October 2026. In progress.** Scoped results and independently checked computations; no major conjecture closure or literature-priority claim.

## Verified results

[Aliquot continuation](aliquot/ALIQUOT_CONTINUATION.md) classifies the complete binary response under prime dilations by parity and the odd squarefree kernel. Adding one actual aliquot step splits a single such class into infinitely many distinguishable classes. Thus no bounded finite-valued sidecar repairs that quotient for the expanded alphabet of all prime dilations and aliquot steps. A fixed finite dilation alphabet has a separate exact finite response count; no infinite lower bound for one fixed finite alphabet is claimed.

The producer's 35,153 explicit gates pass normally and under Python optimization. Independent direct divisor-sum checks cover 30,000 source/multiplier pairs, and 42 ordered Mersenne-family distinguishing pairs. Proof review found and corrected a preliminary exponent typo before publication: sources are 2^p, successors 2^p−1. Scripts and receipts are included under `aliquot/` and `audit/`.

[One-coordinate Keller symmetry](jacobian/README.md) proves that one genuinely mixed-weight semi-invariant coordinate forces every polynomial Keller companion into an explicit triangular automorphism, without bounding the companion's degree. Finite-order certificates and the first aliasing shell are classified. Actual symmetry is a hypothesis; averaging does not manufacture it.

[Parent multiplicity](charts/PARENT_MULTIPLICITY.md) counts and phase-classifies all regular cyclic parents inside a composite wreath product. A singleton selected output still has exponentially many named parent owners. The result is not a lower bound for every compressed representation or unrestricted tournament enumeration.

Both packets passed independent proof review and a separately written [audit](audit/SYMMETRY_AUDIT.md). Producer tests cover 3,977 Keller gates and 38,707 parent gates, with matching normal/optimized runs. The independent checker adds 3,130 algebraic gates and direct permutation universes.

## Session direction

The remaining research lanes examine physical boundary owners in Lonely Runner and periodic Collatz return phases. Their results will be added only after verification.

The shared question is precise: which information must a quotient retain to preserve the target under its allowed continuations? Different problems have different carriers and operations. This is a research principle, not an asserted isomorphism between the conjectures.
