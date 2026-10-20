# Validation receipt

2026-10-20. Proof audit and exact verification completed.

## Mathematical status

- The proof is self-contained after the definitions and the inherited prime-circulant isomorphism dictionary. It proves all-order dyadic-part lower bounds and unbounded necessary degree across actual primes.
- The dyadic source and target orbit counts are proved inside the note, rather than merely read from a finite table.
- The odd-repeat restriction is proved at the level of full antiperiodic words, so no coordinate or stabilizer quotient is implicit.
- A separate adversarial proof audit found no mathematical correction necessary. Its independently written verifier completed 249,374 explicit checks, reproduced the full quadratic census in ordinary binary order, and checked the witnesses and restriction arguments. See [the complete review](review/REVIEW.md). This receipt was updated after that audit; its archived source hash therefore records the earlier pending-review wording, while the mathematical proof source is unchanged.

## Completed exact checks

1. The standalone verifier `verify_degree_boundary.py` completed 268,536 explicit checks. Its normal and `python -O` receipts are byte-identical. It imports no sibling scratch packet or third-party dependency.
2. Full-word verification establishes the degree-two `p=13` witness and the seven-term degree-three `p=17` witness. Degrees are obtained by independent Möbius transforms in both the full and complement-quotient coordinate systems. All source words, equivariance, target landing, type bijectivity, labelled surjectivity, and labelled fiber counts are checked.
3. Independent input regeneration exactly matches the original 29-dimensional `m=8` quadratic basis and all target canonicalization entries.
4. Recompiled C++ exhaustive enumeration visits all 536,870,912 quadratics and exactly reproduces the original histogram. It finds no type bijection, and the optimal target-valued type coverage is fifteen of sixteen. The total number of target-valued quadratics is 22,136,841.
5. Independent direct coefficient enumeration through `m=6` matches the original small-case census, including 432 quadratic type bijections at `p=13`.
6. The elementary Walsh divisibility formula is checked against all Boolean functions in at most four variables. This is a finite control, not the proof of the all-order statement.
7. Exact dyadic target counts and weight sums are checked at lengths four, eight, and sixteen. Odd-repeat embedding checks cover dyadic lengths four and eight with odd repeat factors three and five.
8. Rebuilding the cubic search inputs reproduces their stored bytes. The seeded search's successful witness is independently validated without relying on the search.

## Deliberate scope limits

- The cubic search is not a full cubic census. Exact minimum degree three at `p=17` follows from the explicit witness plus the general lower bound, not from exhaustive cubic search.
- The seven-term presentation is offered for readability; the main theorem makes no sparsity-optimality claim.
- A further exploratory sparse search receipt is present but is not a dependency or promoted claim of this packet.
- The larger prime examples are lower bounds, not asserted minima.
- Composite and prime-power graph-isomorphism generalizations, arbitrary first-order nondefinability, and sharp circuit complexity are not claimed.

## Reproduction

Use the commands in `BOOLEAN_DEGREE_BOUNDARY.md`. The complete mathematical verifier uses only Python's standard library; the optional exhaustive and discovery programs require a compiler supporting C++17 and `__uint128_t`.
