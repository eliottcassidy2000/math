# Validation receipt

2026-10-20. **PASS.**

The complete checker ran normally and under optimization:

    python3 -B verify_lfsr_tournament_clocks.py > verification.out
    python3 -B -O verify_lfsr_tournament_clocks.py > verification.optimized.out
    cmp verification.out verification.optimized.out

The outputs agree byte-for-byte. All 42,860 explicit gates remain active under `-O`; the checker contains no `assert` statements.

Declared universes:

- 2,046 monic F2 polynomials, degrees 1 through 10, including singular transitions.
- 1,398,100 states traversed in that census, plus 524,288 states across the two degree-18 hostile graphs.
- 11,479 gcd/fixed-point cross-checks; 621 graph/profile classes in the census.
- Independent tournament isomorphism canonicalization through order 5, source-adjunction matrices through cutoff 5, full state-graph traversal through cutoff 4, and Pascal/Lucas orbit checks.
- Explicit noncyclic same-characteristic and same-cycle/different-Jordan-partition controls.

Independent proof review passed the factor/jet formulas, nilpotent product graph, profile inversion, degree-18 minimum, dual inverse-limit direction, and faithful Z2 action. The notes retain the pointed-torsor qualification: an orbit conjugate to the adding-one odometer is not thereby an additive subgroup of the F2-state module.

No Lean formalization is supplied. Exhaustive finite checks are distinguished from the all-degree proofs. Generic LFSR cycle structure and procyclic classification are classical and are attributed in the notes; historical novelty of the assembled results has not been established.
