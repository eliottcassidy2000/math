# Independent audit of detector.md and detector.py

Mathematical verdict: PASS. The axial shell bijection, two-orientation Collatz oracle, eight nearest-neighbor arrows, finite shell bound R<=6, carry equation and alphabet bounds, valuation congruences, and 1->4->2->1 convention are correct. The stated interface and communication caveats correctly distinguish digit-strip arithmetic from a one-site bounded-state geometric decoder.

Independent test: `audit_detector_independent.py` constructs coordinates by direct perimeter walking and tests nearest neighbors using the cube-coordinate norm rather than the submitted neighbor routine. It checks 120,601 sites, 40,000 Collatz arrows, boundary labels at radii up to 10^100, and arithmetic certificates up to 4096 bits. Both normal and optimized Python executions passed.

Two changes requested from the author:
1. Replace assert-based checks in the submitted harness with explicit runtime checks so optimized Python cannot print PASS without running the claimed comparisons; also compare the edge list against the literal eight expected arrows.
2. Explicitly require branch tags to agree between neighboring certificate cells (and the oddness condition at the starting marker). Otherwise locally legal mixed-multiplier cells need not certify either advertised affine relation. The scalar mathematics and executable function already hold the multiplier fixed.

Prior-work comparison is limited to the source identifiers and provenance supplied by the parent; the indicated prior NOTE.md is not present in this checkout. No independent historical novelty claim is made.

## Resolution

Both requested changes have now been applied and inspected. The revised submitted harness passed in normal and optimized Python; the independent optimized audit also passed against the revised module. Final verdict: PASS, with no outstanding correction requested by this audit.
