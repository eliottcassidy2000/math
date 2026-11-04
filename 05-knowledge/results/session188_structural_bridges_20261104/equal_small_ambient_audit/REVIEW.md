# Independent review: small equal-profile cyclic cores in arbitrary ambient tournaments

**PASS.** The theorem and its scope are correct: three fixed, nonempty child name sets carry one common complete ordered-path-cover profile of actual order 1 through 7; the target reverses their cyclic quotient and may change child interiors while preserving that profile. The one specified physical two-port query separates the cores in every fixed ambient tournament.

## Proof review

1. Actual tournament profiles lie in the same-parity Stirling cone by the cited odd-cycle expansion, with leading coefficient one and lower-weight sum H−1. The elementary bound H≤n! places them in the exact continuous simplex audited in `../profile_cone_bounded_audit`. This step applies to actual profiles, not to every arbitrary positive sequence.
2. The all-coordinate augmented-core identity from `covered_ambient` and the bounded-cone certificate give nonpositive source-minus-target differences in every coordinate h≥2. The separately audited positive-square identity gives a strictly negative h=1 difference. A scalar nonzero gap alone would not justify this ambient lift; coordinatewise domination is essential and is present here.
3. Giving both ports the old root's incidence to every exterior vertex makes the augmented core a module in both actual graphs. Its collapse is the same ambient quotient Q. The nonnegative run-substitution coefficients are independent of the interior, and c1=H(Q)>0. Thus the actual contrast is at most −H(Q)P<0.
4. Child order one, odd child orders, arbitrary fixed strong ambient quotients, and arbitrary profile-preserving target interiors are included. Different child profiles, orders above seven, changing exterior arcs, and general named-response reconstruction are not claimed.

## Independent executable checks

`verify.py` does not import producer code. It computes child profiles by an unordered set-partition recurrence with a canonical least element, multiplying by k! only afterward. It reconstructs the exact three-variable cover functional independently, verifies all augmented cover coordinates for the nine retained examples, and computes the scalar certificate directly as a sum of squares.

It reconstructs every source and target from their physical child and ambient adjacency rows, checks the actual exterior masks and the common augmented-module condition, and verifies each quantitative lower bound. `count.cpp` uses a backward first-vertex subset recurrence, independently of the producer's forward recurrence, to count all 72 full physical augmented tournaments, of order at most 18. Counts are bounded by 18!, safely within unsigned 64-bit range.

Normal and optimized Python audit runs are byte-identical. Normal and optimized producer replays reproduce both frozen receipts byte for byte. These finite controls supplement the universal substitution argument; they are not used to infer arbitrary ambient behavior from a census.

## Frozen files

Audit allowlist: `REVIEW.md`, `verify.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Source allowlist: `PROOF.md`, `verify.py`, `literal.cpp`, `normal.json`, `optimized.json`, `manifest.sha256` in `../equal_small_ambient`. The audit receipt pins all five substantive source files; `SHA256SUMS` pins the audit artifacts.
