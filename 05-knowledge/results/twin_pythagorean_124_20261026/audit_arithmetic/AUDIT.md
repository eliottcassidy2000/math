# Independent arithmetic proof and implementation audit

**PASS.** The integrating root reviewed the final arithmetic proof
`README.md`, SHA256
`9d274221f4afb60e4ac6b87928231acd5e3d00605b90ffaab0d94c0204a8fb6a`,
and checker SHA256
`6056612a678848899ed0e06805990b418d9abf8664e2d58b0778ce1b0cbc1adb`.
The review imports no result from a reserved theorem stub.

The factorization `rho=n(m+n)` forces n1 when rho is prime. The remaining
radius is m-1 and parity forces m even, giving the twin predicate and its
inverse. The second prime-inradius branch follows from the two factorizations
of a prime as `n(m-n)`; the p2 parity exception is necessary. Its cardinality
is credited to Robbins'2006 Theorem2, whose official two-page paper was read
independently by the reviewer. The explicit branches and marked-exradius
selector are proved in the note rather than attributed to that count alone.

The inverse-fibre bound ensures that {70,142,5329} is the complete fibre of74,
not a search truncated at the target. The s(6p) diagram has its p>=5 guard;
the two-predecessor diagram includes p3 and has the unique perfect target6.
The audit checked these are distinct claims. No prime-inradius or fixed
Gaussian-norm condition alone implies twin radii.

Both abc maps retain their domains. The full centre triple is invertible on
the marked image, but relaxing the domain admits composite neighbours.
Its radical is `p(p+2)rad(p+1)>(p+1)^2`. Squaring the Euclid construction
adds exactly the disjoint support `rad((p+1)^2+1)`, with the stated height
change. The abc-derived inequality keeps its conjectural antecedent and
unknown constant. It is not a prime-counting lower bound.

The F1-compatible monoid-ring calculation retains prime labels. Arithmetic
evaluation kills `X_2+X_p-X_(p+2)`, while augmentation gives1, so augmentation
does not factor through arithmetic evaluation. Conversely equal unweighted
counts do not imply equal aliquot values. This is an exact obstruction for
these maps, not a claim that every F1 framework loses all addition. The
reviewer read Deitmar§1 for the stated monoid/base-extension interpretation.

The independent checker imports no producer code and enumerates the positive
odd/even leg pairs with hypotenuse<=1000, obtaining158 primitive triangles.
Area and semiperimeter fractions independently recover the two radii. It
checks complete prime-inradius fibres through19 (the proven parameter bound
places both branches inside the enumeration), a literal divisor-pair search
through5329 for fibre74, twin abc support identities, and map-loss controls.
It passes522 active gates, with identical normal and optimized receipts.
These finite checks supplement the all-parameter proof review.

The producer's81,771-gate receipt was replayed from its frozen script and
matched exactly. The final packet verifier replays both scripts normally
and under optimization; this report does not count repetitions as new tests.
