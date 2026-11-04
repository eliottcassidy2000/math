# Independent synthesis audit

**PASS, 2026-11-04.** Reviewed the synthesis against all six source packets and final enumeration receipts. No substantive errors found; no edits to `NOTE.md` were required.

## Mathematical consistency

- The projective-response quotient is stated for one exact target mass, preserves only its declared weighted responses, and is not conflated with Hamilton-path enumeration. Its recurrence-class obstruction excludes P-recursiveness, not arbitrary recurrences.
- The one-feedback finite-core bound, asymmetric family, cyclic low spectrum, endpoint-flip identity, orphan correction, and finite-clone obstruction retain the hypotheses and limitations of their source proofs.
- The two infinite strong-witness lanes and carry decomposition are logically sufficient as stated. Neither the exponential strong family nor finite interface witnesses is claimed to settle them.
- Strong components are correctly restricted to nontrivial ones in the clone argument. The source clone proof's initially missing singleton qualifier has now been repaired and independently checked.
- Finite interface counts are not presented as new global attainment, and no component or wrapper obstruction is extended to arbitrary tournaments.

## Source and exclusion checks

Independently opened Busch's primary paper and verified Corollary 1, which gives strong minima 3,5,9,15,25 at orders 3 through 7 and increasing minima thereafter. Recomputed every labeled tournament spectrum through order 6 using the independent recursive endpoint implementation. Neither 7 nor 21 occurs. Combined with the strong minimum and SCC multiplicativity, this supports the stated all-order exclusions. The computation is reproduced by `audit.py` and recorded in the audit receipts.

Independently opened Royle's MathOverflow computation report and verified its strict bound below 80,557. The synthesis correctly says through 80,555 and does not infer nonattainment of 80,557. The Chan–Pak and Grinberg references remain cited external status sources; this review does not claim a new literature search or priority result.

## Numerical and link checks

The final n=200 state counts, transition counts, normal and optimized timing ratios, and n=150 traced allocation values match their source receipts. In particular the speed comparison changes sign, and projected-state compression does not imply a measured memory improvement. The distinction between one-coefficient and whole-prefix APIs is maintained.

All relative links resolve. `audit.py` pins the synthesis hash, checks the metrics and relative links, and recomputes the small spectra. Normal and optimized receipts are compared without timing fields and are byte-identical. This is a scoped proof and source review, not a claim to reprove every inherited theorem.
