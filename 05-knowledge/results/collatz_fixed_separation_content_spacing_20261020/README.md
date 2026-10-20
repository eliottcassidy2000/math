# Fixed source separation: content and spacing

20 October 2026. **PROVED with the cited inputs; independently audited.**

For primitive lower mechanical words with p<q=p+r<2p, consider legal disjoint adjacent swaps moving two source ones whose minimum cyclic label separation is h, with 3<=h<=p/2.

- The full carry gcd satisfies g<=(34*4^h)^(4p^(2/3)), including prime powers and either sign of 2^q-3^p.
- For p>=64^3, an integral output must have h>p^(1/3)/8.
- For fixed h, integral outputs are confined to p<max(100,ceil(720(18+7h)/881)^3). This is a finite-weight reduction, not a claim that every residual is empty.

Start with [the proof](FIXED_SEPARATION_CONTENT_AND_SPACING.md), then [the audit](AUDIT_SUMMARY.md) and [reproduction instructions](REPRODUCE.md). [Sources](SOURCES.md) and [exact dependency hashes](DEPENDENCIES.json) identify every imported input. [Publication lineage](PUBLICATION_LINEAGE.json) records the proof-status promotion and copy transformations.

The independent finite controls cover 650,572 distinct ordered physical paths and 7,474 short-pair resultants, plus three special controls. They retain illegal moves, half-turn orientation duplication, a prime-cube gcd, and a resultant-converse hostile. They do not enumerate the entire fixed-h residual.

Unrestricted rank-four noncancellation, complete fixed-h residual emptiness, and arbitrary two-edit exclusion remain open.
