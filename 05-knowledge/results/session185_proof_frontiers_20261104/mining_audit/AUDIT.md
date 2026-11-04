# Independent audit: negative-constant ordered-context barrier

**ACCEPT after minor corrections, 2026-11-04.** Inspected source hashes are in `inspected.sha256`.

The all-order theorem is valid as stated. In particular:

- The general inequality c(X→B→Y)≤h(X)h(Y)c(B) follows by discarding only nonnegative terms in the exact physical interface law.
- Because c(B)<0, a nontransitive context with h≥3 overshoots below both desired constants. The proof now explicitly invokes insertion monotonicity to justify h≥3 from an induced directed triangle.
- With both contexts transitive, a marked left context subtracts at least F(B); if the left context is unmarked and the right is marked, its subtraction is at least G(B). The strict hypothesis min(F,G)>−2c−2 places both cases strictly below 3c+2. No context-size or marking assumption is hidden.
- Both unmarked contexts preserve c exactly, which is strictly above both targets when c≤−3.
- Empty contexts are covered. Fixed integer clock shifts change the slope equation but never the constant equation, so they cannot evade the exclusion. The result is not an obstruction to arbitrary replacements, internal edits, or inflation of a mixed-port module.

The original checker used a formal empty signature that was not a two-sided identity for all endpoint coordinates. This did not corrupt its final constant checks, but the author corrected it to an explicit empty-context identity. The audited version contains that correction.

Independent `audit.py` uses endpoint-resolved subset DP, not the author's permutation enumerator, to confirm the complete B7 and B8 seven-coordinate signatures. It verifies all11,250 ordered-context inequalities and digit exclusions with actual marked tournaments of orders0 through3 on each side. Eight whole context graphs, including empty contexts and nontransitive contexts on both sides, are independently counted and compared with all seven assembled coordinates.

Both the author's final verifier and the independent verifier pass identically in normal and optimized Python modes. The author's exhaustive module and feedback-arc checks were replayed; the stated order (5,6,4,2,3,0,1) has exactly the back arcs (4,5),(1,4),(1,2).

The finite controls do not prove the general quantifiers; the preceding inequality and exhaustive symbolic case split do. No novelty claim or global spectrum result is certified.
