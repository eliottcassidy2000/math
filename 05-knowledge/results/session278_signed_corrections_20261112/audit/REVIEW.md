# Independent review of new Session278 claims

**Status: PASS WITH MINOR CORRECTIONS (applied).** Review scope is the four
new Session278 mechanisms, not a re-audit of inherited Burnside formulas,
S25 physical join law, or JC/GMC literature. Nothing here promotes an
inherited pre-audit result to independently audited canon.

## Mathematical review

1. **Enumeration.** Independently derived the product defect, size-state
   exponent update, and Bernoulli remainder identity. The two-sided constants
   have the claimed direction. Verified the four-term resummed threshold:
   fewer than four residual cycle vertices cannot have four edges; four
   residual lengths below total24 cannot do so; (5,5,7,7) realizes the boundary.
   The fixed-depth obstruction has consistent singleton cross-exponent,
   factorial denominator and sign. The min-minus-gcd kernel is correctly
   restricted to the occupancy cone, not claimed positive semidefinite.
2. **Physical chains.** The noncommutative matching expansion follows the
   two-step recursion with its stated initialization. K lowers history by at
   least two, so the geometric identity terminates, while physical positivity
   proves the signs that an arbitrary alternating polynomial would not.
   Matching-number degree and strict first-coordinate remainder hold. The
   repeated-10 positive matrix has characteristic polynomial
   t^3-7t^2+4t-4 and discriminant -3376. The individual-layer generating
   functions have the stated poles and positive leading coefficient; hence
   every fixed-depth relative-error obstruction follows as written. This is
   a consequence of the inherited law, not universal coverage.
3. **Cross/Stein.** Product-rule derivation of (1)-(3) has correct signs and
   applies induction to actual polynomials, never multiplies quotient classes.
   Laguerre identity and multi-defect differential system are consistent.
   The signed Taylor brackets require finite displayed moments and do not
   claim convergence of the divergent perturbation expansion. Exact rational
   containment of the reduced intervals was added to the code.
4. **Coupled edits.** The joint inventory formula and partial-binomial identity
   prove the stated finite residual. Shared tails cancel before absolute
   values; triangle bounds and their equality condition are correct. Forced
   paths contract through ordered endpoints, including singleton/one-block
   boundaries. This is a paired certificate application of inherited/classical
   forced-arc expansions, not a new forest identity.

## Corrections made during review

* Enumeration's odd-index D_a notation now defines its even-size endpoint
  using the largest odd index <=n and separately sets U_0=1.
* Cross's extension from polynomials to general nonnegative functions now
  requires finite moments through the displayed Taylor remainder.
* Cross replay explicitly asserts exact containment, not just overlap, for
  the four finite interval comparisons.
* Edits now credits THM-4757's forced-linear-forest expansion in addition to
  the single-arc contraction law and OCF cluster background.

No theorem was refuted. No changes to mathematical numerical outputs were
required by these corrections.

## Controls and receipts

`independent.py` independently tests forced-path contraction on 1,284 cases
of arbitrary binary loopless digraphs of orders1..6, including inadmissible
forced sets, and the Bernoulli remainder with exact rational arithmetic on
420 cases. `mixed_defects.py` independently compares the operator recursion
with the formal generating formula for interacting weights2,3, a=3/2,
through degree6. See INDEPENDENT.json and MIXED.txt.

Author programs were replayed: edits output and physical output agree byte
for byte; cross output agrees byte for byte after repairs; enumeration agrees
semantically after ignoring wall-clock seconds. Author scopes remain 255
edit/DP controls, 1,554 physical word algebraic checks and258 independent
physical subset-DP checks; quartic/sextic polynomial divergence certificates;
and the enumeration controls described in its packet. These finite checks
are controls on the all-order arguments, not substitutes for their proofs.

The exact reviewed file hashes are in HASHES.sha256. Reproduction commands
are the four lane scripts plus the two audit scripts. No conjecture is marked
closed: all-target Chan–Pak realization, universal exact classification, and
planar JC remain outside the proved conclusions.
