# Constructive maximal-prefix signed repair

Session 201. **PROVED, independently audited in the stated scopes.**

The main algorithm retains the longest possible prefix of a supplied physical
source recipe, replaces its terminal segment, and returns a section witness
at the requested target exactly when the target is1 or belongs to C1.
Two infinite append-hostile families have explicit optimal replacements.
Universal C1/C2/C3 target coverage remains open.

The target table is complete, not heuristic. Prefix optimality concerns the
supplied chronological recipe, not arbitrary graph edits or isomorphisms.
The general compiler optimizes the retained prefix, not appended length.
For example, --repair-word ZOOZ retains Z and constructs a valid nine-letter
suffix. The separately proved explicit family suffix ZZZZOJJ has seven
letters and is shortest for that retained prefix.

## Files and separate scopes

- PROOF.md and verify.py: constructive maximal-prefix replacement, target
  table, and the two optimal infinite-family repairs.
- DIGIT-ANSATZ.md and check_digit_ansatz.py: separately audited obstruction
  to coefficient-wise even-scaling rules under separate binary index changes.
- DYADIC-SPAN.md and check_dyadic_span.py: separately audited linear rank3
  of actual homogenized g-index vectors on every dyadic residue-class tail.
  This strengthens the preceding rule to full products of independent classes.
  Coupled or thinner target-dependent index subsets remain unrestricted.
- VALIDATION, DIGIT and DYADIC ordinary/optimized receipts agree in each pair.
- dependencies/: frozen independently audited physical input theorems.
- audit/: independent proofs, implementations, receipts, preserved original
  inputs, terminology correction and final promotion verification.

## Reproduction

From this directory:

    python3 verify.py
    python3 -O verify.py
    python3 verify.py --repair-word ZOOZ
    python3 verify.py --repair-word OZJJ --target 5
    python3 check_digit_ansatz.py
    python3 -O check_digit_ansatz.py
    python3 check_dyadic_span.py
    python3 -O check_dyadic_span.py
    python3 audit/check_independent.py
    python3 audit/check_digit_independent.py
    python3 audit/check_dyadic_independent.py

The independent scripts use their preserved audit inputs. Their reports retain
original scratch-directory names as provenance; the scripts run from these
published paths. See audit/README.md for the original review workflow.

No all-target occupancy assertion is made. A failed root query excludes a C1
section repair, not a possible C2 or C3 realization.
