# Second independent review: repaired consumer and new section 8

Date: 2026-10-21 UTC. Current verdict: PASS within the explicitly stated analytic and finite scopes. This is a separate review, not a retroactive claim that the initial audit covered new section 8.

## Consumer repair verified

Frozen checker SHA256: 44ae88a6d25a922348e668ec537d70f567bdb1b21120635c567fbd920ec4a992.

New copies executed independently under normal Python and -O both pass 10,374 always-active producer checks. All nonvolatile receipt content agrees; actual timing fields differ. All 16 genuine prefix receipts still pass independent strict verification. The exact formerly accepted 17-to-minus-1 alias is now rejected.

An additional independent schema campaign constructs 1,232 malformed receipts, across all 16 genuine representatives: invalid types for each required field, missing keys, nested labels, negative/out-of-range/bool/float/string labels in cycles/words/terminals, invalid direction, extra keys, and malformed/nonbinary rows. Every mutation returns False without an exception in both normal and -O runs. These tests use explicit RuntimeError checks, not assertions disabled by -O. Source inspection confirms validation occurs before indexing, so the observed negative-index and unchecked-direction defects are resolved. This is a bounded adversarial test suite, not a proof against every possible Python object.

The original finding in REVIEW.md is superseded as an outstanding issue but preserved as provenance. Artifacts: repaired-normal/, repaired-optimized/, repaired_audit.py, repaired_audit_results.json, repaired_schema.py, repaired_schema_results.json, current_source_pins.json.

## New section 8 independently reviewed

The score-singleton obstruction is correct in its now-explicit scope: strong self-converse T with Delta<=6, q>=5, m>=5.

A proper four-cycle in a singleton-score shore exists by strong pancyclicity. Its cross-shore arc sum is 4k, whereas universal balance at two per vertex would be 2m, impossible for odd m. Three exits or entries yield at least four black orbits by the already-reviewed half-support lemma. Complement strongness is supplied by the proper-shore deletion argument, and n>=25 is more than sufficient for four-prefix half-support.

The selected words contain four consecutive vertices of one full score. Marked-cycle words have isolated occurrences of this score, since exactly one shore carries it. Automorphisms preserve degree, hence selected orbits are not cut orbits. Self-converseness pairs degree s with n-1-s, preserves class cardinality, and therefore sends the singleton shore's score class to another singleton shore's score class. Reversal/converse mates also have a four-vertex same-score run and cannot be cut orbits. Thus these four black orbits are disjoint from the cut family's at least q-gcd(2,q) black orbits. Their total is at least eight for every q>=5, giving the contradiction. The entry dual preserves these score-pattern arguments.

No score class can have multiplicity one in shore units. Complementary score classes have equal multiplicity; at even n there is no integer self-complementary degree. Consequently:

- q=5: only the single central class of size five is possible.
- q=6: only a complementary pair of size three each is possible.
- q=7: either the single class of size seven, or a size-three central class flanked by a size-two complementary pair.
- q=8: either a size-four complementary pair, or two complementary pairs of size two each.

This is a necessary multiplicity restriction, not a census or an existence assertion. It does not affect the previously proved degree-56 bound and does not extend to m=1 or3. The cut-family black bound remains an explicitly inherited analytic input, as in the surrounding research packet.

Only the review directory was written. No Git or remote writes were performed.
