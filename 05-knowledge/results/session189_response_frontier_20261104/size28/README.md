# Exact 28-name original-module boundary

Producer result, pending independent review: every prime nine-child substitution with eight C3 children and one nontransitive four-child is distinguished from its quotient-only converse. A 14-query child-uniform battery covers all 21,132 endpoint-zero quotient/placement cases; adding endpoint anchors gives 22 queries. All original-module family interiors preserve these child profiles. Together with the reviewed structural restrictions and the exact one-exterior-singleton argument, this gives a source-relative minimum of 29 old names.

Read PROOF.md for all physical-label, family, ambient and finite-census qualifications. This is not unrestricted named-response injectivity.

- bank.cpp implements the exact normalized run-state recurrence modulo eight.
- residues.tsv retains the full 12-query bank; selected_residues.tsv retains the two additional queries on all positions of the 13 affected quotient nodes.
- witnesses.tsv retains one separator for every quotient/placement case.
- physical_witnesses.jsonl.gz retains all three four-child realizations, their exact physical masks, both complete 30-row augmented graphs and residues.
- four_child_records.json retains all 64 labeled four-tournaments, their modules, profiles, and original-module family outputs.
- named_controls.json retains complete actual query values and every named Boolean coefficient for the three canonical four-children and their literal converses, including cubic and quartic terms.
- verify.py replays the producer calculation and all record checks. Run normally and with python3 -O; receipts must agree.
- DEPENDENCIES.json pins the reviewed inherited packets. RELEASE_FILES.json is the publication allowlist.

No software installation or startup helper is needed. The scripts use only the standard library and installed C++17 compiler.
