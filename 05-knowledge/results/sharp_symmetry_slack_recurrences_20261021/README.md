# Sharp symmetry-slack recurrences

Research date: 2026-10-21. Status: PROVED in the stated scopes with independent mathematical audit PASS and FINITE-EXACT controls. Read the [proof](session23-sharp-slack/SHARP_POLE_ORDER.md), [producer overview](session23-sharp-slack/README.md), and [independent audit](session23-sharp-slack-review/AUDIT.md).

For every nonempty exact-slack layer the proof identifies the exact ordinary, self-converse and converse-merged reduced denominators and minimal eventual recurrence orders. The upper bound uses actual exceptional automorphism orbits; injective nested-wrapper families give lower bounds. Sufficient onset bounds are not earliest-onset claims, and the denominator alone does not supply the numerator or all required initial coefficients.

## Reproduction

From this packet directory:

    python3 -B session23-sharp-slack/verify_sharp_poles.py
    python3 -B -O session23-sharp-slack/verify_sharp_poles.py
    python3 -B session23-sharp-slack-review/check_sharp_slack.py
    python3 -B -O session23-sharp-slack-review/check_sharp_slack.py

The producer rewrites sharp_pole_rows.json; use a disposable copy to preserve archived bytes. Its 1,027 checks and the independent review's 17,182 checks are distinct finite universes, not formal verification of the unbounded theorem. All source and review files are copied unchanged; original source pins and historical receipts retain their scope. Publication is established by the enclosing verified checkpoint.
