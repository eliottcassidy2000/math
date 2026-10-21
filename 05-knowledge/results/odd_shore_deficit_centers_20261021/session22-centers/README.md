# Deficit-six centers: odd-shore reduction

The analytic result is **deg R6 <=47**, improving the inherited4096 bound (and Session21's4095 arithmetic trim). For every marked Hamilton cycle, its exact rotation index q and odd shore size m satisfy either q=1 with n<=47, or 2<=q<=8 with m<=7. This is a necessary structural reduction, not a full center classification.

New mechanism: a power-of-two cycle within one rotation shore has strong complement. Another odd shore cannot be exactly half adjacent to every cycle vertex in aggregate; exact edge counting supplies a strict majority attachment and a mate-free prefix clique. Section6 of PROOF.md extends this to a quadratic self-converse bounded-defect order theorem. Section8 adds the five-vertex score-run obstruction: q5/q7 residuals must be regular, q6/q8 with shore size5 or7 are impossible, and the regular49 case is excluded.

Files:
- PROOF.md: complete proof, dependencies, scope, and generalization.
- check.py: dependency-free exact finite controls and source-bound witness verification.
- receipt_normal.json and receipt.json: normal and optimized fresh runs, including actual clocks and literal prefix receipts.
- check.out and check_optimized.out: short run summaries.
- source_pins.json: inherited source hashes and checkout revision.

Reproduce with `python3 -B check.py` and `python3 -B -O check.py` from this directory. Results agree after omitting real clock/elapsed fields. Both runs perform10,374 always-active checks. Small exact controls authenticate C3, R5, and Paley7 only; no previously missing all-class audit is represented as recovered.

The inherited degree-three rational series consequently has safe order18 recurrence onset n>=60. Minimal order18 and the asymptotic are inherited algebra, not new credit here. Literal metagraph incidences still require source/word/endpoint/target data.
