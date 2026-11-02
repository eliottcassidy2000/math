# Session137 — realizable odd-cycle packing carriers

Read `THEORY.md` for proofs, scope, the complete 29-type list through component cycle-budget 15, sharp all-order packing cost, and the 31 resource tradeoff.

## Reproduce

```bash
python3 -B verify.py > normal.json
python3 -B -O verify.py > optimized.json
cmp normal.json optimized.json
```

Standard library only; normal runtime about half a minute per run in the research environment. All checks raise exceptions explicitly; none depend on removable Python assertions. Output is deterministic.

`normal.json` contains all 29 positive source witnesses with literal cycles, the complete small-source census,300 bounded-budget values with their exact minimum costs and source expressions, per-cost value lists, interval-family controls, the 31 Pareto obstruction, and the same-carrier/different-module-context hostile. `SOURCE_PINS.json` pins copied inherited source bytes to Session136 commit `bb7cb2d1f61a14632389814469317131d897b640`.

The complete fixed-budget result is not a universal tournament spectrum theorem. Full spectrum completeness remains OPEN. Inherited SCC/forest classification, Moon's construction and cycle minimum, OCF, tribonacci counts, omissions7/21, and Session136 module obstructions are credited separately from present consequences.

## Cross-lane determinant audit

`DETERMINANT_SIDE_AUDIT.md` independently accepts the root's global nonempty inequality `H(T)>=det(I+A)+|det(A)|`, supplies the exact equality boundary (transitive or directed triangle), and records the indispensable empty-sidecar exception. Run `determinant_audit.py` normally and under `-O`; it checks7,412 complete classes through8 plus470 literal parent/cut recurrences. Both final outputs match; see the side-audit note and `VERIFICATION.md`.
