# Independent audit: PASS

Auditor: separate runner/clock-review worker, session173, 2026-11-03.

## Scope and outcome

Reviewed the full PROOF.md, the exact targeted inherited Monster source, and check.py. No material mathematical issue found. The packet carefully distinguishes the chosen affine clock and chosen extension of the inherited marked cycle from a uniquely determined integer dynamical system. Its limit on tournament-spectrum and Monster conclusions is appropriate.

The full original checker passed normally and with `python -O`; outputs were byte-identical. Independent verification additionally used SymPy1.14.0 Schreier–Sims, rather than the author's exhaustive Cayley enumeration. `independent_audit.py` also passed normally and optimized with identical outputs. The original checker remains dependency-free; the optional independent audit requires SymPy.

## Independent results

- Exact affine CRT reductions: mod8 is -x+4; mod9 is5x+8.
- Affine conjugacy at fixed point34 and cycle type1^2 2^11 6^8 confirmed.
- Functional-composition words C=(A^3 B)^4=(0 5 2) and (A B^3)^3=identity confirmed.
- All five displayed conjugate supports confirmed. The elementary connected3-cycle-support lemma proof is sound, including the one-point-intersection reduction.
- Independently computed generated order181440 on9 states, proving A9 since both generators are even.
- Independently computed order362880 on18 states and on the actual72-state action. This corroborates C2×A9, rather than merely inferring the group from an equal order.
- Actual72-state orbit sizes9,9,18,18,18 confirmed.
- Negative controls: A alone has order6 and B alone order5. B violates congruence preservation mod3, so no integer polynomial can induce it modulo9. The literal six-state cycle fails the observer factor criterion exactly as stated.
- Reviewed the explicit source at `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`: the valuation reading and deliberate13↔4 identification match the inherited note.

## Reproduction

`python check.py`

`python -O check.py`

`python independent_audit.py`

`python -O independent_audit.py`

The independent script does not import or execute the author's script. It constructs the permutations separately and uses Schreier–Sims on9,18,72 points. No forbidden ledger was read or scanned.
