# Validation and reproduction

The algebraic proofs are in RESET_DEBT_COMPOSITION.md. They establish arbitrary finite length; finite checks are controls, not the proof of that quantifier. No external literature theorem or novelty claim is a dependency. Inherited statements are explicitly separated from newly supplied proofs and the guarded examples.

Run from this directory:

    python3 verify_reset_debt.py > /tmp/reset-debt-normal.json
    python3 -O verify_reset_debt.py > /tmp/reset-debt-optimized.json
    cmp /tmp/reset-debt-normal.json /tmp/reset-debt-optimized.json
    cmp /tmp/reset-debt-normal.json receipt.json

Observed result: PASS, 321,690 explicit exception-based checks; normal and optimized outputs identical.

Declared finite universe:

- 625 pairs of rational affine parameters, including signed debts and slopes on both sides of one.
- Every exact odd-valuation word of length 1–4 with each exponent in 1–3: 120 words, independently testing all 41,370 positive odd residue representatives in their full odd-terminal moduli.
- Three source lifts for each word, including complete literal trajectory and inverse-image checks.
- 101 first-reset-bank words, with valuation-one run length 0–100 and the least final exponent satisfying the unit budget.
- Every schedule of 1–6 macros over {A,B,E}, plus four much longer schedules: 1,096 schedules altogether, maximum literal shortcut clock 1,200.
- Full arbitrary-finite-depth identities for repeated valuation 2 and repeated macro A checked through depth 100.

Hostile and positive controls:

- Actual `(1,2)` reset: 11 -> 17 -> 13 grows and fails the unit budget.
- Macro A `(1,2,3)` pays that growing reset; its first literal descent is at step 5, not its final step 6.
- `(1,1,3)`: 55 -> 83 -> 125 -> 47 descends on its entire cylinder while failing the unit budget; barrier 4 repairs the sufficient certificate.
- Formal arbitrary finite prefixes do not imply an infinite positive-integer orbit: repeated A's unique 2-adic source is 23/37.

Scope boundaries:

- The exact cylinder includes the odd-terminal bit, modulo twice the full denominator. Final integrality alone is not silently treated as an exact valuation address.
- ROOT-stopped prefixes are separately enforced where used.
- Reset count, odd-step count, and literal shortcut clock are different quantities.
- The guarded height banks do not prove the native signed-row `7X` estimate or arbitrary-source coverage.
- No claim of independent mathematical peer review is made by this staging artifact.
