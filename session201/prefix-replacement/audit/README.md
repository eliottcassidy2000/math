# Independent maximal-prefix repair audit

PASS in the exact C1-conditional scope. See AUDIT.md for the proof review, computational universe, hashes, and the distinction between maximal prefix retention and shortest appended suffix.

- `frozen-inputs/`: byte-preserved original proof, README, verifier, receipts, and dependencies.
- `check_independent.py`: independent positive-coordinate and literal tournament verifier.
- `independent-normal.json` and `independent-optimized.json`: matching independent receipts.
- `author-replay-normal.json` and `author-replay-optimized.json`: matching replays of the frozen author verifier.

The original proof and implementation were preserved. No GitHub writes or Git index changes were made.

A separate stopping-rule sidecar is reviewed in `DIGIT-AUDIT.md`. Its independent verifier and matching receipts are `check_digit_independent.py` and `digit-independent-*.json`; author replays are `digit-author-replay-*.json`. The original sidecar inputs are also archived under `frozen-inputs/`. This addendum does not change the main audit's scope.

The separate actual dyadic-class strengthening is reviewed in `DYADIC-AUDIT.md`. Its checker is `check_dyadic_independent.py`, with matching `dyadic-independent-*.json` receipts and `dyadic-author-replay-*.json` author replays. This proves rank three of actual homogenized index vectors in every dyadic tail and supports only the independent-product extension of the stopping rule.

Final promoted bytes and portable release-tree execution are verified in `PROMOTION-VERIFICATION.md` and `PROMOTION-VERIFICATION.json`. The self-contained `check_promotion.py` checks exact source deltas and replays all three author and all three independent checkers in both normal and optimized Python. Final sources are archived in `promoted-inputs/`; run `python3 audit/check_promotion.py` from the published packet root.
