# Repair through 24 interleaved expansion episodes

Read INTERLEAVED_HEIGHT_REPAIR.md first. At the first explicit F/H boundary where the cumulative multiplier falls below one, a lawful word with at most 24 expansion episodes either produces actual descent below its starting value within 92 additional Terras steps, or returns to a cycle with a distinct smaller periodic obstruction. The total word length is unrestricted.

The proof includes an exact episode-dependent carry bound, source guards, 12,798 all-lift offset certificates, the K^4F=V normalization, and a constructive reason no episode-independent carry bound can replace it. Arbitrary pointwise bank entry remains open. The 25-episode offset obstruction discussed here concerns only the inherited K^4 prefix; later guard refinement is a separate result.

From this directory run:

    python verify_interleaved_height_repair.py
    python -O verify_interleaved_height_repair.py
    python independent-review/check_interleaved.py
    python -O independent-review/check_interleaved.py

Writer outputs match verification-receipt.json. Independent outputs match independent-review/receipt.json. The full separate PASS review, checker, receipt, and pinned review manifest are included. These are offset and identity controls, not an integer-source census.
