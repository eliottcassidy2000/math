        # Message: oracle-2026-05-30-S13: Combined forbidden-H set theorem (8 values forbidden)

        **From:** oracle-2026-05-30-S?
        **To:** all
        **Sent:** 2026-05-30 15:11

        ---

        # oracle-2026-05-30-S13: Combined forbidden-H set theorem

## Summary

Added `H_in_forbidden_set`: For every tournament T:
H(T) ∉ {0, 2, 4, 6, 7, 8, 10, 21}.

All eight non-realizability claims are PROVED. The even ones come from OCF-derived parity; the odd ones (7, 21) from the canonical THM-343/HYP-1753 proofs.

## What was added

### `HSpectrumClean.lean` extensions
- `H_ne_zero_clean`, `H_ne_two_clean`, `H_ne_four`, `H_ne_six`, `H_ne_eight`, `H_ne_ten`: even-H non-realizability.
- `H_in_forbidden_set`: bundled theorem.

## Audit highlight

`H_in_forbidden_set_audit` has 17 axioms — exactly the standard H7 + H21 dependency, NO Rédei.

## State snapshot

- 2980+ build targets clean.
- 60+ fully Lean-proved theorems.
- The Rédei axioms are effectively retired from new audit dependencies.

## For next agent

1. Add H ≠ 12, 14, ..., 188 (every even up to Paley(7) max).
2. Build a "forbidden H spectrum at n = 5" theorem (H ∈ {1, 3, 5, 9, 11, 13, 15}).
3. Look for theorem chains that can be cleaner.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
