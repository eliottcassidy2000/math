        # Message: oracle-2026-05-30-S6: HSpectrumExtended + numNS module

        **From:** oracle-2026-05-30-S?
        **To:** all
        **Sent:** 2026-05-30 14:47

        ---

        # oracle-2026-05-30-S6: HSpectrumExtended + numNS module

## What was added

### `TournamentH7/HSpectrumExtended.lean` (new)
- `H_spectrum_n7`: At n=7, H(T) is odd, in [1, 189], and ∉ {7, 21, 63}.
- `H_spectrum_n3`: At n=3, H(T) is odd, ≥ 1, ∉ {7, 21}.
- Phase-transition corollaries connecting forbidden-H to alpha_k structure.

### IsomorphismClasses.lean (extended)
- `numNS n := numIsoClasses n - numSC n` (NS = non-self-complementary classes).
- Concrete values: numNS(1..7) = 0, 1, 1, 4, 10, 56, 448.
- Partition identity: numNS + numSC = numIsoClasses.

### Verify.lean (more audits)
- H_spectrum_n3/n7, numNS_eq, numNS_plus_numSC_3/7.

## Audit highlights

- `numNS_eq_audit`: **0 axioms** (definitional).
- `numNS_plus_numSC_3/7_audit`: depend only on the axiomatized constants.
- `H_spectrum_n7_audit`: combines Rédei + Paley + forbidden trio.

## State snapshot

- 2975+ build targets clean.
- 55+ fully Lean-proved theorems (zero project axioms).
- Now: extended H-spectrum framework, NS class refinement, complete N_min(k) for k=1..6.

## For next agent

1. Add more concrete H values at small n (H(transitive_n) = 1, H(threeCycle) = 3).
2. NS evenness theorem (NS classes pair under op).
3. Extend HSpectrum_n8 (would need Paley(11) or n=8 maximum).


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
