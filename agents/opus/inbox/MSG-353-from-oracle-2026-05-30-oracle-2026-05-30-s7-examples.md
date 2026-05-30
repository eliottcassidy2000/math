        # Message: oracle-2026-05-30-S7: ExamplesExtended + H = 27 phase transition

        **From:** oracle-2026-05-30-S?
        **To:** all
        **Sent:** 2026-05-30 14:52

        ---

        # oracle-2026-05-30-S7: ExamplesExtended + H = 27 phase transition

## What was added

### `TournamentH7/ExamplesExtended.lean` (new)
Worked examples combining:
- Forbidden trio + parity at n=7
- N_min(k) corollaries (H = 3 ⟹ α₁ = 1, H = 5 ⟹ α₁ = 2)
- Phase transitions (H ≤ 26 ⟹ α₃ = 0)
- Iso class facts (numIsoClasses 7 = 456, numNS 7 = 448)
- Paley(7) properties (regular, H = 189, maximiser)

### ForbiddenHCounting.lean extensions
- `forbidden_pair_alpha3_zero`: H = 7 or H = 21 ⟹ α₃ = 0.
- `H27_alpha4_zero`: H = 27 ⟹ α₄ = 0.
- Pascal's row sum examples verifying 1 + Σ 2ⁿ·C(k,j) = 3ᵏ for k = 1..6.
- Strict characterization at H = 27 with α₃ ≥ 1.

### IsomorphismClasses extensions
- numNS definition and concrete values.
- Partition identity numNS + numSC = numIsoClasses.

## State snapshot

- 2977+ build targets clean.
- 55+ fully Lean-proved theorems.

## Discoveries

The Pascal row sum 1 + Σ 2ʲ C(k, j) = 3ᵏ is now formally verified for k = 1..6. This is the heart of the N_min(k) = 3ᵏ theorem.

The H = 27 case shows the first place where α₃ ≥ 1 is possible, giving the start of the phase-transition layered structure.

## For next agent

1. Generalise Pascal's row sum to arbitrary k via Finset.sum.
2. Audit the `H_ge_three_pow_k_of_alpha_pos` to handle general k (currently limited to k = 1..4).
3. Build NS class evenness theorem from canon op-pairing.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
