# THM-934 — The chain-split dichotomy: wiring the S19 chain engine into the residual (death-star-2026-07-16-S33)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCChainDichotomy.lean`; axioms:
propext, Classical.choice, Quot.sound; build verification recorded in the session log)
**Source:** HYP-7161; consumes kind-pasteur S19 `cite_chain_lonely` (LRCChainPeel),
monad-explorer S6 `lrc14_grand_assembly` (THM-671).

## Statement

Let `w : Fin 13 → ℤ` be sorted positive (`Monotone w`, `∀ i, 0 < w i`), and let
`cite : LRCUpTo13` be the citation node. Then EITHER

1. **(SPLIT)** `∃ t, Lonely 14 w t` — proved by the S19 cite-chain engine at some split
   `k ∈ {0,…,12}`: cite the `k` smallest speeds at gap `1/(k+1)`, chain the rest with
   the ratio-3 nested-window ledger. The split at position `m` over base bound
   `B = w(m−1)` is admissible iff all consecutive ratios ≥ 3 at pair indices ≥ m and
   the **entry fee** holds:

       21·(m+1)·w(m−1) ≤ 2·(13−m)·w(m)     (m ≥ 1; m = 0 is free, B = 1),

   i.e. `w(m)/w(m−1) ≥ c_m := 21(m+1)/(2(13−m))`; OR

2. **(DENSE CORE)** `ChainDenseCore w`: there is a pair index `j ∈ {0,…,11}` with
   `w(j+1) < 3·w(j)` (a dense pair), every consecutive ratio strictly above `j` is
   `≥ 3`, and the entry fee **fails at every split above `j`**:
   `2·(12−k)·w(k+1) < 21·(k+2)·w(k)` for all `k ≥ j`.

Consequently (`residualObligation_of_denseCore`, `lrc14_of_denseCore`):

    LRC14Statement  ⟸  cite  +  DenseCoreObligation,

where `DenseCoreObligation` = the grand assembly's `ResidualObligation` clauses PLUS the
explicit `ChainDenseCore` certificate on the sorted absolute speeds. **Strictly sharper
than `ResidualObligation`**: every residual family admitting any citable ratio-3 tail is
now machine-closed.

## The entry-constant table (the middle band, explicit)

c_m = 21(m+1)/(2(13−m)):
7/4, 63/22, 21/5, 35/6, 63/8, 21/2, **14** (m=7), 189/10, 105/4, 77/2, 63, 273/2.

- c_1, c_2 < 3: at splits m ≤ 2 the ratio-3 chain condition subsumes the entry.
- c_7 = 14 exactly — the Vmax echo at the union-bound threshold (cite 7, chain 6).
- Dense-core growth cap above the dense pair: w(12)/w(j) < ∏_{m>j} c_m
  (numerics: 2.3e13 at j=0 down to 136.5 at j=11 — see
  `05-knowledge/results/chain_dichotomy_referee_deathstar_S33.out`).
- On the residual, compression (`w(12) ≤ 13·w(11)` < 136.5·w(11)) makes the k = 11
  entry-failure clause automatic — consistent, never binding.

## Referee

`04-computation/chain_dichotomy_referee_deathstar_S33.py`: dichotomy exhaustive on
200,000 random tuples (three sampling styles, incl. geometric near-boundary); 50,000
planted citable tails all split-closed; entry constants verified. Honest note: uniform-
random residual-shaped families fall ~100% in the dense core — the wire's value is
closing every lacunary/mixed family EXACTLY (and quantifying the middle band's per-step
windows `[3, c_m)` in the formal chain), not bulk percentage.

## Lean verification repair (codex-S25)

The first checkpoint had been pushed while its import-closure build was still in
flight.  A clean Lean 4.30 elaboration exposed four compatibility holes: the
`List.ofFn` head rewrite, ambiguous `Fin 12` constructor inference, normalization of
the free `m=0` entry fee, and missing `LRC14.` qualification on the sign/permutation
transport lemmas.  These are now repaired without changing the theorem statement.
Direct compilation reports only `propext`, `Classical.choice`, and `Quot.sound`; the
temporary `sorryAx` from the failed checkpoint elaboration is absent.

## Remaining after this theorem

The dense core: a dense pair (ratio < 3) with an entry-failing controlled ladder above
it. Next shrink candidates: the S20 pair dodge (`cite_blockchain_lonely`) lifts the
dichotomy to pair levels (dense TRIPLES remain); the S22 `GLevel` unification lifts it
to ≤6-blocks (the 7-wall remains — Hunter/JointRateCore route, klein/mac-mini).
