# THM-938 — The pair lift of the chain-split dichotomy (death-star-2026-07-16-S34)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCPairLiftDichotomy.lean`;
verify the build report in the session log). Lifts THM-937 through kps-S20's pair
dodge (`LRCPairBlock.cite_blockchain_lonely`). Source: HYP-7171.

## Statement

`lonely_or_tripleCore`: every sorted positive 13-family is lonely — via a THM-937
singles split, or via the PAIR-AT-THE-DENSE-POSITION split (cite the j smallest,
dodge the dense pair `(w(j), w(j+1))` as one BLevel pair, chain the rest as
singles) — or lies in **TripleDenseCore**: the THM-937 certificate PLUS

    w(j+2) < 21·w(j+1)                    (the compressed triple), or
    (13−j)·w(j) < 13·(j+1)·w(j−1)         (pair-entry failure, j ≥ 1).

Wire: `lrc14_of_tripleCore (cite) (htriple : TripleCoreObligation) : LRC14Statement`
— strictly sharper than THM-937's `DenseCoreObligation`.

## Key constants

- Pair entry over cited base B: `w(j)·2δ ≥ 13/7 ⟺ 13(j+1)·B ≤ (13−j)·w(j)` —
  entry constant `13(j+1)/(13−j)`, CHEAPER than a single's `21(j+1)/(2(13−j))` by
  the factor 21/13; but the pair's output window `1/(14·w')` is 7× tighter than a
  single's `1/(2·w')`: transition fee pair→single = 21 (vs single→single = 3).
- At j = 0 the pair entry is FREE (`w₀ ≥ 1`); at j = 11 the singles tail is EMPTY,
  so every top-pair dense core with entry closes outright.
- `bChainOK_of_chainOK`: the S19 singles ledger embeds in the S20 block ledger —
  the S33 arithmetic reused verbatim.

## Referee

`04-computation/chain_pairlift_trap_referee_deathstar_S34.py` (out in
05-knowledge/results/): lifted-dichotomy exhaustiveness PASS on 200k tuples;
20k planted dodge families all close. Honest measure: in a uniform-random mix the
dodge closes ~0.2% of the former dense core — the value is structural (every
21×-jump-after-the-pair family and every entry-passing top-pair family, exactly),
plus the automaton constants {3, 26/7, 21, 26} now formal.

## Remaining after this theorem

TripleDenseCore: a dense pair whose successor stays within 21× (or whose entry
fails). Next lift: full BLevel search (pairs anywhere + ≤6-blocks, S22 GLevel) —
the 7-wall; and the trap-supplied B5 route (THM-939) on what remains.
