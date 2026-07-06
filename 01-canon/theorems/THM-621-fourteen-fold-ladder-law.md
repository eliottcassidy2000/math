# THM-621: The fourteen-fold ladder law — hdich's surviving lifts, their closed form, and the self-referential witness

**Status:** VERIFIED (exact, all rows, witness-certified) + proof mechanism identified; corrects the set-identification in HYP-4098's prose
**Author:** mac-mini-2026-07-05-S52 (HYP-4099)

## The law

For the hdich lift families `W_k(r) = ({1..12}\{r}) ∪ {r + 13k}` at prime modulus 13:

1. **Survivor characterization.** For `r ≥ 7` (where `r` is the base's unique multiple of `r`), sieve survival forces `13k ≡ 0 (mod r)`, hence `k = r` in range: **the unique surviving lift is the fourteen-fold `r + 13r = 14r`** — the n = 14 shadow inside the n = 13 problem.
2. **Closed form.** `M({1..12}\{r} ∪ {14r}) = 14/(13(r+1))` for `r = 7..12` (exact: 7/52, 14/117, 7/65, 14/143, 7/78, 14/169), with rigidity gap `(13−r)/(13(r+1))`, minimized at `r = 12`: gap `1/169`.
3. **Witness law.** The binding witness lives on the `13(r+1)`-grid (verified: 15/104, 44/117, 29/130, 43/143, 71/156, 14/169). At `r = 12` the witness POSITION equals the VALUE: `t* = M = 14/169` — the same self-referential signature as the master construction `{1..12, 182}` at `t* = M = 14/183` one level up. The ladder is the deep-well family of the descended problem, with `13² = 169` playing `Φ₆(14) = 183`'s role.

## Correction to HYP-4098 (S51 prose)

S51's sweep DATA were correct (minima per r, zero violations, gap 1/169) but the prose misidentified the extremal SET as `{1..11, 25}` (the k = 1 lift — which is sieve-exposed with `M = 1/12` exactly). The extremal is `{1..11, 168}` (the 14·12 lift). All downstream constants (β = 1/169-scale for the window-margin leg) are unaffected; consumers should cite THIS file for the sets.

## Why (mechanism, for the assembly's per-r certificates)

Each row is one witness verification (the grid point above, 12 clearance checks, exact rationals) plus the covering direction (no better `t`: the binding-pair grid check — finite). The general-`r` proof: at `t = a/(13(r+1))` the residues `v·a mod 13(r+1)` place the twelve runners at distance ≥ 14 grid units with the binding pair achieving exactly 14 — the CRT interleaving of the 13-grid and the `(r+1)`-grid; the `14r` runner sits at `14r·a ≡ −14·(...)` closing the system. Formalization shape: six `decide` rows (the certificates above) — no general lemma needed for the assembly, though the CRT structure suggests one exists at all `n` (parked: `M({1..n−2}\{r} ∪ {(n+1)r}) = (n+1)/(n(r+1))`, the descended-deep-well family).

-> HYP-4097 (opus, mechanism), HYP-4098 (S51 sweep + this correction), THM-610 (hiding at 13·(r+1) ≥ 2n ✓ consistent), klein lonely_of_window_margin (consumer, β = 1/169), OPEN-Q-108.
