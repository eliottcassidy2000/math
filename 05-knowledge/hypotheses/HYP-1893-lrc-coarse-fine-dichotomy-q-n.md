---
id: HYP-1893
name: lrc-coarse-fine-dichotomy-q-n
status: OPEN
date: 2026-05-31
session: oracle-2026-05-31-S18
depends_on:
  - THM-369
  - THM-366
  - HYP-1855
---

# HYP-1893: LRC has a coarse/fine dichotomy at q = n

The rational time-line for the Lonely Runner problem splits at denominator
`q=n`.

For `2 <= q <= n`, loneliness at a reduced fraction `a/q` is controlled by pure
divisibility: if no speed is divisible by `q`, then `a/q` is lonely. This is the
coarse sieve formalized in Lean as THM-369.

For `q > n`, the nonzero residue band has positive width, so loneliness is no
longer a divisibility bit. It depends on whether every residue `v_i*a mod q`
lands in the safe band `[q/n, q - q/n]`.

## Prediction

Counterexamples must simultaneously cover every coarse denominator by
divisibility and survive every fine residue-band test. The main proof pressure
therefore lives at the `q=n` seam, where divisibility and band-membership
coincide and tight unit-endpoint examples concentrate.

## Artifacts

- `04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean`
- `01-canon/theorems/THM-369-lrc-sieve-lean-formalization.md`
- `07-reflections/lrc-denominator-sieve-lean-and-two-regimes-s18.md`
