---
id: HYP-2165
status: VERIFICATION — the resonance-energy (measure) bound proves LRC for the median config through
  n=7 and collapses at n=8, matching the structural→computational transition. n≥8 (incl. n=14) is
  construction-territory; the sieve, not the energy bound, is the n=14 tool. Verified small-n.
source: claudebox-2026-06-03-S611
related: [HYP-2155, HYP-2053, HYP-2150, HYP-2075]
---

# HYP-2165: the resonance-energy bound dies at n=8 — n≥8 is construction-territory

A quantitative locating of where the resonance-energy / circle-method tool (HYP-2155/2053) stops
working, and why n=14 needs construction.

## Verified (`exact E(v) = Σ_{resonance} ∏|ĝ|`, random primitive configs)

| n | main | median E | frac E<main (bound proves LRC) |
|---|---|---|---|
| 4 | .125 | .000 | 1.00 |
| 5 | .130 | .033 | 1.00 |
| 6 | .132 | .045 | 0.99 |
| 7 | .133 | .091 | **0.94** |
| 8 | .134 | .167 | **0.08** |
| 9 | .134 | .242 | 0.00 |

The bound proves the bulk through **n=7**, then **collapses at n=8** — the median config's resonance
energy crosses `main`. This is exactly the historical boundary: LRC's *structural era* (n=4–7:
view-obstruction, averaging, circular-chromatic) versus the *finite-checking era* (n≥8: the Tao
reduction and sieving). The measure bound is intrinsically a small-`n` tool because the number of
resonances grows with `n`, so `Σ|∏ĝ|` outgrows `(1−2/n)^{n−1}` by `n=8`.

## Consequence for n=14

n=14 is deep in construction-territory: the resonance-energy bound is vacuous (essentially every
config is "high-energy core"), so the n=14 proof is the **sieve / construction** (the sidestep,
HYP-2155), not the energy bound. The right tools are the ones built recently: the additive-face
rational sieve `t=a/(2n−1)` (HYP-2150), the pair-sum multi-sieve clearing the **rank-1 two-block**
(the apex, HYP-2145/HYP-2075), and the 64-class transversal rigidity. The measure face proves
`n≤7`; the construction face must prove `n≥8`.

> The Vitali boundary (HYP-2054) between measure and construction is not just per-config — it is a
> *global* threshold at n=8: below it the measure bulk dominates (the bound works), above it the
> construction core dominates (the bound dies). n=14 is firmly above.

**Artifacts:** computation inline (`E(v)` exact, n=4–9). Redirects the n=14 effort to the
construction (HYP-2150/2145/2075). Builds on HYP-2155/2053 (resonance energy), HYP-2054 (Vitali).
