---
source: klein-2026-07-23-S410
status: cross-cutting adjudication of the owner's gap-vs-surplus fork for the snippet's ">1/25". Computational,
  robust. Two parallel per-reading investigations (subagents) run alongside; this is the reading-independent
  spine. Conclusion: the GAP reading (13 speeds, wider-gap, incomplete) is correct; the SURPLUS reading
  (24 speeds, clears the tight floor) is untenable.
tags: [lrc14, snippet, gap-vs-surplus, variational, wider-gap, lossy-factor, fork]
---

# Gap vs Surplus — a universal lossy factor decides it (and both readings share one wall)

**klein-2026-07-23-S410.** Owner: pursue ">1/25" as a GAP (reading a: wider-gap `>1/(2·13−1)`, incomplete)
AND as a SURPLUS (reading b: clears the tight floor `1/25 = ML({1..24})` for 24 speeds, could be full) in
parallel. Here is the reading-independent spine; the two detailed investigations are separate.

## The universal lossy factor (the key fact)
The snippet is `X = ∫g dμ` (kps-S132): a single-amplitude Riesz-product variational lower bound on the gap.
For the tight AP `{1..n}` at `τ*=1/(n+1)`, `best_a ∫g dμ / gap` is a **near-universal constant ≈ 0.66**:
`n=6: 0.668, 9: 0.666, 12: 0.661, 13: 0.660, 16: 0.657, 20: 0.654, 24: 0.652`. So
`∫g dμ ≈ 0.66/(n+1)` — a uniform **fraction ~2/3 of the tight floor**. Consequences:
- It **clears `1/(2n−1)` for all `n≥6`** (`0.66/(n+1) > 1/(2n−1) ⟺ n>5.5`): a genuine wider-gap improvement over
  the union bound `1/(2n)`, at value `~0.66/(n+1)`.
- It **never reaches `1/(n+1)`** (the conjecture): the ~34% loss is structural (a band-limited `μ` is not a point
  mass, so `∫g dμ < max g = gap` strictly). So the method is **INCOMPLETE for the full conjecture under EITHER
  reading** — the same wall, whichever floor you name.

## The fork is broken — GAP wins, SURPLUS is untenable
The two readings are NOT symmetric under this bound:
- **GAP (13 speeds):** `∫g dμ({1..13}) ≈ 0.047 > 1/25 = 1/(2·13−1)`. The snippet genuinely clears the wider-gap
  threshold. Reading (a) is **coherent and correct** as a (lossy, incomplete) wider-gap bound.
- **SURPLUS (24 speeds):** the tight `{1..24}` has `gap = 1/25` exactly, but `∫g dμ({1..24}) ≈ 0.026 < 1/25`.
  The variational method **cannot even reach** the 24-speed tight floor — there is no surplus over it to certify.
  A "surplus clearing `1/25`" can only exist for LOOSE 24-configs (gap `≫1/25`), which are not the hard cases,
  and it collapses on the tight ones. Reading (b) is **untenable** for the actual crux.

## Arithmetic corroboration (13, not 24)
The snippet's integers carry `13²` (`2974400 = 2⁶·5²·11·13²`) and `3·Σ₁¹³k² = 2457`, with **zero** 24-signals
(no factor 24/25/49, no `Σ₁²⁴k² = 4900` / `9800` / `14700`). And `1/25 = 1/(2·13−1)` is already the natural
13-speed wider-gap threshold, so the coincidence `1/25 = ML({1..24})` is not needed to explain it. Reconstruction
agrees: `∫g dμ({1..13}) ≈ X`, `∫g dμ({1..24}) ≈ 0.026 ≠ X`.

## Verdict
`">1/25"` is a **GAP** (wider-gap bound, 13 speeds), not a surplus. As a proof step it is SOUND but LOSSY and
therefore **INCOMPLETE for the full conjecture** — it reaches `~0.66/(n+1)`, above `1/(2n−1)` but below `1/(n+1)`.
If the outside LLM claimed a *full* LRC proof and this is a load-bearing step read as a surplus that closes the
conjecture, that is the **flaw** (the ~34% variational loss cannot be recovered by this functional). The escape —
if any — is a HIGH-DEGREE concentrator (Fejér/Jackson climbs toward the gap, kps-S132); the snippet's low-amplitude
Riesz instance is structurally a wider-gap bound, not a conjecture-closing one. → kps-S132 (∫g dμ), klein-S407/S409,
opus-S4b (24-floor), THM-518 (exact-value wall).
