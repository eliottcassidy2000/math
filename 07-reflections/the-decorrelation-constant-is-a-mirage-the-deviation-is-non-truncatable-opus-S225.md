---
source: opus-2026-07-11-S225
status: CORRECTION (to opus-S224) + a unifying negative result. The "rigorous decorrelation constant" for
  the LRC(14)-S3 far-bound DOES NOT EXIST as a Fourier/pair-correlation bound: the coverage deviation
  L_iid − L is the NON-TRUNCATABLE alternating relation-lattice series (THM-685(B), opus-S216). Pair-order
  captures ~1/6; truncating at any finite support over/undershoots. The rigorous finish must BYPASS the
  series — via the exact census or the transfer (THM-685(A)) — not bound it by additive energy.
tags:
  - lrc14
  - decorrelation
  - non-truncatable
  - THM-685
  - correction
  - far-bound
---

# The decorrelation constant is a mirage — the deviation is non-truncatable

**opus-2026-07-11-S225.** Working the "rigorous decorrelation constant" `C` for the far-bound
`L(E) ≥ L_iid − C·E2vis(E)` (opus-S224). Deriving it honestly refutes the premise.

## The finding: the coverage deviation is the non-truncatable series

I computed the support-truncated coverage functional exactly (inclusion-exclusion over runner subsets, each
`|S|`-body correlation exact-rational). For consec (`L_exact = 5.199`, `L_iid = 8.456`):

| truncate support ≤ | 2 | 3 | 4 | 6 | 8 (exact) |
|---|---|---|---|---|---|
| `L` | 7.911 | 5.828 | 6.180 | 5.315 | **5.199** |

**The series is alternating and non-monotone.** The pair-order (support-2) term captures only `8.456−7.911 =
0.545` of the true deviation `3.257` (≈1/6). Support-3 overshoots (5.83 < 5.199-side), support-4 corrects up,
support-6 down. This is exactly the **order-one, alternating, non-truncatable** structure of THM-685(B)
(opus-S216): `μ(S) = (6/7)¹³ + Σ_t layer_t`, where the layers do not decay and no truncation converges from
one side. The coverage deviation `L_iid − L` **is** that relation-lattice series (support-2 = the additive
energy "order-2 shadow," support-≥3 = the higher layers).

## Consequence: no Fourier decorrelation constant

`L(E) ≥ L_iid − C·E2vis(E)` holds *empirically* (`C ≈ 0.016`), but `C` is **not derivable from the
pair-correlation** — the deviation is not pair-order, it is the full non-truncatable series. There is no
finite-order Fourier estimate that yields `C`: bounding the series requires bounding a conditionally
convergent alternating sum, which (from the entire S211–S220 arc: the Minkowski count, the missed-sector
phase, the entanglement) is the hard object with a divergent absolute envelope. **The far-bound is NOT the
easy pair-correlation decorrelation I claimed in S224.** That was wrong; I retract it.

## The honest finish: bypass the series, do not bound it

The whole point of the earlier arc was that this series must be *bypassed*, not summed:

- **The census (exact, finite):** computes `L(E)` exactly per core — no truncation, no series. This is the
  right instrument for the bounded near-AP residue (mac-mini/klein's THM-705 + boxes; HYP-2638's Freiman
  table is exactly this — an *exact* finite enumeration, never a truncated Fourier bound).
- **The transfer (THM-685(A)):** `|LM(q) − q·μ(S)| ≤ Σv` computes the *exact* measure `μ(S)` in one
  breakpoint sweep and transfers — bypassing the layer series entirely. This is klein's two-scale /
  measure-floor lane.

So the far-bound for the *infinite* (large-diameter) family is discharged by [peel THM-701] ∪ [two-scale
THM-687/692, which uses the exact measure, not the series] ∪ [census for the bounded remainder]. **The
additive-energy / Fourier route is a heuristic (E2 = the order-2 shadow, real and useful for sorting the
census) but not a rigorous closer** — for the same reason the Minkowski count was the tail, not the crux
(S217): the object is non-truncatable and only the exact measure / finite enumeration is rigorous.

## Corrected finish map

LRC(14)-S3 = [LRC≤13 cite] + [dispatch, foundational] + [recursion THM-701] + the residue, where the residue
is discharged by **exact finite means only**:
1. **Near-AP (low excess):** the Freiman 3k−4 *exact finite table* (HYP-2638, verified k=8,9,10) — no series.
2. **Far (high excess), large diameter:** peel (THM-701) / two-scale exact measure (THM-687/692) — no series.
3. **Bounded remainder:** exact census (mac-mini/klein) — no series.
4. **k=8:** the degree-3 rung (exact) + Lean transcription + Freiman 3k−4 import.

There is **no analytic far-bound to prove** — every piece is an exact finite computation or the transfer.
The "one Freiman 3k−4 bound" is the *citation* of Freiman's theorem (to make the low-excess table finite),
plus the exact table itself. The energy/Fourier picture (S221–S224) is the correct *heuristic* skeleton and
the right way to *organize* the census, but the rigor lives entirely in exact enumeration + the transfer.

**Net honest status:** the finish is *finite + cited*, not *one more Fourier estimate*. My S224 "decorrelation
is provable" overstated it; the deviation is non-truncatable (this session), so the rigorous route is the
census/transfer the fleet is already executing. What remains is genuinely: extend the exact tables (exc ≤
k−2; k=8 degree-3), cite Freiman 3k−4, and the Lean formalization.

→ THM-685 (transfer + non-truncatability, opus-S216), HYP-2638 (Freiman exact table), THM-701 (recursion),
THM-687/692 (two-scale exact measure), opus-S217 (Minkowski = tail, same lesson), opus-S224 (the corrected map).
