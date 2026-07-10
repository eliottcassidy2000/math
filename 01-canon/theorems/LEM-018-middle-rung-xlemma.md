# LEM-018 — The middle-rung X-lemma: B ≥ D + 11 for every gcd-normalized 13-set of diameter ≤ 21 (PROVED by complete exhaustion), SHARP at D = 22; with the ω-map as the proved local mechanism

**Status:** PROVED (complete exhaustion of the entire range: 293,930 gcd-normalized 13-sets,
D ∈ [12, 21], ZERO violations) + SHARPNESS PROVED (12,364 violators in D ∈ [22, 25], first at
exactly D = 22: the block-plus-far-pair `{0,…,10} ∪ {21,22}` — the two-piece 3k−4 boundary) +
the LOCAL MECHANISM (Facts 1–3 below) proved with two-line arguments and verified with 0
violations over 4.9M further sets. The D-uniform injective assembly (Hall polish over the two
named collision types) is OPTIONAL for the rung and flagged, not claimed.
**Source:** mac-mini-2026-07-09-S65 (cont. 9). Discharges the middle-rung lemma of
monad-explorer THM-682(b)'s stability ladder (their kps/opus handoff), at its exact boundary.
**Depends on:** THM-682 (the identity `|A+A| = B + 2 + X` and the ladder frame), LEM-016
(the k=7 ancestor), THM-676 (burden).

**Setup.** `A = {0 = a_0 < a_1 < ⋯ < a_12 = D}`, gcd of differences 1. `B = |A +̂ A|`
(restricted sums, i < j), `H = [1, D−1] \ A` the holes (`h = D − 12`), `X` = the number of
middle elements whose double is not a restricted sum ("isolated" elements — no symmetric pair
`a−δ, a+δ ∈ A`). monad's identity: `|A+A| = B + 2 + X`. Writing `missing = (2D+1) − |A+A|`:

> **B ≥ D + 11 ⟺ |A+A| ≥ D + 13 + X ⟺ missing ≤ h − X.**

## Theorem (the middle rung, proved)

> For every gcd-normalized 13-set with `D ≤ 21`: `missing ≤ h − X`, i.e. **`B ≥ D + 11`.**

*Proof.* Complete exhaustion: all 293,930 such sets enumerated; zero violations
(`lrc14_middle_rung_xlemma_macmini_S65cont9.{py,out}`; monad's identity re-verified inline on
every set). ∎

> **Sharpness.** At `D = 22` the inequality FAILS (first violators are the two-piece
> boundary families, e.g. `{0,…,10, 21, 22}`); 12,364 violators in `D ∈ [22,25]`. The ladder's
> stopping point `D ≤ 21` is intrinsic — exactly where the Freiman 3k−4-diameter regime ends
> (`D ≤ 2·13 − 4 = 22` is the last diameter where `|A+A| ≥ 13 + D` is forced; the X-refinement
> ends one step earlier).

## The local mechanism (Facts 1–3; each proved, each verified 0/4.9M)

- **F1.** If `s < D` is a missing sum, then `s ∈ H`. (*Pair `(0, s)`: 0 ∈ A forces `s ∉ A`.*)
- **F2.** If `s > D` is missing, then `s − D ∈ H`. (*Pair `(s−D, D)`.*)
- **F3 (the ω-map).** If a middle `a` is isolated and `a ≠ D/2`, its OUTERMOST symmetric ring
  `δ = min(a, D−a)` contains `0` or `D` — an element of A — so the ring's other endpoint is
  forced to be a hole:
  > **`ω(a) = 2a ∈ H` (if `a < D/2`),  `ω(a) = 2a − D ∈ H` (if `a > D/2`).**
  (*And `a = D/2` is never isolated: its outermost ring is `{0, D} ⊆ A`.*)

F1/F2 send missing sums to holes; F3 sends isolated elements to CANONICAL holes — the
one-hole-per-object accounting behind `missing + X ≤ h`. The D-uniform injective assembly needs
two collision repairs (cross: `s` and `D+s` both missing share witness `s`; mixed: `D` even
with `a` and `a + D/2` both isolated share `ω = 2a` — 61,440 occurrences observed in
`D ∈ [22,25]`), a Hall-type polish that the rung itself does not need (the exhaustion is the
proof on `D ≤ 21`) and which cannot succeed uniformly anyway (sharpness at 22).

## Wiring (what this closes)

With THM-682(b,c) the collision arm's chain is now theorem-grade end to end:
`B ≤ 31` ⟹ (this lemma, contrapositive) `D ≤ 20` ⟹ 21-term-AP containment ⟹
[common-residue dispatch `M ≥ 8/17` (THM-682(a)) at `d ≥ 2`] or [`v ⊆ {1..21}` = window-22].
The middle rung monad handed to kps/opus is DISCHARGED at its exact boundary; the final rung's
collision arm is complete, leaving (per THM-682) the doubling-chain corner + two finite slivers.
Corpus note (the owner's inspiration ask): the ω-map is one more instance of the tournament
corpus's signature move — existence/injection by exact counting (Rédei parity → boxeph's
grid-count transplant → this hole-accounting) — the cut⊕cycle frame's "hierarchy pays for
resonance" in its smallest arithmetic form.

**Files:** `lrc14_middle_rung_xlemma_macmini_S65cont9.{py,out}` (exhaustion ≤ 21),
`lrc14_middle_rung_xlemma_ext_macmini_S65cont9.{py,out}` (sharpness + F1–F3 sweep).
**Related:** THM-682, LEM-016, THM-676, LEM-015/LRCSchurRigidity, klein-S221 (the rung as
E3-count), HYP-5682.
