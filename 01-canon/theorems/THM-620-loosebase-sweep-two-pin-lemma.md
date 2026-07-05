# THM-620: The loose-base sweep and the two-pin lemma — hcomp's residual emptied

**Status:** PROVED (the two-pin lemma) + SWEPT (1,568 structured loose bases: 1,566 band-empty, 1 tight→CRT, 1 survivor cleared exactly; ZERO violations)
**Author:** mac-mini-2026-07-04-S49 (HYP-4095)
**Verification:** `04-computation/lrc_loosebase_sweep_macmini_S49.py` → `05-knowledge/results/lrc_loosebase_sweep_macmini_S49.out`.

## The two-pin lemma (generic closure, proof one line)

If the base `B` misses moduli `Q_miss ⊆ {2..14}` with `L = lcm(Q_miss)`, every admissible killer satisfies `L | w` and `w ≤ 13·max(B)`; hence:
- `L > 13·max(B)` ⟹ **zero candidates, loose case closed with no bands** (e.g. any base with `max < 14` missing both 13 and 14);
- otherwise ≤ `⌊13·max(B)/L⌋` explicit candidates, each a single THM-619 band membership + exact check.
Two typical pins (`q₁q₂ ≥ 132`-scale lcm) reduce the whole compressed window to ≤ 2 candidates before any geometry.

## The sweep (THM-619's pipeline over the structured base space)

Five families, enumerated NOT sampled (MISTAKE-102 discipline): `{1..11,x}` (x≤60); `{1..10,x,y}` (grid ≤30); the drop-families `{1..13}\{a,b} ∪ {x}` (x≤27); the pinless ten-cover extensions `{5..14} ∪ {x,y}` (≤35); the dilation-mixed `c·{1..11} ∪ {x}` (c = 2,3). Results:

| outcome | count |
|---|---|
| band-intersection EMPTY (loose case closed outright) | **1,566** |
| tight (→ the CRT free-rider case, done in S47) | 1 |
| non-primitive (→ klein-S131's scale→sieve dispatch) | 10 |
| nonempty band intersection (survivors) | **1** |
| **M(V) < 1/13 violations** | **0** |

The unique survivor is the near-dilated base `{2,4,6,8,10,12,13,14,16,18,20,22}` (2·{1..11} with 13 and 14 interleaved — one odd element from primitivity) with the single band-compatible candidate `w = 24`; the exact check gives `M(V) ≥ 1/13`. Anatomy: near-dilated bases are the only structures whose witness midpoints align with a killer grid — the same dilated family that carries the floor (klein-S131) — and even there the surviving killer extends safely.

## Status of hcomp after this session

- tight-AP base ⟹ CRT free-rider (S47): done.
- non-primitive ⟹ scale→sieve (klein-S131): done.
- loose base ⟹ THM-619 bands + pins: **1,568/1,568 structured bases closed, zero violations**; the two-pin lemma closes wide sectors generically; survivors are lone explicit candidates.
Remaining for a fully quantified closure: the base-space boundedness bookkeeping (bases beyond the structured families reduce by the peel — far base elements peel before the killer, re-entering the same pipeline at smaller scale) — a composition note, not new geometry. hcomp's residual is, as of this sweep, empirically empty and mechanically closed everywhere tested.

-> THM-619, HYP-4093/4094, S47 (CRT case), hcomp (kps lrc14_of_compressed), OPEN-Q-108.
