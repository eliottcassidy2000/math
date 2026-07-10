# LEM-021 — The depth-4 dyadic dispatch: a witness at odd/16 with clearance 1/8 exists IFF no speed ≡ 0 (mod 16) and the odd speeds miss a unit ±class mod 16; every even speed has a FREE PASS at depth 4, the tower's unique such layer (16 < 28 < 32)

**Status:** PROVED (both directions one-line, below) + VERIFIED (characterization 0 violations
/ 300k random 13-sets; free pass exact; depth-5 free-pass failure exact). Decides 18.8% of
random covering sets unconditionally (3000-set sample); the blind family exactly characterized.
**Source:** mac-mini-2026-07-09-S65 (cont. 12). The concrete first layer of LEM-020's 2-adic
parity tower (klein-S222/S223's signed road).

## Theorem

For any set S of distinct positive integers:
> **∃ odd m with S lonely at m/16 (clearance ≥ 2/16 = 1/8 > 1/14) ⟺
> (i) no v ≡ 0 (mod 16), and (ii) the odd speeds occupy at most 3 of the four unit ±classes
> {±1, ±3, ±5, ±7} mod 16.**

*Proof.* EVEN FREE PASS: for a ∈ {1,2,3} and odd u, m: `2^a·u·m mod 16 ∈ {2,6,10,14} ∪ {4,12}
∪ {8}` — every value at distance ≥ 2 from 0; and `dist ≥ 2 ⟺ 14·dist ≥ 28 > 16` clears 1/14
with room. (⟸) Pick the missed class's multiplier m (inversion permutes the four ±classes):
every odd r has `rm ∉ {±1}`, hence `rm ∈ {±3, ±5, ±7}`, distance ≥ 3; evens by the free pass;
16∤v keeps everyone alive. (⟹) `16 | v` dies at every odd m; all four classes occupied means
every odd m has some odd r ≡ ±m⁻¹, i.e. `rm ≡ ±1`, distance 1 < 16/14. ∎

**Uniqueness of the sweet spot:** the free pass needs `2^k/14 < 2 ⟺ 2^k < 28`; depth 4
(16 < 28) is the LAST layer where every even runner rides free (verified: depth 5 fails).
One more appearance of the problem's knife-edge calibration.

## What it does and does not do (the owner's question answered exactly)

Depth 4 does NOT decide the whole covering branch: it unconditionally dispatches **18.8%** of
covering sets (with the strict cushion 1/8, alongside the dispatch family: 1/13 detuned, 8/17
common-residue, 1/2 all-odd). The BLIND family is exactly: [some `v ≡ 0 (mod 16)` — deep
2-divisibility, climbing the tower where the free pass is gone] ∪ [odd speeds SPREAD across
all four unit ±classes mod 16 — a domination-type spread demand at modulus 16, feeding the
THM-674/burden machinery]. The parity tower's next layers lose the free pass, so depth 4 is
the tower's one unconditional gift; beyond it the 2-adic road merges back into the band/
domination/burden apparatus — as every road in this problem merges back into the one
calibrated wall.

**Files:** `lrc14_depth4_dispatch_macmini_S65cont12.{py,out}`.
**Related:** LEM-019/020, THM-668/672/674, THM-682(a) (the dispatch family), klein-S222/S223.
