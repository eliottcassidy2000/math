# The core-length monotonicity is anchored by `1/13 > 14/183`; the elementary witness closes only single-killer

*kind-pasteur-2026-07-11-S127 cont.59. Owner: "prove the core-length monotonicity as a global inequality." I
proved the rigorous global lower bound the balance gives, found it closes the single-killer case tightly, and
proved — with an explicit undershoot — that it CANNOT close genuine multi-killer. The honest outcome: the
monotonicity is real, but its engine is the inequality `1/13 > 14/183`, and the tight multi-killer bound reduces
to LRC(13)-escape + finite check, not to an elementary global witness.*

---

## The rigorous global lower bound (proved)

For a primitive covering family `F` with maximal interval core `{1,…,k}`, perturb `t = 1/(k+1) − δ` (`δ>0`):

- **Core:** for `1 ≤ i ≤ k`, `‖i·t‖ ≥ 1/(k+1) − δ`, with runner 1 binding — all other core runners stay above
  it for `δ ≤ 1/(k+1)` (runners `i > (k+1)/2` rise; runners `2 ≤ i ≤ (k+1)/2` start higher and fall slower
  relative to their head start). So **core-min `= 1/(k+1) − δ`.**
- **Resonant killers** (a multiple of `k+1`, so `κ/(k+1) ∈ ℤ`): `‖κt‖ = κδ` — they *rise*.
- **Non-resonant killers** `κ'`: `‖κ't‖ ≥ ‖κ'/(k+1)‖ − κ'δ ≥ 1/(k+1) − κ'δ` — they fall at rate `κ'`.

Balancing the smallest resonant killer `a` (rising) against the fastest-falling term (runner 1 at rate 1, or a
non-resonant killer at rate `R = max κ'`):

> **`M(F) ≥ (1/(k+1)) · a/(a + max(1,R))`,   `a` = smallest resonant killer, `R` = fastest non-resonant killer.**

Verified `M(F) ≥ B(F)` on the extremizers. This is a genuine, proved **global inequality**.

## It closes single-killer — tightly — and nothing more

| family | `k` | `a` | `R` | `B(F)` | true `M` | closes? |
|---|---|---|---|---|---|---|
| deep well `{1..12,182}` | 12 | 182 | 1 | **14/183** | 14/183 | ✔ tight |
| `{1..11,13,84}` | 11 | 84 | 13 | 7/97 ≈ .0722 | 7/89 ≈ .0787 | ✘ |
| `{1..10,13,22,84}` | 10 | 22 | 84 | .019 | 2/23 ≈ .0870 | ✘ |

Single-killer has `R=1` (no non-resonant killer), so `B = (1/13)·182/183 = 14/183` — the bound *is* the floor.
But a **genuine multi-killer family always has a fast non-resonant killer**: it must cover a divisor `d ≠ k+1`
(e.g. `d = 13` via the outlier `13`), and that covering killer is not a multiple of `k+1`, so it is
non-resonant at `t_core` and falls at its own (large) rate `R`. That makes `B(F) < 14/183`. The elementary
witness provably undershoots.

## Why no elementary witness works — the family-specific modulus

The undershoot is not slack in the bound; the **core-resonance time itself is the wrong place to look.** For
`{1..11,13,84}`, the best time near `t_core = 1/12` gives only `7/95 ≈ 0.0737 < 14/183`, because lifting the
resonant killer `84` off `0` drags the non-resonant killer `13` down at rate 13. The **true** optimum `7/89`
lives at `t* = 37/89` — a *different, larger, family-specific* modulus (89) where `13` and `84` are placed
simultaneously. Locating that modulus is exactly the band-packing / clearing problem, i.e. the covering-min
itself. **There is no uniform elementary witness for multi-killer.**

## The real engine: `1/13 > 14/183`, and the honest reduction

So why is the monotonicity nonetheless true (cont.58)? Because of one inequality:

> **`1/13 = 0.076923 > 14/183 = 0.076503`** — the LRC(13) floor sits *above* the LRC(14) covering-min.

Every multi-killer family contains a **12-speed sub-family** (drop the largest killer) which, by *settled*
LRC(13), has `M ≥ 1/13 > 14/183`. The full family drops below only through the one dropped killer `L`:

- **`L` large:** decorrelation (the cont.48 atom) gives `M(F) ≥ 1/13 − B/L ≥ 14/183` — the family stays
  anchored to its `1/13` sub-family. **Escape.**
- **`L` bounded:** `M(F)` is a finite check (klein's ILP, `≤ 182`) — and the binding multi-killer minimizers
  (`{1..11,13,84}`, outliers `≤ 84`) live here.

So the tight core-length monotonicity `multi-killer ⟹ M ≥ 14/183` reduces to **LRC(13)-escape + finite check**,
anchored by `1/13 > 14/183`. This is the same escape+ILP structure as the covering-min proper — the
monotonicity is **not easier than the main problem**, and the razor-thin single-killer drop (`1/13 → 14/183`,
a gap of `1/2379`) is the *only* place a covering family pierces the `1/13` anchor to reach the floor.

## Net

Answering the owner's request honestly: the **elementary global inequality** `M(F) ≥ (1/(k+1))·a/(a+max(1,R))`
is proved and closes the **single-killer** case tightly at `14/183`. The **tight** monotonicity for genuine
multi-killer is **not** an elementary global inequality — the core-resonance witness provably undershoots
(shown: `7/95 < 14/183`), the optimum is family-specific, and the correct proof is the LRC(13)-escape + finite
check anchored by `1/13 > 14/183`. The positive residue is a clean structural understanding: the covering-min
sits just below the LRC(13) floor, multi-killer families are pinned to that floor by their sub-families, and
only the single lcm-killer reaches through. The genuinely open computation is unchanged — the bounded finite
check (ILP) — now with the multi-killer part reduced to it explicitly via settled LRC(13).

*Files: lrc14_corelength_monotonicity_kps_S127.py, lrc14_corelength_balance_bound_kps_S127.out. Builds on
kps cont.58 (multi-killer enumeration), opus-S253 (balance), the cont.48 decorrelation atom, LRC(≤13) citation;
rests on THM-366, klein-S267 (14/183). HYP-6228.*
