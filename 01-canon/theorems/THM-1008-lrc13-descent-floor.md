# THM-1008 — The LRC(13)-descent floor: `M ≥ v_max/(13(v_max+v_2nd))`, and the compact residual is exactly `ρ<13` (boxeph-2026-07-18-S83)

**Status:** PROVED (elementary, cites LRC(13) SETTLED) **and the dominant regime `ρ≥13` FORMALIZED
in Lean, kernel-pure** (`TournamentH7/LRCDescentFloor.lean`, `LonelyRunner.descent_dominant`, axioms
`[propext, Classical.choice, Quot.sound]`, no `sorry`; the `1/13`-lonely time is the hypothesis, so
citation-free at the Lean level; companion to `sieve_frac`/`fill1_perturbation`). Verified exactly on
the crux families, 12 random families, and the dominant regime (`lrc_descent_floor_boxeph_S83.py` + `.out`). Gives a
universal floor, proves LRC(14) on the dominant far-element regime by a pure ratio condition
(cleaner than [[THM-1003-fill1-perturbation-base-case]]), and **localizes the compact-core residual
exactly** as `ρ := v_max/v_2nd ∈ [1,13)`.

## Statement

Let `V = {v_1 < … < v_13} ⊂ ℤ_{>0}`, `v_max = v_13`, `v_2nd = v_12`. Then

> **`M(V) ≥ v_max / (13·(v_max + v_2nd)) = ρ/(13(ρ+1))`**,  `ρ := v_max/v_2nd`.

Corollaries:
- **Universal floor:** `ρ ≥ 1 ⟹ M(V) ≥ 1/26` (every 13-family).
- **Dominant regime ⟹ LRC(14):** `ρ ≥ 13` (i.e. `v_max ≥ 13·v_2nd`) `⟹ M(V) ≥ 1/14`.
- **The compact-core residual is exactly `1 ≤ ρ < 13`** — where the floor drops below `1/14`.

## Proof (LRC(13) + one perturbation)

Remove the largest speed: `V' = V \ {v_max}`, twelve distinct speeds. **LRC(13) is SETTLED**, so
`V'` has a `1/13`-lonely time `t₀`: `‖v t₀‖ ≥ 1/13` for all `v ∈ V'` (`‖x‖ =` dist to `ℤ`). The kept
runners carry a **margin** `1/13 − 1/14 = 1/182` above threshold.

Let `d₀ = ‖v_max t₀‖`. If `d₀ ≥ 1/14`, then `t₀` already witnesses `M ≥ 1/14 ≥` the floor. Otherwise
kick by the **minimal amount that lifts `v_max` to the band edge**: with `f = v_max t₀ − round(v_max t₀)`
(signed fractional part), set
`` `s = sign(f)·(1/(14·v_max))`,  `t = t₀ + s`. ``
Then:
- **`v_max`:** `v_max t = round(v_max t₀) + f + sign(f)/14`, so `‖v_max t‖ = ‖f + sign(f)/14‖ ≥ 1/14`
  (for `|f| ≤ 1/2` the sign choice moves *away* from the nearest integer, reaching exactly `1/14`
  when `f=0` and more otherwise). ✓ (`= 1/14` at the extreme.)
- **kept `v ∈ V'`:** `‖v t‖ ≥ ‖v t₀‖ − |v s| ≥ 1/13 − v/(14 v_max) ≥ 1/13 − v_2nd/(14 v_max)`.
  In the dominant regime `v_max ≥ 13 v_2nd` this is `≥ 1/13 − 1/(14·13) = 1/14`. ✓

So `ρ ≥ 13 ⟹ min_V ‖v t‖ ≥ 1/14`. For general `ρ`, run the same kick to threshold `1/T`: the window
`|s| ≤ (1/13−1/T)/v_2nd` sweeps `v_max` by `ρ(1/13−1/T)`, which reaches `1/T` iff
`ρ(1/13−1/T) ≥ 1/T`, i.e. `T ≥ 13(1+1/ρ)`; the best (smallest) `T` gives
`M ≥ 1/T = ρ/(13(ρ+1)) = v_max/(13(v_max+v_2nd))`. ∎

*(Removing `v_max` is optimal: removing any other runner leaves `v_max` in the kept set, so the
window shrinks to `1/v_max` and the fixed runner's ratio is `< 1` — a worse floor `< 1/26`.)*

## Where it sits (verified)

| family | `ρ` | floor `ρ/(13(ρ+1))` | `M` | descent `ρ≥13`? | fill-1 (THM-1003)? |
|---|---|---|---|---|---|
| deep well `{1..12,182}` | 15.2 | 7/97 | 14/183 | **✓** | ✓ |
| `{1..12,161}` | 13.4 | — | 1/13 | **✓** | ✗ |
| residue `{1..11,13,84}` | 6.5 | 84/1261 | 7/89 | ✗ | **✓** |
| AP `{1..13}` / kps-floor-min | ≈1.0 | ≈1/25 | 1/14 / 1/9 | ✗ | ✗ |

**Descent (`ρ≥13`) and fill-1 (THM-1003) are complementary far-element tools** — each certifies
families the other misses. Together they close the far-element regime; the residual is the compact
core `ρ<13` with no fill-1 dominant circle.

## Significance — the residual is the "factor-2 gap"

- LRC(13) gives `1/13`; the naive single-runner descent loses **exactly a factor 2** to `1/26` at
  `ρ=1`. LRC(14)'s target `1/14` sits just above `1/26`, so the compact residual is precisely the gap
  between `1/26` (elementary) and `1/14` (target) for `ρ ∈ [1,13)`. Closing it **requires beating
  single-runner perturbation** — i.e. the measure/harmonic content of Route A (`μ_0>0`), or the
  kind-pasteur covering floor (empirically `M ≥ 1/9`, [[THM-995-trapped-cut-excludes-tight-locus]] X).
- This *quantifies* the finish-map's "bounded-`Vmax` compact core": it is the `ρ<13` band, and the
  analytic content is exactly the sub-factor-2 improvement.
- Uses only **LRC(13) [SETTLED]** as a citation. Formalizable as a conditional perturbation lemma
  (the `1/13`-lonely time is the input), a companion to `LonelyRunner.sieve_frac` and
  `LonelyRunner.fill1_perturbation`.

Related: [[THM-1003-fill1-perturbation-base-case]], [[THM-999]] (pair-sum denominator),
[[THM-995-trapped-cut-excludes-tight-locus]] (sieve-margin IX, covering floor X), HYP-7335.
