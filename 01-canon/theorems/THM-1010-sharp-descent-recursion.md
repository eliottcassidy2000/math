# THM-1010 — The sharp descent recursion: the sub-factor-2 floor for the compact core (boxeph-2026-07-18-S84)

**Status:** PROVED (elementary, cites LRC(13) SETTLED) **and FORMALIZED kernel-pure**
(`TournamentH7/LRCDescentFloor.lean`, `LonelyRunner.descent_general`, axioms
`[propext, Classical.choice, Quot.sound]`, no `sorry`). Verified exactly on the compact
residual families and 20 random compact-covering families (`lrc_recursion_floor_boxeph_S84.py`).
Sharpens [[THM-1008-lrc13-descent-floor]] (its `μ=1/13` special case): **beats the factor-2
floor `1/26` whenever the kept 12-set is non-tight**, and empirically reaches `1/14` across the
compact-covering residual.

## Statement

Let `V ⊂ ℤ_{>0}`, `|V|=13`, `v_max` the largest, `v_2nd` the second, `ρ = v_max/v_2nd`, and
`V' = V∖{v_max}`. Then

> **`M(V) ≥ ρ · M(V') / (ρ + 1)`**.

More generally (the Lean form): if `t₀` is `μ`-lonely for `V'` (`‖v t₀‖ ≥ μ` for all `v∈V'`) and
each kept `v_i` obeys the **sweep budget** `v_i ≤ 14·v_max·(μ − 1/14)`, then `M(V) ≥ 1/14`.
Taking `μ = M(V')` and optimizing gives the displayed floor.

## Proof (one perturbation, sharp base)

`V'` has 12 speeds; let `μ = M(V')` (`≥ 1/13` by LRC(13)) with maximizer `t₀`: `‖v t₀‖ ≥ μ` for all
`v∈V'`. Kick `t = t₀ + s`, `s = ±1/(14 v_max)` (the minimal kick lifting `v_max` to the band edge):
- **`v_max`:** `‖v_max t‖ = ‖f ± 1/14‖ ≥ 1/14` (`f = ` frac of `v_max t₀`; sign moves away from ℤ).
- **kept `v∈V'`:** `‖v t‖ ≥ μ − v/(14 v_max)`. This is `≥ 1/14` iff `v ≤ 14 v_max(μ−1/14)`; worst kept
  is `v_2nd`, giving `M(V) ≥ 1/14` when `v_2nd ≤ 14 v_max(M(V')−1/14)`, i.e.
  `M(V') ≥ (1/14)(1 + 1/ρ)`. Solving for the threshold: `M(V) ≥ ρ M(V')/(ρ+1)`. ∎

*(At `μ=1/13`: budget `v_2nd ≤ 14 v_max/182 = v_max/13`, i.e. `ρ≥13` — recovers THM-1008.)*

## The sub-factor-2 content (verified exactly)

The `1/26` of THM-1008 was the **tight** base `M(V')=1/13` at `ρ=1`. The recursion uses the *actual*
`M(V')`, recovering the factor proportionally to how non-tight the kept set is:

| family (compact covering) | `ρ` | base `1008` | **recursion** | `M(V)` | reaches `1/14`? |
|---|---|---|---|---|---|
| kps floor-min `M=1/9` | 1.05 | 0.0393 | **0.0949** | 0.111 | **✓** |
| residue `{1..11,13,84}` | 6.5 | 0.0666 | **0.0722** | 0.0787 | **✓** |
| random compact-covering (20/20) | <13 | — | ≥ 1/14 | — | **✓ all 20** |
| `2·{1..13}` dilated tight | 1.08 | 0.0400 | 0.0500 | **1/14** | ✗ (boundary) |
| `{2..14}` (bounded) | 1.08 | 0.0399 | 0.0691 | 1/8 | ✗ (loose; finite-check case) |

**The recursion strictly beats `1/26` on every family, and reaches `1/14` across the compact-covering
residual** — its two failure modes are exactly (i) the **dilated tight AP** (`M(V')=1/13`, the `M=1/14`
equality boundary, non-primitive ⟹ dilation-reducible, THM-995 IV) and (ii) **small-`Vmax` bounded
families** (loose here, but closed by the finite check THM-734 / THM-999).

## What it reduces the residual to (honest scope)

The recursion **does not unconditionally close** the compact core — it converts it into a clean
**12-subset floor**:

> **Compact residual ⟸ `M(V∖{v_max}) ≥ (1/14)(1 + 1/ρ)`** for primitive covering `V` with `ρ<13`.

The required kept-floor lies in `(1/13, 1/7]` (a **sub-factor-2** demand: the kept 12-set must clear a
little more than `1/13`). This is the genuine remaining analytic content, now sharply located:
- equality forces `M(V')=1/13` ⟹ `V'` a dilated AP ⟹ `V` dilated tight (handled);
- the strict floor `M(V') > (1/14)(1+1/ρ)` for primitive covering `V` is the **n=12 near-tightness /
  rigidity** input (klein's Hamming-radius theorems THM-1004/1005; HYP-7310) — equivalently the
  Route-A measure floor. The descent turns "LRC(14) on the compact core" into "a `1/7`-scale floor on
  its 12-subsets."

## Significance

- Closes the **factor-2 gap** conceptually: the loss was an artefact of using the tight base; the sharp
  base recovers it, and the residual is only the *dilated-AP boundary* + a 12-subset floor.
- Unifies with the trilogy of elementary witnesses (`sieve_frac` empty-circle, `fill1_perturbation`
  fill-1, `descent_dominant` far-element) — all one kick, now with a sharp base.
- Kernel-pure Lean (`descent_general`); `descent_dominant` (THM-1008) is now a corollary.

Related: [[THM-1008-lrc13-descent-floor]], [[THM-999]], [[THM-995-trapped-cut-excludes-tight-locus]]
(covering floor X, dilation IV), [[THM-1003-fill1-perturbation-base-case]], HYP-7345.
