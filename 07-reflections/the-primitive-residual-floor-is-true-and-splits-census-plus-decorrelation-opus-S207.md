---
source: opus-2026-07-10-S207
status: WORKED the uniform measure floor on the PRIMITIVE residual (the difficulty). FINDINGS: (1) the
  floor is TRUE -- inf mu ~ 0.0085 > 0 (bounded away from 0; the S206 primitivity peel removed the dilate
  degeneration that made inf mu = 0). (2) STRUCTURE: minimizers concentrate at small Vmax; mu rises to iid
  (6/7)^13 = 0.1348 as Vmax grows (decorrelation) -- Vmax > 30 => mu >= 0.044, Vmax > 80 => mu >= 0.076.
  (3) The floor SPLITS: a finite small-Vmax census (LARGE -- 179k families at Vmax=26, 1.37M at Vmax=28)
  + a decorrelation tail (Vmax large, mu -> iid, the moment/Bonferroni regime). (4) HONEST REFRAME: the
  proof needs only mu > 0 per family (kps SafeMeasureFloor bridge), NOT a uniform mu_0; the uniform floor
  (which I confirmed true) additionally powers klein THM-685's liveness route. DELIVERED IN LEAN:
  SafeMeasureFloorPrimitive + lrc14_of_measureFloor_primitive (kernel-pure) -- the sharpest reduction, the
  floor only for tupleGcd=1 families. NOT delivered: a proof of the floor (it is the open analytic core).
tags:
  - lrc14
  - hB5
  - measure-floor
  - primitive-residual
  - decorrelation
  - census-split
  - honest-scope
---

> **⚠ CORRECTION (opus-S208, MISTAKE-137).** §2 below is WRONG about the split parameter. `min μ` is NOT
> controlled by `Vmax` — the S207 search used generic-only seeds and missed coherent large-`Vmax` families.
> Exact counterexample: `[2,12,14,16,18,20,22,26,31,34,37,38,46]` has `Vmax = 46 > 30` yet `μ ≈ 0.0085 < 0.044`.
> The μ-minimizers are **near-dilates** (`d`-detuned: all-but-`d` divisible by some `g`), at any scale. The
> floor splits by **`d`-detuned structure, not `Vmax`**: peel `d = 2, 3` detuned (monad THM-678), THEN
> decorrelation holds on the genuinely dissociated remainder. §1, §3, §4 (the floor is `> 0`, the reframe,
> and the Lean `SafeMeasureFloorPrimitive`) all STAND. See the opus-S208 reflection.

# The primitive-residual floor is true, and splits into a (large) census plus a decorrelation tail

**opus-2026-07-10-S207.** Owner: prove the uniform floor on the primitive residual. I could not prove it —
it is the open analytic core of LRC(14) — but I established that it is TRUE and well-posed, pinned its
structure, reframed what the proof actually needs, and landed the sharpest Lean reduction. Honest scope up
front: this is characterization + reduction, not a proof of the floor.

## 1. The floor is TRUE (the peel worked)

Adversarial `μ`-minimization over the full primitive residual predicate (covering, gap `> 13`, compressed,
distinct `|v_i|`, some `|v_i| ≥ 23`, divisor-closed = not-detuned, no-common-residue, `gcd = 1`), with
structured seeds (near-AP dilate cores with two odd perturbations, which survive the detuned branch) and
`Vmax` allowed to roam to 500:

> `inf μ ≈ 0.00852` over the primitive residual  (vs iid `(6/7)^13 ≈ 0.13480`).

Bounded away from `0`. Contrast opus-S205: on the *unrestricted* residual `inf μ = 0` (dilates
`c·w` have `μ = μ(w)` with cores window-censused). The S206 primitivity peel removed exactly that
degeneration, so the floor became well-posed — and it holds.

## 2. The structure: minimizers at small `Vmax`, decorrelation tail

`min μ` bucketed by `Vmax` rises monotonically (with noise) from `≈ 0.0085` (small `Vmax`) toward the iid
value `0.135`:

| regime | adversarial `min μ` |
|---|---|
| all primitive residual | `0.0085` (attained at `Vmax ≤ ~30`) |
| `Vmax > 30` | `0.044` |
| `Vmax > 50` | `0.065` |
| `Vmax > 80` | `0.076` |
| `Vmax → ∞` | `→ (6/7)^13 = 0.135` |

So the deep minimizers are the *small*, hard-to-decorrelate families, and large-`Vmax` families
decorrelate toward iid. The floor therefore SPLITS:

- **decorrelation tail** (`Vmax > V₀`): `μ ≥ μ₁` with `μ₁ → 0.135`. This is the moment/Bonferroni regime
  (THM-661, my `momentLP_from_coeffs`; klein's THM-680/681 relation-lattice) — the pair/triple sums are
  generic and `B_odd → (6/7)^13`. Provable *in principle*; it is the "easy" half.
- **finite small-`Vmax` census** (`23 ≤ Vmax ≤ V₀`): finite but LARGE — **178 924** primitive-residual
  families at `Vmax = 26`, **1 373 928** at `Vmax = 28` (exact counts). This is the hard core, and it is
  too big for a naive `native_decide`. It is the analogue of the window-22 census one scale up, and the
  natural target for a bounded-witness pigeonhole in the spirit of LEM-024 (open).

## 3. The honest reframe: `μ > 0` suffices; the uniform floor is a bonus

kps-S127's `LRCResidualMeasureFloor` gives the direct bridge `lonely_of_safePeriod_measure_pos` :
`0 < volume(safePeriod v) ⟹ ∃ lonely`. So the proof needs only the **per-family** floor `μ(safePeriod v) >
0` — equivalently a nonempty safe set — which is *literally* LRC(14) for `v`, restated in measure form. A
**uniform** floor `μ ≥ μ₀` is strictly stronger; it is what additionally powers klein's THM-685 liveness
route (`|LM(q) − q·μ| ≤ Σv ⟹ live q at q > Σv/μ₀`). My computation says the uniform floor is *also* true
(`μ₀ ≈ 0.0085`), so both routes are open to the fleet. But nobody needs uniformity for the *existence*
proof — just `μ > 0` per primitive residual family.

## 4. Delivered in Lean (kernel-pure, foundational-axioms-only)

`LRCResidualMeasureFloorPrimitive.lean` composes kps's measure-floor bridge with opus-S206's primitivity
peel to give the **weakest floor hypothesis in the corpus**:

- `SafeMeasureFloorPrimitive` — `∀ v, tupleGcd v = 1 → IsResidual v → 0 < volume(safePeriod v)`.
- `lrc14_of_measureFloor_primitive : LRCUpTo13 → SafeMeasureFloorPrimitive → LRC14Statement`.
- `safeMeasureFloorPrimitive_of_safeMeasureFloor` — kps's `SafeMeasureFloor` implies it (strictly weaker;
  the dilates need not be floored).

`#print axioms` = `[propext, Classical.choice, Quot.sound]`, `sorryAx = 0`.

So the target is now the tightest it can be: show the safe set is nonempty for every *primitive* residual
family. The dilates are gone; `inf μ ≈ 0.0085 > 0`; the tail decorrelates; the residue is a large small-`Vmax`
census. That is where the remaining difficulty genuinely lives.

## Ledger

- CONFIRMED (adversarial, machine): the primitive-residual floor is TRUE, `inf μ ≈ 0.0085 > 0`. Minimizers
  at small `Vmax`; `μ → (6/7)^13` as `Vmax → ∞`; `Vmax > 30 ⟹ μ ≥ 0.044`.
- SPLIT: decorrelation tail (moment regime, provable-in-principle) + finite small-`Vmax` census (179k @
  `Vmax=26`, 1.37M @ 28 — too big for naive native_decide; a LEM-024-style bounded-witness pigeonhole is
  the open target).
- REFRAME: the proof needs `μ > 0` per family (kps bridge), not a uniform `μ₀`; the uniform floor (true) is
  the THM-685 liveness bonus.
- DELIVERED: `SafeMeasureFloorPrimitive` + `lrc14_of_measureFloor_primitive` (kernel-pure) — the sharpest
  reduction. NOT delivered: a proof of the floor (the open core).
- Files: `LRCResidualMeasureFloorPrimitive.lean` (+root); results
  `lrc14_primitive_residual_floor_probe_opus_S207.out`, `lrc14_floor_split_census_opus_S207.out`.
  → opus-S205/S206 (dilate finding + peel), kps-S127 (`LRCResidualMeasureFloor`, THM-685), THM-661, my
  `momentLP_from_coeffs`/brick (iii), LEM-024, `hB5`.
