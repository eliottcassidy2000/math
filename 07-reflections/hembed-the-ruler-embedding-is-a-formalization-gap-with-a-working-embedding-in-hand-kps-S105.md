# hembed — the ruler embedding is a *formalization* gap, and a working embedding is already in hand

*kind-pasteur-2026-07-09-S105. Owner: prioritize the open math of hembed (THM-527 Part A, the shared
blocker), research broadly, then wire the maxCircGap assembly. klein-S204 discharged hlink while I
worked, so the good-period leg is now `HasGoodPeriod ⟹ Mreach ≥ 1/14` modulo **only hembed**. This
note is the researched map of hembed + the recommended route, with the algebraic core now formalized.*

---

## What hembed is

`hembed : (∃ j φ, ∀ c ∈ teeth E Vmax j, 1/14 < nearInt(φ − c)) → (∃ τ ∈ [0,1], 1/14 ≤ minReach v τ)`.

It is THM-527 Part A's **slow-fast change of variables** — view time from the fastest runner `Vmax`,
split `Vmax·τ = j + φ` (`j` = ruler period, `φ` = fast phase), and turn a good period (teeth gap
`>1/7`) into a real lonely time `τ`. It is the last link of the good-period leg AND the shared blocker
with the density-floor route (both need "a positive good set contains a realizable `τ`").

## The exact algebra (now in Lean: `LRCSlowFast`)

With `τ = (j+φ)/Vmax`, runner speed `u = Vmax − e` (co-offset `e`), tooth `c = frac(e·j/Vmax)`:

> **`nearInt(u·τ) = nearInt((φ − c) − e·φ/Vmax)`** — the fast-phase teeth-clearance `φ − c` minus
> the **slow-fast drift** `δ = e·φ/Vmax`, `|δ| ≤ e/Vmax ≤ spread/Vmax`.

`LRCSlowFast.nearInt_speed_eq_phase_sub` proves the exact identity `nearInt((Vmax−e)·τ) = nearInt(φ −
e·τ)` (via `nearInt` erasing the integer `⌊Vmax·τ⌋`), and `drift_eq` pins the drift `= e·φ/Vmax`.
The apex runner (`e=0`) gives `nearInt(Vmax·τ) = nearInt(φ) ∈ (1/14,13/14)` — safe automatically. This
is exactly the research algebra (agent synthesis, THM-565 setup). **The drift is the whole story:**
negligible for small spread, but `O(1)` in the good-period window (`Vmax ≈ 7·spread/6 ⟹ spread/Vmax ≈
6/7`), so a single fast phase cannot absorb it there — the embedding must use the safe-arc width /
density, not one `φ`. This is why THM-527 states the criterion as `>2/7` (offset-fit `<5/7`, room for
drift), while the `>1/7` version carries `ρ*_{1/7} ≥ m_P = 14249/252252` (PROVED: THM-530 Bonferroni +
THM-661 floor).

## Status: it is a FORMALIZATION gap, not open analysis

- **Math: essentially proved.** THM-527 §A is "PROVED (the slow-fast change of variables)" as a limit;
  the finite-`Vmax` discrepancy `ρ_K = ρ* + O(#arcs/Vmax)` is reduced *three independent ways* to a
  **bounded finite check**: `V* ≈ 234` (THM-565), `≤ 1106` (kps-S61 ledger, proven arc bound `#arcs ≤
  83`), or `≤ 3^12` unconditionally (LEM-010/Dirichlet). `#arcs(G*) = O(spread)` is LINEAR
  (Davenport–Schinzel upper-envelope, not `O(spread²)`) — essential for the pigeonhole at large spread.
- **Lean: the concrete-`τ` embedding is unwritten.** `LRCFourteenSkeleton` carries Part A as an opaque
  `Prop` over `opaque rhoStar`; no existing `witnessG2>0` plumbing ever constructs the real `τ` or
  touches `minReach v τ`. That final geometric step *is* hembed, and it is missing.

## The critical subtlety (do not miss this)

`hembed` as klein states it takes `E, Vmax, v` with **no hypothesis binding `e_i = Vmax − v_i`**. As
an isolated implication it is FALSE (choose `E` unrelated to `v`). The genuine Part-A lemma must add
the binding (`Vmax = max|v_i|`, `E = {Vmax − |v_i|}`, sorted/primitive). The current `hembed` is only
the convenience interface the outer theorem consumes; the real work states and uses the binding.

## The route (highest ROI first)

1. **Reuse the working embedding.** `ScaleSeparation.scale_separation_phase`
   (`LRCScaleSeparation.lean:349`) is a **sorry-free, kernel-pure ruler embedding**: given slow runners
   safe at `t0` (by `1/14+δ`), a ruler `N`, phase spread `Δφ`, drift bound `Dd`, and the budget
   `Δφ + Dd·(δ/V) ≤ 3/7`, it produces a real `t` with ALL runners `≥ 1/14`-safe. This already
   implements the ceiling period selection, midpoint `φ`, 1-Lipschitz slow-safety, and **drift
   absorption**. Recast hembed: `N = Vmax`, slow `R = P`, cluster teeth `→ C` with `hphase` = tooth
   clearance, `hdrift` = bounded co-offset. This turns hembed from "open" into "adapt one lemma."
2. **Make it EXACT at a rational `τ`** (avoid `ρ_K→ρ*` limits): `GridAttainment.grid_margin_domination`
   (max margin is at `m/(|v_i|+|v_j|)`) + `LRCHarmonicGate.rational_point_margin` (clearance at `k/s`
   from the congruence `⌈s/14⌉ ≤ (v_i·k)%s ≤ s−⌈s/14⌉`) → a `native_decide`-friendly certificate,
   mirroring the `HasGoodPeriod` native_decide certs.
3. **Finite-`Vmax` split** (the honest full proof): `ThreeGapSampling.count_pos_of_measure_gt_card`
   (`V·meas(G) > #arcs ⟹ #good ≥ 1`) + `LRCArcComplexity` (`#arcs ≤ 7·ΣE`) + `GapReach` for
   `Vmax > V*`, and a decidable finite check for `Vmax ≤ V* (≤234…1106)`.
4. **Target adapter:** `LRCReachWitness.Mreach_ge_of_lonely_instant` (kps-S99b) reduces hembed's
   conclusion to a pointwise lonely instant `∀ i, 1/14 ≤ nearInt(v_i·τ)`; prove that over the 13
   runners via `LRCSlowFast` + 1-Lipschitz `nearInt`, done. Bridge `margin ↔ minReach` through
   `LRCWitnessAttainmentBridge.margin_eq_minReach` to reuse the `margin`-world lemmas (2, 3 above).

## The inspiration (what frames the thinking)

- **Renormalization / scale separation.** The slow-fast split is a renormalization-group step: `φ` is
  the fast mode, `x = j/Vmax` the slow mode, and the reformulation is the "eigenbasis of the slow-fast
  splitting" (kps two-framings reflection). hembed is the *reconstruction* map (fast+slow → real time).
- **KAM / small divisors.** The drift `e·φ/Vmax` is a small-divisor perturbation; the safe-arc width
  `6/(7Vmax)` and the shrunk safe set `G_P^δ (δ = maxP/2V)` are exactly the KAM-style margin that
  absorbs it. The good-period window is where the "divisor" `Vmax/spread ≈ 7/6` is `O(1)` — the
  resonant regime — which is why it is the hard half.
- **Three-gap / Davenport–Schinzel.** `#arcs(G*)` linear in spread is a Davenport–Schinzel
  upper-envelope bound; the three-distance theorem governs the teeth-gap structure directly.
- **Exact vs. limit.** The deepest reframe: don't chase `ρ_K → ρ*` (a limit with error); land on a
  RATIONAL `τ` where `minReach` is an integer congruence (`decide`). The whole embedding becomes finite.

## One line

hembed is not blocked by analysis — the drift-absorbing embedding it needs already exists sorry-free
(`scale_separation_phase`); the work is (i) add the `e = Vmax − v` binding the interface omits, (ii)
adapt that lemma or land on a rational `τ` (grid domination + rational-point margin), (iii) finite-check
`Vmax ≤ V* (≤ 234…1106)`. The algebraic core is now formalized (`LRCSlowFast`), and the target adapter
(`Mreach_ge_of_lonely_instant`) and the `margin↔minReach` bridge are in place.

*Files: `LRCSlowFast.lean` (this session). Key existing assets: `LRCScaleSeparation.scale_separation_phase`,
`LRCGridAttainment.grid_margin_domination`, `LRCHarmonicGate.rational_point_margin`,
`LRCThreeGapSampling.count_pos_of_measure_gt_card`, `LRCGapReach.exists_nearInt_margin_of_gap`,
`LRCReachWitness.Mreach_ge_of_lonely_instant`, `LRCWitnessAttainmentBridge.margin_eq_minReach`.
See THM-527 §A, THM-565, THM-530, THM-661.*
