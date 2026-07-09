---
source: opus-2026-07-09-S184
status: MAPPED the LRC(14) Lean endgame precisely and supplied the conceptual bridge for its lone
  hfloor crux. The proof is a SORRY-FREE CONDITIONAL (no sorry/admit/sorryAx; native_decide in censuses):
  LRC14Statement reduces (skeleton lrc14_from_witness_floor) to exactly hfloor + hpartA. hpartA reduces
  (proven sorry-free lemmas) to the single REALIZATION step hrefl/hembed, discharged EXCEPT the binding
  good-period window Vmax~7*spread/6 (THM-608 covers the drift-controlled regime). hfloor reduces
  (LRCWitnessBonferroni) to hbonf (PROVED, kps-S30 toReal_bonferroni) + Lemma B (proved rational) +
  LEMMA A nu(E)>=nuConsec(k) -- the LONE genuinely-open hfloor node. So the WHOLE remaining surface is:
  {Lemma A} and {hrefl binding window}. BRIDGE: Lemma A = "the consecutive cluster MINIMIZES the good-set
  measure nu" = the MEASURE-level version of my LEM-015 (interval MAXIMIZES Schur triples E3). VERIFIED:
  consecutive achieves min nu=0.4425 (=all dilated intervals; nu is scale+rotation invariant), and
  corr(E3, nu) = -0.911 -- the AP minimizes nu BECAUSE it maximizes E3 (additive resonance =>
  equidistributed orbit => small maxgap => small good-set). Converges kps-S113/mac-mini-THM-668 (E3 = the
  hardness coordinate).
tags:
  - lrc14
  - formalization
  - endgame
  - lemma-a
  - schur-triples
  - density-floor
  - realization
---

# The LRC(14) formalization endgame, and the Schur bridge for Lemma A

**opus-2026-07-09-S184.** Owner: finish all LRC(14) formalization. Mapping it precisely (via a dedicated
scan) shows the proof is a **sorry-free conditional** whose entire remaining surface is TWO analytic
cruxes — and my session's Schur/`E₃` work is the conceptual key to one of them.

## The endgame map (precise)

`LRC14Statement : ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∃ t, Lonely 14 v t`. No `sorry`, `admit`, or
`sorryAx` anywhere in the LRC files; the only trust beyond Lean foundations is `native_decide` in the
finite censuses. `lrc14_from_witness_floor` reduces it to exactly two parameters:

- **`hpartA`** (`0 < witnessG2 → Mreach ≥ 1/14`). The whole reach chain is proved sorry-free
  (`LRCReachWitness`, `LRCCriterionC`, `LRCHembedIdentity`, `LRCSmoothBridge`, `LRCGoodPeriodReach`,
  `LRCGapReach`), leaving one step: the **realization** `hrefl`/`hembed` — realize a real `τ` whose fast
  phase `frac(Vmax·τ)` lands in the slow teeth's good-period gap. Discharged in the drift-controlled regime
  (`LRCHembedScaleSep`, THM-608, `Δφ + Dd·δ/Vmax ≤ 3/7`); **open only in the binding window** `Vmax ≈
  7·spread/6` that saturates the `3/7` budget. Fleet-active (kps-S112 continuum, monad/mac-mini THM-666/668
  grid-port + drift-embed, boxeph/death-star realization).
- **`hfloor`** (`witnessG2 ≥ m_P`). Reduced (`LRCWitnessBonferroni`, THM-527 Part G) to three nodes:
  - `hbonf` (Bonferroni `μ(A∩B) ≥ μA+μB−1`): **PROVED** — kps-S30, `toReal_bonferroni`, sorry-free.
  - Lemma B (`measGP ≥ cap_k`): proved finite rational.
  - **Lemma A** (`ν(E) ≥ νConsec(k)`): the **lone genuinely-open hfloor node**, a scale-invariant
    three-distance measure inequality framed as a compactness crux.

So the entire remaining mathematical surface of LRC(14) is **{Lemma A} ∪ {hrefl binding window}**.

## The Schur bridge for Lemma A

Lemma A says the **consecutive co-offset cluster `{0,1,…,k−1}` minimizes** `ν(E) = meas{x : maxgap{frac(eᵢx)}
> 1/7}` (the good-set measure). This is precisely the **measure-level version of LEM-015** (the interval
*maximizes* the Schur-triple count `E₃`): the consecutive cluster = the AP speed set = my interval
extremal. Verified (`lrc14_lemmaA_consecutive_minimizes_nu`):

| cluster | `ν(E)` | `E₃` |
|---|---|---|
| consecutive `{0..12}` / interval / `2·AP` | **0.4425 (min)** | 78–91 |
| near-2-AP | 0.4711 | 81 |
| GAP 4×4 | 0.7709 | 79 |
| Sidon-ish | 0.9949 | 25 |

The consecutive cluster attains the minimum `ν` (tied across all dilated intervals — `ν` is scale- AND
rotation-invariant, so translates/dilates of the AP agree), and **`corr(E₃, ν) = −0.911`**. The mechanism:

> **The AP minimizes `ν` BECAUSE it maximizes `E₃`.** Maximal additive resonance (Schur triples) ⟹ the
> orbit `{frac(eᵢx)}` is maximally equidistributed for most `x` ⟹ smallest maxgap ⟹ smallest good-set ⟹
> minimal `ν`. Sidon/dissociated clusters (few Schur triples) have clumpy orbits, big gaps, `ν ≈ 1` (easy).

This is the resonance-side "why" for the compactness crux: the extremal of Lemma A is forced by the same
additive-coherence coordinate `E₃` that LEM-015 makes rigorous discretely, and that kps-S113 and
mac-mini's THM-668 ("the Schur-triple kill rule mechanizes the `E₃` law; the AP has exactly one live ruler
`q=14`") independently identified as the hardness axis. Three routes (my resonance sum, kps's density-floor
cert, mac-mini's ruler witnesses) now agree: **`E₃` is the LRC(14) hardness coordinate, and the AP is its
maximizer.**

## Honest status

The two cruxes are genuine analysis (a compactness minimization and an equidistribution realization),
actively worked by the fleet; neither closes in a single solo session. What this session settled: `hbonf`
is already done, so **Lemma A is the *entire* open content of `hfloor`**, and it is an interval-extremal
statement whose extremal and mechanism are now pinned (consecutive minimizes `ν`; `E₃` is why). The
realization crux `hrefl` is likewise reduced to one regime. The scaffold is otherwise sorry-free.

## Ledger

- LRC(14) Lean = sorry-free conditional; entire remaining surface = **Lemma A** (lone open hfloor node) +
  **hrefl binding window** (`Vmax≈7spread/6`). hbonf PROVED (kps-S30), Lemma B proved, reach chain proved,
  hpartA reduced to hrefl, hrefl discharged outside the binding window (THM-608).
- Lemma A = measure-level LEM-015: consecutive MINIMIZES `ν` (verified, min 0.4425), `corr(E₃,ν)=−0.911`;
  the AP minimizes `ν` because it maximizes `E₃` (resonance ⇒ equidistribution ⇒ small good-set).
- 3-way convergence on `E₃` as the hardness coordinate: opus (resonance sum) + kps-S113 (density cert) +
  mac-mini THM-668 (ruler witnesses).
- File: `lrc14_lemmaA_consecutive_minimizes_nu_opus_S184` (+out). -> LEM-015, opus-S182/S184 (step 2),
  THM-661 (moment-LP), LRCWitnessBonferroni (Lemma A/hbonf), THM-608/THM-666/THM-668 (realization),
  kps-S30 (toReal_bonferroni), kps-S113, mac-mini-S65.
