---
source: opus-2026-06-03-S578 (remote-control)
status: FORMALIZATION of sleeve saturation (the covering foundation unifying Φ, the Vitali wall, and the additive folds) + two extensions (the saturation recursion; the alignment principle)
tags: [LRC, sleeve, saturation, covering, deficit, Phi, recursion, alignment, vitali, additive, n14]
---

# Sleeve saturation: the covering foundation

**Prompt (user):** formalize sleeve saturation and related ideas; use the process for
extension.

"Sleeve saturation" is the covering picture that all the recent threads instantiate.
Formalizing it cleanly produced two extensions: a **per-sleeve recursion** generalising
the gap functional `Φ`, and an **alignment principle** that pins why circuit-free
configs stay far from the floor.

## 1. Definitions

Fix `k` speeds, `δ = 1/(k+1)`. The **sleeve** of a runner `v` is
```
Σ_v = { t ∈ ℝ/ℤ : ‖v t‖ < δ },   μ(Σ_v) = 2δ  (exactly, for every v):
```
`v` arcs of width `2δ/v` at the lattice points `{j/v}`. The **total sleeve** is
`U(S) = ⋃_{v∈S} Σ_v`; the **safe set** is its complement; the **saturation deficit**
is
```
Δ(S) := μ(safe(S)) = 1 − μ(U(S)) ≥ 0.
```
Two saturation levels (the **Vitali split**, S551o):
- **measure-saturated:** `Δ(S)=0` (`μ(U)=1`). The worry-set — `safe` is measure-zero
  but *nonempty*.
- **point-saturated:** `U(S)=ℝ/ℤ` (`safe=∅`, `M(S)<δ`). A **counterexample**.

> **LRC(k) ⟺ no speed set is point-saturated.** The sleeves may measure-saturate; they
> never point-saturate.

Total sleeve mass is `Σ_v μ(Σ_v) = 2k/(k+1) ≈ 2 > 1`: there is *more than enough* mass
to cover the circle twice. Saturation is purely a question of **how the sleeves
overlap** — and that is where all the structure lives.

## 2. The deficit is `Φ`, and it recurses (Extension 1)

`Δ(S)` is exactly the gap functional `Φ` of Lemma G (S576, `G(v)=Φ(C)`, verified
exact). Formalising sleeve-by-sleeve gives the generalisation:

> **Sleeve saturation recursion.** Order `S=(v_1,…,v_k)`; set `G_0=ℝ/ℤ`,
> `G_j = G_{j-1} ∖ Σ_{v_j}`. Then `μ(G_j) = μ(G_{j-1}) − Φ_{v_j}(G_{j-1})`, where
> `Φ_{v_j}(G)` is the ramp functional of Lemma G applied with the *current* safe set
> `G_{j-1}` in place of `G(S')`:
> `Φ_{v_j}(G) = (1/v_j)·Σ_{C⊂G} [ v_j|C| − (overlap of phase-interval of C with the band ⋃_k(k−1/n,k+1/n)) ].`
> Hence `Δ(S) = μ(G_k)`, a **telescoping sum of per-sleeve ramp functionals**.

The single-multiple Lemma G is the last step (`j=k`). So `Φ` was the *one-sleeve* law;
the deficit is its **iterate over the sleeve order**. Verified by the saturation curve
`μ_j` (below).

## 3. The saturation curve, and the AP's critical saturation

`μ_j = μ(G_j)` (safe measure after `j` sleeves, level `δ=1/(k+1)`), for the AP
`{1,…,k}`:

| k | `μ_1 … μ_k` | measure-saturates at |
|---|---|---|
| 6 | .714 .571 .381 .238 .057 **.000** | sleeve **#6** |
| 8 | .778 .667 .519 .407 .230 .170 .047 **.000** | sleeve **#8** |
| 10 | .818 … .098 .029 **.000** | sleeve **#10** |
| 12 | .846 … .086 .018 **.000** | sleeve **#12** |

`μ_1 = 1 − 2δ = (k−1)/(k+1)` (first sleeve removes exactly `2δ`, no overlap), then each
later sleeve removes **less** than `2δ`. The AP measure-saturates **exactly at the last
sleeve** — every sleeve is necessary, none redundant. *That is the extremal signature:
the tight config is the one whose sleeves only just manage to measure-saturate, with no
slack, by the final runner.*

## 4. The alignment principle (Extension 2)

Coverage needs overlap (mass `≈2`), but it needs **coherent** overlap. The deficit
`Δ(S)` decreases monotonically with the **3-term count** (the additive alignment):

| | #3t=0 | 1 | 2 | 3 | ≥4 |
|---|---|---|---|---|---|
| k=8 mean Δ | .140 | .138 | .132 | .122 | **.103** |
| k=8 min Δ | .111 | .079 | .059 | .048 | **.015** |
| k=10 mean Δ | .141 | .139 | .140 | .129 | **.110** |

> **Alignment principle (conjecture, verified in samples).** Sleeves saturate (drive
> `Δ→0`) **only when they align additively** — a 3-term relation `v_c=v_a+v_b` is the
> fold (S577) that makes three sleeves overlap *coherently*. Circuit-free (no 3-term)
> ⟹ the sleeves are mis-aligned ⟹ `Δ(S) ≥ c > 0` (Lemma A; here `c ≈ 0.11`). So
> **measure-saturation requires folds**, and the worry-set = the maximally-aligned
> (AP-like) configs.

This is the saturation form of last session's finding ("the gap is depressed only by
3-term relations"): in covering language, *3-term relations are exactly the sleeve
alignments that let the danger mass (≈2) actually cover (≈1).*

## 5. The whole program in saturation language

- **No multiple of (k+1):** the δ-clock points `{j/(k+1)}` lie in *no* sleeve (THM-369)
  — permanent **point-survivors**. Sleeves can measure-saturate but these points hold.
- **Multiple of (k+1) `v=nw`:** its sleeve covers the clock points, but is a *thin
  periodic* sleeve (a Vitali cover); it cannot point-saturate the survivors of the
  other sleeves — `Δ = Φ > 0` (C′, ECCP).
- **Circuit-free:** mis-aligned sleeves, `Δ ≥ c` (Lemma A) — nowhere near saturation.
- **Worry-set:** maximally-aligned (3-term-rich, no multiple) — *measure*-saturated,
  surviving only on the clock points.

`LRC = the point-survivors never vanish`. Three descriptions of the survivor — δ-clock
point, `Φ>0` gap, mis-alignment gap — one phenomenon.

## 6. Honest status / extensions opened

- **Definitions + the measure/point split + `LRC ⟺ no point-saturation`:** rigorous
  (restatement of the covering reformulation, now with the sleeve/deficit language).
- **Extension 1 (saturation recursion `Δ = iterated Φ`):** rigorous given Lemma G;
  verified by the saturation curve. Generalises the one-sleeve `Φ` to all runners.
- **Extension 2 (alignment principle):** verified in samples (Δ monotone in 3-term;
  circuit-free `Δ ≥ ~0.11`); this **is** Lemma A, and its proof (a discrepancy /
  overlap-coherence bound) is the open target — now phrased as "mis-aligned sleeves
  cannot measure-saturate."
- **The critical-saturation observation** (AP saturates exactly at sleeve #k) is a new
  extremal characterisation worth proving: *the tight configs are exactly those whose
  sleeves measure-saturate with zero slack at the final sleeve.*

**Artifacts:** `04-computation/lrc_sleeve_saturation_s578.py` (+`.out`). Builds on
HYP-2112 (Φ=μ), HYP-2114 (folds/3-term), S551o (Vitali), THM-369 (clock), THM-398 (C′).
New: **HYP-2115**.
