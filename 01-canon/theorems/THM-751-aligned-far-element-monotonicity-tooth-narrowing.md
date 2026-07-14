---
id: THM-751
title: Aligned far-element monotonicity (the tooth-narrowing lemma, general form) — if t₀ is a tight point of a speed set P (M(P)=μ₀) and an outlier w is ALIGNED (w·t₀ ∈ ℤ), then M(P ∪ {w·m}) ≥ μ₀·wm/(wm+p_max) for every m ≥ 1, a bound strictly increasing in m up to μ₀=M(P). This is the RIGOROUS proof of THM-726's Step 1 (far-element monotonicity) for aligned lcm-carriers — the extremal multi-killers — replacing the verified-not-proved THM-720/THM-717 with a 4-line explicit-witness argument; non-aligned outliers leave M=M(P) constant (looser). Closes the multi-killer tail: M(P∪{wm}) ≥ 1/13 for all m ≥ an explicit m₀ whenever M(P) > 1/13.
status: PROVED (elementary: triangle inequality + alignment + one balance point; explicit rational witness). VERIFIED (mac-mini-S101): the bound is ≤ actual M and monotone on {1..12}+182m and {1..11,13}+84m; for {1..11,13}+84m (M(core)=1/12) the bound is ≥1/13 for all m≥2 (m=1 by finite check, M=7/89). Non-aligned control ({1..10,13,14}+45m, {1..10,13,22}+84m): M=M(core) exactly constant.
source: mac-mini-2026-07-14-S101 (generalizing the S87 single-killer tooth-narrowing)
depends_on:
  - THM-726   # multi-killer covering-min rigidity — this proves its Step 1 (far-element monotonicity) for the aligned case
  - THM-720   # looseness dichotomy (verified) — this is the rigorous mechanism behind it
  - THM-724   # single-killer covering-min (the μ₀=1/13 boundary case, separate)
related:
  - THM-736   # mac-mini Farey-disc single-killer (the m→∞ limit object)
  - HYP-6660  # the complete tiling map (this closes the loose/spread tile's monotonicity)
external: LRC(≤13) SETTLED.
---

# THM-751 — Aligned far-element monotonicity (the tooth-narrowing lemma)

**One line.** Scaling an *aligned* outlier out **cannot lower M below an explicit bound that rises to
`M(core)`** — the far element's danger tooth, pinned at the core's tight point, only narrows. This is
the rigorous form of THM-726's far-element monotonicity (Step 1), by a 4-line witness, for the aligned
lcm-carriers that realize the extremal multi-killers.

## Statement

Let `P` be a finite set of positive speeds with `M(P) = μ₀ := max_t min_{p∈P} ‖pt‖`, achieved at a
**tight point** `t₀` (`min_{p∈P} ‖pt₀‖ = μ₀`). Let `p_max = max(P)`. Let `w` be an outlier that is
**aligned at `t₀`**: `w·t₀ ∈ ℤ` (equivalently, if `t₀ = a/q` reduced, `q ∣ w`). Then for every integer
`m ≥ 1`,

> **`M(P ∪ {wm}) ≥ μ₀ · \dfrac{wm}{wm + p_max}`.**

The right side is **strictly increasing in `m`** and **`→ μ₀ = M(P)`** as `m → ∞`. In particular, if
`μ₀ > 1/13`, then `M(P ∪ {wm}) ≥ 1/13` for every `m ≥ m₀`, where
`m₀ = ⌈ p_max / (w(13μ₀ − 1)) ⌉` (explicit); the finitely many `m < m₀` are a bounded-outlier check.

## Proof

Take the witness `t = t₀ + s` with `s = μ₀/(p_max + wm) > 0`.

- **Core stays high.** For each `p ∈ P`, by the triangle inequality for the distance-to-ℤ norm
  (`‖x+y‖ ≥ ‖x‖ − ‖y‖ ≥ ‖x‖ − |y|`),
  `‖p t‖ = ‖p t₀ + p s‖ ≥ ‖p t₀‖ − p s ≥ μ₀ − p_max·s`.
- **The tooth is narrow.** Alignment gives `w m t₀ ∈ ℤ`, so
  `‖wm·t‖ = ‖wm t₀ + wm s‖ = ‖wm s‖ = wm·s`, valid because
  `wm·s = μ₀ wm/(p_max+wm) < μ₀ ≤ 1/2`.
- **Balance.** At `s = μ₀/(p_max+wm)` the two bounds coincide:
  `μ₀ − p_max s = wm s = μ₀·wm/(p_max+wm)`.

Hence `min_{v ∈ P∪{wm}} ‖v t‖ ≥ μ₀·wm/(p_max+wm)`, so `M(P∪{wm}) ≥ μ₀·wm/(p_max+wm)`.
Monotonicity and the limit are immediate (`wm/(wm+p_max) = 1/(1+p_max/(wm)) ↑ 1`). ∎

## Why this is THM-726 Step 1 (and what it replaces)

THM-726 (multi-killer, `M ≥ 1/13`) is certified by a **finite check** (bounded outliers) plus **Step 1
— far-element monotonicity** ("scaling an outlier out does not decrease `M`"), which was **verified,
not proved** (via THM-717/THM-720). THM-751 **proves** the mechanism for the aligned case with an
explicit rising bound. It applies exactly to the extremal multi-killers, whose outliers are aligned
lcm-carriers:

| family | core `P` | `μ₀=M(P)` | outlier `w` | aligned? | bound `≥1/13` for |
|---|---|---|---|---|---|
| `{1..11,13}+84m` | `{1..11,13}` | `1/12` at `t₀=1/12` | `84=lcm(12,14)`, `12∣84` | yes | `m≥2` (m=1: `7/89>1/13`) |
| `{1..12}+182m` (single-killer) | `{1..12}` | `1/13` at `t₀=1/13` | `182`, `13∣182` | yes | (`μ₀=1/13`; single-killer, THM-724) |

The single-killer row has `μ₀ = 1/13` exactly, so the bound rises *to* `1/13` but not past it — correct,
since single-killer sits at the covering-min `14/183 < 1/13` and is THM-724's domain, not the multi-killer
`≥1/13` claim.

**Non-aligned outliers are looser (no work needed).** If `w` is *not* aligned (`q ∤ w`) but is safe at
`t₀` (`‖w t₀‖ ≥ μ₀`), then `t₀` itself witnesses `M(P∪{wm}) = μ₀` — **constant** in `m` (verified:
`{1..10,13,14}+45m ≡ 1/11`, `{1..10,13,22}+84m ≡ 2/23`). Either way `M(P∪{wm}) ≥ μ₀·wm/(wm+p_max)`, so
the far-element **minimum over outlier values lives at bounded outliers** — the finite check — which is
THM-726 Step 1.

## Honest scope

- **Proved:** the explicit rising lower bound for aligned outliers, and the constant-`M` fact for
  non-aligned-but-safe outliers. Together: scaling one outlier out never drops `M` below the aligned
  bound, so the multi-killer min is at bounded outliers (Step 1, rigorous, for these strata).
- **Remaining for full THM-726:** (i) the induction when several outliers are scaled (apply the lemma
  one outlier at a time, each step reducing to `M(rest)` of a one-fewer-outlier family); (ii) the
  non-aligned-*unsafe* case (`‖w t₀‖ < μ₀`), where the tight point shifts — empirically still looser,
  but not covered by this witness; (iii) the finite check itself (THM-726 Step 2). These are the
  bounded-diameter residue of the loose tile (HYP-6660), not a monotonicity question.

*Artifacts:* `04-computation/lrc14_far_monotonicity_tooth_macmini_S101.py`,
`lrc14_aligned_tooth_lemma_macmini_S101.py` (+outs). Generalizes the S87 single-killer tooth-narrowing
`M({1..12,182m}) = 14m/(182m+1)`. Credits: THM-726 (Step 1 target), THM-720/717 (the verified mechanism
now proved for the aligned case), THM-724 (single-killer boundary).
