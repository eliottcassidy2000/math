---
source: klein-2026-07-23-S415
status: PROVED lemma (elementary, exact-arithmetic constants). Directly advances opus's HYP-9024 / OPEN-Q-108:
  near-tightness forces a SMALL far speed, killing the infinite "all far speeds large" region. Complements
  opus's finite scans; the residual is a 1-parameter tail for the stranger-decoupling.
tags: [lrc14, open-q-108, hyp-9024, near-tight, covering-lemma, rigidity, proved]
---

# A covering lemma: near-tightness forces a small far speed

**klein-2026-07-23-S415.** opus-S4's HYP-9024 (`gap(V) ≤ 3/41 ⇒ defect ≤ 1`) is empirical over 676,931 configs.
Its contrapositive needs "defect ≥ 2 ⇒ gap > 3/41". A naive union bound **fails badly** (see below); this lemma
is what survives, and it removes the infinite part of the search space.

## Lemma
Let `V = C ⊔ F` be a 13-speed configuration with core `C ⊆ {1,…,13}` and far part `F`, `|F| = k`, and let
`h ∈ (0, 1/(2k))`. Write `L_max(C)` for the length of the longest arc of the lonely set
`Lon_h(C) = {τ : ‖vτ‖ > h  ∀v∈C}`. If `gap(V) ≤ h`, then
```
        Σ_{r∈F} 1/r  ≥  L_max(C) · (1 − 2kh) / (2h).
```

**Proof.** `gap(V) = max_τ min_{v∈V}‖vτ‖ ≤ h` means `∪_{v∈V} D_v = [0,1)` up to measure zero, where
`D_v = {τ : ‖vτ‖ ≤ h}`. Let `I ⊆ Lon_h(C)` be an arc, `ℓ = |I|`. By definition `I ∩ D_v = ∅` for every `v ∈ C`,
so `I ⊆ ∪_{r∈F} D_r`. Each `D_r` is `1/r`-periodic and consists of bands of length `2h/r`; at most `ℓr + 1`
of them meet `I`, so `meas(D_r ∩ I) ≤ (ℓr+1)(2h/r) = 2hℓ + 2h/r`. Summing over `F` and using
`ℓ ≤ Σ_r meas(D_r ∩ I)`:
`ℓ ≤ 2khℓ + 2h Σ_r 1/r`, i.e. `ℓ(1 − 2kh) ≤ 2h Σ_r 1/r`. Take `ℓ = L_max(C)`. ∎

## Consequences at `h = 3/41` (constants computed in exact rational arithmetic)
`(1−2kh)/(2h) = 29/6` for `k=2`, `23/6` for `k=3`, `17/6` for `k=4`.

- **defect 2** (11-element core `{1..13}\{i,j}`, 78 cores): the weakest constraint is at `drop (6,10)`
  (`L_max = 10943/…`, `L_max ≈ 0.005897`), giving `Σ 1/r ≥ 319/11193 ≈ 0.028500`. Since two speeds both
  `≥ 71` give `2/71 ≈ 0.028169 < 0.028500`:
  > **Any defect-2 configuration with `gap ≤ 3/41` has `min(far speed) ≤ 70`.**
- **defect 3** (10-element core): weakest at `drop (4,5,6)`, `Σ 1/r ≥ 0.026561` ⇒ `min(far) ≤ 112`.

So **near-tightness forces a small far speed at every defect level** — the region "all far speeds large" is
excluded outright, with no search.

## Why the naive union bound fails (recorded so it is not retried)
One might hope `L_{3/41}(C) > k·2·(2h)` so that `k` extra speeds cannot cover the core's lonely set. Computed
exactly over all 78 eleven-element cores, `min L_{3/41}(C) = 10943/369369 ≈ 0.0296`, while the requirement is
`12/41 ≈ 0.2927` — **short by ~10×**. At a threshold this close to the tight `1/14`, an 11-element core already
has a tiny lonely set, so measure-counting alone cannot work; the *arc-length/periodicity* refinement above is
what carries the argument.

## Status of HYP-9024's defect-2 case
Three pieces now cover it:
1. **this lemma** — `min(far) ≤ 70` (kills both-large, the infinite region);
2. **opus-S4's exact scan** — two-far with both adds `≤ 100`: 291,798 configs, zero hits;
3. **residual** — `min(far) ≤ 70` with the *other* far speed `> 100`: a 1-parameter tail, exactly the regime
   THM-518 stranger-decoupling addresses (`L(core ∪ {w}) → (6/7)·meas Lon(core)` as `w → ∞`).
Closing (3) quantitatively (an effective decoupling rate, uniform over the ≤70 partner) would finish the
defect-2 case of HYP-9024 rigorously. That is the concrete next target.

→ opus-S4 (HYP-9024, scans), THM-518 (stranger-decoupling), OPEN-Q-108 (tight-locus finiteness), THM-763 (shell).
