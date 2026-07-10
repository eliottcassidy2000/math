---
id: THM-677
title: The D_m sharpening (boxeph-S3's gap) — the exact harmonic identity D_m = Σ_{h≥1} 2k(h)·Re Z_m(h) (k(h) = sin(πh/7)/(πh), fast convergence by h ≤ 14–28), the L²-SHIFT-AVERAGE clears the needed 5.7% bound by 3.2–3.7× (measured 0.21·√N vs needed 0.057N, TIGHT to the true D_m ≈ 2–4), and the PARSEVAL–BERNSTEIN assembly D_m ≤ √(2H_eff+1)·ℓ₂ closes a-priori for N > N₀ ≈ 1100 (V ≳ 5000) MODULO one one-sided ingredient: the pair-correlation upper bound of the Good set under m-dilation (an E2-type near-coincidence count — divisor/gcd territory, the most tractable residue yet); per-cell DK and plain Erdős–Turán both honestly REFUTED as closers (cells are singletons; ET's harmonic log factor is 4× lossy)
status: MIXED, decisive progress. PROVED: the harmonic identity (the danger window is symmetric ⟹ real coefficients k(h); exact, absolutely convergent); the Parseval identity ℓ₂² = Σ_h 2k(h)²|Z_m(h)|² = mean-over-shifts of the squared deviation; the Bernstein-type assembly shape. MEASURED (boxeph-bank reproductions, exact pipeline): true D_m = 0.4–2.2% of |Good| (4–14× below needed — sharper than the 4% first reported); signed partials converge by h ≤ 14–28; ℓ₂ = 2.8–3.2 vs needed 10.4 (bank2), tight to truth; spectrum square-root-cancelling (max|Z| ≤ 2.6√N, no resonant spikes on covering instances). REFUTED: (i) the per-cell Denjoy–Koksma route — the widest-gap cells of the Good set are ~95% SINGLETONS (Good is j-fragmented; no AP runs to exploit); (ii) plain ET assembly — 40–71 vs needed 10.4–19.4 (the Σ|Z(h)|/h log factor eats the margin; truth has a second cancellation layer across harmonics). OPEN (the named ingredient): the one-sided pair-correlation bound #{(j,j′) ∈ Good² : ‖m(τ_j − τ_j′)‖ ≤ 1/(2H+1)} ≤ (1+ε)·N²/(2H+1) + c·N — an additive-energy/near-coincidence count, upper-bound-only.
source: klein-2026-07-09-S214 (owner-directed "sharpen boxeph's D_m gap"; builds directly on boxeph-S3/HYP-5760)
depends_on:
  - boxeph-S3 (HYP-5760)  # the D_m frame, the Good/τ_j pipeline, the 5.7% budget
  - THM-667   # the mid-band localization this realizes
  - monad-S2 phi-interval  # the drift-valid realization inside good periods
related:
  - THM-676 (the modular-route mirror of the same cancellation wall)
  - LEM-011 (exact Ŵ machinery — the pair-correlation's natural tool)
  - HYP-5692 (one node — this is its sharpest quantitative form on the τ-line)
---

# THM-677 — the D_m Parseval–Bernstein sharpening

**Setting** (boxeph-S3): Good = realizable good periods j (widest cluster gap +
monad's drift-valid φ-interval + exact clearance), τ_j the realized time, N =
|Good|. Killer m ∈ P̃ kills j iff ‖m·τ_j‖ < 1/14. Survival at the 5-killer
budget needs each D_m := |kill_m − N/7| ≤ (1/5 − 1/7)N ≈ 0.0571·N.

## 1. The exact harmonic identity (proved)

The danger window (−1/14, 1/14) is symmetric, so with k(h) = sin(πh/7)/(πh)
(k(h) = 0 at 7 | h) and Z_m(h) = Σ_{j∈Good} e(h·m·τ_j):

> kill_m − N/7 = Σ_{h≥1} 2·k(h)·Re Z_m(h) — exact, absolutely convergent.

Measured: the signed partial sums reach the true D_m by h ≤ 14–28 (the first
two sign-blocks of k dominate); truth D_m = ±2..4 on N = 182–340.

## 2. What fails (measured, honest)

- **Per-cell Denjoy–Koksma**: the widest-gap pair identity fragments Good into
  ~95% singleton cells (166/174, 290/315) — no AP runs; the route is dead.
- **Plain Erdős–Turán**: with the TRUE spectrum (max|Z| ≤ 2.6√N — genuine
  square-root cancellation, no resonant spikes on covering instances), the ET
  assembly still gives 40–71 vs needed 10.4–19.4: the Σ|Z(h)|/h log factor is
  4× lossy. The truth cancels ACROSS harmonics too (signed k(h) blocks).

## 3. The L²-shift-average (proved identity + measured clearance)

ℓ₂² := Σ_{h≥1} 2k(h)²|Z_m(h)|² = the mean over window shifts θ of the squared
deviation (Parseval). Measured: ℓ₂ = 2.8–3.2 on bank2 (N = 182, needed 10.4) —
**3.2–3.7× clearance, and tight to the true D_m**. Since Σ 2k(h)² = 6/49
exactly, ℓ₂² ≈ (6/49)·(mean |Z|²) and mean |Z|² ≈ N·(1 + pair-correlation
deviation): measured ℓ₂ ≈ 0.21·√N.

## 4. The assembly and the one named ingredient

Bernstein-type sup-vs-L² for the shifted-deviation function (effective degree
H_eff ≈ 14 from the k-decay): D_m ≤ √(2H_eff+1)·ℓ₂ ≈ 5.4·0.35·√N ≈ 1.9·√N,
which beats 0.0571·N for **N > N₀ ≈ 1100, i.e. V ≳ 5000** (N ≈ 0.2–0.35·V).
The single a-priori ingredient this rests on:

> **(PC) the pair-correlation upper bound**: #{(j,j′) ∈ Good² :
> ‖m(τ_j − τ_j′)‖ ≤ 1/(2H+1)} ≤ (1+ε)N²/(2H+1) + c·N.

One-sided, a counting statement about near-coincidences of the dilated Good
differences — divisor/gcd/additive-energy territory (LEM-011's exact Ŵ and
the project's E2 machinery are the natural tools), NOT a signed-cancellation
statement. This is the most tractable form the τ-line node has ever taken.

## 5. The state of the finish (τ-line side)

[V ≳ 5000: D_m closes by (PC) + this assembly] + [V < 5000: per-instance
exact certificates (boxeph's pipeline; already realizing lonely τ through the
full 5-killer budget)] + [coherent branch: LEM-012] — with (PC) the one open
lemma. Mirror status on the modular side: THM-676's high-tail bound. Both
routes now name ONE one-sided counting ingredient each.

## Files

`04-computation/lrc14_Dm_cell_sharpening_klein_S214.py` (+.out; cell refutation
+ true-margin ledger), `05-knowledge/results/lrc14_Dm_spectrum_klein_S214.out`
(spectrum + ET loss), `05-knowledge/results/lrc14_Dm_L2_nearmiss_klein_S214.out`
(signed partials + ℓ₂ clearance).
