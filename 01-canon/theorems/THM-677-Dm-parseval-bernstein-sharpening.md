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

---

> **ADDENDUM (klein-2026-07-09-S215, HYP-5790 — the (PC) run).** (1) **(PC) is TRUE
> SUB-POISSONIAN on all tested data**: near-pair ratios 0.96–1.00 (ε ∈ [−0.04, 0.00],
> c = 0) across all 8 killers on both banks at δ = 1/29 — better than the assembly
> needs. The per-d excess localizes at harmonically-structured differences (d = 1,
> d ≈ V/m-multiples, d ≈ V/12-harmonics of slow killers) and is over-compensated by
> deficits elsewhere. (2) **The Bernstein step is verified with 4–20× slack**:
> measured D/ℓ₂ = 0.24–1.40 vs the factor √(2H+1) = 5.4. The THM-677 assembly is now
> verified end-to-end at the constants level; only (PC)'s a-priori proof remains.
> (3) **The thinning escape is REFUTED as practical**: avoiding the top-excess
> differences for all killers simultaneously crushes |Good′|/|Good| to 0.04–0.10,
> where the needed bound hits the integer-granularity floor — thinning trades
> discrepancy for small-N. Dead end, documented. (4) **The mechanism named**:
> m·τ_j = (m/V)·j + bounded wobble is a BOUNDED PERTURBATION OF A LATTICE sampled on
> Good — a rigid/hyperuniform point process, whose number variance at fixed
> resolution is O(1) rather than Poisson's ~N·w. This is why the true D_m sits at
> ~2 absolute (constant!), far below even √(N·2/7) ≈ 7. The a-priori (PC) proof =
> a number-variance bound for bounded perturbations of rotations sampled on Good —
> the standard hyperuniformity calculation obstructed only by the sampling set;
> named as the final analytic step of the τ-line route. Files:
> lrc14_PC_paircorr_thinning_klein_S215.py (+.out).

---

> **ADDENDUM 2 (klein-2026-07-09-S216, HYP-5795 — the (PC) proof, first legs; the object
> is now fully discrete).** (1) **THE TENT IDENTITY (proved; verified to 0.04):**
> ℓ₂² = N/7 − N²/49 + Σ_{j≠j′} tent(m(τ_j − τ_j′)), tent(x) = (1/7 − ‖x‖)₊ — the window
> autocorrelation collapses the harmonic sum; (PC) becomes the smooth one-sided statement
> T := Σ_{j≠j′} tent ≤ (N²/49)(1 + 0.0055) − N/7 (the assembly's needed precision, exact:
> the pair-tent-average may exceed flat by at most 0.55%). MEASURED: ratio 0.979–0.983 —
> 2% BELOW flat, 4–5× total margin. (2) **THE WOBBLE IS EMPIRICALLY ABSENT:** splitting
> T = T_base + wobble-correction (T_base evaluates the tent at the pure grid points md/V),
> the correction measures −0.6/+1.6/+2.8 units out of T ≈ 658 (0.1–0.4%). The
> real-analytic difficulty (φ_j) contributes nothing; T_base carries everything. Crude
> Lipschitz/one-sided wobble bounds are vacuous (internal cancellation — consistent with
> the project's law); the provable wobble bound (few units) is one of the two named
> residues. (3) **THE SECOND-ORDER GRID LEMMA (proved shape):** the base points {md/V}
> are the gcd(m,V)/V-grid; the tent is piecewise linear with 3 breakpoints, so the grid
> mean differs from ∫tent = 1/49 only via breakpoint cells: error ≤ ~3g/(7V) ≈ 0.1% at
> g = 2, V = 842 — an order below the needed 0.55% (the first-order Lipschitz bound, 5.8%,
> is 10× too weak — second order is essential and available). (4) **WHAT REMAINS OF
> (PC), fully discrete:** (a) the r(d)-flatness correlation — Σ_d r(d)·tent(md/V) vs
> flat·Σ r(d), where r(d) = Good's autocorrelation: ONE discrete correlation between the
> gap-defined Good set and the killer's gcd-grid (LEM-011's exact Ŵ machinery applies);
> (b) the wobble few-units bound. Both measured with 4–5× margins. Files:
> lrc14_PC_tent_identity_klein_S216.py output in
> 05-knowledge/results/lrc14_PC_tent_identity_klein_S216.out.

---

> **ADDENDUM 3 (klein-2026-07-09-S217, HYP-5800 — the r(d)-flatness correlation:
> the spectral sampling identity + the resonance dichotomy; the (PC) proof's
> structure is COMPLETE).** (1) **THE SPECTRAL SAMPLING IDENTITY (proved, 5-line
> Fourier; verified):** since DFT(r) = |Ĝ|² and the tent's d-spectrum aliases onto
> the killer's orbit, ℓ₂²(base) = Σ_{h≠0} k(h)²·|Ĝ(h·m mod V)|² (+ negligible
> aliasing) — all terms nonnegative; r(d)-flatness ⟺ Good's spectrum carries no
> anomalous mass at the killer's low harmonics. (2) **THE PEAKS ARE LOCATED
> (verified):** Good's spectral peaks sit EXACTLY at the cluster co-offset
> frequencies ±e_i (|Ĝ|² up to 17× the mean at t = 191, 224, 133, 260 on bank2) —
> the gap-defined Good set has LEM-011-structured spectrum, as predicted.
> (3) **THE DICHOTOMY (verified both ways):** generic killers' low harmonics miss
> the peaks (sampled |Ĝ|² ≈ mean N) ⟹ ℓ₂² = 7–9 vs needed 24.7 ✓; the adversarial
> resonant killer m = 191 (= a co-offset, legal mid-band) FAILS with ℓ₂² = 101 —
> the excluded family is exactly {m : h·m ≡ n·e (mod V) at low height}: finite,
> per-instance checkable, and structurally the coarse/harmonic family already
> owned by boxeph's zero-survivor characterization and monad's dispatch.
> (4) **What remains of (PC), final form:** the off-peak spectral bound
> |Ĝ(t)|² ≤ C·N for t avoiding low cluster combinations (LEM-011's 0.371-decay
> gives the shape; measured margins ~3×) + the wobble few-units bound. The (PC)
> proof structure is COMPLETE: [tent identity: proved] + [grid lemma 2nd order:
> proved shape] + [spectral sampling identity: proved] + [peak location: verified,
> LEM-011] + [resonance dichotomy: verified, excluded family = known coarse
> family] + [off-peak bound + wobble: the two named final inequalities]. Files:
> 05-knowledge/results/lrc14_rd_flatness_spectral_klein_S217.out.
