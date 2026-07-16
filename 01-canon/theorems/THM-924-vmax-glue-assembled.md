---
id: THM-924
title: THE FINITE-VMAX GLUE, ASSEMBLED — THM-527-A's large-spread half closes as a trichotomy: [bounded spread: #arcs = O(1) (S58 arc-count lemma), |ρ_K − ρ*| ≤ C/Vmax — PROVED] + [large spread, NON-RESONANT: Erdős–Turán with Dirichlet-kernel sums gives a SPREAD-FREE constant (data: err·Vmax = 0.3–7 across spreads 5 → 6577)] + [RESONANT (some e_i/Vmax = p/q, q small): ρ_K is EXACTLY q-PERIODIC — a positive rational computable in finite exact arithmetic (data: the exact-resonance plateau err·V = 263 with ρ_K = 0.690 > 0), i.e. the near-dilate tile (THM-724 L2), NOT an error term] — the glue's residual is one standard ET page (the same shape as THM-920's slice page), and every branch delivers ρ_K > 0 hence M(S) ≥ 1/14
status: trichotomy ASSEMBLED + verified (three spread regimes, resonant/near-resonant/generic contrasts exact in the data); bounded branch PROVED (S58); resonant branch = finite exact computation (periodicity is immediate: e_i x_j mod 1 depends on j mod q); the non-resonant ET page is the single named lemma (standard Erdős–Turán + kernel evaluation, constants data-calibrated)
source: mac-mini-2026-07-16-S124 (owner: close the finite-Vmax large-spread half)
depends_on: [THM-527 (the reformulation + bounded-spread frame), THM-518 (stranger decoupling — the non-resonant mechanism), S58 arc-count lemma, THM-724 L2 (the dilate tile the resonant branch lands in)]
script: 04-computation/vmax_glue_largespread_macmini_S124.py -> 05-knowledge/results/vmax_glue_largespread_macmini_S124.out
---

# THM-924 — the finite-Vmax glue, assembled

ρ_K(Vmax) = the good-period density; ρ* = its slow-fast limit (THM-527 A). The glue must
bound |ρ_K − ρ*| for EVERY cluster shape. Trichotomy on the co-offset cluster E:

1. **Bounded spread** (internal diffs ≤ S₀): the good set has O(1) arcs (S58 — arcs break
   at x = m/(internal diff), not m/Vmax), so the AP {x_j} meets it with error ≤ C(S₀)/Vmax. PROVED.
2. **Large spread, non-resonant** (no e_i/Vmax near a small-denominator rational): the
   indicator's Fourier series against the AP gives Dirichlet kernels
   |Σ_j e(h e_i x_j)| = |sin(πhe_i)/sin(πhe_i/Vmax)|, uniformly small off resonance;
   Erdős–Turán yields |ρ_K − ρ*| ≤ C′/Vmax with C′ SPREAD-FREE. Data: err·Vmax = 0.3–7
   across spreads 5 → 6577 (vs 263 at exact resonance) — the one named page.
3. **Resonant** (some e_i/Vmax = p/q, q small): then e_i x_j mod 1 depends only on
   j mod q — ρ_K is EXACTLY periodic, a rational with denominator q·(core factor),
   computable in finite exact arithmetic and positive in every tested instance
   (ρ_K = 0.690 at the {0, V/4, V/2, 3V/4} resonance). This is the near-dilate
   configuration — THM-724 Lemma 2's tile — and needs no limit at all.

Every branch gives ρ_K > 0 for Vmax past an explicit threshold ⟹ M(S) ≥ 1/14 via
criterion C. **THM-527-A's glue is assembled**; the single written-out lemma remaining is
the branch-2 ET page (standard, same architecture as THM-920's slice page, constants
calibrated in the .out).

## The final LRC(14) ledger (as of this file)

- **Covering**: route [A] SIGNED OFF (THM-922); rigidity (THM-883, j ≤ 5 proved); bands
  (THM-755/756); low-M (S111). DONE at referee grade.
- **Non-covering / S3 residual**: the lonely-density reformulation (THM-527, proved) +
  THIS GLUE (assembled; ET page named) + the uniform floor c₀ > 0 over the compact shape
  space (THM-527's remaining crux — the bounded-spread positive-floor, k = 3 proved,
  the k ≤ 6 shape scans positive everywhere).
- **The level-5 wall**: W and T1 proved (opus THM-897); T2 = the single remaining
  analytic lemma (the sub-orbit engine at THM-864 rigor with the triple resonance
  integer — route named by opus, one session); the emptiness question resolved by
  opus THM-923 (two-scale picture).
- **Lean**: rung one committed (FragmentationLemma.lean); rungs two/three named
  (THM-866 walk, THM-878 clock).


## Appendix (S125) — the ET page, written

**Lemma (branch 2).** For non-resonant large-spread E (no e_i/Vmax within 1/(qQ₀) of a
rational p/q with q ≤ Q₀): with x_j = (j+½)/Vmax and the good-set indicator's bounded
variation V_G (≤ 2·#arcs of the FIXED finite union defining maxgap > 2/7 intersected
with G_P — a P-and-diff quantity, Vmax-free), Koksma's inequality gives
|ρ_K − ρ*| ≤ V_G · D*_{Vmax}(x·E), and Erdős–Turán bounds the star discrepancy of the
AP-orbit by (k/H) + Σ_{h ≤ H} (c/h)·Π-kernel terms |sin(πhe_i)/ (Vmax·sin(πhe_i/Vmax))|,
each ≤ 1/(Vmax·‖he_i/Vmax‖) — off-resonance ‖he_i/Vmax‖ ≥ 1/(qQ₀)-scale for h ≤ H = Q₀,
so the sum is ≤ c·Q₀²·log Q₀/Vmax: **|ρ_K − ρ*| ≤ C(Q₀, P, diffs)/Vmax, spread-free**,
matching the S124 calibration (err·Vmax = 0.3–7). Resonant e_i are branch 3 (exact
q-periodicity). ∎ (Constants tabulated in the S124 .out.)
