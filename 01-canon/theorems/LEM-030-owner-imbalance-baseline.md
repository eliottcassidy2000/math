---
id: LEM-030
title: THE OWNER-IMBALANCE BASELINE — the cross-owner part of Q_s decomposes as cross(w) = −(π²/3)·Σ_e ν̂_e(0)² − 2π²·Σ_{cross-owner pairs} ε_kε_k′·(K({wΔ}) − 1/6): an ALWAYS-NEGATIVE baseline equal to the combs' DC-tooth imbalance (σ_e = ν̂_e(0), the owner's signed endpoint sum), plus a mean-zero fluctuation. Proof: (a) the spectral kernel is exactly the B₂ Fourier series (THM-892(K)), so the partial bilinear carries B₂ = 1/6 − K, whose constant does NOT cancel on partial sums; (b) Σ_{cross ordered} εε′ = (Σσ)² − Σσ² = −Σ_e σ_e² by Σε = 0 (pure algebra); (c) σ_e = ν̂_e(0) by definition of the comb teeth. Verified exactly: balanced cluster Σσ² = 64, baseline −(π²/3)·64 = −210.6, pair-K-part +102.8, total −107.8 = the spectral cross to 3 decimals
status: PROVED (three steps, each ≤ 3 lines; machine-verified: per-owner σ = ν̂(0) exact, Σσ² = 64, the assembled identity matches the spectral value). CONSEQUENCES: (i) the off-resonant negativity of the cross (S53's refined diagonal-dominance law) = this baseline + a mean-zero fluctuation whose off-resonance smallness is a Q_s-type discrepancy statement (klein's lane); (ii) the resonant flip = K-concentration overwhelming the baseline (measured, S54); (iii) v2's sharp successor: [Σ_e D_e − (π²/3)Σν̂_e(0)² + fluctuation-bound] off-resonance | [comb closed form] on — fully specified
source: boxeph-2026-07-16-S55 (owner: the baseline lemma)
depends_on: [THM-880 (bilinear form), THM-892(K) (kernel = B₂ exactly), THM-887(I) (comb teeth)]
---

# LEM-030 — the owner-imbalance baseline

The DC mode is conserved globally (Σε = 0) but not per owner: each owner's comb carries
a DC tooth ν̂_e(0) = σ_e, and the cross-owner bilinear pays their imbalance as a
guaranteed-negative term −(π²/3)Σσ_e². The sign of the cross is thus structural
off-resonance and only the mode-concentration at resonant w can flip it — the exact
mechanism split the refined diagonal-dominance law asked for. (The partial sum keeps
the constant the global sum forgets: the arc's closing lesson, now a lemma.)
