# T_11 from-scratch Betti verification — monad-compute progress notes

**Script:** `04-computation/verify_t11_betti_s_monad.py`
**Live output:** `05-knowledge/results/verify_t11_betti_s_monad.out`
**Task:** INV-143 next step — recompute full T_11 GLMY Betti from scratch, verify β_5=5, β_6=15.

## Confirmed from scratch (use_cache=False)
- raw |A_m| (directed path counts) = `[1, 5, 25, 110, 430, 1430, 3970, 8735, 14395, 15745, 8645]`
- **Ω dims = `[1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]` — MATCHES cached, χ = 1.** ✓
- root-field prime = 23 (no int64 overflow risk: 23² · 15745 ≈ 8.3e6 ≪ 9.2e18).

## Eigenspace boundary ranks rank(d_m^(k)), m=0..11
- **k=0:** `[0, 0, 5, 15, 55, 150, 305, 390, 300, 150, 30, 0]`  (388.3 s)
- **k=1:** `[0, 1, 4, 16, 54, 151, 309, 390, 300, 150, 30, 0]`  (381.8 s)

### Per-eigenspace Betti contributions β_m^(k) = (Ω_m − rank_m) − rank_{m+1}
- **k=0:** `[1, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0]`
- **k=1:** `[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]`  ← contributes **+1 to β_6 only**

## FINAL RESULT — COMPLETE (S2, all 11 eigenspaces, from scratch)
- **All k=0..10 boundary ranks computed from scratch** (use_cache=False). The 10
  non-principal eigenspaces k=1..10 ALL have IDENTICAL ranks
  `[0,1,4,16,54,151,309,390,300,150,30,0]`; k=0 is `[0,0,5,15,55,150,305,390,300,150,30,0]`.
- **β = `[1,0,0,0,0,5,15,0,0,0,0]` — MATCH `KNOWN_BETTI[11]`. β_5=5, β_6=15. ✓ CONFIRMED.**
- Per-eigenspace β_6 contributions = `[5,1,1,1,1,1,1,1,1,1,1]` → 5 + 10·1 = **15**.
- Per-eigenspace β_5 contributions = `[5,0,0,0,0,0,0,0,0,0,0]` → **5** (k=0 only).
- Timings per non-principal eigenspace: 384–413 s (cache-bloat fixed: clearing
  `_omega_basis_cache` per eigenspace held timing flat; no slowdown observed).

### Structural finding → predicts the cached Betti exactly (refines HYP-453)
- β_5 = 5 comes *entirely* from k=0; non-principal eigenspaces add 0 to β_5.
- β_6: k=0 contributes **5**, and each non-principal eigenspace contributes **+1**.
  All 10 eigenspaces k=1..10 verified identical → **β_6 = 5 + 10·1 = 15**.
- This **refines HYP-453** ("all T_11 homology at k=0"): correct for β_5, but β_6 is
  distributed — k=0 carries 5 and every eigenspace carries 1 more.

### Euler-characteristic cross-check clarification (script "MISMATCH" was a FALSE ALARM)
- The first OVERALL print said MISMATCH because it compared χ_Betti=Σ(−1)^m β_m=**11**
  against the single-copy χ_Omega=Σ(−1)^m Ω_m=**1**. These are NOT supposed to be equal.
- **By THM-125 each of the n=11 eigenspaces carries a FULL copy of Ω_m** (dim Ω_m^(k)=omega[m]
  for all k; see `circulant_homology.betti_numbers` docstring). So the total chain dim at
  degree m is n·Ω_m and the correct Euler identity is **χ_Betti = n·χ_Omega = 11·1 = 11**. ✓
- Fixed the script's final check to compare χ_Betti against n·χ_Omega. The Betti numbers
  themselves were never in doubt — they reproduce the library's official `betti_numbers()`
  formula exactly and match `KNOWN_BETTI`.

## Performance / handoff
- One eigenspace ≈ 384–413 s (pure-Python null-basis on n_junk × |A_m| matrices, |A_9|=15745).
- Full 11-eigenspace from-scratch run ≈ **70 min** on this node — exceeds the 30-min
  compute-node budget, so this was run in background with per-eigenspace JSON persist +
  auto-commit (`verify_t11_betti_s_monad_ranks.json`, resumable). Recommend a C/LinBox
  reimplementation for routine re-verification (same conclusion as INV-141: degree-9+ needs C/C++).
- Cache-bloat fix (clear `_omega_basis_cache` per eigenspace) resolved the S1 slowdown.
