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

### Structural finding → predicts the cached Betti exactly
- β_5 = 5 comes *entirely* from k=0; non-principal eigenspaces add 0 to β_5.
- β_6: k=0 contributes **5**, and each non-principal eigenspace contributes **+1**.
  With 10 eigenspaces k=1..10 (k=1 verified, k=2..10 pending) this predicts
  **β_6 = 5 + 10·1 = 15** — matching cached `KNOWN_BETTI[11]`. ✓ (pending k=2..10 to
  confirm every non-principal eigenspace follows the +1 pattern).
- This **refines HYP-453** ("all T_11 homology at k=0"): correct for β_5, but β_6 is
  distributed — k=0 carries 5 and every eigenspace carries 1 more.

## Performance / handoff
- One eigenspace ≈ 388 s (pure-Python null-basis on n_junk × |A_m| matrices, |A_9|=15745).
- Full 11-eigenspace from-scratch run ≈ **70 min** on this node — exceeds the 30-min
  compute-node budget. Recommend a C/LinBox reimplementation for routine re-verification
  (same conclusion as INV-141: degree-9+ needs C/C++).
- Background run continues; later eigenspaces append to the .out and are committed as they land.
