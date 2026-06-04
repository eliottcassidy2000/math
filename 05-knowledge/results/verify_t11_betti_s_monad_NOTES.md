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

### Betti contribution of k=0 alone
β_m^(k=0) = (Ω_m − rank_m) − rank_{m+1}:
`[1, 0, 0, 0, 0, 5, 5, 0, 0, 0, 0]`

**Finding:** β_5 = 5 comes *entirely* from k=0. But k=0 contributes only **5** to β_6,
while the cached total is β_6 = **15**. So the remaining **+10** must come from the
non-principal eigenspaces k=1..10 (≈ +1 each). This refines HYP-453 ("all T_11 homology
at k=0"): true for β_5, but β_6 receives contributions from every eigenspace.

## Performance / handoff
- One eigenspace ≈ 388 s (pure-Python null-basis on n_junk × |A_m| matrices, |A_9|=15745).
- Full 11-eigenspace from-scratch run ≈ **70 min** on this node — exceeds the 30-min
  compute-node budget. Recommend a C/LinBox reimplementation for routine re-verification
  (same conclusion as INV-141: degree-9+ needs C/C++).
- Background run continues; later eigenspaces append to the .out and are committed as they land.
