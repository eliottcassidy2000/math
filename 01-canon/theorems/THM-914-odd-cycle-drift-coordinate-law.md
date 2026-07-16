---
id: THM-914
title: THE ODD-CYCLE DRIFT COORDINATE LAW (n=4..7, machine-exact over ALL tournaments) — for the uniform one-arc-flip walk, the Hamiltonian-path drift E[ΔH | T] is a class function of the ODD-CYCLE CENSUS (H; c₅, c₇): (n=4) 6·E[ΔH] = 12 − 4H, affine in H alone; (n=5) 5·E[ΔH] = 12 − H − 6c₅, EXACT affine — c₅ enters with weight 6 (explains the split H=15 fiber: c₅=3 vs 2); (n=6) NOT affine in (1,H,c₅), but drift = f(X) EXACTLY under the UNIVERSAL COORDINATE X = H − 6c₅ — the 29 (H,c₅) fibers COMPRESS onto X-lines with equal drift: (1,0),(13,2) → 52/15; (3,0),(9,1),(15,2) → 4; (5,0),(11,1) → 68/15; (9,0) → 12/5 (the n=5 gradient (−1,−6) survives as the LEVEL-SET direction while the response goes nonlinear); (n=7) (H,c₅,c₇) determines the drift EXACTLY — 203 tuples, ZERO split fibers over all 2^15 tilings — and the triple is MINIMAL: every sub-coordinate fails ((H,c₅): 16 collisions, (H,c₇): 19, (c₅,c₇): 62, H alone: 88, X=H−6c₅: 166), so c₇ GENUINELY ENTERS at n=7; no linear compression X = H − 6c₅ − w₇c₇ fibers (w₇ over the natural candidate lattice). COVERAGE IS TOTAL: H, c₅, c₇, and the drift are isomorphism invariants and every tournament has a Hamiltonian path (Rédei), so exhaustive-over-tilings = exhaustive over ALL tournaments at each n. READING: the OU law (THM-833) made the c₃-drift exactly affine (fluctuation–dissipation); the H-drift is affine only through n=5, but its ENTIRE state-dependence is channeled through the odd-cycle counts — the odd-cycle census is a SUFFICIENT STATISTIC for the H-drift. This is the drift-side shadow of THM-466 (H ≡ 1 + 2c_odd mod 4) and the OCF: odd cycles are the only coordinates the flip walk feels
status: PROVED at the census level for n = 4..7 (machine-exact, exhaustive over all tournaments via the tiling fibration; scripts cited below). The all-n sufficiency — E[ΔH|T] = f_n(H, c₅, c₇, …, c_{odd ≤ n}) for every n — is HYP-7095 (open)
source: kind-pasteur-2026-07-16-S128 (cont.30–32; owner: find the H-drift's coordinate; verify the odd-cycle functional at n=6,7 and pin whether c₇ enters; harvest the n=7 verdict and canonize)
depends_on:
  - THM-833   # the OU flip atom Δc₃ = d_u − d_v − 1 and the c₃ fluctuation–dissipation law
  - THM-466   # H ≡ 1 + 2c_odd (mod 4) — the parity face of the same census
related:
  - THM-854 (F3 exchange walk — the same flip atom walking c₃ down to transitive)
  - THM-852/853 (self-line atlas — the class functions this drift law lives beside)
  - HYP-6985 (the discovery thread, now CONFIRMED n≤7 and canonized here)
  - HYP-7095 (the all-n sufficiency conjecture)
  - OCF / Rédei (H odd; the census as the fundamental parity object)
---

# THM-914 — the odd-cycle drift coordinate law

**Setting.** Fix n and the base path n → n−1 → … → 1. For a tournament T let H = #Hamiltonian
paths, c_l = #directed l-cycles. The walk flips one uniformly random arc (C(n,2) choices);
E[ΔH | T] is the conditional one-step drift of H.

**Law (n = 4..7, exhaustive).**

| n | statement | form |
|---|-----------|------|
| 4 | 6·E[ΔH] = 12 − 4H | affine (H) |
| 5 | 5·E[ΔH] = 12 − H − 6c₅ | affine (H, c₅) |
| 6 | E[ΔH] = f(X), X = H − 6c₅ | exact, nonlinear |
| 7 | E[ΔH] = f(H, c₅, c₇), minimal triple | exact, nonlinear |

n=6 compression (all values verbatim, drift × 1): X=1 → 52/15, X=3 → 4, X=5 → 68/15,
X=9 → 12/5 — tuples (1,0),(13,2) / (3,0),(9,1),(15,2) / (5,0),(11,1) / (9,0) respectively.
n=7: 203 (H,c₅,c₇) tuples, 0 split fibers; sub-coordinate collision counts (H,c₅):16,
(H,c₇):19, (c₅,c₇):62, H:88, X:166 — the triple is minimal and c₇ enters.

**Why coverage is total.** H, c₅, c₇, E[ΔH] are iso-invariants; Rédei gives every tournament a
Hamiltonian path, hence a representative among the 2^{C(n−1,2)} tilings; the sweeps enumerate all
tilings. So the law is verified on every tournament of each order, not a sample.

**Interpretation.** THM-833's c₃-drift is exactly affine — perfect Ornstein–Uhlenbeck. The
H-drift loses affinity at n=6 but NOT coordinate-freeness: everything the walk feels about T is
carried by the odd-cycle census. The n=5 gradient direction (1, 6) persists at n=6 as the
level-set (fiber) direction of X = H − 6c₅; at n=7 the gradient rotates out of the (H, c₅) plane
into c₇. Conjecture (HYP-7095): at every n the census (H, c₅, …, c_{odd≤n}) remains a sufficient
statistic — the flip walk on tournaments is an odd-cycle-census Markov chain in disguise.

## Evidence log
- [x] n=4,5 affine laws exact (hdrift_coordinate_kps_S128c30.py, cont.30)
- [x] n=6: 29 tuples, 0 splits, no affine law, X-compression exact
      (oddcycle_functional_n67_kps_S128c31.py + inline 3-unknown refit, cont.31)
- [x] n=7: 203 tuples, 0 splits; minimality table; linear-X refutation
      (xcoord_n7_test_kps_S128c32.py + .out; tuples saved to n7_drift_tuples.json)
- [ ] HYP-7095: all-n sufficiency; n=8 census (2^21 tilings — needs the bit-packed engine)
- [ ] structure of f at n=6,7: candidate — drift as census-polynomial with Δc-atom weights
      (THM-833 gives Δc₃ exactly; the analogous Δc₅, Δc₇ atoms are the natural next objects)
