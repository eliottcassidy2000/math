---
id: HYP-3611
title: APEX CORE ATLAS -- fine patterns of the Z_p core landscape: (1) the gap g(O) is a CONCENTRATION INDEX ordering cores by spread, binding = the LEAST-spread doublet (opposite of the Paley/random core); (2) FLAT spectrum <=> INTEGER gap {0,1,2} (difference-set/singleton/full) vs IRRATIONAL binding gaps {0.198,0.308} (non-flat); (3) cores ARE circulant graphs (doublet=C_7 cycle, difference-set=K_7, singleton=I, full=7J, m-arc=Dirichlet/FEJER kernel); (4) the ARC cores are exactly the Fejer-Bochner minorant cores, doublet=2-arc with gap 4sin^2(pi/2p)=the Fejer kernel minimum; (5) cross-prime: doublet gap ~ pi^2/p^2, p=3mod4 has Paley flat-QR (gap (p+1)/4 = max), bracket [4sin^2(pi/2p), (p+1)/4]; (6) NEW TRACKABLE the descent gap-profile, near-universal tail [...,0.308,0.198,1.0], head varies (cusp 0.0 vs co-singleton 1.0)
status: VERIFIED (exhaustive over 127 Z_7 cores; cross-prime p=3..19; arc=Fejer verified p=7,11,13; descent profiles over standard+binding coverings). Reference artifact: 05-knowledge/reference/apex-core-atlas.md.
source: klein-2026-06-29-S21
depends_on:
  - THM-590    # the apex gap (g(O) = lambda_min of the autocorr circulant)
  - HYP-3604   # the least-eigenvalue certificate (doublet = signless Laplacian of C_7)
related:
  - HYP-3606   # mac-mini: non-bipartiteness form
  - HYP-3598   # the finite families (all 127 cores arise)
  - HYP-3612   # the gap localizes to level 0 (the descent gap-profile head)
  - HYP-3593   # the Paley / nu_2 = 3 mod 4 theme
results:
  - 04-computation/apex_core_atlas_klein.py
  - 04-computation/apex_cross_prime_family_klein.py
  - 04-computation/apex_core_graph_map_and_profiles_klein.py
  - 05-knowledge/results/apex_core_atlas_klein.out
  - 05-knowledge/results/apex_cross_prime_family_klein.out
  - 05-knowledge/results/apex_core_graph_map_and_profiles_klein.out
  - 05-knowledge/results/apex_arc_fejer_kernel_klein.out
---

# HYP-3611 — the apex core atlas: fine patterns & new trackables

Full reference: `05-knowledge/reference/apex-core-atlas.md`. The new fine patterns:

## 1. The gap is a CONCENTRATION INDEX
`g(O)` orders the cores by Fourier SPREAD: `doublet (0.198) < arc-3 (0.308) < singleton (1) < diff-set (2)`,
full-`Z_7` = 0. **LOW gap = concentrated; HIGH gap = spread.** The binding apex obstruction is the
LEAST-spread core (the doublet), the OPPOSITE of the random/Paley core. The hardest case is maximal
concentration, not equidistribution.

## 2. FLAT spectrum <=> INTEGER gap
Cores with a flat nonzero spectrum (difference sets, singletons, full set) have INTEGER gaps `{0,1,2}`; the
non-flat cores (doublets, generic triples) have the IRRATIONAL binding gaps `{0.198, 0.308}` in
`Q(cos 2pi/7)`. Integer/flat = difference-set-like; irrational = the binding obstruction.

## 3. Cores ARE circulant graphs
`C(O) = [a(i-j)]` is a circulant/Cayley graph; `g(O)` is its least eigenvalue. doublet -> `C_7` (cycle,
`Q(C_7)`); difference set -> `K_7` (`I+J`); singleton -> `I`; full -> `7J` (rank-1, cusp). The gap classes
ARE named graph types.

## 4. ARC cores = the Fejér–Bochner minorant cores
The `m`-arc `{0,..,m-1}` has `|Ohat(k)|^2 = (sin(pi m k/p)/sin(pi k/p))^2` = the **Fejér kernel** (verified
p=7,11,13). So the arc cores are exactly the floor owners' Fejér–Bochner minorant (S75e). The binding
**doublet IS the 2-arc**, and its gap `4 sin^2(pi/2p)` is the Fejér kernel's smallest value. This unifies
the least-eigenvalue certificate with the analytic minorant: same object, the 2-arc Fejér minimum.

## 5. Cross-prime family
Doublet gap `= 2-2cos(pi/p) = 4 sin^2(pi/2p) ~ pi^2/p^2` (smaller odd prime = larger floor; p=7 -> 0.198,
3rd-largest). For `p = 3 mod 4` (3,7,11,19) the QR is a **Paley difference set** (flat, gap `(p+1)/4` = the
MAX); for `p = 1 mod 4` (5,13,17) it is not. `7 = 3 mod 4` is the Paley case. The gap spectrum is bracketed
`[4 sin^2(pi/2p), (p+1)/4]`. # distinct gap values: `1,2,4,14,35` for `p=3,5,7,11,13`.

## 6. NEW TRACKABLES
- **gap-profile** `[g(O_j)]_j` of a covering (descent signature): near-universal TAIL
  `[..., 0.308, 0.198, 1.0]` (arc, doublet, singleton); HEAD (level 0) distinguishes coverings (cusp `0.0`
  for dense-consecutive / fills `Z_7`, vs `1.0` co-singleton). Binding atom universally the doublet `0.198`.
- **concentration/spread index** `g(O)` (normalized `g/((p+1)/4)` in `[0,1]`).
- **core shape type** (doublet/arc/V/difference-set/singleton/full).
- **flat-vs-nonflat** (integer vs irrational gap).
- **level-0 apex type** (cusp vs proper).

## Net
A single concentration index `g` grades the apex problemscape, bracketed by the binding doublet (2-arc,
Fejér minimum, signless-Laplacian of the odd cycle, positive iff p odd) and the Paley difference set
(`K_p`, flat, p=3 mod 4). The descent threads coverings through it with a universal gap-profile tail. All
fine features now mapped in the atlas.
