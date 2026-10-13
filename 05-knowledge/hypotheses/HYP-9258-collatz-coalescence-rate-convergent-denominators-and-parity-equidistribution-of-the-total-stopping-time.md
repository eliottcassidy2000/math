---
id: HYP-9258
title: "Collatz coalescence rates and signed resonance: exact pair survival through depth48; ROOT-supported parity mean zero from coalescence; half-integer transfer spectrum tested through2^34; quantitative mode expansion and marked-source entry remain open"
status: >
  OPEN: quantitative coalescence rate, signed modal expansion, amplitude law, and universal fixed-source entry.
  PROVED from THM-4581's almost-sure coalescence input: the ROOT-supported parity function, set to zero off
  the ROOT basin, has ordinary mean zero. The new bounded-antifunction lemma proves this from vanishing
  adjacent variation using incompatible signs under dilation by2 and3; no sketch-level rate is used.
  FINITE-EXACT: two full parity censuses through2^34 agree on every1024-bin/octave sum; an independent
  jump16 run at4096 bins/octave verifies aggregation and exposes high-frequency aliases. Pair-chain
  unmerged mass is computed exactly through48 by forward/backward composition, with independent
  integer-cylinder and edge-multiplicity controls. NUMERICAL: the signed extension of HYP-9230 predicts
  low frequency2.454923470 and decay0.212159672, compared with2.45485663 and0.21247561 on log2(n)=24..34;
  half drift2.40942084 is its linear approximation and5/2 its unshifted label. The strongest resolved
  fine-bin mode is near26.5036, outside the original scan. The all-frequency expansion remains unproved.
  The earlier exact-half-drift/two-mode conjecture and low-band claim of a globally dominant2.455 mode
  are SUPERSEDED as working models, not refuted by a universal exact spectral theorem.
  See05-knowledge/results/collatz_signed_resonance_pair_mass_20261013b.md for proof, scope and data.
session: opus-2026-10-13-S25
source: 05-knowledge/results/five_manuscripts_chowla_collatz_residue_bsd_20261013.md (sections 3.3-3.4)
scripts:
  - 04-computation/experiments/collatz_tau_runs_20261013.py (+ .out)
  - 04-computation/experiments/collatz_merge_lags_20261013.py (+ .out)
  - 04-computation/experiments/collatz_tau_diffs_20261013.py, collatz_tau_oscillation_20261013.py, collatz_tau_spectrum_20261013.py, collatz_tau_drift_20261013.py (+ .out)
related:
  - THM-4627 (the density-one coalescence and the correlation limits)
  - THM-4581 (P(no merge by T) between c T^(-1/2) and C T^(-1/2) (log T)^2, both bounds at sketch level)
  - HYP-9217 (debt recurrence law, exponent 1/2)
  - HYP-9230 (the signed model has exactly the same characteristic root on the half-integer lattice; actual modal expansion OPEN)
  - Garner, On heights in the Collatz 3n+1 problem, Discrete Math. 55 (1985) 57-64 (consecutive integers with equal total stopping time; CITED via the search of 2026-10-13); LaDue, Clusters of integers with equal total stopping times in the 3x+1 problem, Fibonacci Quart. 56 (2018), arXiv:1709.02979 (merge criterion for consecutive even-odd pairs; the identity sigma_inf(8i+4) = sigma_inf(8i+5), which is the depth-3 mass 1/8 of THM-4627); Elia and Tucker, Consecutive integers and the Collatz conjecture, arXiv:1511.09141 (counterexamples to Garner's conjectural family); Lagarias's survey attributes the clustering to coalescence of trajectories; no source found gives the density of n with tau(n) = tau(n+1), nor a parity bias
---

# HYP-9258 — rates and residues of Collatz coalescence along lags

## Statements

(i) **Rate.** `beta(N) := (1/N) #{n <= N : tau(n) != tau(n+1)}` tends to `0` (THM-4627) and `beta(N) = N^(o(1))`-slowly: conjecturally `beta(N) ~ c (log N)^(-gamma)` with `gamma` between `1/2` (THM-4581's merge law with `T` of order `log N`) and the observed slower decay; the mean length of a maximal run of consecutive integers with equal `tau` is `1/beta(N)`.
(ii) **Residues.** Among `n` with `tau(n+h) != tau(n)`, the difference is `19a + 8b` with small integers `a, b` (both signs), the `19`-direction favoured over the `8`-direction and the frequency falling with the size of the coefficients: the ten commonest absolute values at every lag `h <= 4` are `19, 38, 27, 46, 8, 57, 65, 84, 11, 30`, i.e. `(a,b) = (1,0), (2,0), (1,1), (2,1), (0,1), (3,0), (3,1), (4,1), (1,-1), (2,-1)`; below the Frobenius number `125` of the semigroup `<8, 19>`, `84.2%` of the differences lie in the semigroup `a, b >= 0`, and the commonest exceptions `11, 30, 3, 49, 14, 22, 68` are `19a - 8b` with `a <= 4`, `b <= 3`. The general law should be that unequal-time coincidences occur at time offsets that are signed combinations of the convergent and semiconvergent denominators of `log_2 3`, the shortest first.
(iii) **Signed resonance, revised2026-10-13.** The [new note](../results/collatz_signed_resonance_pair_mass_20261013b.md) proves a qualitative implication: a bounded Terras antifunction with vanishing adjacent variation has ordinary mean zero. Applied using THM-4581, the difference between even- and odd-total-stopping-time ROOT counts is `o(N)`; this does not assert ROOT density1. Conditional on Collatz the two densities are1/2. The quantitative modal expansion remains OPEN. Its model equation is `exp(-lambda)+exp((log_2(3)-1)lambda)=-2`, giving HYP-9230's root at half-integer base frequencies. At base5/2 the predicted frequency is2.454923470 and decay0.212159672; finer censuses through2^34 favour these over both2.40942084 and2.5. The half-drift value is the linear approximation to this same root, so the old two-independent-mode explanation is superseded. The strongest resolved fine-bin mode is near26.5036; high-frequency coefficients must retain their bandwidth and aliasing sidecars.

## Decisive tests and remaining work

- **Completed FINITE-EXACT:** unmerged masses for every depth through48, extending the requested26..40 range. At48, `q_48=150619085192045/2^48`. The finite effective exponent is about0.13, not an asymptotic-rate determination.
- **Completed NUMERICAL comparison:** 1024 and4096 bins/octave through2^34, two exact jump algorithms, the signed characteristic equation and phase/decay fits. A scale holdout and the finer-bandwidth hostile are retained in the note.
- **OPEN:** a rate for the signed residue discrepancy that survives refinement and implies the proposed modal expansion; rigorous renewal tail bounds; a source-specific decreasing certificate. The existing difference-count data in(ii) remain evidence, not a new theorem.

## Why

THM-4627 makes `(-1)^tau` locally constant along lags on a density-one set; the three statements quantify the complement: how fast the exceptional density decays, what the exceptional differences are, and how the global average, which is not locally random, oscillates. The coarse parity bias excludes a simple independent-fluctuation explanation. The refined signed-scaling mechanism now separates a proved qualitative mean implication from the unproved rate and pointwise entry questions.
