---
id: HYP-2711
title: LRC14 phase-carrier analogy filter for Gibbs, propagator, Hopfield, road-coloring, crossing, and Clifford lenses
status: OPEN; synthesis and falsification ledger
source: codex-2026-06-20-S67
tangent: T942
depends_on:
  - HYP-2705
  - HYP-2707
  - HYP-2706
  - HYP-2704
  - HYP-2701
  - HYP-2702
  - HYP-2698
  - HYP-2684
  - HYP-2675
  - OPEN-Q-108
related:
  - THM-557
  - THM-558
  - THM-555
  - THM-554
  - THM-071
  - HYP-2694
  - HYP-2695
  - HYP-2696
  - HYP-2697
---

# HYP-2711 - Phase-Carrier Analogy Filter

## Claim

The user's Gibbs/cat-map/Fubini-Study/road-coloring/Hebbian/propagator/
Clifford/crossing prompts are useful only after being filtered through the
repo's exact carriers.  The surviving architecture is:

```text
exact carrier identity
  -> generated phase/mask profile
  -> signed deviation from the free/decorrelated quotient
  -> finite low-rank atlas or analytic projective-angle bound.
```

The dangerous mistake is to treat an analogy as a proof.  HYP-2711 instead
classifies each analogy as one of:

```text
EXACT:       a theorem/identity already present in the carrier,
COORDINATE:  a useful metric/basis but not a bound,
SIGN-ONLY:   predicts direction but not magnitude,
DEAD-END:    refuted as a literal proof route.
```

The LRC payoff is a sharper version of HYP-2705: the remaining true-wide
deviation should be bounded by phase degree / magic rank / Fubini-Study angle
of the mod-7 relation profile, while low-rank defects route to finite AP,
cube-root, Freiman, squarefree, and synchronizer atlases.

## Evidence Ledger

### Exact or Useful

1. **LRC path integral is exact at the indicator level.**
   `lrc14_pathintegral_amplitude_macmini_0620s8.py` verifies that the
   all-seven-colors indicator is exactly an alternating missing-set sum and
   also an additive-character expansion on `Z/7`.  This is a genuine
   phase/interference identity, not metaphor.

2. **Death-chain Gibbs transfer matrix is exact.**
   `gibbs_variational_extremality_macmini_S7.py` writes the HYP-2704
   death-chain recurrence as a lower-triangular transfer operator with
   eigenvalues `1,6/7,5/7,...,1/7`.  The spectral gap is `1/7`, and the
   decorrelated cover is the absorbing covered coordinate.

3. **Hopfield energy is exact after signed incidence correction.**
   `hopfield_ising_signed_incidence_opus_s_hop.py` proves that score variance
   is a two-body Ising Hamiltonian on the line graph only with signed
   incidence couplings `+/-1/2`, not with the naive uniform line-graph
   coupling.  This is the exact cut-layer/Hopfield weight matrix.

4. **Clifford/magic is exact for tournaments through HYP-2707.**
   Incoming HYP-2707 proves the core instance: `c3 mod 2` is a GF(2)
   quadratic Clifford Gauss sum, and higher OCF/H layers are magic.  HYP-2711
   imports that lesson into LRC as a proposed phase-degree measure for mod-7
   relation corrections, not as raw relation count.

5. **Even-page crossing carrier is real but not decisive.**
   The crossing scripts verify that the 2-page complete-graph optimum is
   attained by an even-graph page through `n=8`; this ties crossing carriers
   to the cycle/even-graph side of the repo.  It does not by itself bound LRC.

### Dead Ends or Guardrails

1. **Literal Arnold cat map is wrong.**
   `gibbs_catmap_partition_function_macmini_S7.py` finds the mod-7 multiplier
   action is finite-order/torsion, while the band map is an expanding base-7
   endomorphism, not a unimodular hyperbolic toral automorphism.  Use
   transfer operators and Markov partitions, not literal cat-map dynamics.

2. **Strict road coloring is wrong for the tiling hypercube.**
   The LRC road is a coverage/reachability process, not Cerny-style state
   collapse by an invertible map.  The useful automaton is the non-invertible
   residual deletion automaton in HYP-2698/HYP-2702/HYP-2705.

3. **Gibbs convexity does not prove extremality.**
   `gibbs_variational_extremality_macmini_S7.py` records that
   `log sum exp(beta E)` convexity is tautological.  Meanwhile `measS7` is not
   concave along natural top-slide deformations; consec is a boundary
   extremizer, not a convex variational maximum.

4. **Raw path coherence is sign-level only.**
   `lrc14_pathintegral_decisive_macmini_0620s8.py` shows the raw clock
   coherence and `measS7` share the consecutive peak in some boxes, but their
   top non-consecutive lists are disjoint and the interference/off-diagonal
   margin ratios vary from about `11.6` to `77.0`.  It predicts direction, not
   magnitude.

5. **Low-order log-linear free energy does not capture `measS7`.**
   `gibbs_aggregate_order_test_macmini_S7.py` uses overdetermined train/test
   fits.  Pairwise and triple models do not generalize to `R^2=1`; the
   aggregate obstruction is really the log of a sum over clock states.

## External Source Mapping

The arXiv paper
[`1803.07039`](https://arxiv.org/pdf/1803.07039) supplies the quantum
Hebbian/QSE analogy: batches of quantum data simulate a density matrix, and
the controlled partial-swap operation is decomposed into Clifford+T gates.
HYP-2711 uses only the structural lesson: density-matrix/Hebbian weights are
good for the cut layer, while the non-Clifford correction is the signed
relation/cycle layer.

The
[propagator-as-weights article](https://medium.com/@denizaskin/feynman-propagators-as-weights-a-quantum-neural-network-architecture-b65415ff41cc)
supplies the double-slit dictionary: slits as classes and propagator
amplitudes as weights.  HYP-2711 accepts this only at the exact
character-expansion level.  The LRC slits are mod-7 residual classes;
interference signs exist, but raw coherence is not the proof metric.

## Proof Program

1. **Define the LRC phase degree.**  Build a normalized mod-7 phase profile
   `psi_E` whose coordinates are generated residual masks / relation packets.
   The stabilizer layer is the death-chain/decorrelated profile; the magic
   layer is the signed deviation.

2. **Prove a projective deviation bound.**  Use HYP-2705's inequality
   `|L(psi)-L(phi)| <= 2||L|| sin(d_FS(psi,phi)/2)` with a number-theoretic
   estimate for `d_FS(psi_E, psi_0)` outside low-rank relation atlases.

3. **Classify low magic rank.**  Route small phase-degree defects to the
   existing finite atlases: AP-prefix, cube-root/AP-triple, Freiman/GAP,
   squarefree crossing carrier, and residual synchronizer.

4. **Keep the generated profile until the last comparison.**  HYP-2704 and
   HYP-2706 show that scalar/per-band/FOSD quotients fail.  The proof should
   compare signed generated profiles against the HYP-2701 boundary margin,
   then scalarize.

## Status

No LRC14 proof is claimed.  HYP-2711 is a research filter and proof-program
refinement: it turns the new analogy batch into exact identities, dead-end
warnings, and one actionable next target, namely a phase-degree/FS-angle bound
for the signed true-wide deviation.
