---
id: HYP-3127
title: The uniform multi-far floor R′≥c is an ASANO contraction of single-far Lee-Yang factors — the coverage is multilinear (inclusion-exclusion) so Asano applies; each far element is a "tip" whose single-far zero-free region is HYP-2829; contracting the tips against the bounded-core "tail" preserves the floor; Gaussian = the decoupled free-field limit, Elliott-Halberstam = the level-of-distribution bounding SPEC, tip/tail recursion (HYP-3124) = the contraction order
status: VERIFIED structural evidence (multilinearity; R′_cov floor ≈0.98 distinct, reproduces R′∈[0.81,1.0]; Gaussian decoupling) + REDUCTION (multi-far floor ⟸ single-far Lee-Yang region + Asano + EH SPEC bound). Not a proof.
source: mac-mini-2026-06-27-S68
extends:
  - HYP-3122   # the cap is a φ⁴ field theory (S67) — single-tip is the φ⁴ measure
  - HYP-3124   # codex edge-witness tip/tail recursion = the Asano contraction order
related:
  - HYP-2829   # single-far closure = the tip factor's zero-free region (the load-bearing lemma)
  - HYP-2900   # Node-3 equidistribution = the tip peel
  - HYP-2692   # the two-far residual R_2
  - HYP-2797   # the genuine-wide doublet (the correlated r=2 hard case)
  - THM-563    # single-far periodicity
  - OPEN-Q-108
external: Asano contraction (Lee-Yang); Lieb-Sokal; Elliott-Halberstam; (φ⁴)₂ Euclidean QFT
---

# HYP-3127 — The multi-far floor is an Asano contraction of single-far factors

## The reduction
The covering bound's open core = the **multi-far floor** `R′ = 1 + SPEC ≥ c` over `r=2..6` far placements
(`R′ = meas(GOOD∩G_P)/(meas(GOOD)·meas(G_P))`, `SPEC = Σ_{n≠0} ĉ(n)ĝ(n)*`). Bold reduction (first use of
Asano/EH/Gaussian in the LRC context):
- **Coverage is MULTILINEAR in the runners** (`p0 = E_x[∏_e (1 − 1_{miss,e})]`, degree 1 per runner) = exactly
  the hypothesis of **Asano's contraction lemma** (multi-affine + zero-free polydisk ⟹ the contracted polynomial
  is zero-free in the disk).
- So **each far element is a "tip"** (a single-variable factor), its **single-far Lee-Yang region IS HYP-2829**,
  and **Asano-contracting the `r` tips against the bounded-core "tail" preserves the zero-free region ⟹ R′≥c**.
- This is the multi-variable extension of S67: the single tip is the φ⁴ measure; **Asano is the Lee-Yang-
  preserving coupling of φ⁴ measures** (Lieb–Sokal).

## Verified evidence (`lrc_multifar_asano_floor_macmini_S68.py`)
`R′_cov(F) = p0(B∪F)·p0(B)^{r−1}/∏_f p0(B∪{f})` (=1 = Asano-factorized):
- **`R′_cov ∈ [0.87, 1.05]`, floor `≈0.98` for distinct far** (0.87 only at a degenerate repeated-speed config)
  — reproduces the kps/codex `R′∈[0.81,1.0]`. The multi-far coverage **quasi-factorizes over the tips** with a
  **positive floor** = the Asano floor.
- **Single-far factor `d(f)=p0(B∪{f})/p0(B)≈1.10` STABILIZES for large f** = the **Gaussian / free-field decoupled
  limit** (`R′→1`). Tight doublets `R′_cov>1` (the HYP-2797 correlated hard case); separated `R′_cov<1`.

## The three tools, placed
- **Asano** = the engine (multi-affine coverage ⟹ contraction preserves the zero-free region).
- **Gaussian** = the `λ→0` free-field limit of the φ⁴ single-tip (large/separated tips decouple, `R′→1`,
  `SPEC=0` baseline).
- **Elliott–Halberstam** = the level of distribution of the far tips; `SPEC` small when the far cluster
  equidistributes; `EH (θ→1) ⟹ R′→1`; unconditionally a Bombieri–Vinogradov-level input + signed cancellation
  bounds `SPEC`.
- **tip/tail recursion (HYP-3124)** = the Asano contraction order: peel one tip (Node-3 pull-back), contract,
  recurse on the tail (push-forward); the edge-witness compatibility = the Asano zero-free preservation.

## Proof obligations (what remains, named)
1. **Single-far Lee-Yang region** (Asano applicability): state HYP-2829 as a zero-free polydisk for the single-far
   factor — *the load-bearing lemma*.
2. **The SPEC bound / the constant `c`:** `|SPEC| ≤ 1−c` via EH-level far-equidistribution + signed cancellation
   (verified `c≈0.98` distinct).
3. **Recursion termination + monotonicity:** `r↑ ⟹ R′↑` (contraction improves), bottoming at `r=1` (closed) + the
   bounded core (finite check).

## Next
1. Compute the single-far factor's zeros in the far variable (test the Lee-Yang polydisk = obligation 1).
2. Bound `SPEC` for the worst (resonant) `r=2..6` placements via the φ⁴/Gaussian spectral decay + signed sum.
3. Verify the `r↑ ⟹ R′↑` monotonicity rigorously (the decorrelation-improves-with-contraction step).
