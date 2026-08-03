# External snippet decode: threshold/SINR/warm-start/RMT fragment

**Session:** boxeph-multifront (subagent decode task).
**Status:** OPEN (ranked provenance hypotheses; no byte-level match found in public record) + CITED (all anchor claims below carry web sources).
**Companion precedent:** `amm12592-snippet-context-confirmation-deathstar-coinC2.md` (previous fragment, confirmed home = AMM 12592 proposer notes). This fragment has a **disjoint vocabulary fingerprint** and a different home.

## 1. The fragment (13 thinking-summary labels, in given order)

1. Calibrating threshold analysis with negligible quadratic terms.
2. Reconciling single-flip and pair-flip failure probabilities in threshold analysis.
3. Reconciling pair-flip probabilities with empirical observations.
4. Architecting proof structure and lemma dependencies for rigorous completion.
5. Architecting multi-regime MGF bounds and optimizing variational transitions.
6. Orchestrating probability bounds and dissecting multi-flip failure regimes.
7. Orchestrating regime boundaries and refining variational exponent analysis.
8. Architecting SINR bounds and warm-start error analysis rigorously.
9. Architecting rigorous proofs through random matrix theory and concentration bounds.
10. Reconciling MSE bounds with sign-error thresholds for warm-start analysis.
11. Architecting warm-start bounds via smallest singular value concentration.
12. Rigorously bounding small eigenvalue counts for Gaussian matrices.
13. Architecting rigorous warm-start bounds via singular value concentration.

The style ("Calibrating/Reconciling/Architecting/Orchestrating ...") matches machine-generated reasoning-summary headlines of a long LLM proof session, read chronologically: lines 1-3 = threshold pinning via Hamming-1/2 flips, line 4 = proof skeleton, lines 5-7 = union bound over all Hamming distances via Chernoff/MGF with regime-dependent variational exponents, lines 8-13 = a second, separate lemma: a *warm start* whose error is controlled by SINR/MSE/sign-error counts through smallest-singular-value and small-eigenvalue-count bounds **for Gaussian matrices**.

## 2. Vocabulary fingerprint and what it forces

| Term cluster | Forced structural reading |
|---|---|
| single-flip / pair-flip / multi-flip **failure probabilities** | The estimand is a vector `x* in {±1}^n`; "failure" = a Hamming-1/2/d flip of `x*` beats it in the (likelihood) objective; union-bound achievability + single-flip converse. Same grammar as SBM/Z2-sync exact-recovery proofs. |
| threshold analysis, negligible quadratic terms | Sharp (constant-level) threshold; the flip-objective difference has a deterministic quadratic term (`4\|a_k\|^2`-type) that concentrates, plus cross terms between flipped columns (`8 x_i x_j <a_i,a_j>`) shown negligible. |
| multi-regime MGF, variational exponents, regime boundaries | `P(fail at Hamming distance d)` bounded by `inf_λ E exp(-λ Δ_d)`; the optimizer and the exponent shape change across regimes of `d` (linear exponent for small `d`, logarithmic saturation for `d` comparable to `σ^2`), traded against binomial entropy `d log(n/d)`. |
| SINR + warm start + MSE + sign-error thresholds | Stage-1 estimator is a **linear** detector (ZF/decorrelator or LMMSE) of a *linear Gaussian model* `y = A x* + σ z`; per-coordinate post-detection SINR (`SINR_k = 1/(σ^2 [(A^T A)^{-1}]_{kk})` for ZF) gives per-coordinate sign-flip probability `Q(sqrt(SINR_k))`; number of initial sign errors ≤ MSE via Markov over coordinates with error ≥ 1. |
| smallest singular value concentration, **small eigenvalue counts** for Gaussian matrices | Warm-start error `e = σ (A^T A)^{-1} A^T z` requires control of `Tr((A^T A)^{-1}) = Σ 1/λ_i`; for aspect ratio near 1 this is dominated by the Wishart hard edge, needing (i) `σ_min(A)` lower bounds (Rudelson–Vershynin / Davidson–Szarek type) and (ii) a *counting* bound `#{λ_i ≤ t} ≲ n sqrt(t)` summed over dyadic shells. This lemma is only needed when `m/n → 1` (square or near-square design) — a strong structural tell. |
| reconciling with **empirical observations** | The author cross-checks derived flip-failure formulas against simulations — characteristic both of an LLM validating itself numerically and of a literature where near-ML behavior is an *empirical* fact awaiting proof (see H1). |

Immediate exclusions: no graph/Bernoulli vocabulary (excludes SBM proper), no eigenvector-perturbation vocabulary (weakens spectral-initialization models), no sparsity vocabulary (weakens LASSO/1-bit CS/group testing), no dual-certificate vocabulary (weakens SDP-tightness as the *method*, though SDP-tightness literature still frames the open problem).

## 3. Reconstructed theorem being attempted

> **(T)** Let `y = A x* + σ z`, `A` an `m×n` (m = δn, δ ≈ 1; iid Gaussian) matrix, `x* ∈ {±1}^n`, `z` Gaussian. Then
> **(i)** [sharp threshold] the ML detector `argmin_{x∈{±1}^n} ||y − Ax||^2` recovers `x*` exactly w.h.p. iff `m/(2σ^2) ≥ (1+o(1)) log n` (single-flip-dominated threshold, i.e. `σ*^2 = m/((2+o(1)) log n)`); and
> **(ii)** [efficient achievability] a two-stage algorithm — linear warm start (ZF/LMMSE) followed by greedy 1-flip/2-flip local search (LAS-type likelihood ascent) — attains the same sharp threshold in polynomial time.

Line-by-line mapping of the fragment to (T)'s proof:

| Fragment line | Proof step of (T) |
|---|---|
| 1 | Single-flip difference `Δ_k = 4||a_k||^2 − 4σ<a_k, z>`; failure `≈ Q(||a_k||/σ)`; calibrate `n·exp(−m/(2σ^2)) → 0` vs `→ ∞`; cross/quadratic corrections negligible at scale. |
| 2 | Check `C(n,2)·P(pair fails) = o(n·P(single fails))` near threshold (exponent doubles), so the constant is set by single flips; converse from second-moment on single flips. |
| 3 | Monte-Carlo of pair-flip failure frequency vs the derived exponent/prefactor (validation loop). |
| 4 | Skeleton: Lemma A = `x*` is the unique minimizer within (and beyond) a Hamming ball of radius ρn, w.h.p.; Lemma B = warm start lands inside radius ρn with monotone-descent path; Thm = A ∧ B. |
| 5-7 | Lemma A for all `d = |S|`: Chernoff on `P(||A d_S||^2 ≤ 2σ<A d_S, z>)` with `||A d_S||^2 ~ 4d·χ^2_m`; variational parameter optimum shifts across `d ≪ σ^2 / d ≍ σ^2 / d ≫ σ^2`; regime boundaries against entropy `d log(n/d)`. |
| 8 | Lemma B, part 1: per-coordinate SINR of ZF/LMMSE warm start; initial per-bit error `Q(sqrt(SINR_k))`. |
| 9 | Toolbox declaration: nonasymptotic RMT + concentration (Hanson–Wright for quadratic forms, `σ_min` bounds). |
| 10 | Lemma B, part 2: MSE `||e||^2 = σ^2 z^T A (A^T A)^{-2} A^T z` vs the *sign-error threshold* `|e_k| ≥ 1`; #sign errors ≤ MSE ⇒ within basin radius ρn. |
| 11-13 | Lemma B, part 3 (the hard case δ ≈ 1): `Tr((A^T A)^{-1})` via `σ_min` concentration + hard-edge eigenvalue *counting* `N(t) ≲ n·sqrt(t)` summed dyadically. Two "Architecting ... singular value concentration" passes = iterated repair of this lemma. |

Every one of the 13 lines is consumed by this single theorem; no line is left over. That tightness of fit is the main internal evidence.

## 4. Ranked hypotheses

### H1 (best): Large-MIMO detection — rigorous near-ML optimality of flip-based (LAS-type) detection at the sharp exact-recovery threshold. Confidence ~0.40.

**Named open problem.** Since Vardhan–Mohammed–Chockalingam–Rajan (2008), likelihood-ascent-search detectors (MMSE/ZF warm start + 1-flip/2-flip ascent) have been *empirically* observed to approach ML bit-error performance in large square MIMO (e.g. 64×64 V-BLAST, ~1 dB off the single-input matched-filter bound); the companion analysis paper [arXiv:0806.2533](https://arxiv.org/pdf/0806.2533) explicitly leaves the rigorous BER/ML analysis as an **open problem** (its 4-QAM "convergence to ML" argument is heuristic, and the extension is stated as a conjecture). The SDP-relaxation line (Lu–Liu et al., [arXiv:1710.02048](https://arxiv.org/abs/1710.02048), SIAM J. Optim.; Jiang–Liu–Bao–Jiang, [arXiv:2102.04586](https://arxiv.org/pdf/2102.04586)) proves polynomial-time exact recovery only under `σ ≲ sqrt(n/log n)` with the **sharp constant open**. A proof of (T) resolves both: it pins the ML constant and shows cheap local search achieves it — precisely "reconciling ... with empirical observations".

**Evidence.** "SINR" is native, load-bearing vocabulary here (decorrelator/LMMSE SINR, Tse–Hanly); "warm start + flips" is literally the LAS architecture; square Gaussian channel forces the small-eigenvalue-count lemma (line 12); the empirical near-ML literature supplies the observations of line 3. Surveys: [Fifty Years of MIMO Detection](https://arxiv.org/pdf/1507.05138); AMP-side optimality with its uncovered gap regime: [LAMA, Jeon–Ghods–Maleki–Studer](https://arxiv.org/pdf/1510.06095).

### H2: Same theorem in statistics dress — sharp exact-recovery threshold for dense ±1 (binary/Rademacher) linear regression, with efficient achievability. Confidence ~0.22.

The identical mathematical object phrased as "linear regression with ±1 coefficients": adjacent to the Gamarnik–Zadik sparse-binary-regression program (thresholds + local search + OGP; [arXiv:1711.04952](https://arxiv.org/pdf/1711.04952), [Zadik thesis](https://iliaszadik.github.io/files/PhDThesisZadikFinal.pdf)) and to the binary-compressed-sensing "relaxation gap" (region where the signal is determined but convex methods provably fail; [arXiv:2606.00806](https://arxiv.org/html/2606.00806), June 2026, with *empirical* quantum-annealing recovery inside the gap). A dense-±1 sharp-threshold + local-search theorem would be a significant entry in this program. Distinguishing H2 from H1 is mostly sociological (venue/notation); "SINR" is the one term that leans H1 — a statistician would write "effective SNR".

### H3: Z2/phase synchronization — sharpness of the `sqrt(n/log n)` tightness threshold for SDP/GPM. Confidence ~0.12.

The threshold scale and the flip grammar match (Z2-sync exact recovery at `σ* = sqrt(n/(2 log n))`, single-flip-dominated; spectral methods optimal: [arXiv:2209.04962](https://arxiv.org/pdf/2209.04962)); the sharpness of `σ ≲ sqrt(n/log n)` for SDP/Burer–Monteiro tightness is explicitly open ([arXiv:2510.01522](https://arxiv.org/pdf/2510.01522)). **Against:** synchronization proofs need eigenvector perturbation, not `σ_min` of an inverted design, no SINR, no least-squares warm start; lines 8, 11-13 have no natural home. Ranked third only because the problem is famous and the flip-side (lines 1-7) fits well.

### H4: Two-component Gaussian mixture / label recovery sharp threshold (Ndaoud-style iterative sign updates). Confidence ~0.08.

Exact recovery of `η ∈ {±1}^n` from `X = η μ^T + E` with sharp threshold and an iterative flip algorithm ([arXiv:1812.08078](https://arxiv.org/pdf/1812.08078)). Flips + threshold + warm start fit; but the sharp-threshold result is already *solved*, no SINR, no design-matrix inversion ⇒ no small-eigenvalue counting. Plausible only as an extension (e.g. growing dimension regimes).

### H5: Box-relaxation / regularized-LS BPSK error analysis (CGMT, Thrampoulidis–Stojnic school). Confidence ~0.05.

SINR-style effective quantities and BPSK sign errors are native ([Thrampoulidis–Abbasi–Xu–Hassibi ICASSP 2016](https://ieeexplore.ieee.org/document/7472383/), [regularized LS BER, ICASSP 2017](https://dl.acm.org/doi/10.1109/ICASSP.2017.7952960), [Stojnic fl-RDT ML line, arXiv:2410.07651](https://arxiv.org/pdf/2410.07651)); but CGMT/RDT proofs contain no flip union bounds and no warm start — lines 1-7 would be homeless. More likely these works are *inputs* to H1 (the warm-start MSE/SINR sub-lemma) than the home.

### H6: Perceptron / discrepancy / 1-bit CS / matrix completion local algorithms. Confidence ~0.03 combined.

Each misses ≥ 3 fingerprint clusters (no noise ⇒ no SINR/MSE for perceptron and discrepancy; no ±1 flip structure for 1-bit CS target vector; no Gaussian-inverse RMT for matrix completion). Retained only as residual mass.

Residual "other/unknown" ~0.10 (e.g. an unpublished problem set posed privately to the LLM; note the previous fragment's home was proposer-side notes, i.e. private material is a live possibility for this corpus).

## 5. Decisive tests (what would confirm, as done for the AMM 12592 fragment)

1. **Model-equation test.** Any further context containing `y = Hx + n` / `y = Ax + σz` with `x ∈ {±1}^n` (or `S = {±1}` symbol set), Gaussian channel/design ⇒ confirms H1/H2 family and kills H3/H4. Complex field + QAM/BPSK words ⇒ H1; real field + "design matrix"/"regression" ⇒ H2.
2. **SINR-formula test.** Appearance of `1/(σ^2 [(H^*H)^{-1}]_{kk})` (ZF) or the LMMSE analogue ⇒ H1 essentially confirmed (this exact expression is the decorrelator/MMSE per-stream SINR of multiuser-detection literature).
3. **Threshold-constant test.** Extracted sharp constant of form `σ*^2 = m/((2+o(1)) log n)` (equivalently `SNR_min ~ (2/m) log n`; for m=n: `σ* = sqrt(n/(2 log n))`) ⇒ H1/H2. Constant of form `λ ≥ (1+o(1)) sqrt(n/(2 log n))` attached to a GOE/Hermitian observation matrix ⇒ H3. Separation-type constant `Δ^2/σ^2 ≍ (1+sqrt(1+c p/(n log n))) log n` ⇒ H4.
4. **Hard-edge-counting test.** A stated lemma of shape `#{λ_i(A^T A) ≤ t} ≲ n sqrt(t)` (or dyadic-shell summation of `Σ 1/λ_i`) ⇒ the design is square/near-square ⇒ H1 (large square MIMO) strongly over H2-with-δ-large; this lemma is pointless for δ bounded away from 1.
5. **Empirics test.** If the fragment's reconciled numbers are BER-vs-SNR curves matching MMSE-LAS simulation folklore (~1 dB off MFB at 64×64, 4-QAM) ⇒ H1 confirmed at the provenance level (the "empirical observations" would be the published LAS curves, not the author's own).
6. **Negative test.** Any appearance of edge probabilities `p = a log n / n` (SBM), "communities", or eigenvector `ℓ_∞` perturbation ⇒ demote H1/H2, promote H3/SBM-adjacent readings.

## 6. Contrast with the previous fragment (provenance bookkeeping)

The AMM 12592 fragment was elementary-analysis flavored (artanh sandwiches, exact rationals, coin-flip deadlines) and resolved to *proposer-side* notes. The present fragment is high-dimensional-probability flavored (MGF regimes, Wishart hard edge, SINR) with an *achievability-proof* arc — a different author profile (an LLM run at a rigorous-proof task, plausibly against an open-problem benchmark or a researcher's private conjecture). No public 2025-2026 arXiv item was found that already states (T) with this two-lemma architecture; if the fragment's session succeeded, the result would be new. The closest public statements of the open problem are [arXiv:0806.2533](https://arxiv.org/pdf/0806.2533) (rigorous LAS/ML analysis open) and the sharp-constant gaps in [arXiv:1710.02048](https://arxiv.org/abs/1710.02048)/[arXiv:2102.04586](https://arxiv.org/pdf/2102.04586)/[arXiv:2510.01522](https://arxiv.org/pdf/2510.01522).

## 7. Sources

- LAS asymptotic analysis, open-problem statement: https://arxiv.org/pdf/0806.2533
- LAS near-ML empirics (2008 ICC): https://ieeexplore.ieee.org/document/4595342/
- Fifty Years of MIMO Detection (survey): https://arxiv.org/pdf/1507.05138
- LAMA / AMP optimality for large MIMO: https://arxiv.org/pdf/1510.06095
- Enhanced SDR tightness for MIMO detection: https://arxiv.org/abs/1710.02048
- Tightness and equivalence of SDRs for MIMO detection: https://arxiv.org/pdf/2102.04586
- SDP/Burer–Monteiro tightness for phase synchronization, high-noise regime (sharpness open): https://arxiv.org/pdf/2510.01522
- Exact minimax optimality of spectral methods in phase/orthogonal synchronization: https://arxiv.org/pdf/2209.04962
- Sharp optimal recovery in two-component Gaussian mixture: https://arxiv.org/pdf/1812.08078
- Sparse binary regression thresholds + local search (Gamarnik–Zadik program): https://arxiv.org/pdf/1711.04952 , https://iliaszadik.github.io/files/PhDThesisZadikFinal.pdf
- Binary compressed sensing relaxation gap (empirical recovery inside the gap): https://arxiv.org/html/2606.00806
- Box relaxation BER for BPSK (CGMT): https://ieeexplore.ieee.org/document/7472383/
- Regularized LS BER for BPSK: https://dl.acm.org/doi/10.1109/ICASSP.2017.7952960
- Stojnic fl-RDT descending-ℓ0 ML line: https://arxiv.org/pdf/2410.07651
- Smallest singular value of rectangular random matrices (Rudelson–Vershynin exposition): https://pages.cs.wisc.edu/~brecht/cs838docs/rv_rectangular.pdf
- High-dimensional statistics open-problems survey (checked, no verbatim match): https://arxiv.org/pdf/2605.05076
