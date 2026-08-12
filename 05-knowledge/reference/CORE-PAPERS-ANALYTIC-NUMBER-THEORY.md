# Analytic number theory source sidecar

## Claude — *More than two thirds of the zeros of the Riemann zeta function lie on the critical line*

- **Primary / freshness:** [author PDF, 2026-08-10](https://www-cdn.anthropic.com/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf),
  35 pages. **EXTERNAL PREPRINT CLAIM + FORMAL ARTIFACT; PARTIAL LOCAL
  AUDIT.** The companion [Lean repository](https://github.com/anthropics/zeta-23-lean)
  publishes a tagged, sorry-free formalization and audit instructions.
- **Claimed results:** lower density `2/3` for simple critical-line zeros,
  `5/6` for distinct zeros, optimized constants `0.6725` and `0.83625`, and
  an analogous primitive Dirichlet `L`-function theorem. Treat these as the
  paper's claims pending ordinary external validation; they are not repo
  canon.
- **Imported mechanism:** compress Weil's Hermitian form to a finite Gabor
  family. On-line zeros contribute positive rank-one fixed blocks, while an
  off-line functional-equation pair contributes a signature-`(1,1)` block.
  Prime-side trace and Frobenius-square asymptotics, a Hermitian rank--trace
  inequality, and Sylvester inertia then force positive fixed blocks. The
  reusable lesson is to retain the indefinite block-signature sidecar through
  compression; trace and second moment alone lose the distinction.
- **Audit note:** [upstream issue 7](https://github.com/anthropics/zeta-23-lean/issues/7)
  identifies a scale typo in the equality sentence of Lemma 3.2: the equality
  model should read `P=(c/2)Pi_1`. The inequality and downstream applications
  are unaffected. We inspected the paper, rendered its proof core, and audited
  the formal source tree. A local Lean build was not completed because the
  toolchain download reset/timed out; do not report local kernel reproduction.
- **Repo consumers:** [THM-3337](../../01-canon/theorems/THM-3337-cross-shell-compression-attains-the-T4-floor.md)
  and [THM-3338](../../01-canon/theorems/THM-3338-horizon16-cross-shell-surgery-closes-the-first-fifteen-pointwise-AMM-values.md)
  use the finite-compression lesson motivationally; [THM-3340](../../01-canon/theorems/THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md)
  is the resulting AMM continuation with an independent cyclic-rotation proof.
  No analytic statement from the paper is a dependency.
- **Does not prove:** any AMM 12592, LRC(14), FC, GMC, or Jacobian claim. The
  proposed `C_13` owner-phase transfer to LRC remains conditional on native
  address, ancestry, outer co-support, and tail-norm sidecars; raw Gram rank is
  provably blind to the AMM integer-capacity obstruction.
