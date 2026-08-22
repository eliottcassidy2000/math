# Analytic number theory source sidecar

## Integer power-series rigidity used by THM-3342

- **Polya--Carlson:** F. Carlson, *Uber Potenzreihen mit ganzzahligen
  Koeffizienten*, Math. Z. 9 (1921), 1--13. An integer-coefficient power
  series of radius one is rational or has the unit circle as a natural
  boundary. THM-3342 uses continuation through an open boundary arc to force
  the rational branch.
- **Fatou--Kronecker pole step:** P. Fatou, *Series trigonometriques et series
  de Taylor*, Acta Math. 30 (1906), 335--400; L. Kronecker, *Zwei Satze uber
  Gleichungen mit ganzzahligen Coefficienten*, J. reine angew. Math. 53
  (1857), 173--175. Applied to the rational integer series with
  subexponential coefficients, these put every pole at a root of unity.
- **Finite-alphabet endpoint alternative:** G. Szego, *Uber Potenzreihen mit
  endlich vielen verschiedenen Koeffizienten* (1922). Bell--Chen,
  [*Power Series with Coefficients from a Finite Set*](https://arxiv.org/abs/1606.04986),
  gives a modern route to the classical finite-alphabet theorem. Szego alone
  handles bounded additive slack; THM-3342 needs Polya--Carlson for the full
  `o(n)` range.
- **Scope:** these classical theorems are cited dependencies, not formalized
  locally. The repo's exact companion checks the algebraic reductions and
  hostile controls, not the analytic dichotomies.

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
  proposed `C_13` owner-phase transfer to LRC fails its first native hostile:
  the support involution has inertia `(91,78,338)`, but exhaustive tests in
  both certified fields find no finite-Gabor amplitude covariance and the
  outer semantic co-support is empty.  Raw Gram rank is likewise blind to the
  AMM integer-capacity obstruction.
