# A native odd-depth boundary for target coverage

Session 39, 23 October 2026. **PROVED elementary counting bound.** This is a separate extension; earlier frozen packets are unchanged. No novelty or convergence claim.

For either signed accelerated map Uε, a positive odd target y, real X≥1 and integer r≥0, let Hε,y(X,r) be the number of positive odd starts n≤X that reach y within at most r accelerated odd steps. Count a start once even if it revisits y.

**Proposition.** With empty product equal to 1,

Hε,y(X,r) ≤ sum_(j=0)^r product_(i=0)^(j−1) floor(log_2(3·2^i X+1))
             ≤ sum_(j=0)^r (floor(log_2 X)+r+3)^j.                 (1)

The bound is uniform in y. For a finite target set Y, multiply the right side by |Y|; overlapping paths only improve the upper bound.

**Proof.** At a positive odd integer x, both maps obey

Uε(x) ≤ (3x+1)/2 ≤ 2x.

Thus along any candidate start n≤X, the i-th odd state xi is at most 2^i X. Its valuation ai satisfies 2^ai≤3xi+ε≤3·2^i X+1, giving the individual letter bound in (1). A word of length j and a fixed endpoint y determine at most one start, because its affine branch has nonzero slope 3^j/2^A; explicitly n=(2^A y−εB)/3^j. Multiplying the possible letter counts and summing over j proves the first inequality; it only overcounts inadmissible words and multiple visits. The second is an intentionally loose uniform bound on every letter count. The j=0 term covers the possible start y itself. ∎

**Asymptotics.** For fixed r, Hε,y(X,r)=O_r((log X)^r), hence has natural density zero among positive odd starts or all starts. More generally, if r=r(X)=o(log X/log log X), put L=floor(log_2 X). From (1),

Hε,y(X,r) ≤ (r+1)(L+r+3)^r = exp(o(log X))=X^o(1)=o(X).

Indeed r=o(log X/log log X) implies r=o(log X), log(L+r+3)=O(log log X), and log(r+1)=o(log X). The same conclusion holds for a fixed finite Y. A fixed positive proportion of starts cannot therefore be covered at such a uniformly small actual odd-depth horizon. This is only a necessary growth boundary, not a sufficient depth estimate or a typical stopping-time law.

**Why this matters, and what it does not say.** The closing-step construction in CLOSED_TARGET_CYLINDERS.md puts infinitely many points of each known minus basin inside every finite valuation cylinder, but topological ubiquity is not positive density at bounded arrival depth. Exact arrival is different from a local modular resemblance. Reaching any of the finitely many odd vertices in the known cycles is covered by the finite-Y version.

The clock is actual accelerated odd iteration depth. This proposition does NOT bound proof-search running time, symbolic certificate size, shortcut-step depth, or total arithmetic cost. It does not refute a linear 153N annulus budget using unbounded path depth, nor does it prove any such budget. It also does not control reaching a moving target set whose size grows with X without accounting for that growth.

**Repository comparison.** Read-only searches at snapshot e5d18e55 located `05-knowledge/results/collatz_depth_layers_zeta5_preprints_20260927.md`, whose §1 already retains the distinction between odd-step depth and division/shortcut depth, and `01-canon/theorems/THM-4517-no-root-uniform-positive-density-of-collatz-basins.md`, whose disjoint-root density obstruction varies the target rather than bounding native arrival depth to a fixed target. Neither result is needed as a premise. The present note is an elementary scoped boundary consolidating those clock distinctions, not a claim that an exhaustive search established originality.

**Controls.** `verify_coverage_clock.py` independently enumerates starts for X=1,2,16,64,256,1024, depths 0 through 4, both signs, and targets 1,5,17: 180 coverage cases and 40,920 step-height/valuation controls. It also checks that no two starts with the same endpoint and valuation word disagree. Normal and optimized Python match `coverage-clock-receipt.json` and `coverage-clock-receipt.optimized.json`. These are intentionally small validation cases; the displayed bound is loose and can exceed X at these sizes. Its density conclusion is asymptotic, proved above, not read from the finite table.
