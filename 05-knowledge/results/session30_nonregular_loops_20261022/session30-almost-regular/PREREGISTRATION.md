# Session30 independent almost-regular experiment

Frozen before computing outcomes, 2026-10-22 15:05 UTC.

Definitions: M=S/n, delta=M1, kappa=2^(n-1)H/n!, rho=exp(1/2 logdet(I-M²)-1/2 logdet(I+M²)). Evaluate audited regular two-loop formula and its unaudited path-tree extension y0+t, t=n delta^T(I+M²)^(-1)delta. No coefficients fitted. These are finite predictions, not an asserted asymptotic theorem.

Controls: all labelled almost-regular tournaments n=2,4,6; all permutations verify H at n<=6 against subset DP on representatives. Deterministic families n=8,10,12,14,16,18,20: delete vertex zero from rotational regular n+1; delete zero from a fixed pseudorandom circulant (seed 30001+n); delete zero from directed triangle-flip chain starting rotational (seed 30002+n, 1000 accepted flips); also vertex deletions from Paley prime orders 7,11,19,23 (largest n22 if feasible). Nonregular hostile controls transitive n4,6,8,10,12, and regular controls rotational n5,7,9,11. Chain outputs are not asserted uniformly random.

Record exact adjacency bits, exact H, all predictor ingredients, y, residual, n*residual, and n²*residual. Independently implement formulas without importing producer code. If a revised formula arrives after outcomes, label it prospective only on newly generated cases, or retrospective on frozen cases.

15:10 UTC supplement, before held-out lollipop outcomes: theory lane supplies L=−2diag(H−G)^T G1, y_new=y_path−nL. All original cases evaluated retrospectively. Fresh held-out family: first regular7 chain sample (seed30500, batches10flips) with H171; lexicographic product with C3 on21vertices, delete vertices0,3,6,9,12,15,18. No coefficients fitted; finite first-order diagnostic only. Nonvertextransitive base is chosen structurally by known H171, not predictor outcome.
