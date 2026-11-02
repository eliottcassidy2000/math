# Verification and review

Status: PROVED for the written stability theorem, supported by exact arithmetic and separate combinatorial checks. No independent agent was delegated or claimed.

Reproduce from repository root:

```sh
python3 05-knowledge/results/stability_748_guarded_encodings_20261102/verify_748.py
python3 05-knowledge/results/stability_748_guarded_encodings_20261102/audit_groups.py
```

The first verifier checks all rational endpoint constants, square completion, the exact exponential endpoint,1253 finite corroborating orders, and the ideal two-bin ceiling. Monotonicity in the written proof extends the endpoints to all n>=748; the finite sweep is not the theorem.

The second implementation builds greedy groups and counts incident edge pairs directly. It exhausts729 ternary-weighted K4 graphs at two cuts, plus800 seeded rational controls at orders5..12. It checks disjoint group energy, mass conservation and both Cauchy bounds. Its finite samples supplement the all-size algebraic proof.

Proof review checked: both rounding inequalities have the SAME independent terminal distribution; the retained heavy matching is nonempty; its denominator is safely replaced by R-L only with positive retained mass; the retained matching has an unmatched vertex; strict inequalities handle equality candidates; the heavy/light wedge pairs are disjoint; the last scalar minimum lies in[0,1]; the inherited H=1 equality classification applies only after B=0. No termwise Pfaffian upper bound is assumed.

The inherited low-order hostile C4 disjoint-union K3 still obstructs a generic unsigned graph shortcut. The new proof does not supply that shortcut or an all-order Hamiltonian spectrum theorem.

## Current n>=514 continuation

```sh
python3 05-knowledge/results/stability_748_guarded_encodings_20261102/verify_514.py
python3 -O 05-knowledge/results/stability_748_guarded_encodings_20261102/verify_514.py
python3 05-knowledge/results/stability_748_guarded_encodings_20261102/audit_dual.py
python3 -O 05-knowledge/results/stability_748_guarded_encodings_20261102/audit_dual.py
```

The exact endpoint verifier checks the bootstrap, heavy count17, both convex tail endpoints, the paid-small-m tail, and final strict surplus187/41120000. It additionally checks141556 unsplit rational(R,m) controls. Normal and optimized output must agree. Removing U fails the declared hostile; the same constants fail at512.

The separate dual checker minimizes the quadratic directly in S, retaining the nonnegative residual w(t-w)/4. It checks5151 rational points, all729 three-weight K4 graphs and800 rational graphs at orders5..12 by literal incident-edge accounting. This is an independent implementation, not a delegated independent mathematical reviewer.

The all-n proof depends on monotonicity and convexity, not the finite sample. Audit review checked that U enters both the retained matching mass and the SAME local-energy budget; in the m<=31 branch its lower bound comes from M_h<=m/4, while the m>=32 branch pays the tail without U. The lower bound may be negative outside its branch and is not used there. All inequalities use the same B; D=0 is included in the strict contradiction. No order below514 is promoted.
