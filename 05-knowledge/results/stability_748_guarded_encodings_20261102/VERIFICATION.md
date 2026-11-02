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
