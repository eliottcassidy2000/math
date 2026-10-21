# Kepler boundary audit and current concordance

**Status: independently audited scoped results, with concurrent correction preserved, 2026-10-21.** No new Lonely Runner case is proved. The [frozen audit](AUDIT.md) concerns the [original inspected THM-4737](original/THM-4737.md), not the subsequently corrected live theorem. Its statements that the old K2 criterion is refuted and the old K3 finiteness argument is circular must be read against that preserved source.

## Current concordance

Concurrent [commit cfcb930](https://github.com/eliottcassidy2000/math/commit/cfcb930344470c60cf13ffc8694e990eaa074209) corrected the same endpoint and existence issues before this audit packet was published. Its [THM-4737](../../../01-canon/theorems/THM-4737-kepler-polygon-reformulation-of-the-lonely-runner-conjecture-observer-angular-momentum-vertex-times-and-the-relation-torus.md) and [boundary proof](../runner_area_clock_20261020j/boundary_and_colouring.md) remain intact, with their original attribution. This packet is add-only: it neither overwrites that correction nor adds a duplicate MISTAKES entry.

The two treatments agree on the substantive repairs:

- For positive speeds and n>=3, closed loneliness means `max(A_i^-(t),A_i^+(t))>0` separately for every moving runner. Far endpoints of incident sides are allowed; a collision at the observer is not. For D={1,2} at t=1/3 the runners require opposite one-sided germs.
- For integer speeds, define the least lonely rational denominator N_min as infinity when no witness exists. It is finite exactly when the row reaches its LRC threshold. Rational attainment of the maximum alone does not prove this.
- The real-speed closed-orbit-closure equivalence is a supremum statement, not by itself proof of actual boundary attainment. The rank-at-most-one actual-time conclusion has a separate valid argument: interior hits for n>=4 and a closed commensurable orbit in the two-moving-runner case.

The frozen audit also proves the conditional bound `N_min<=2 max(D)` once the integer row reaches its threshold. It makes the density convention explicit: for noncommensurable real speeds, Haar measure denotes asymptotic time density, not the measure in an arbitrary unit window. The old broad language about a standard real-to-integer reduction is not certified merely by proving the closed-torus equivalence; this packet does not supply a general attainment/reduction proof. None of these cautions is a claimed real-speed LRC counterexample.

The concurrent packet additionally contains its own attributed distance-graph colouring counterexample; this scoped boundary audit does not claim authorship of that result or broadly re-audit all claims in the live theorem.

## Exact lineage and reproduction

The corrected THM-4737 as read at cfcb930 has SHA-256 `859453551375340557d52de9f9cfa19c92722f0fa2c650832f7f59ef76b93f8a`. The original inspected source has SHA-256 `0ea008a5ab25cfde79fc11b5cce785104af2205e047022a613e4de85d89cf769`; it is retained only as superseded lineage. `inspected-source.sha256` names the pre-correction repository path and matches `original/THM-4737.md`, not the live theorem.

From this directory run `python3 -B check_kepler_boundaries.py` and `python3 -B -O check_kepler_boundaries.py`; both outputs must equal `check_kepler_boundaries.json`. The original `manifest.json` is preserved unchanged and pins the frozen audit's four deliverable files. The bounded rational checks corroborate the displayed proofs; they are not a general LRC search or theorem proof by exhaustion.
