# Dynamics and labeling: entry, finite memory, and bounded release

Research checkpoint, 19 October 2026. These are proofs and finite-exact controls, not resolutions of the Collatz or Graceful Tree conjectures. No historical novelty is claimed for classical coding, matching-switch, or parity-density mechanisms.

## Read in this order

1. BERGGREN_POSITIVITY_NONREGULAR.md. The positivity language for Pythagorean AP square-sum triangles is nonregular despite the universal four-step cone entry rule. More generally, finite-jump predicates have regular word languages exactly when their genuine slope boundaries are rational or quadratic. The necessity proof is elementary; the quadratic sufficiency invokes the attributed Romik Lagrange theorem.
2. GRACEFUL_RELEASE_AND_MATCHING_SUPPORT.md. An exact component-reconnection test identifies the rank gained by releasing a fixed edge. A one-edge repair criterion gives a connected completion fiber while retaining every other fixed edge. Tree-supported fibers reduce exactly to perfect matchings and alternating-cycle switches. The prior 13-label obstruction admits eight successful one-edge releases and one ineffective forced-edge release.
3. COLLATZ_DENSITY_ENTRY_AND_POINTWISE_GAP.md. The classical parity-word bijection gives density-one eventual entry into the explicit all-lift descent cylinder KF, with an effective power-saving bound for the exceptional set. Pointwise entry, infinite recurrence, and descent below the original source are explicitly not established.

## Verification

Run Python 3 on verify_dynamics_session2.py. It prints the included verification-receipt.json. Running Python with -O gives byte-identical output: all checks use explicit exceptions, not assertions.

The replay passes 191,138 checks. Controls include 64 finite-word residuals distinguished in all 2,016 pairs; 29,524 ternary Pythagorean nodes; the complete fixed-edge release census; the shortest three-move repair; the finite parity bijection through length 14; 256 lawful KF lifts; 513 independently counted pattern-avoidance lengths; and 256 increasing all-odd runs.

The proofs establish the unbounded statements. The finite controls check composition conventions, exact formulas, and concrete witnesses.

## Files

- Three proof notes listed above
- verify_dynamics_session2.py
- verification-receipt.json
- MANIFEST.json and SHA256SUMS

The exploratory release search is intentionally omitted. This package contains no Git changes or publication claims; it is ready to be copied to the repository-relative directory shown by its path.
