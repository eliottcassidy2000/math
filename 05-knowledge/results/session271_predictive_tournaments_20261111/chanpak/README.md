# Predict large Hamiltonian counts through occupied-slot phases

The main advance is a **proved, finite-support classification of exponential
growth** for the fixed-core/transitive-spine carrier. It replaces an uncalibrated
random-matrix forecast with an explicit integer invariant r(C,S).

* For iid colors supported on S, almost surely log H/m tends to log r.
* Every positive moment has pressure s log r; Var(log H) stays bounded.
* Full support gives r=k+1 for every k-core, regardless of its orientations.
* All 15 nonempty two-core supports fall into exactly three phases: bases1,2,3.
* With bit marginals alpha_j in (0,1), the leading mean coefficient is
  k! product alpha_j(1-alpha_j), even for correlated bits in each letter.

These are fixed-core carrier theorems, not a classification of all recursive
tournaments. The cone is the positive occupied-slot state space; the17-dimensional
rational quotient remains useful for exact algebra but hides this geometry.

For a two-core target of scale10^100, the language-specific order lower bound
is211 and the ensemble mean-scale proxy is214. The latter is **not** an existence
or witness-size upper bound. The scaling table extends to10^(1,000,000) without
trying to construct or enumerate those objects.

## Honest comparison of target heuristics

All methods get1,200 exact transfer evaluations per target, at a fixed predicted
length. Relative errors below are distance from the target, not exact solutions.

| Target | iid | local absolute-distance mutation | residue tie-break |
|---|---:|---:|---:|
| 10^12+1 | .003699 | .0007285 | .00007017 |
| 10^30+1 | .001408 | .000002051 | .000006486 |
| 10^100+1 | .0005364 | .00003670 | .01509 |

There were **zero exact hits**. Residue tie-breaking helps one case and damages
another. Different random streams and one seed do not establish a statistical
ranking. The productive lesson is that scale prediction is easy while exact
integer inversion remains difficult. A residue sieve should be integrated
with compatible prefix/suffix states, not used as a superficial scalar score.

## Files and reproduction

* [PROOF.md](PROOF.md): four proved statements, cone/tropical interpretation,
  negative-moment variational bounds, q-signal decay, symmetry and SCC sidecars,
  exact-residue/returning-cone target mechanisms, and explicit limits.
* [TABLES.md](TABLES.md), [predictions.json](predictions.json): bounded growth
  probes, all support phases, random residue observations, huge-target scales.
* `predict.py`: seed27120261111, five iid laws, 1,200 words per law through
  length512; plus10,000 small random support probes. Requires NumPy.
* `target_probe.py`, `target_probes.json`: exact integer heuristic scores and
  best words for nine runs. Seed27191. No claim of optimality or absence.
* `controls.py`, `controls.json`: new sandwich controls on525 support cases
  for three cores; not an audit of the inherited transfer or prior atlas.

From the repository root:

    OPENBLAS_NUM_THREADS=1 python3 05-knowledge/results/session271_predictive_tournaments_20261111/chanpak/predict.py
    python3 05-knowledge/results/session271_predictive_tournaments_20261111/chanpak/target_probe.py
    python3 05-knowledge/results/session271_predictive_tournaments_20261111/chanpak/controls.py

The programs import the inherited Session270 positive transfer. They do not
regenerate its million-target atlas. Floating normalized matrices are used
only for distribution statistics; target scores and sandwich checks use integers.
The displayed integer order scales use60-digit Decimal logarithms and are far
from rounding thresholds, not interval-certified logarithm computations.

Chan–Pak remains open. Large class entropy, broad residues, and positive mass
cannot establish full scalar coverage. The next concrete prize is an integral
returning cone with residue-covering inverse macros and a decreasing height.
