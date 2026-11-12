# Tournament compatibility-fiber lane

[PROOF.md](PROOF.md) contains the mathematics and provenance.

**Inherited:** quadratic attachment law, seven-coordinate shared-root fusion,
full-central criterion, C5's41/43 gate, Paley unordered-pair symmetry.

**New scoped results/application:** exact five-coordinate fixed-deletion affine
span; transitive-spine compatibility filter and first indefinite curvature;
vertex-transitive response fibers as an orbital energy problem; actual
order10 central14207/14209 pair with deletion3357; controlled cyclic-family
experiments and separately labelled conjectures.

No tournament-stress isomorphism or Chan–Pak completion is claimed.

## Reproduce

From the repository root (Python3 standard library only):

```sh
python3 05-knowledge/results/session285_five_web_20261112/tournaments/verify.py
python3 -O 05-knowledge/results/session285_five_web_20261112/tournaments/verify.py
python3 05-knowledge/results/session285_five_web_20261112/tournaments/scan_orbits.py
python3 -O 05-knowledge/results/session285_five_web_20261112/tournaments/scan_orbits.py
```

The first pair emits `replay.json`; the second emits `orbit_scan.json`.
Scripts write to stdout only. All guards use explicit exceptions and survive
optimization. `explore.py` supplies the explicit path-coefficient compiler;
its optional main block records the early constructive discovery search.
The verifier separately computes full pointed path data, using literal
permutations for every labelled core through order4 and an independent
subset recurrence for larger controls. The separate audit lane independently
rechecks the main examples and new pair.

The scan covers all2^((n−1)/2) cyclic step sets for n=9,11,13 and all fixed
cardinalities2..(n−1)/2, plus the specified carousel cores5..17. It does not
claim an isomorphism-class census or all vertex-transitive cores.
