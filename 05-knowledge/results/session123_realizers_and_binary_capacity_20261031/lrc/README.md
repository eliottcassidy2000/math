# Session123 LRC: witness extraction without losing the clock

* [Constructive interval conductor](CONSTRUCTIVE_CONDUCTOR.md): an exact
  polynomial-bit common-time selector, including the equality boundary,
  using rational measure oracles and at most one overlap repair.
* [Same-clock categories](SAME_CLOCK_CATEGORIES.md): a precise closed symmetric
  monoidal predicate category, the failure of existential truth to preserve
  its tensor, and a continuous-selection obstruction for phase families.
* [Full family selector](FAMILY_SELECTOR.md): constructive extraction for every
  certified eight-core/five-tail row, including all larger bounded cores.
  Universal completion rules out point-only unfinished prefixes.
* `family_selector.py` implements that dispatcher without expanding unbounded
  speed combs. `check_family.py` tests all 1,287 eight-cores and all 1,093
  larger cores, including 4,096-bit tail controls.
* `witness.py` implements the rational interval / integer-speed contract.
* `check_witness.py` supplies literal geometry, equality, hostile-domain and
  large binary-input controls. `check_categories.py` supplies finite algebraic
  controls and exact phase-fibre examples.

Reproduce from this directory:

    python3 -B check_witness.py > normal.json
    python3 -B -O check_witness.py > optimized.json
    cmp normal.json optimized.json
    python3 -B check_categories.py > categories.normal.json
    python3 -B -O check_categories.py > categories.optimized.json
    cmp categories.normal.json categories.optimized.json
    python3 -B check_family.py > family.normal.json
    python3 -B -O check_family.py > family.optimized.json
    cmp family.normal.json family.optimized.json

The companion independent audit is in `../audit_lrc/`. The mathematical proofs
are not inferred from the finite tests. No new LRC family, unrestricted LRC
proof, realizability-topos construction, continuous uniform selector, or
literature-priority claim is made.

Independent proof and implementation audits passed for all three scoped
results. The full-family selector inherits the Session122 universal-family
theorem; its independent event-cell controls cover all 1,287 eight-cores,
all 1,093 larger cores, repeated bounded tails, a point-only completed cover,
exact one-tail equality, invalid inputs, and huge mixed tails. See
`../audit_lrc/AUDIT.md` for counts, limitations and reproducible receipts.
