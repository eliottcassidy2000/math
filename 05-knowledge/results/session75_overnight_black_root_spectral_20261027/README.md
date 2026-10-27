# Session75 — overnight black-core, ROOT and spectral programme

**Checkpoint 1: audited spectral exception budgets and recovered inheritance.**
Research source: `96e35642110104cc547516a63f37eebcf930f707`.
This is a continuing overnight programme, not a final session close-out.

## What is established here

- **PROVED, audited:** an arbitrary-depth, disjoint substitution parse with
  quotient arities at most a fixed `b>=2` has skew nuclear norm at most
  `C_b n log₂ n`, where `C_b=2√(α(1−α))/h₂(α)` and `α=1/(b+1)`.
  The proof covers unbalanced combs, not just balanced trees.
- **PROVED, audited:** a flat target `KKᵀ=nI−J` incurs at least
  `max(0,ceil((n(n−1)−C_b n^(3/2) log₂ n)/4))` named arc exceptions.
  For fixed `b`, this is at least `(1/2−o(1))` of all arcs, with sharp leading
  constant when transitive baselines are available. Arc edits, additive
  rank-two corrections and fixed-size legal path complements have distinct
  budgets; arbitrary block reversals are not silently priced as rank two.
- **FINITE-EXACT, independently recomputed:** join/C3 arc-edit distances are
  `0,3,11` for Paley orders `3,7,11`. Full subset optimum vectors and attaining
  labelled parses are checked, not merely asserted terminal numbers.
- **PROVED lower bounds, not attained distances:** Paley versus regular C3
  towers needs at least `54` arc edits at order `27` and `9234` at order `243`;
  their fixed-size native-complement lower bounds are `3` and `39`.
- **VERIFIED inherited source recovery:** 21 source snapshots and 37501 literal
  primitive ROOT edges through seed `100001`. This verifies the affine carry
  interface, not an all-depth energy or entry theorem.

Read the [spectral theorem and scope](session75-spectral/RESULT.md),
[independent audit](session75-spectral-audit/REPORT.md),
[inheritance report](session75-inheritance/INHERITANCE_REPORT.md), and
[overnight plan](session75-root/OVERNIGHT_PLAN.md).
The producer's reviewed candidate, status-only promotion diff and metadata
are retained in its `status-promotion/` lineage; the separate independent audit
and final-byte acceptance bind the promoted result without changing mathematics.

## Boundaries retained

The fixed-`b` qualification is essential for the asymptotics. Source/sink or
score information does not substitute for the matrix hypothesis. A field of
order `27`, not integers modulo `27`, is used for Paley controls. Session74's
order-eleven native distances `6` to forest frames and `15` to transitive
frames are not the join/C3 arc-edit distance `11` above.

This checkpoint claims no native all-word reachability closure, ROOT positivity,
all-depth source-energy theorem, universal dyadic debt control, or LRC(14)
closure. Black-core, native ROOT transport and dyadic pairing-debt audit rounds
are not included as established results here. No historical novelty claim is made.
Concurrent mainline RESERVED stubs are preserved but are not dependencies.

## Reproduction and provenance

See [verification scope](VERIFICATION.md). From this directory:

```sh
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
python3 -I -B -O verify_provenance.py
python3 -I -B test_packet_corruption.py .
python3 -I -B test_provenance_corruption.py .
```

The frozen manifest selects current replay scripts; inherited source programs
are inert unless explicitly approved. Native reports, ledgers, inputs, modes
and normal/optimized receipts are preserved byte-for-byte.
