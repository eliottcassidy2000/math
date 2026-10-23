# Second frozen audit: compiled shared-state equations and owner-indexed enumeration

Status: PROVED theorem review; FINITE-EXACT independent replays PASS.
This packet supplements, and does not rewrite, the earlier first audit pins.

The current residual README retains the previously audited bounded conformal
atom theorem, all-size finite-domain compiler and frame comparison criterion.
Its reference internal isomorphisms explicitly preserve retained marks. The
added executable compiler validates complete frame domains and rejects duplicate
whole-cell frames. Independently replayed check_compiler.py normal and -O:
32,768 assignment-row comparisons, 4,096 encode/decode inversions, and 324
actual six-vertex block/frame comparisons all pass.

OWNER_INDEXED_ENUMERATORS.md is sound under its declared complete finite-state
carrier. The forest recurrence keeps the shared parent state fixed until every
incident child factor has been multiplied, then sums over the state. Separator
contraction is exact with the running-intersection property and each factor
assigned once. Scalar min/max messages compute true extrema; equality in every
coordinate is equivalent to a nonempty constant vector observation. Zero
polynomial means infeasible and is not called a uniquely realized wiring.

The distinction between 1+zw and (1+z)(1+w) is exact. Identifying polynomial
variables after discarding a shared owner does not recover the missing equality
constraint. The zero-slot elimination example is also exact: a forced-absent
zero slot induces the constant-one factor, while a freely selectable zero slot
induces multiplicities (2,1,1,1); neither equals reconnecting its two neighbors.
The coupled cancellation example X00+Y01=1 on the two-switch fibre is correct.

A wording overreach was corrected during audit: variable stabilizers alone do
not invalidate group-order division. For a complete finite invariant inventory,
sum_x |Stab_G(x)| / |G| is exactly its orbit count. The revised text requires a
specified common group action, a complete inventory and the appropriate weights.
Frame-comparison counts by themselves are not unlabelled tournament counts.

Independent replay of check_enumerators.py in normal and optimized Python
matches every coefficient for all 4,096 binary four-vertex-star relation systems,
with 2,727 nonempty min/max certificates. All new finite checks use explicit
exceptions rather than optimization-removable assertions. Reviewed source
hashes are recorded alongside this report.

The Fibonacci main theorem file remains unchanged at e43d0b155982b8f40b8fef5038bf2dd392bb552218b0d27ea026332cf2daf541 and covered by the first audit. Native order31 attachments and the newly found indexed-versus-unindexed arc hostile are excluded from this second packet and require their own audit.
