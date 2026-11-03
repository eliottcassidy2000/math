# Verification scope

**PASS, with explicit mathematical boundaries.** This addendum has two frozen,
separately implemented finite computations plus their independent comparison.
It does not promote any unproved general MSOS statement or identify the user's
lambda referent.

## Exact universe and independent paths

- All 25 pairs `(a,b)` in F5, exactly 20 satisfying the short-Weierstrass
  nonsingularity condition. All 20 are counted over F5, F25 and F125.
- Producer F125 arithmetic uses `z^3+z+1`; the auditor uses `z^3+2z+1`.
  Both count points directly, including infinity, rather than manufacturing
  the extension counts from a recurrence. Field inverses and finite-field
  Frobenius identities are checked.
- The 12 scaling-isomorphism classes, nine trace states, twist action, and
  every common point-count tuple agree between implementations.
- Magic arrangements use an exhaustive centered parameterization in the
  producer and all `8!` center-fixed permutations in the auditor. There are
  eight for the cubic-extension counts, zero for the literal and symmetric
  cubes on the same nine traces. F25 is separately ruled out by its mean.
- All 1024 labeled five-vertex tournaments are classified under 120 vertex
  permutations. Characteristic polynomials use a Leibniz determinant versus
  Newton identities. The unequal quotient fibers and involution fixed-point
  counts are therefore actual obstructions, not just matching numerology.

## Normal, optimized and hostile replay

`replay.py` runs both checkers, both mutation suites and the auditor's producer
review in normal and optimized Python. It compares each stdout byte for byte
with its frozen certificate and directly compares the common mathematical
outputs. The producer's 10 mutation cases (five in each mode) are rejected and
two unchanged runs pass. The auditor rejects 22 negative cases and accepts two
unchanged certificates per suite invocation. No correctness gate relies on
Python `assert`.

The top-level manifest inventories all packet payload files except itself;
its SHA-256 and byte counts are independently checked before and after replay.
The producer and auditor keep their own original manifests. The top-level
source-pin verifier checks all six copied snapshots offline using both Git
blob identity and SHA-256. With `--repo /path/to/math`, it additionally checks
all seven contextual source identities against their exact published commit
and path. No Git repository or network is required for mathematical replay.

```sh
python3 replay.py > /tmp/frobenius-replay.json
cmp /tmp/frobenius-replay.json REPLAY.json
python3 -O replay.py > /tmp/frobenius-replay-optimized.json
cmp /tmp/frobenius-replay-optimized.json REPLAY.json
python3 verify_source_pins.py
python3 -O verify_source_pins.py
python3 verify_manifest.py
```

## Publication and preservation

The addendum is prepared against published mainline
`fda1d845ebbf0f0d632cfeeac96afb517bece1db`. Its mathematical repository inputs
are pinned to the already published Session143 final checkpoint
`74e81ed37c9175d94e8fad57fc5809cf78cee572`, not an unpublished local commit.
The original 183-file packet has tree
`b3284c991abb0abcb6738506ee4641fe31b542d5` at both commits and remains unchanged.
The only existing file extended is `00-navigation/S143.md`, by a scoped suffix.
All concurrent tree entries and startup-document budgets are preserved.

Documentation checking is deliberately scoped to these owned files, their
curated publication links, conflict markers, source identities, and unchanged
startup surfaces. No repository-wide `check_docs.py`, mistakes-ledger read,
Lean build, or external-PDF redistribution is claimed. The historical source
snapshots retain their original link context rather than being rewritten.

## Cited versus independently checked

Frobenius characteristic polynomials, extension point-count formulas,
short-Weierstrass normal forms/isomorphisms and the stated Hecke normalization
are cited standard inputs with primary locators in the two source-pin files.
The cubic identities and finite consequences are checked here. The complete
all-real cubic classification is inherited from the original audited packet,
not re-proved by a finite q-search. The formal determinant corollary is a direct
substitution into that classification; it does not assert curve realization
for arbitrary real determinant.

The [final wrapper audit](WRAPPER_AUDIT.json) separately confirms the older
lambda identities, determinant-sign distinction, isogeny boundary, and the
short obstruction to a common base/extension magic placement.

The older golden-field lambda annotation and the Frobenius eigenvalue pair
have opposite quadratic norm/determinant signs and are not identified.
No ordering of curve classes, canonical representatives, squarehood map,
curve-level cycle successor, or all-extension theorem is produced. The bounded
follow-up stops once the finite bridge, its controls and its published bytes
are verified.
