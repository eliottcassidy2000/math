# Verification scope and replay

The session supplies sixteen actual finite ordered FRACTRAN lists. Each lane
records its own input convention, all-order proof, finite execution limits,
and counterexamples to invalid shortcuts. In particular native evaluation
is not reported as a literal FRACTRAN run.

Run `python3 verify_session.py` to validate all manifests, reject corrupted
artifacts, and relocate every lane to independent temporary directories for
normal and optimized Python replays. Only the standard library is required.
The resulting `SESSION_REPLAY.json` records exact receipts and coverage.

- `fractran/verify.py`: eight compiled and two direct counters,602 primitive
  checks and14 hostile controls; separate source/IR/exponent/integer ranges.
- `distribution/check.py`: generic positive matrix programs, three
  multiplicity streams, literal-integer and valuation checks, exact costs.
- `core_reduction/check_core.py`:33868 labelled tournaments,2131019
  induced-equality tests, exact partition/budget controls and module hostiles.
- `spectrum/run_checks.py`: all odd thresholds1..21, complete56-class basis,
  separate fixed-Hamilton-cycle cover and434 literal checks, normalized
  polynomial controls and exact basis-set plateaus.
- `path_counter/check.py`: general evaluator, one-apex and square-shell
  programs, proof-backed inverse fibres; million-bound counts are native
  arithmetic results, not a million literal FRACTRAN executions.
- `independent/replay.py`: seven independent audits, each in both modes,
  including full integer first-eligible scans and independent graph censuses.

The first checkpoint's polynomial-padding hostile error is explicitly repaired
in [CORRECTION](spectrum/CORRECTION.md). No counts, strong inventories or
obstruction theorem depend on that erroneous distinction.

Classical source theorems remain CITED dependencies, not machine-formalized
proofs. All-order statements are proved structurally; finite tests corroborate
them. Universal realization of all other odd H values remains OPEN.

## Final root replay

PASS:14 isolated principal replays (seven verifiers, each normal and
optimized),9 manifests with208 checked entries, and3 artifact-corruption
negative controls. All normal/optimized receipt hashes match. Every replay
ran in its own relocated temporary packet; no original-workspace import or
network source was required. The independent driver also validates69 input
source/artifact pins and executes its seven normal/optimized audit pairs.
