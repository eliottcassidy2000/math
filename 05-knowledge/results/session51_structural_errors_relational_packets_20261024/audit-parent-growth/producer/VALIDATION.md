# Validation

## Source and input contract

The order-eight input is inherited from the complete ordinary-class Session50
atlas, with actual orientation-preserving vertex permutation actions. The source
JSON and its earlier completeness receipt are archived. The verifier checks
its byte hash against that receipt and compares every parent and every action
to `parents-eight.txt`. Thus the production input is not authenticated merely
by count, group closure, or mass. The inherited completeness proof uses actual
pairwise nonisomorphism, complete groups, and exact labelled mass2^28.

The independent structural correction proof and receipt are also archived;
SOURCE_PINS.json records their exact hashes and the frozen correction packet's
manifest hash. The input used for the predictor is5508, not a target count.
The order-nine error758 is a validation-only crosscheck after the atlas has
been generated; it never decides which states are selected.

## Builder replay

Normal-O3 and-DNDEBUG builds produce exactly identical receipt JSON and exactly
identical order-nine atlas bytes. In each mode:

- All1716608 pointed states are covered through actual parent subset actions.
-191480 classes are selected by an invariant singleton.
-232 singleton-free pointed states give56 more classes via minimum exact
  parent/subset keys, using2088 order-eight canonicalizations and7638 candidate
  eight-vertex labellings.
- Every recovered child action is checked as a permutation and against every
  directed arc; identities, distinctness and group closure are verified.
- Labelled mass is2^36, error is758, and the prior upper-trace identity agrees.

No ordinary nine-vertex canonicalizer occurs in the production source.

## Trace replay

The normal trace has a safe early exit precisely when the root becomes a
singleton. The-DNDEBUG -DFULL_REFINEMENT build disables that exit, completing
ordinary stable refinement for every pointed extension. Their191536 per-parent
summaries agree byte-for-byte. Both visit97213472 states and obtain upper
trace97336068, correction5508, and predicted count9733056.

Only refinement-round and singleton-shortcut counters may differ. The
arithmetic result, cell histogram, subset-orbit coverage, input action checks,
and per-parent outputs are compared exactly.

## Boundaries and portability

No count/trace output literal appears in either C++ source. No assertion is
used for validity; all guards throw exceptions and remain active with NDEBUG.
Header extraction is checked before allocating input groups, so empty or
truncated headers fail instead of using uninitialized sizes.

`reproduce.py` compiles both builder and both trace modes in a separate output
directory, validates pinned input bytes, checks exact output equality, and
writes a timed replay receipt. It requires Python3 and a C++17 compiler, but
no network or third-party Python packages. `verify_packet.py` passes in ordinary
and Python-O modes. Binaries and caches are excluded from all published payload
manifests; use-B or PYTHONDONTWRITEBYTECODE=1 for Python subprocesses.

Separate reviewer work is not substituted for the production proof. An
independent audit may canonicalize order-nine objects to check this generated
atlas, but any such audit-only calls must not be reported as production calls.
