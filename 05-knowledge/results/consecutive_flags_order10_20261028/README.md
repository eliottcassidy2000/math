# Consecutive tournament flags, order ten, and closed observation carriers

**Status:** FINITE-EXACT computation; elementary proofs in the companion notes;
same-agent alternate-path audit, not an independent human review. No claim to
solve LRC(14), the graceful tree conjecture, or unrestricted carrier minimization.

The refined parent-only generator now constructs **all 9,733,056 order-ten
tournament classes**, from 97,213,472 pointed extension orbits. It makes only
15,107 competing-root canonical queries, all at order nine. The exported bank
is a dictionary of nine-vertex parents and delta-coded attachment addresses:
829,899 compressed bytes. A separate full-child canonicalizer checks every
representative, all automorphism orders, uniqueness, and total labelled mass
2^45. The generation run took 10.88 seconds and the separate audit 15.46 seconds
on this local machine; these are receipts, not cross-machine benchmarks.

The square-matrix suggestion has an exact structural interpretation. An aligned
pair consisting of an n-tournament and a rooted (n+1)-tournament corresponds,
up to simultaneous relabelling, to a self-converse (2n+1)-tournament. The
construction records a normalized reversing involution. Forgetting that
involution preserves the unlabelled class bijection but can enlarge the
automorphism group: already 1 becomes 3 at order three. This is an attachment
boundary, not a reason to discard the correspondence.

For the explicitly weaker observer that sees only the child class and its
deletion-parent class, we compute the smallest strongly lumpable carrier under
native fixed-path complement and a uniform fresh re-marking preserving the
distinguished vertex. At order seven it has **207,070 states instead of
229,376 literal path/root states**. This protocol does not expose named
Hamilton guards or arbitrary attachment operations. Those richer languages
have the inherited faithful-reconstruction obstructions.

## Packet

- [Generation proof and bank format](GENERATION.md), [paired operations and carrier proofs](THEORY.md),
  [historical recovery and ETH scope](HISTORY_AND_COMPLEXITY.md), and [audit](AUDIT.md).
- `verify.py`, `source_pins.json`: lightweight validation or `--full` rebuild
  of every producer and the complete alternate bank audit.
- `check_structure.py`, `structure_audit.json`: literal permutation/Hamilton
  reset audit, deletion forests, signed scores, native nonclosure and81/27.
- `generate.cpp`, `generation.json`: complete constructive census through ten.
- `verify_bank.cpp`, `verification10.json`: full alternate canonicalization audit.
- `publish.py`, `bank_manifest.json`, `parents9.u64.gz`, `flags10.varint.gz`:
  reproducible bank and checksums. Run `python3 publish.py` from this directory;
  a C++17 compiler named `g++` is required. Temporary files stay outside the repo.
- `flag_states.cpp`, `closed_flags.py`, `closed_flags.json`,
  `closed_flags2.json.gz` through `closed_flags7.json.gz`: complete finite
  classifiers for the declared reset/complement protocol.
- `pair_bridge.cpp`, `pair_bridge.json`: all square-bit inputs through n=4,
  checking the doubled class correspondence and the full-group/centralizer gap.

## Inheritance and portfolio

**Anchor:** exact order-ten generation and native complement/re-marking closure.
The closest proved mechanism is stable ordered-cell root selection with exact
parent transports in [the root-port packet](../root_ports_closed_carriers_20261028/CARRIERS.md).
The canonical hostile is the twelve-vertex equal-cyclic-triple 81/27 pair; its
full automorphism action is retained. The corrected near miss is replacing
exact root keys by refinement colours alone. The least-used sidecar is the
ordered path forest left after deleting an interior marked vertex, from
[the strong-atom packet](../strong_atom_continuation_vitali_20261026/README.md).

**Niche:** aligned consecutive pairs and reversing-involution doubles.
**Wildcard:** what an ETH failure changes for compact decision/counting
certificates, and what output size and erased information still forbid.
The live concept board is: canonical root; relative alignment; reversing
involution; ordered path forest; future-observation quotient; output cost.

Earlier order-ten work already supplied a nauty census and a parent-nine scalar
trace. The advance here is this refined generator's complete bank, compact
attachment-address export, full alternate audit, and the paired-operation
analysis. It is not the first order-ten tournament census.

Eplett's 1979 representation theorem is the classical predecessor of the odd
self-converse bridge. Concurrent Session89 recovered a paired codec and exact
double-coset/stream formulas; this packet integrates them and adds the even and
odd exhaustive centralizer checks. The signed B/D root operator connects the
paired score map to graceful differences and square-sum equations, while
retaining their unsolved label-selection constraints.

The closing integration also routes Session89's all-size rectangle-code future
kernel, giving a64-state arithmetic decision carrier at n=3. Its named-toggle
and paired-relabel protocol is distinct from this packet's Hamilton-reset
protocol; THEORY records the exact connection and the F64 kernel obstruction.
