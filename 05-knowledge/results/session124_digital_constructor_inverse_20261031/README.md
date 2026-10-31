# Digital normal forms and exact inverse tournament constructors

**Proved within a specified actual-seed construction language.** This packet
provides exact digital identities, infinite frozen-fiber descriptions, and a
complete finite decision/reconstruction algorithm for non-3-divisible targets
in the singleton/triangle language. “Complete inverse” does **not** mean that
every target has a positive answer. No infinite tournament-spectrum coverage,
new forbidden unrestricted count, or square-sum conjecture proof is claimed.

## Results

* For chronological singleton bits b_i, set `B=sum_i(1-b_i)2^i`. Then
  `F_word=F0+qB+g(B mod 2)+h(f(B)-1)`, with the exact recursive digital f.
* An exact binary/run formula gives
  `f(n)=n(1-rho_2(n))+O(log(n+1))`. Here **rho_2 is the explicitly defined
  base-two radical inverse, not a 2-adic valuation**.
* Unique backward q ancestry gives the invariant `F-Phi_h(q)`. Actual source
  nodes remain attached to their offsets; arbitrary positive states are never
  substituted for actual tournament seeds.
* Triangle blocks are two-bit macros in normalized q, with explicit local
  defects in the remaining coordinates.
* For non-3 outputs, modulo-three identities permit only zero, one, or two
  adjacent nonconstant triangles. After removing trailing singleton 1s, a
  nonempty witness satisfies `q_final<4F`. This yields the exact inverse
  normal form without unsound arbitrary-block pruning or a depth cutoff.

See [PROOF.md](PROOF.md) for all proofs, conventions, exact scalar inverse,
the global *language-specific* exclusion of 13, fiber tables, complexity
limitations, and the remaining universal-existence obligation.
[Independent audit](audit/AUDIT.md) checks the proof quantifiers separately.
The straightforward inverse is pseudo-polynomial in the numeric target, not
polynomial in its bit length.

## Actual inputs and independent checks

`actual_seed_roots.json` retains all 273 distinct seed states with actual
(order, edge-code, cut-mask) witnesses. The independent auditor reconstructs
all **33,866 labeled seed/cut profiles**, every retained witness, and compares
5,733 inverse decisions through 63 with an independently saturated forward
closure having **no word-depth limit**. It checks 111 returned witnesses by
literal adjacency DP. The independent fiber audit reconstructs the complete
finite boundary data rather than assuming the listed offsets are available.

`frontier_repairs.jsonl` contains the four actual block recipes for 103061,
594793, 774073 and 949777. The independent inverse audit reconstructs their
known-Q macro decodings and actual graph orders. Their full adjacency matrices
are retained in `audit/`; a separate checked-arithmetic C++ endpoint DP directly
counts all four graphs. These controls do not assert an enlarged coverage range.

Producer and independent algorithms are distinct. The independent auditors
import the producer only for comparisons against its output, not for their
own seed profiles, forward closure, graph construction, DP, or fiber closure.
All validation remains active under Python `-O`.

## Reproduce everything, including relocation

Requirements: Python 3.9+ standard library, g++ with C++17, and enough memory
for the retained graphs' subset DPs (largest order 20). No Python packages or
network access are needed. Use a fresh working directory outside this package:

```sh
python3 run_checks.py --work-dir /tmp/digital-inverse-check
```

The driver:

1. runs all producer checks, scalar/macro inverse checks, and fiber generation
   normally and with `-O`;
2. runs both independent auditors normally and with `-O`;
3. recompiles and runs the four full-adjacency controls at `-O0` and `-O2`, then
   checks that edge mutations are rejected;
4. copies this package to an unrelated directory, repeats the complete suite
   from a different working directory, and compares every substantive output.

Every data path resolves from the source file. Generated fibers, executables,
and fresh receipts remain in the selected work directory. Retained inputs are
not overwritten. `receipts/verification.json` records the completed replay and
source/input hashes. `receipts/relocation_verification.json` records the second
complete run.

Individual commands from the package root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 check.py
PYTHONDONTWRITEBYTECODE=1 python3 check_macro_inverse.py
PYTHONDONTWRITEBYTECODE=1 python3 fibers.py --output /tmp/fiber_boundaries.json
PYTHONDONTWRITEBYTECODE=1 python3 audit/independent_check.py
PYTHONDONTWRITEBYTECODE=1 python3 audit/structure_check.py
```

The `digital` module exposes `singleton`, `triangle`, `run_blocks`,
`inverse_addresses`, and `inverse_non3`. The last function takes a target not
divisible by three and one **actual** seed state `(h,F,q,g)`, returning an actual
block word or `None`. Stored seed JSON uses `(h,F,q,t)`; convert `g=h-t` first.
The proof explains the seed and adjacency encodings in full.

## Curation and integrity

`CURATION.json` records original audited source hashes and the minimal curation:
removal of one unused legacy helper with an inaccurate docstring, explicit
radical-inverse notation, portable audit inputs, and an explicit fiber-output
path. The actual forward actions and exact inverse algorithms are unchanged.
Historical independent receipts remain under `provenance/`; current receipts
are under `receipts/` and match the historical mathematical results, apart from
the intentionally changed producer-source hash.

`SHA256SUMS` covers proof, code, actual data and receipts. `WHITELIST.txt` is the
exact intended repository tree. `payloads/` holds local-only Git-tree JSON
transport batches, each at most 125,000 UTF-8 bytes, targeting
`05-knowledge/results/session124_digital_constructor_inverse_20261031/`.
No executable, cache, support bitmap, external discovery path or unused scratch
program is included. Preparation did not modify the repository worktree or
publish anything.
