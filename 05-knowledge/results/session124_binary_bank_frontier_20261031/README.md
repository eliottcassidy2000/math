# Finite binary-bank frontier: independently verified

For every tournament of order at most nine with `H>=9`, arbitrary one-vertex extensions cannot form a depth-two binary tree realizing both targets `2H-1,2H+1` at each nonterminal node. Independently generated isomorphism classes, exact automorphism orbit masses, and actual-child subset DP verify this **finite** statement.

The included 15-node positive certificate instead uses level offsets `(3,-3,27)` on one `H=5` seed, reaching leaves `73,75,...,87`. It is not an infinite compiler or an unrestricted spectrum theorem.

The response identity used here is documented in the companion [canonical cut-response theorem](../session124_tournament_cut_response_20261031/README.md); the present package independently implements and tests it.

Read **[FINITE_BANK_PROOF.md](FINITE_BANK_PROOF.md)** for the precise theorem, proof of finite coverage, arithmetic bounds, independent methods and exclusions.

Requirements: Python 3 standard library and a C++17 `g++` compiler. No network, third-party packages, external representative lists, or prebuilt executable is required.

```bash
python3 run_audit.py
```

The complete run independently generates all tournament isomorphism classes through order nine twice, computes automorphism orbit masses, audits the balanced depth-two bank predicate, and compares every potentially decisive grandchild count against direct subset DP. It also runs the literal-permutation phase certificate, the submitted polynomial cross-check through child order eleven, and the literal automorphism test. All Python checks run normally and with `-O`; C++ runs use `-O2` and `-O3 -DNDEBUG`. Explicit checks do not disappear under optimization. Paired stdout receipts must match byte for byte.

The audit takes several minutes. Its binaries stay in local `_build/`; its receipts are compact JSON. All subprocesses run from `_build/`, not from the source directory.

To additionally compare the producer's scratch representative lists with the independent universe:

```bash
python3 run_audit.py --producer-input /path/to/original/scratch/bandwidth
```

This optional comparison was performed during the original audit and is retained separately in `PRODUCER_BIJECTION_RECEIPT.json`. The publication and relocation runs deliberately omit this argument, so their `AUDIT_RECEIPT.json` correctly records `producer_representative_bijection_checked: false`. It requires `reps1.txt` through `reps9.txt` but does not use them to generate the independent universe or prove coverage. Do not publish the large representative lists or binaries.

`AUDIT_RECEIPT.json` records final source hashes, compiler, flags, output hashes and execution results. `SUBMITTED_SOURCE_HASHES.json` identifies the producer inputs that were reviewed. `PUBLIC_FILES.txt` is the compact publication allowlist; `MANIFEST.json` and `SHA256SUMS` bind the delivered bytes. `PORTABILITY.json` records the actual full relocated replay. The copied producer polynomial and helper are snapshots solely for an explicit cross-implementation check; the main census does not import them.

## What may be stated

* Independently verified finite nonexistence: arbitrary one-vertex cuts, standard targets `2H-1,2H+1`, both branches complete through depth two, parent order at most nine, parent `H>=9`.
* Independently verified explicit positive tree: seed `H=5`, three offsets `(3,-3,27)`, 15 nodes, leaf values `73,75,...,87`.

No all-order bank obstruction or global forbidden-value theorem follows. Other offset/phase negative searches are not promoted by this audit.

A supplementary `SANITIZER_RECEIPT.json` records an undefined-behavior/bounds-sanitized build on the complete order-seven census and three order-ten full-cut cases. Its exact flags and test arguments are recorded there; it is an additional safety check, not a claim of a sanitized full order-nine run.
