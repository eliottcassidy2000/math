# Every allowed odd tournament path count through 999,999

**Finite, verified construction theorem.** Every positive odd integer
`H <= 999999`, other than `7` and `21`, is the number of directed Hamiltonian
paths of an explicitly constructible finite tournament.

This does **not** prove infinite completeness, that this grammar is universal,
or the square-sum conjecture. It does not independently re-prove the global
exclusions of 7 and 21. No novelty claim is made.

## Evidence and proof structure

* Exact singleton-spine transitions, with all-finite-word pruning proof:
  [PROOF.md](PROOF.md).
* All labelled seeds of orders 1–5 and every external cut give 273 states.
  Actual `(order, edge-code, cut-mask)` seed witnesses are retained in
  `constructors/seed_roots.json`.
* The primary subset-DP/state-deduplication implementation visits 10,621,974
  states and obtains 499,980 different counts below 1,000,001.
* Independent literal-permutation seed enumeration and a word-tree DFS
  **without global state deduplication** visits 37,079,906 nodes and obtains
  exactly the same support size and 18 allowed misses. The independent audit
  verifies that these are the only misses, not an assumed interval.
* Four explicit triangle-block witnesses have counts 103061, 594793, 774073,
  949777, on 16, 20, 19, 20 vertices. Their complete adjacency matrices and
  constructor recipes are retained. Independent forward and backward endpoint
  DPs count paths directly from the actual edges, not transfer coordinates.
* Fourteen ordinal products of directly attained factors fill all other misses.
  Every factor is checked against the independently generated support bitmap.
  The product for 618723 is `49*9*23*61`, never a fictitious factor with count 7.
* The 499,980 direct values plus 18 repairs give the 499,998 allowed positive
  odd values at most 999,999. Witness existence is inductive over exact
  constructors; this packet does not store 499,998 separate adjacency files.

See [audit/AUDIT.md](audit/AUDIT.md) for the independent proof audit.
The modulo-three and general p-adic grammar restrictions in `PROOF.md` are
separate proved statements; finite hostile tests are not substitutes for their
proofs. No arbitrary-block completeness claim relies on singleton pruning.

## Reproduce all checks

Requirements: Python 3.9 or later, C++17 compiler `g++`, no third-party Python
packages. Allow several GiB of RAM for the primary 10.6-million-state traversal
and several minutes of CPU time. All input paths resolve from each source file;
the caller's current directory is irrelevant. Outputs go only to the explicitly
selected work directory. Run from the package root:

```sh
python3 run_checks.py --work-dir /tmp/spectrum-check
```

The driver runs the primary traversal and all Python proof/DP checks both
normally and with `-O`, compares outputs, compares regenerated seed data with the
retained bytes, compiles both independent C++ programs at `-O0` and `-O2`,
compares support bitmaps, verifies all factor repairs and all four matrices,
and checks edge-flip and malformed-input rejection. It also exercises the
portable witness constructor in normal and optimized modes and independently
counts the five generated tournaments, including H=999999 on 22 vertices. Fresh execution receipts are retained
under `receipts/`; regenerated receipts appear in `WORK/receipts/`.

Individual main commands (run from the package root):

```sh
mkdir -p /tmp/spectrum-check/independent
python3 constructors/seeded_words.py 5 1000000 --output-roots /tmp/spectrum-check/seed_roots.json
python3 -O constructors/seeded_words.py 5 1000000 --output-roots /tmp/spectrum-check/seed_roots.optimized.json
python3 constructors/audit_seeded_words.py
python3 constructors/check_target_witnesses.py
python3 constructors/check_block_targets.py
python3 constructors/check_mod3_and_products.py
python3 constructors/check_padic_filtration.py
g++ -O2 -std=c++17 audit/independent.cpp -o /tmp/spectrum-check/independent-audit
/tmp/spectrum-check/independent-audit /tmp/spectrum-check/independent
python3 audit/check.py --work-dir /tmp/spectrum-check/independent
g++ -O2 -std=c++17 audit/adjacency_dp.cpp -o /tmp/spectrum-check/adjacency-audit
/tmp/spectrum-check/adjacency-audit audit/adjacency_*.txt
python3 audit/mod3.py
python3 audit/mutation.py /tmp/spectrum-check/adjacency-audit
```

Generated executables, `support.bin`, temporary mutations and Python caches are
not package inputs or publication files. The independent audit imports no
primary constructor or transfer code. Its 64-bit DP accepts orders
at most 22 and explicitly checks every addition for overflow. The four main
repair witnesses have at most 20 vertices, so `20! < 2^64` already bounds
every partial count there; the five reconstruction samples also pass without
overflow (including the two 22-vertex examples).

## Construct a requested witness

```sh
python3 construct_witness.py 103061 --adjacency /tmp/witness.txt --recipe /tmp/recipe.json
/tmp/spectrum-check/adjacency-audit /tmp/witness.txt
```

`construct_witness.py` accepts every positive odd target at most 999999 except
7 and 21. It uses all-word singleton traversal with exact-self-loop elimination,
not an unproved depth bound, plus the frozen block repairs and fourteen product
recipes. It writes the full actual tournament and a reproducible recipe.
Search may take time and the witness is not claimed to have minimum order.
Five checked sample adjacencies and recipes are retained under `examples/`,
including the endpoint H=999999 and the product repair H=34869, each on
22 vertices. Independent checked-arithmetic DP counts both directly.
Product witnesses or some large targets may exceed the independent DP's
22-vertex limit; their existence and count then follow from the proved
constructor/product identities, not from a purported completed direct DP.

Adjacency format: first line `n target`; next `n` lines are binary rows, with
row `i`, column `j` equal to 1 precisely when `i -> j`. Diagonal entries are 0.
Seed edge codes enumerate `i<j` in lexicographic order, bit 1 meaning `i -> j`.
A cut-mask bit 1 means apex `x -> i`. Singleton blocks are `0` or `1` with
that same apex convention. `Cj` is a directed triangle `0->1->2->0` appended
as a right ordinal block, with x beating local vertices `0,...,j-1`.
An ordinal product puts all edges from each earlier component to each later one.

## Integrity and publication boundary

`SHA256SUMS` covers research sources, witness data and receipts. `WHITELIST.txt`
is the exact intended repository file list, including these two manifests.
`publication/` contains local-only Git-tree JSON payload batches plus their
manifest; those transport envelopes are not themselves research tree entries.
Every batch is at most 140,000 UTF-8 bytes and every entry is a regular text blob
under `05-knowledge/results/session124_million_path_spectrum_20261031/`.
No executable, support bitmap, cache, external discovery path or unverified
scratch note is included. This package was prepared without modifying the
source worktree or publishing a commit.
