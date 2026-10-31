# Packing obstructions for tournament odd-cycle realization

**Verified restricted-family results, not a proof of the unrestricted
Hamiltonian-path spectrum conjecture.**

1. Across all tournament orders, the Hamiltonian-path counts at most 73 with
   directed-triangle packing number at most one are exactly

   ```
   1,3,5,9,11,13,15,17,19,23,25,27,31,33,37,39,
   45,47,49,51,53,61,63,65,69,73.
   ```

   In particular the ordinarily attainable values
   **29,35,41,43,55,57,59,67,71** require two disjoint triangles.
2. Triangle-core presentations with packing at most one have an exact
   forbidden-subsequence basis: **240 words of length three and 35 of length
   four**. They admit at most 84 binary free phases separated by at most 83
   exceptional letters, giving `O(m^83 2^m)` presentations of spine length m.
3. For fixed k, the corresponding bounded-packing isomorphism-class count is
   at most `O(n^P_k 2^(kn))`, with an explicit effective constant P_k. This is
   **not** a bound on the density of attainable H-values.
4. **K14 is not the literal odd-cycle conflict graph of any tournament**,
   although `I(K14,2)=29` is attainable. Avoiding the factor `1+3x` and having
   an allowed numerical value do not certify conflict-graph realizability.

The all-order cutoff uses Busch's cited strong-tournament bound, not an
assumption that a finite census is complete at every order. The word encoding
is surjective, not a canonical labeling. Literal cycles sharing a vertex
support are kept distinct.

## Read and reproduce

- `RESULTS.md`: proof, scope boundaries, nine positive witnesses, and historical
  poison-factor reconciliation.
- `AUDIT.md`: independent proof and computational audit.
- `CITATION_CORRECTION.md`: source-number correction and lineage.
- `PROVENANCE.json`: exact original input identities and portability changes.

Python standard library only. No downloaded census, external executable, or
repository checkout is needed. From any working directory run:

```sh
python /path/to/package/run_checks.py
```

This runs all four programs in both normal and optimized mode, compares every
output with frozen audited baseline hashes, and verifies mode agreement.
The default writes temporary receipts and prints the summary. To retain them:

```sh
python /path/to/package/run_checks.py --output-dir /path/to/new-receipts
```

Individual programs are also directly runnable:

```sh
python /path/to/package/check.py
python /path/to/package/pattern_check.py
python /path/to/package/check_witnesses.py
python /path/to/package/independent_check.py
```

The first three checkers are exact copies of the audited original programs.
The independent checker changes only its source-data directory to its own
directory; its algorithm and output are unchanged. It independently rebuilds
all 37,449 presentations, counts paths on all 3,624 accepted presentations,
regenerates the forbidden basis, and validates all nine ordinary witnesses.
It performs **38,207 explicit checks per mode**. All programs use explicit
exception checks rather than optimization-sensitive assertions.

`receipts/` contains fresh normal/optimized outputs, their portable aggregate
verification, and the successful relocated-package test. `MANIFEST.json` and
`SHA256SUMS` freeze the package. Publishing transport metadata, payloads, and
the strict target whitelist live separately in the preparation directory's
`publication/` and are not part of the mathematical package.

**Boundary:** no claim of novelty priority, unrestricted spectrum completeness,
or H-spectrum density follows. The global claim that only 7 and 21 are missing
is not proved here.
