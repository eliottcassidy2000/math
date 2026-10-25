# Cyclotomic support dynamics — Session 59

## Scope

An explicit pointed-monoid projective model over `A_m=C_m∪{0}`, with complete proofs of its support count, power dynamics, unit-group orbits, and a fixed-point compiler for simultaneous diagonal twists, coordinate permutations, and power maps. A finite-field descent corollary shows that identical projective totals can conceal different rational support strata. The Session-58 squarefree gap index is exactly the `m=2` specialization; higher Gaussian exponents remain a sign quotient, not a larger scalar quotient.

The ordinary frame/count construction is established literature. This packet offers a reusable proved synthesis and implementation, not a literature-priority claim. No new hourglass existence/nonexistence result, no F1 descent of the additive equations, and no actual destroyed Schur triple are claimed.

## Files

| File | Purpose |
| --- | --- |
| `CYCLOTOMIC_SUPPORT_DYNAMICS.md` | Definitions, category/base guardrail, full proofs, arithmetic boundaries, primary references |
| `verify_cyclotomic.py` | Standalone standard-library exact verifier; every gate is explicit, not a Python assertion |
| `verification.json` | Machine-readable PASS certificate, exact finite ranges, and example outputs |
| `normal.log` | Verifier stdout under normal Python |
| `optimized.log` | Identical stdout under `python3 -O` |
| `SOURCE_SHA256SUMS.txt` | Hashes and repository-relative paths of the two inherited Session-58 source notes |
| `SHA256SUMS.txt` | Relative-path hashes of this frozen proof/test packet |

## Reproduction

From this directory, using Python 3.9 or newer:

```sh
python3 verify_cyclotomic.py > normal.log
python3 -O verify_cyclotomic.py > optimized.log
cmp normal.log optimized.log
sha256sum -c SHA256SUMS.txt
```

The verifier overwrites only `verification.json` alongside itself and prints the same JSON to stdout. It reads no repository source or network data. The proofs, not bounded computation, establish universal claims. The stored run passed 6,192 semilinear fixed-count cases, 440 ordinary power fixed counts, 468 complete finite functional graphs, 3,200 exact tail-period cells, 2,009 composition point checks, 21 cyclic Burnside cases, 32 unit-group orbit cases, 144 coefficient-map cases, and 36 descent-formula cases (22 also checked by raw point enumeration).

## Source basis

The exact primary URLs and sections are recorded in the proof note. Thas supplies the cyclic pointed-monoid frame convention and known point count; Deitmar's original unpointed convention is explicitly distinguished. The arithmetic gap bijection and valuation/deletion results are inherited from the two source-hashed Session-58 notes. Neither source was modified. No git operation was performed by this worker.
