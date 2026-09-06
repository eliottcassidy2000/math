# Complete n8 stress: Booleanization lifts twenty-four weighted zero modes

Status: **FINITE-EXACT; full source and consequence independently audited.**
This is a complete finite universe, not an all-order theorem. The general
native nullity, free-pair transversality, and Boolean Laplacian gap remain
**OPEN**.

## 1. Question and inheritance

The closest proved mechanism is the conditional
[full Fourier repair](eulerian-boolean-fourier-repair-overnight-hexagon-sep05.md).
Write the actual unweighted triangle-toggle graph on Eulerian isomorphism
classes as `B=[[0,C^T],[C,0]]`, with C mapping even-edge functions to
odd-edge functions. Let S be the self-complementary switching-class
Fourier columns, restricted to even classes, and P the sums of columns in
each free complementary pair. The complete even Fourier basis is `[S P]`.
If `A=CP` is invertible, the independent kernel columns are

```text
K=S-P A^(-1)CS.                                         (1)
```

This requires the **actual Boolean C**. The canonical hostile is already
n=4: a self-complementary multiplicity-weighted Fourier mode is not
Boolean-null. The all-order four-entry minor in the
[native-kernel note](eulerian-boolean-kernel-overnight-hexagon-sep05.md)
also rules out repairing the operators by diagonal row/column scalings.
The corrected near miss is to infer native transversality from a weighted
rank or index. The least-used sidecar is the *whole free-pair block*, rather
than one putative eigenvector.

The complete n3..7 bank happened to have weighted nullity equal to the
parity index. The n8 stress is the first order in this bank where that
coincidence fails, so it tests a genuinely different obstruction.

## 2. Exact result and its actual consequence

There are `2^21=2,097,152` labelled Eulerian graphs on eight vertices and
243 isomorphism classes: 131 even-edge and 112 odd-edge classes. There
are 19 self-complementary dual classes and 112 free complementary pairs.
The weighted triangle matrix M and Boolean adjacency B have:

| Quantity | Exact value |
|---|---:|
| `nullity M` over Q | 43 |
| Number of self-complementary Fourier columns | 19 |
| `rank(CP)` over Q | 112 |
| `rank B` over Q | 224 |
| `nullity B` over Q | 19 |
| `rank(CS)` over Q | 19 |

The square block `A=CP` has determinants

```text
det A = 68     modulo 101,
det A = 946    modulo 1009,
det A = 57183  modulo 1000003.                           (2)
```

Any one nonzero residue proves its integer determinant is nonzero; no
rational-zero inference is made from a zero residue. Thus (1) applies at
n8, and all native kernel modes vanish on odd-edge classes. This establishes
rank and nullity without printing a large rational inverse. Since CS has
full column rank, **no nonzero combination of the nineteen original
self-complementary modes is already Boolean-null**.

Meanwhile M has 43 independent zero modes, certified using its complete
Fourier diagonalization. Of these, 24 occur in twelve free complementary
pairs. Booleanization removes the excess nullity while retaining the
forced index. This is a finite rank comparison, not a spectral-gap bound,
an eigenvalue-by-eigenvalue correspondence, or a general monotonicity law.

The next useful operation is native residual elimination, as developed
in the [isolated-cycle and complement-character note](overnight_hexagon_sep05_boolean_native_structure.md).
Its unitriangular minor is unconditional at every order; its residual rank
is not supplied by this n8 computation. The two methods retain different
coordinates and should be compared on the residual matrix, not identified.

## 3. Complete universe and independent paths

The compiled producer uses the anchored triangle basis of the entire
cycle space and the corresponding root-zero gauge of the dual cut quotient.
Adjacent transpositions generate full S_n orbits in both spaces. Nothing
is sampled; there are no inherited pruning filters.

Controls explicitly include:

- All labelled states, including **every triangle-neighbor row**, not just
  chosen representatives. Every full row agrees with the quotient row.
- A second, literal full-permutation path: for every representative and all
  `8!=40320` permutations, verify canonical minimality and the stabilizer
  count; each orbit size equals `8!/|Aut|`. Distinct minima plus the complete
  size sum certify coverage independently of generator traversal.
- An independent Burnside partition count; raw-matrix reversibility using
  orbit sizes; exact edge-parity bipartition; dual complement involution.
- The full identity `M Psi=Psi diag(lambda)` and separate nonsingularity
  certificates for both parity Fourier bases. These justify the exact
  weighted zero multiplicity rather than presuming completeness.
- Bit-for-bit comparison of both full matrices, all representatives, orbit
  sizes and complement maps against the previous Python implementation
  for n3..7. At n8, the full-permutation representative/size check is the
  independent native path; a second independent full n8 matrix census is
  **not** claimed.
- Positive, negative and singular determinant controls, including the
  nonzero integer101 that vanishes modulo101 but not modulo1009.

All arithmetic is integer, rational, or finite-field arithmetic. The
native matrices are never replaced with the weighted ones after the
Fourier identity has been checked. The C++ program uses bounded32-bit
states (at most28 edges), and64-bit reversibility products. The Python
script uses exact rational-domain rank for the small CS matrix. Semantic
checks remain active under Python optimization.

## 4. Reproduction and audit

```bash
python3 04-computation/eulerian_n8_fourier_stress_overnight_hexagon_sep05.py
python3 -O 04-computation/eulerian_n8_fourier_stress_overnight_hexagon_sep05.py
```

The driver compiles its adjacent C++ source with `g++ -O3 -std=c++17` in
a temporary directory. Normal and optimized outputs are byte-identical.
The independent source audit checked the coordinate actions, all chunk
bounds, permutation and orbit controls, the modular determinant algorithm,
the weighted/native distinction, and every rank consequence above; no
repair was requested. This is a source audit, not a claim of an additional
independent n8 enumeration by the referee.

Frozen files:

```text
C++ source abe81cc4f161aad73578cb3c75b1faf0606497bebe418d66655a9478ff4af733
Python source c35ed4320142e4f1b710fb8b408fc628fd83cfb7bf68401dc4ac64386a3e28cd
output a6f0c04274237b69447fd7c24545a7c102371bda45bb2df3823e01be5a9e470a
semantic db558b4a0c625e815dd6f02a4fabbbdfe4f29f5c27b1256f2a280e6b896b039d
```

The full-data hash at n8 is
`ea2ebdecd17ce7aa026b00edf5a71a1c24da6ddb9db2f1dc8a5c7dace59d55da`;
the CP matrix hash is
`8ac4fd229448580a09feddb644bea193199351fc07eda10eff7c86cdaa9cdba7`.
The [frozen output](eulerian_n8_fourier_stress_overnight_hexagon_sep05.out)
records the complete n3..8 summary, not a large floating-point spectrum.
