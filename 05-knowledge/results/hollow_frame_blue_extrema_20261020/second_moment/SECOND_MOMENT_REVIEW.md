# Second-moment proof and validation review

Status: **PROVED; INDEPENDENT REVIEW PASSED**, 2026-10-20.

## Proof review

The author derived the exact signed-prefix overlap expression and its
formal-series reduction. A separate mathematical review independently checked:

- the normalization when the two final pair indices differ or agree;
- the constant relative direction on each shared adjacency interval;
- the expansion into edge-disjoint dominoes of weight 1/4;
- the additional terminal monomer of weight 1/2 when final indices agree;
- the orientations and factorials for contracted vertex blocks;
- the coefficient extraction giving the stated second moment.

The independent reviewer approved the theorem after reviewing these steps.
The complete proof is in `HOLLOW_FRAME_SECOND_MOMENT.md`.

## Exact computational checks

`verify_second_moment.py` does not import the frame producer or any of its
algorithms. It directly sums signed-prefix overlap factors over every
permutation and difference-sign vector, using exact rational arithmetic.
The coefficient formula and this independent overlap sum agree for k=1..7.

For k=1..5, their second moments also agree with the complete external-frame
census files `exact_k1.json` through `exact_k5.json`. The recorded JSON gives
SHA-256 hashes of those inputs and distinguishes these five census checks
from the two higher-order formula/overlap checks.

The verifier was rerun with and without Python optimization. Both modes
produce identical recorded results. All validation failures use explicit
exceptions, not assertions. A standalone run without census files also
passes and correctly marks every external-census check as false. A deliberately
tampered census second moment was rejected under `python3 -O`, confirming
that optimized execution does not disable the checks.

## Preservation and scope

This supplement adds new files only. It does not change the main maximum
census, source, witnesses, proof packet, or existing manifest.

The variances 89685/2 at k=6 and 1875825 at k=7 are exact consequences of
the proved formula and are independently checked by overlap sums. They are
not exhaustive external-frame results. No new maximum, higher-order frame
census, asymptotic expansion, or priority claim is made here.
