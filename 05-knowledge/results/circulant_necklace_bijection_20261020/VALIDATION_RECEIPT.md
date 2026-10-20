# Validation receipt

2026-10-20. All listed checks completed successfully. The mathematics is proved in the companion note; this receipt concerns the implementation.

## Reproduction

    python3 verify_necklace_bijection.py --output verification.json > verification.out
    python3 -O verify_necklace_bijection.py > verification.optimized.out
    cmp verification.out verification.optimized.out

Normal and optimized receipts are byte-identical: **1,757,041 checks**. The checker uses explicit exceptions rather than Python assertions, so optimization does not remove its validation conditions. It requires only the Python standard library.

## Complete finite universes

- Every binary necklace for m=1,...,16. Its indegree is independently obtained by constructing the full derivative graph and compared with the integration solver and primitive-root classification.
- Every odd-total source and every even-primitive-root target in that range, with forward/inverse identities, surjectivity, odd-period preservation, all cyclic shifts, stable-subspace test, and at-most-2^v2(m) path bound.
- The alternative generic lexicographic cut rule, with its own forward/inverse identity on every source in the same range.
- Every labelled prime circulant tournament and every symmetric connection word at p=3,5,7,11,13,17,19,23,29,31. All multipliers are tested as cyclic shifts; translations fix these circulants. The labelled forward map is a surjection, induces a type bijection, and has the proved dyadic fiber sizes.
- Independent full graph-automorphism/sign checks for every graph connection word at p=3,5,7,11,13. These use adjacency-preserving backtracking, without assuming the affine criterion. Empty and complete graphs have direct certificates.

## Selected results

| p | Labelled tournaments | Labelled sign-even graphs | Types |
|---|---:|---:|---:|
| 3 | 2 | 1 | 1 |
| 5 | 4 | 1 | 1 |
| 7 | 8 | 4 | 2 |
| 11 | 32 | 16 | 4 |
| 13 | 64 | 28 | 6 |
| 17 | 256 | 117 | 16 |
| 19 | 512 | 256 | 30 |
| 23 | 2,048 | 1,024 | 94 |
| 29 | 16,384 | 8,128 | 586 |
| 31 | 32,768 | 16,384 | 1,096 |

For m=1,...,16 the observed maximum ordinary derivative path lengths are

    1,2,1,4,1,2,1,8,1,2,1,4,1,2,1,16.

Each equals the proved upper bound 2^v2(m) in this tested range. The one initial derivative from the antiperiodic tournament word is additional.

## Example type table

For p=13 and primitive root g=2, the first six bits of the least antiperiodic tournament representative map as follows. Graph words are displayed up to rotation; k counts ordinary derivative path edges after the initial tournament derivative.

| Tournament half-word | Graph necklace | k |
|---|---|---:|
| 000000 | 000101 | 2 |
| 000010 | 011011 | 2 |
| 000100 | 010111 | 1 |
| 000110 | 001111 | 2 |
| 001010 | 000011 | 1 |
| 001100 | 000000 | 2 |

The table is only illustrative; the algorithm does not look it up.

## Independent coverage

The root's separately implemented `check_bounded_routing.py` does not import this producer implementation and independently verifies the necklace routing for every m≤14. An independent mathematical review of the root's general Fitting theorem is supplied in `REVIEW_BOUNDED_FITTING_ROUTING.md`. A separate reviewer should audit the prime-circulant translation before promoting its proof status to independently audited.
