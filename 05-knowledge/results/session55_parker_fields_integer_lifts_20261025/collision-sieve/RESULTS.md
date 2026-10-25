# Collision sieve: bounded results

Center roots 1 through 100000 inclusive; arbitrary content retained.

- Exact candidate pairs: 347561
- Remaining: 474
- Rejection percentage: 99.86362106220203
- Every baseline pair fails the difference-gap test; no full pattern survives.

## Sequential sieve, ascending primes

| Prime | Rejected here | Remaining | Standalone rejected |
|---|---:|---:|---:|
| 3 | 0 | 347561 | 0 |
| 5 | 0 | 347561 | 0 |
| 7 | 0 | 347561 | 0 |
| 11 | 0 | 347561 | 0 |
| 13 | 41649 | 305912 | 41649 |
| 17 | 76861 | 229051 | 88127 |
| 19 | 88822 | 140229 | 134345 |
| 23 | 51324 | 88905 | 130672 |
| 31 | 58111 | 30794 | 219414 |
| 43 | 20939 | 9855 | 237537 |
| 47 | 7600 | 2255 | 264971 |
| 67 | 1781 | 474 | 284311 |

## Measured filter-only timings

Three samples each, seconds; medians below. No guarantee beyond this observed run.

- bank_signature_bitsets_then_gaps: 0.473569
- difference_first_gap_baseline: 0.028578
- difference_then_mask_then_sum: 0.030260
- full_gap_baseline: 0.086815
- mask_then_gaps: 0.376488
- polynomial_then_gaps: 1.425101
- reduced_descending_mask_then_gaps: 0.203665

No tested sieve implementation beats difference-first exact gap lookup. Difference-first plus sieve reaches no sieve calls, so its small timing difference from the difference-first baseline is noise, not modular benefit. The modular sieve rejects 99.864% of the structured candidate set, whereas uniform CRT density is about 0.00207% survival; those populations are not interchangeable.

## Files

- PROOF.md: theorem, exact density, content conditions, completeness and scope.
- verify.py / certificate.json: exact replay including local controls and candidate counts.
- survivors.csv: all 474 passing pairs with actual opposite-entry roots, not claimed full solutions.
- benchmark.json: three actual samples per implementation.
- test_controls.py: direct-discriminant and pair-collision checks, upstream hash equality, tamper controls.
- SHA256SUMS: manifest.

Run `python verify.py --verify` and `python test_controls.py`. Input dependency is the Parker classification and exact bank certificate in sibling session55 directories; main verifier itself imports no producer code.
