# Pre-registration (S44): transfer orders of the degree-N generators of Vir_(6,7) and Vir_(5,7)

Written and committed BEFORE running `vir_transfer_orders_20261018.py` on these models.

**Hypothesis.** HYP-9315(b): in degree N, the minimal generator of grI at Li degree 2k+1 is a transfer of order exactly
k (k = 1, ..., p-2). The generator weights are from kind-pasteur's table (singular_support... note, section 3.2).

**Predictions.**
- P1. Vir_(6,7), N = 15. The degree-15 generators at weights 33, 35, 37, 39 (Li 3, 5, 7, 9) are transfers of orders
  1, 2, 3, 4 respectively: the first appearance of order k at length 15 is at weight 31 + 2k. No order >= 5 appears at
  length 15 through weight 40.
- P2. Vir_(5,7), N = 12. The degree-12 generators at weights 27, 29, 31 (Li 3, 5, 7) are transfers of orders 1, 2, 3:
  the first appearance of order k at length 12 is at weight 25 + 2k.
- P3. At both models the recorded orders sum to dim ker pi at every bidegree (no mismatch).

**Runs.** `vir_transfer_orders_20261018.py 6 7 40` and `... 5 7 34` (single prime 2^31 - 1).

## Outcome (recorded after the runs): ALL PASSED

- The pre-registration was committed in 94fbe3b30b at 13:37:51. The outputs `orders_57_D34.out` and `orders_67_D40.out`
  were written at 13:37:58 and 13:38:05.
- **P1 PASSED.** In Vir_(6,7), N = 15, the first appearances at length 15 are: order 1 at weight 33 (Li 3), order 2 at
  35 (Li 5), order 3 at 37 (Li 7), order 4 at 39 (Li 9). There is no order >= 5 at length 15 through weight 40.
- **P2 PASSED.** In Vir_(5,7), N = 12, the first appearances at length 12 are: order 1 at 27, order 2 at 29, order 3 at 31.
- **P3 PASSED.** There are 0 MISMATCH lines. The totals are {1: 24, 2: 14, 3: 5} = 43 = dim ker π for (5,7) through 34,
  and {1: 24, 2: 14, 3: 5, 4: 2} = 45 = dim ker π for (6,7) through 40.

**Audit note (N4).** All three are satisfied and were re-verified at two primes. Not every clause was informative:
- P3 is automatic for a correct implementation.
- The "order <= k" halves and "no order >= 5 through weight 40" are forced by weight counting.
- k = 2 is forced by HYP-9315(a) + (d).

The informative content was k = 3, 4.
