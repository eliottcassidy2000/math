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
