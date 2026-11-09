# Pre-registration addendum: the connectivity floors at order 11

kind-pasteur-2026-11-08-S17. Frozen once pushed. It follows [PREREGISTRATION.md](PREREGISTRATION.md) (af0112b596).

## Data seen since the first registration (orders <= 10, exploratory)

**The order-10 generation census.**
- R6 and R9 hold. R7 and R8 fail: 344 order-10 classes lie outside FC, EXT, SUM, SUB and the pendant extensions.
- 30,649 of the 30,869 REST classes have a vertex whose deletion is not strong. Such a class is an ordinal sum plus one
  vertex.
- The other 220 have arc-strong connectivity 2.

**The connectivity floors** (`lam.c`). Here lambda(T) is the arc-strong connectivity: the fewest arcs whose reversal
makes T non-strong. f_k(n) is the least H over n-tournaments with lambda >= k.

| n | f_1 (Moon-Busch) | f_2 | f_3 | f_4 | max H(n - 1) |
| --- | --- | --- | --- | --- | --- |
| 5 | 9 | 15 | | | 5 |
| 6 | 15 | 41 | | | 15 |
| 7 | 25 | 99 | 171 | | 45 |
| 8 | 45 | 195 | 621 | | 189 |
| 9 | 75 | 385 | 2,055 | 3,159 | 661 |
| 10 | 125 | 727 | 5,753 | 14,877 | 3,357 |

- At orders 4-10, every class with lambda >= 3 is a first creator, because f_3(n) > max H(n - 1).
- Non-first-creators with lambda = 2 first occur at order 9: 124 classes at order 9 and 40,401 at order 10.
- f_2(9) = 385 is the H of the exceptional order-9 class of the first registration.

## Predictions for order 11

| ID | Prediction | Expected |
| --- | --- | --- |
| R10 | f_3(11) <= max H(10) = 15,745. Some 11-tournament with lambda >= 3 has H <= 15,745. | CONFIRMED |
| R11 | The implication "lambda >= 3 implies first creator" fails at order 11: some 11-tournament with lambda >= 3 has H in S(10). | CONFIRMED |

The basis is the ratios f_3(n + 1)/f_3(n) = 3.63, 3.31, 2.80, which suggest f_3(11) is about 13,000-15,000.
f_3(11)/max H(10) would then fall below 1.

**Test.**
- Search: extend every order-10 class with lambda >= 2 by every pattern whose new vertex has in- and out-degree at
  least 3. Compute H through the attachment form. Test lambda >= 3 directly whenever H <= 15,745, and keep the least H
  and every H in S(10).
- A witness confirms each prediction. If no witness is found, the search reports its exact scope; it does not compute
  f_3(11) exactly.
