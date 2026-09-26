---
id: HYP-9141
title: "Is 5n+1 in the Haar closure of the bounded-lookahead-provable sign strategies? By THM-4481 this happens only if provable strategies concentrate their invariant measure on their flips; and is class (i) ever nonempty for 7 <= q <= 21?"
status: >
  OPEN QUESTION (both directions live). PROVED (THM-4481): the Haar
  distance of any class-(i) 5n+-1 strategy is > 0.01391 / max_(R u R*) f,
  where f = dpi/dU is its invariant density. It is also at least
  N_k(5)/2^(k-1) ~ 2/k (necklaces, THM-4479's Theorem 2, which does not
  depend on the multiplier).
  FINITE-EXACT: delta_7(5) = 29 and 44 <= delta_8(5) <= 50. VERIFIED:
  certified sets with Haar fraction 3.09/k to 3.41/k for 10 <= k <= 19,
  whose fits cannot separate 1/k (rising constant), k^(-0.84), and
  0.0355 + 2.73/k. The same sets show no concentration (max f <= 2.4 for
  k <= 14).
  Second question: for 7n+1, 9n+1 and 11n+1 the least rho_max is 3/7, above
  their critical densities, through k = 9, 8, 8. Does rho*(q,k) ever drop
  below log_q 2 for some 7 <= q <= 21?
  UPDATE 2026-09-26 (THM-4486): the "3/7 floor" is REFUTED. For 5n+1,
  rho*(5,k) = 2/5 exactly for k >= 15. For 7n+-1 the certified values
  reach 14/37 = 0.3784 at k = 22, still above log_7 2 = 0.3562; the
  proved all-level floor is 1/3 (Theorem N). For q = 9..21, no provable
  strategy exists at k <= 18. The second question stays OPEN, and
  entropy-type arguments provably cannot settle q <= 7 (cap 1/3).
source: collatz-procgen-20260922 session, cube-distance lane (DRIFT, 2026-09-25) and drift lane (2026-09-26)
depends_on:
  - 01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md
---

# HYP-9141 -- 5n+1 and the provable sign strategies

**Question A.** Does `delta_k(5)/2^(k-1) -> 0`?
* **A "yes"** needs a *funnel*: provable strategies whose orbits revisit a Haar-vanishing flip set at a positive stationary frequency. This follows from THM-4481's constant flip mass `pi(R ∪ R*) > 0.01391`. The obvious funnels are too expensive, because forcing an orbit into a flip `j` steps later requires flipping all `2^(j+1)` completions.
* **A "no"** would follow from a uniform bound on invariant densities of provable strategies on their flips.

**Question B.** Is class (i) ever nonempty for `q = 7, 9, ..., 21`? THM-4481 empties it for `q >= 23` at every level. The computed min-max densities stay at `3/7` from `k = 7`, observed to be the same for `q = 5, 7, 9, 11` (not explained).

**Why it matters.** The strategy-cube picture is:
* exponential distance for negative drift (`q = 3`, THM-4479);
* impossibility for large drift (`q >= 23`, THM-4481);
* the middle band of positive drift, which is where "the drift changes the rate but not the limit" could first fail.

For arbitrary (height-aware) edits, THM-4480 shows the price tends to 0 for every `q`. So a "no" here would make periodicity, the inability to see height, the precise obstruction.
