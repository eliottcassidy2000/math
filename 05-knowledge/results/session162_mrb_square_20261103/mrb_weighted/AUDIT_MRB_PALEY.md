# Independent audit of the MRB wrap and Paley certificate

Audited `RESULT.md`. Result: the specified claims pass. The code `audit_mrb_paley.py` uses subset/end-vertex dynamic programming for all eight Paley flip states, independently of the original permutation-based script, and exact rational arithmetic for the weighted wraps. Floating-point checks of the MRB finite identities are only numerical controls, not the convergence proof.

## Wrap proof and exceptional path

In A->b->c->A, a Hamilton path containing b not last must continue immediately to c. If b is last, c must be first, since its only possible predecessor is b. Hence paths are precisely `(A-path on S), b, c, (A-path on A\S)` and `c, (A Hamilton path), b`. With A transitive the two A paths in each split are uniquely ordered, so the weighted sum is `prod(1+w_i)+prod(w_i)`; the second summand is essential. More generally, the first family counts ordered two-path covers, allowing empties, and the second counts P(A). For a transitive composition H1->H2, restricting either of the two ordered paths to each block gives a bijection to the product of blockwise two-path covers, so Q2 is multiplicative.

For `w_i=(-1)^i (i^(1/i)-1)`, one has `|w_i|<1/2`: `(3/2)^i>i` for i>=2, by induction. Thus `1+w_i>0` and the real logarithm is well-defined. The inequality `|w-log(1+w)|<=w^2` holds for |w|<=1/2, while `w_i=O(log i/i)` and `sum_i w_i^2` converges. Also `sum_i w_i` converges by eventual monotonicity of `i^(1/i)-1` and the alternating-series test. Therefore the sum of logarithms converges, the product of `1+w_i` has a positive finite limit, and the exceptional product tends to zero geometrically. The finite logarithmic identity and its limit follow; these arguments need no interchange of non-absolutely-convergent sums.

## Exact DP output

For the Paley orientation on F_7 with positive differences {1,2,4}, independently counting path endpoints by `D[S,v]=sum_{u in S\{v}:u->v} D[S\{v},u]` gives state totals for bit order `(1->2,2->3,5->6)`:

    state 0..7: 189, 159, 159, 159, 159, 141, 141, 145.

Each equals `189-30(x+y+z)+30xy+12xz+12yz-8xyz`. In state 7 the out-neighbor masks are `[22,40,82,53,98,5,43]`, endpoint counts `[21,24,20,17,18,28,17]`, and the directed spanning cycle `0,1,3,2,4,6,5,0` verifies strong connectivity. The wrap count `2^n+1` and rational weighted identity were checked for n=1,...,6. No mathematical discrepancy identified.
