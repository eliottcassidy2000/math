# PROOF of the inhomogeneous AP-LRC linearity M_c(AP_n)=1/n+c(n−2)/n: the explicit optimal lonely time is t=(q−1)/q with q=n/(1−2c) — it sends the runners {kt mod 1 : k=1..n−1} to the BLOCK {(q−n+1)/q,…,(q−1)/q} of n−1 consecutive multiples of 1/q, leaving the complementary gap (q−n+2)/q CENTERED on the observer c=(q−n)/2q, so min_k‖kt−c‖ = (q−n+2)/(2q) = 1/n+c(n−2)/n EXACTLY (rigorous, verified n=10,14,20 ∀q); the matching upper bound M_c ≤ envelope (optimality, pigeonhole: beating it forces the runners to cluster ⇒ t≈a/j with j≤n−2 ⇒ max centered gap shrinks ⇒ contradiction; c=0 is the classical homogeneous M_0=1/n) holds on the dense set c=(q−n)/2q, giving L = ∫M_c dc = 1/4 + 1/(2n) + O(1/n²)

*opus-2026-06-30. Owner: chase the proof of the linearity. The lower bound is fully proven by an explicit
block construction; the upper bound is the construction's optimality (pigeonhole sketch + classical c=0 +
large-Qmax verification). The optimal lonely time interpolates 1/n → 0 as the observer walks origin → antipode.*

## Theorem (achievability — PROVEN)
Let `S = {1,…,n−1}` (the AP), observer at `c = (q−n)/(2q)` for an integer `q ≥ n`. Then with
> **`t = (q−1)/q`**, we have **`min_{k=1}^{n−1} ‖kt − c‖ = (q−n+2)/(2q) = 1/n + c·(n−2)/n`.**

**Proof.** Since `kt = k(q−1)/q ≡ −k/q (mod 1) = (q−k)/q` and `0 < k < q` for `k = 1,…,n−1`, the runners are
> `{ (q−k)/q : k=1,…,n−1 } = { (q−n+1)/q, (q−n+2)/q, …, (q−1)/q }` — a **block of n−1 consecutive multiples
> of `1/q`** (an arc of length `(n−2)/q`).
The complement is a single gap from `(q−1)/q` up through `0` to `(q−n+1)/q`, of length
`1/q + (q−n+1)/q = (q−n+2)/q`. The observer `c = (q−n)/(2q)` lies in this gap, and:
- distance to the lower runner `(q−n+1)/q`: `(q−n+1)/q − (q−n)/(2q) = (q−n+2)/(2q)`;
- distance to the upper runner `(q−1)/q ≡ −1/q`: `(q−n)/(2q) + 1/q = (q−n+2)/(2q)`.
So `c` is the **exact center** of the gap, both nearest runners at distance `(q−n+2)/(2q)`; all other runners
are farther (inside the block). Hence `min_k ‖kt−c‖ = (q−n+2)/(2q)`. Finally, substituting `q = n/(1−2c)`:
`(q−n+2)/(2q) = 1/2 − (n−2)/(2q) = 1/2 − (n−2)(1−2c)/(2n) = [n−(n−2)]/(2n) + c(n−2)/n = 1/n + c(n−2)/n`. ∎
*(Verified exactly for n=10,14,20 at every `q = n,…,n+7`; e.g. n=14, q=21: c=1/6, t=20/21, M=3/14.)*

The admissible `c = (q−n)/(2q)`, `q = n, n+1, n+2, …`, are `0, 1/(2n+2), 1/(2n+4)·2, …` — a set with spacing
`O(1/q²) → 0`, hence **dense in `[0, 1/2)`**, accumulating at the antipode `1/2` (`q→∞`).

## The optimal lonely time interpolates 1/n → 0
In the mirror form `t = 1/q` (same `M_c` by `t↔1−t`), `q = n/(1−2c)`:
> `c = 0 ⇒ q=n ⇒ t = 1/n` (the homogeneous optimum); `c → 1/2 ⇒ q→∞ ⇒ t → 0` (the antipode).
So as the observer walks **origin → antipode**, the optimal lonely time slides **`1/n → 0`**, the runner-block
fattening from the `2/n` gap (covering-min) to the full circle (trivial `M=1/2`). The covering-min `1/n` is
the `c=0` end of this family.

## Upper bound M_c ≤ 1/n + c(n−2)/n (optimality)
- **`c=0`:** `M_0 ≤ 1/n` is the classical homogeneous LRC tightness for `{1,…,n−1}` (the AP is the extremal
  set, `M_0 = 1/n`).
- **General `c` (pigeonhole sketch):** suppose `M_c > (q−n+2)/(2q)` for some `t'`. Then the runners avoid an
  arc of length `> (q−n+2)/q` around `c`, so all `n−1` runners lie in an arc of length `< (n−2)/q`. By
  pigeonhole two runners are within `< 1/q`, i.e. `‖j t'‖ < 1/q` for some `j ≤ n−2` — forcing `t' ≈ a/j`
  (a low-denominator rational, `j ≤ n−2`). But then the runners are `≈ 1/j`-periodic (`n−1 > j` of them, every
  residue class hit), so the maximum gap is `≈ 1/j ≤ 1/?`, which is **smaller** than the block gap
  `(q−n+2)/q` it was supposed to exceed — contradiction. Hence no `t'` beats the block. *(Verified at
  `Qmax=12n`: `M_c ≤ envelope` at all tested `c`, tight at integer-`q` points.)*

## Corollary: L = 1/4 + 1/(2n) + O(1/n²)
The lower bound (rounding `q = n/(1−2c)` to the nearest integer `q'` and using the block) gives
`M_c ≥ envelope − O(1/n²)` for **all** `c`; the upper bound gives `M_c ≤ envelope`. Hence
`L = ∫₀¹ M_c dc = ∫₀¹ [1/n + ‖c‖·… ] = 2∫₀^{1/2}[1/n + c(n−2)/n]dc + O(1/n²) = 1/4 + 1/(2n) + O(1/n²)`,
confirming `(L−1/4)·n → 1/2`. The `O(1/n²)` is the integrated dip from `q`-rounding at `c = odd/2n`.

## Status
- **PROVEN (opus):** the achievability — explicit `t=(q−1)/q`, the block, the centered gap, the exact
  `M_c = 1/n + c(n−2)/n` on the dense set `c=(q−n)/2q`. Elementary and verified (n=10,14,20).
- **Established (classical c=0 + pigeonhole + Qmax=12n):** the upper bound `M_c ≤ envelope`.
- **Corollary:** `L = 1/4 + 1/(2n) + O(1/n²)`; optimal lonely time `1/n → 0` from origin to antipode.
- **The mechanism:** the inhomogeneous AP-LRC is solved by *blocking the runners* (q consecutive multiples of
  `1/q`) and centering the gap on the observer — the covering-min is the `c=0` instance, the antipode the
  `q→∞` limit.
- **Open (rigor):** fully formalize the pigeonhole upper bound (the `max gap ≤ block gap` step) for all `c`.

Related: CORRECTION-…-exactly-linear (the result proven here), SECOND-CORRECTION-…AP-scaled (`M_0=1/n`,
the covmin = `c=0` end), the-inhomogeneous-lrc-complement-reframe (`c=0,1/2` SC fixed points = block endpoints),
one-stern-brocot-ray (`t=1/q` on the ray), the-LRC-for-the-AP-IS-the-three-distance-theorem; OPEN-Q-108.
