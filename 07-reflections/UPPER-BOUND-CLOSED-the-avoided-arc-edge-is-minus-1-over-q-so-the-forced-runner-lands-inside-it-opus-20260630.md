# The inhomogeneous AP-LRC upper bound is now FULLY RIGOROUS and ELEMENTARY: the avoided arc (c−env, c+env) has LEFT EDGE exactly at −1/q (because c−env = −(1−2c)/n = −1/q, q=n/(1−2c)), so the min-gap pigeonhole — which forces a runner j≤n−2 with ‖jt‖<1/q whenever the n−1 runners are squeezed into the length-(n−2)/q complementary arc — puts that runner INSIDE the avoided arc, an immediate contradiction; this replaces the clumping inequality's structural step, closing THM-591 in BOTH directions with elementary arguments (and it recovers the homogeneous M_0≤1/n as the c=0 case)

*opus-2026-06-30. Owner: work the remaining part of the LRC proof, make it rigorous; search old repo work.
Done — the search surfaced the q-witness (THM-523), the first-window (THM-405), and far_iff_fract (THM-386),
and combining them with my nearest-neighbor exploration gave a two-line upper bound. The whole law is proven.*

## The theorem
For `S = {1,…,n−1}`, every `c ∈ [0,1/2]`, every `t`:
> **`min_{v=1}^{n−1} ‖vt − c‖ ≤ env := 1/n + c(n−2)/n`.**
Together with the dip-upgrade lower bound (`M_c ≥ env`, rational-`q` achievability + 1-Lipschitz), this gives
`M_c(AP_n) = 1/n + c(n−2)/n` **exactly**, with BOTH directions elementary and rigorous.

## The one identity that does all the work
Let `q := n/(1−2c)` (so `q ≥ n`, `1/q = (1−2c)/n`). Then
> **`c − env = c − 1/n − c(n−2)/n = (2c−1)/n = −1/q`**  and  **`1 − 2·env = (n−2)/q`**.
So the observer's **avoided arc** `[c−env, c+env]` has **left endpoint exactly `−1/q`**; reduced mod 1 it is
`[0, c+env] ∪ [1−1/q, 1)`, which **contains `(−1/q, 1/q)`** (as `1/q ≤ c+env`, i.e. `0 ≤ 2cn`). Verified
`c−env = −1/q` for n=9,14,20, all tested c.

## The proof (two lines)
Suppose `μ := min_v ‖vt − c‖ > env`. Then no runner lies in `[c−env, c+env]`, so the `n−1` runners lie in the
complementary **open** arc of length `1 − 2·env = (n−2)/q`. Their `n−2` consecutive gaps sum to *less than*
`(n−2)/q`, so the **smallest gap is `< 1/q`**; being the arc-distance between two consecutive runners `vt, v't`
it equals `‖(v−v')t‖ = ‖jt‖` with `j = |v−v'| ∈ {1,…,n−2}`. Hence `‖jt‖ < 1/q`, i.e. `jt mod 1 ∈ (−1/q, 1/q)`
— which lies **inside the avoided arc**. So the runner `jt` satisfies `‖jt − c‖ ≤ env`, contradicting
`μ > env`. ∎

*(The `c=0` case recovers the classical homogeneous tightness `M_0 = 1/n` (THM-405): the avoided arc becomes
`[−1/n, 1/n]`, and the forced runner with `‖jt‖ < 1/n` lands in it. So the AP is tight and this is its
optimality proof.)*

## Why it's elementary (no three-gap theorem)
The clumping inequality (THM-591 §B, old) needed a structural step — "clumping ⇒ gap `≤1/(2j)` centered at
`c≈(2m+1)/2j`" — with an `O(1/n²)` cluster-width slop. **That step is gone.** The new argument uses only:
- the **min-gap pigeonhole** (`n−1` points in a length-`(n−2)/q` arc ⇒ a gap `< 1/q`) — Dirichlet, exact;
- the **identity `c−env = −1/q`** pinning the avoided-arc edge to the forced runner.
No three-distance/Slater structure, no counts, no `Qmax`. It is the **dual of the q-witness** (THM-523: "no
multiple of `q` ⇒ `M ≥ 1/q`"); here the *forced* near-`0` runner at `‖jt‖<1/q` is exactly what the observer's
`1/q`-aligned danger arc catches.

## Provenance (old repo work that made it click)
- **THM-405** (bounded-ratio first window): the AP `{1..n−1}` is the tight extremal `b=(n−1)a`, window
  collapses to `t=1/n`. My upper bound is its inhomogeneous completion and reproves `M_0≤1/n`.
- **THM-523 / THM-360** (q-witness / divisibility): "no multiple of `q` ⇒ gap `≥1/q`." The dual gives my
  contradiction (the forced runner IS a near-multiple).
- **THM-386** (`far_iff_fract`, Lean-checked): `‖x−c‖ ≥ μ ⟺ frac(x) ∈ [c+μ, c+1−μ]` — the avoided-arc
  characterization I used.
- The nearest-neighbor / three-gap rigidity work (HYP-3741/3717) pointed at the `0`-neighborhood structure;
  the final proof needed only the min-gap edge, not the full three-gap counts.

## Status
- **PROVED (opus, elementary, both directions):** `M_c(AP_n) = 1/n + c(n−2)/n` for all `c∈[0,1/2]`. Lower
  bound = rational-`q` block + 1-Lipschitz (dip upgrade); **upper bound = the avoided-arc-edge argument
  above** (replaces the clumping structural step). `L = 1/4 + 1/(2n)` exactly.
- **THM-591 §B updated** to this proof; status upgraded to fully PROVED (no structural gap remaining).
- **Scope (honest):** this completely solves the inhomogeneous LRC for the AP tight locus; it is NOT LRC(14)
  (OPEN-Q-108 bounds the lonely MEASURE under SPEED perturbation — a different quantity). But the
  avoided-arc-edge trick is a clean new optimality tool that may transfer.

Related: THM-591 (the law), DIP-UPGRADE-… (the lower bound), PROOF-…-block-construction (the old §B),
THM-405/523/360/386 (the tools); OPEN-Q-108.
