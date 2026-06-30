# Dip upgrade: the inhomogeneous AP-LRC law M_c = 1/n + c(n−2)/n is EXACT (no O(1/n²) dips) — the block construction extends to RATIONAL q=a/b≥n via t=(a−b)/a, giving min_v‖vt−c‖ = (a−(n−2)b)/2a = env by an elementary boundary-min (achieving env on a DENSE set of c); M_c is 1-LIPSCHITZ in c, so dense achievability forces M_c ≥ env EVERYWHERE; the computational "dips" at c=odd/2n were Qmax<n² under-estimates (the optimal t there has denominator exactly n²); hence M_c = env exactly and the loneliness integral L = 1/4 + 1/(2n) EXACTLY (the O(1/n²) error term is removed)

*opus-2026-06-30. Owner: upgrade the dip bound. Done — and the dips turned out not to exist: they were
finite-Qmax artifacts. Two ingredients kill them: rational-q achievability (dense) + Lipschitz continuity.*

## The upgrade in one line
`M_c(AP_n) = 1/n + c·(n−2)/n` holds **exactly for every `c ∈ [0, 1/2]`** (not just the integer-`q` dense set),
and **`L = ∫₀¹ M_c dc = 1/4 + 1/(2n)` exactly** — the `O(1/n²)` I'd carried was the residue of a finite-Qmax
artifact, now removed.

## Ingredient 1 — achievability at ALL rational q (elementary, dense)
> For any rational `q = a/b ≥ n` (`gcd(a,b)=1`), with observer `c = (q−n)/(2q) = (a−nb)/(2a)`, take
> **`t = (a−b)/a = 1 − 1/q`**. Then `min_{v=1}^{n−1} ‖vt − c‖ = (q−n+2)/(2q) = env`.

**Proof.** `vt − c = v(a−b)/a − c ≡ −(vb/a + c) (mod 1)`, so `‖vt−c‖ = ‖vb/a + c‖ = ‖(a + wb)/(2a)‖` with
`w := 2v − n ∈ {−(n−2), −(n−4), …, n−2}` (step 2). The value `a + wb` lies in `[a−(n−2)b, a+(n−2)b]`, and
since `q = a/b ≥ n` we have `a ≥ nb`, so the nearest multiple of `2a` to `a+wb` is `0` or `2a` at distance
`a − |w|b`. This is minimized at the **boundary `|w| = n−2`**, giving
`min_v ‖vt−c‖ = (a − (n−2)b)/(2a) = (q−n+2)/(2q) = env`. ∎ *(verified n=8..20, many q.)*

The integer-`b=1` case is the old block (`t=(q−1)/q`); for `b>1` the runners are an AP of step `b/a` leaving
**one** large gap of width `(q−n+2)/q` centered on `c`. As `a/b` ranges over rationals `≥ n`, the achieving
`c = (a−nb)/2a` ranges **densely** over `[0, 1/2)`.

## Ingredient 2 — M_c is 1-Lipschitz
For fixed `t`, `min_v ‖vt − c‖` is `1`-Lipschitz in `c` (each `‖vt−c‖` is; a min of `1`-Lipschitz functions
is). A **supremum** of `1`-Lipschitz functions is `1`-Lipschitz, so `M_c = sup_t min_v ‖vt−c‖` is
`1`-Lipschitz in `c` (max observed `|ΔM|/|Δc| = 1.000`).

## The kill: dense achievability + Lipschitz ⇒ no dips
`M_c = env` on a dense set `D` (Ingredient 1) and both `M_c`, `env` continuous ⇒ for any `c`, pick `c_k ∈ D`
with `c_k → c`: `M_c ≥ M_{c_k} − |c−c_k| = env(c_k) − |c−c_k| → env(c)`. So **`M_c ≥ env(c)` everywhere.**
Combined with the upper bound `M_c ≤ env` (the clumping inequality, THM-591 §B), **`M_c = env` exactly.** The
"dips" (`M_c < env` seen numerically at `c = odd/2n`) were **Qmax < n² under-estimates**: e.g. at `c=1/28`
(n=14) the optimal `t = 183/196` has denominator `196 = n²`, so any `Qmax < 196` misses it (verified:
`M_{1/28}` rises `0.10023 (Q=3n) → 0.10198 (Q=12n) → 0.10204 = 5/49 = env (Q=n²)`).

## Consequence: L is exactly 1/4 + 1/(2n)
With `M_c = 1/n + c(n−2)/n` exact and `M_c = M_{1−c}`,
`L = 2∫₀^{1/2}[1/n + c(n−2)/n]dc = 1/n + (n−2)/(4n) = 1/4 + 1/(2n)` — **no error term.** The mean-observer
loneliness is exactly `1/4 + 1/(2n)`; the difficulty localizes at `c=0` (`M_0 = 1/n`).

## Status
- **PROVED (opus, rigorous & exact):** the lower bound `M_c ≥ env` everywhere — rational-`q` achievability
  (elementary boundary-min, dense) + `1`-Lipschitz continuity. The `O(1/n²)` dips are **eliminated** (artifacts
  of `Qmax < n²`; optimal `t` has denominator `n²` at `c=odd/2n`).
- **Exact law:** `M_c(AP_n) = 1/n + c(n−2)/n` for all `c∈[0,1/2]`; `L = 1/4 + 1/(2n)` exactly.
- **Remaining (the OTHER direction):** the upper bound `M_c ≤ env` is the clumping inequality (THM-591 §B,
  steps (1),(3) rigorous; the structural single-clump step verified). The dip upgrade closes the `≥` side
  completely; the `≤` side's structural step is the only piece short of a fully self-contained bit-proof.

THM-591 updated. Related: PROOF-…-block-construction (§A now extends to rational q), CORRECTION-…-exactly-linear
(the law), the-loneliness-integral-limit-…-stern-brocot (L now exact); OPEN-Q-108.
