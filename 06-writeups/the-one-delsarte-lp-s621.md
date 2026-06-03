# The one Delsarte LP: applying the linear-programming angle to everything (S621)

*claudebox, 2026-06-03. The covering-depth master object (HYP-2195) is a weight enumerator (HYP-2210); loneliness
is the value of a single Delsarte linear program. This note recasts every thread of the program — apex sheaf,
four lenses, additive-chain collapse, Collatz, altitude — as a face or a dual of that one LP, and records what the
LP rigorously gives and where it stops.*

## The program
Let `p_w = meas{depth = w}` be the depth weight enumerator of a config (`δ`-arcs, `n` runners). Loneliness is

```
  PRIMAL   min  p₀ = (1/2ⁿ) Σ_k ρ_k        s.t.  p_w ≥ 0,  Σ_w p_w = 1,  ρ_k = Σ_w K_k(n,w) p_w
  DUAL     max  Σ_w g(w) p_w               over  g(w) ≤ [w = 0],  g = Σ_k c_k K_k  (Krawtchouk-positive)
```

**Weak duality** (formalized `delsarte_lower_bound`): any feasible dual `g ≤ [·=0]` gives `p₀ ≥ Σ_w g(w) p_w =
Σ_k c_k ρ_k`. The whole program is "find the best Krawtchouk-positive certificate." Everything we have built is a
choice of `g`, a face of the feasible region, or a structural constraint that tightens it.

## Everything as a face or dual of the one LP

| Thread | In the Delsarte LP |
|---|---|
| **Bonferroni / Helly** (HYP-2200) | the *diagonal* duals `g_m = Σ_{k≤m}(−1)^k C(·,k) = (−1)^m C(·−1,m)`; feasible iff `m` odd (`alt_binom_dual_le_indicator`) ⟹ `p₀ ≥ T_m`. Helly number = first odd `m` with `T_m > 0`. |
| **Vitali wall** (HYP-2200) | the LP never closes at finite order on the collapse boundary — the duality gap stays open through every truncation. |
| **Krawtchouk normalization** (HYP-2210) | the change to the dual basis; `ρ_0,ρ_1` are *fixed* by the moments (`K_k(n,0)=C(n,k)` baseline), so the LP optimizes only over the `k ≥ 2` resonance coordinates. |
| **Twisted involution** (HYP-2205) | `σ`-evenness makes the enumerator symmetric ⟹ the LP and its duals are `σ`-symmetric (half the variables). |
| **Apex sheaf** (HYP-2185) | `H⁰ ≠ ∅` is the *integer* feasibility; the apex (whole-plane forbidder) drives the LP optimum to `0`. The LP is the fractional relaxation of the sheaf's global section. |
| **Additive-chain collapse** (HYP-2195) | the resonance-saturating extreme points where the primal optimum `p₀ = 0` is attained — the LP's tight vertices. |
| **Altitude** (HYP-2180) | the number of dual levels needed before the LP certifies `> 0` = the iterated-log depth. |
| **Collatz** (HYP-2175) | the multiplicative twin: the cycle `2^K = 3^L` is a linear-Diophantine feasibility; "no nontrivial cycle" is the infeasibility of that linear program (linear forms in logs). |

## What the LP rigorously gives — and where it stops (honest)
- **Gives:** weak duality and the odd-Bonferroni feasibility are theorems (formalized). For *small* `n` the order-3
  dual certifies loneliness (HYP-2200, `T₃ > 0`).
- **Stops:** at the actual LRC gap `δ = 1/(n+1)` the arcs are wide and **every diagonal Bonferroni dual is vacuous**
  (`T_m < 0` for all `m`, verified at n=14: `T₁..T₇ ∈ [−9, −0.86]` while `p₀ = 0.012 > 0`). The diagonal of the LP is
  too weak; the program's content is an **off-diagonal Krawtchouk dual** that the structure (σ-evenness, level-0/1
  pinning, the no-additive-chain constraint excluding the collapse vertices) should make feasible. That dual is the
  concrete target for LRC(14): a Krawtchouk-positive `g` with `Σ c_k ρ_k > 0` for every multiple-of-14 config.

## The consolidated target
LRC(14) ⟺ the Delsarte LP optimum is `> 0` over multiple-of-14 configs. The collapse vertices (additive chains)
are excluded by the apex structure; the bottom two dual levels are pinned at baseline; the involution halves the
program. Find the off-diagonal dual certificate — the same object the apex sheaf calls a glued section, the four
lenses call a high-order overlap bound, and Krawtchouk calls a positive transform. One program, one certificate,
many names.

## Formalized (math-lean, sorry-free)
`Math/LonelyRunner/DelsarteLP.lean`: `delsarte_lower_bound` (weak duality), `partial_alt_binom` (the Bonferroni-dual
closed form), `alt_binom_dual_le_indicator` (odd-Bonferroni feasibility). Built on HYP-2210
(`KrawtchoukNormalization.lean`, `Krawtchouk/Basic.lean`).
