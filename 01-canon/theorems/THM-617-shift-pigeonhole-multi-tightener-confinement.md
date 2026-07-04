---
id: THM-617
title: The shift-pigeonhole closes multi-tightener confinement for large scale — for S = mU ∪ {w_1,…,w_f} (m∤w_i, M(U)>1/14), if f·m/7 + Σ_i gcd(w_i,m) < m then M(S) > 1/14. Since mU is (+1/m)-periodic, M(U) is hit at all m shifts of U's argmax, each tightener is unsafe on ≤ m/7+gcd(w_i,m) of them, so the f tighteners cannot cover all m shifts — a fully-safe shift survives. Extends opus THM-616 (f=1, all m) to f≥2: closes f-tightener confinement for m > 7f/(7−f) (coprime tighteners: f=2⇒m≥3, f=3⇒m≥6, …), leaving only bounded small-m residuals (f=2: exactly m=2, the folding case).
status: PROVED (elementary shift-pigeonhole). VERIFIED: count bound 0/1394 violations; pigeonhole implication 0/1054 failures.
source: mac-mini-2026-07-04-S41
depends_on:
  - THM-616   # opus: f=1, all m (one tightener useless); this is the f≥2 companion
  - THM-612   # tight families decompose as mU ∪ F (Lemma A: binding pair divisible by m)
  - LRCUpTo13 # M(U) ≥ 1/(e+1) > 1/14 for the ≤12-runner even part U
related:
  - THM-615   # opus folding identity (m=2, f=2) — handles the m=2 residual this theorem leaves open
  - HYP-4080  # klein: 11-runner even-parts have a spectral gap (M(U) ≥ 2/23 for loose U)
results:
  - 04-computation/f2_pigeonhole_macmini_20260704.py
  - 05-knowledge/results/f2_pigeonhole_macmini_20260704.out
external: Lonely Runner Conjecture n=14.
---

# THM-617 — the shift-pigeonhole for multi-tightener confinement

**Setup.** A tight-candidate covering family at hiding denominator `q*=14m` decomposes (THM-612 Lemma A) as
`S = mU ∪ F`, where `E=mU` are the `m`-divisible runners (`U` an `e`-set, `e=13−f`), and `F={w_1,…,w_f}` the
`f` non-`m`-divisible **tighteners** (`m∤w_i`). By LRC(≤13), `M(U) ≥ 1/(e+1) > 1/14` (`U` is loose).

## The theorem
> **If `f·m/7 + Σ_{i=1}^{f} gcd(w_i,m) < m`, then `M(S) > 1/14`** (so `S` is not tight — confinement holds).

**Proof.** `g_{mU}(t)=min_{u∈U}‖mu·t‖=g_U(mt)`, so `M(mU)=M(U)`. Let `τ*` be an argmax of `U`
(`g_U(τ*)=M(U)`) and set `t_j=(τ*+j)/m`, `j=0,…,m−1`. Then `g_{mU}(t_j)=g_U(mt_j)=g_U(τ*+j)=g_U(τ*)=M(U)`
for **every** shift `j` — the `m`-divisible part is at its maximum `M(U)>1/14` on the whole shift-orbit.

Call shift `j` *unsafe for `w_i`* if `‖w_i t_j‖ ≤ 1/14`. Now `w_i t_j = w_i τ*/m + w_i j/m`, and
`{w_i j/m mod 1 : j}` is the subgroup of order `m/g_i` (`g_i=gcd(w_i,m)`), each value hit `g_i` times; the
danger arc `‖·‖≤1/14` has width `1/7`, so it contains `≤ ⌊(1/7)(m/g_i)⌋+1` of those values, whence
```
   #{unsafe shifts for w_i}  ≤  (⌊m/(7g_i)⌋ + 1)·g_i  ≤  m/7 + g_i.
```
Summing, the total number of shifts unsafe for *some* tightener is `≤ f·m/7 + Σ g_i < m` by hypothesis. So
some shift `j*` is safe for **all** tighteners: `‖w_i t_{j*}‖ > 1/14` for every `i`. There,
`min_{v∈S}‖v t_{j*}‖ = min(M(U),\ min_i‖w_i t_{j*}‖) > 1/14`, hence `M(S) ≥` that `> 1/14`. ∎

## Consequences
- **`f=1` (recovers/agrees with opus THM-616):** `m/7 + gcd(w,m) < m` holds for all `m≥2` (`gcd<m`,
  `m/7<6m/7`), so one tightener never makes `M≤1/14` — the single-tightener confinement, all `m`.
- **`f=2`, coprime tighteners (`g_i=1`):** `2m/7+2<m ⟺ m>14/5 ⟺ m≥3`. **Confinement for m≥3, f=2** — the
  only residual is `m=2` (opus's THM-615 folding / the argmax conspiracy on the 2-shift orbit).
- **General `f` (coprime):** closes for `m > 7f/(7−f)` (`f=3⇒m≥6`, `f=4⇒m≥10`, …, `f≤6`); tighteners
  sharing factors with `m` need `Σgcd(w_i,m) < m(7−f)/7`.
- **Net:** multi-tightener confinement reduces to a **bounded, per-`f` set of small-`m` cases** (the
  large-scale regime is closed by one pigeonhole). Combined with THM-616 (`f=1`) and klein's HYP-4080
  spectral gap, the residual is the small-`m` argmax problem — the folding endgame.

## Verified (`f2_pigeonhole_macmini_20260704.py`)
Count bound `#unsafe ≤ m/7+gcd`: **0/1394 violations**. Implication (condition ⟹ a fully-safe shift exists):
**0/1054 failures**. The `min(a_j,b_j)>1/14` shift survives at every `m≥3` sampled; `m=2` is where two
tighteners, unsafe at the two *different* shifts (`a_0+a_1=1/2`), can leave no both-safe shift — exactly the
folding residual.

## Not claimed
`M > 1/14`, not the sharp covering-min `M ≥ 14/183`. This closes the *confinement* (no tight covering family)
in the large-scale regime, not the full covering-min. The `m=2` (and other small-`m`) residuals stay open.
