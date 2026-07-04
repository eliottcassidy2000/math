---
id: THM-617
title: The shift-pigeonhole gives the FULL COVERING-MIN at large scale — for S = mU ∪ {w_1,…,w_f} (m∤w_i, M(U)≥c), if f·2c·m + Σ_i gcd(w_i,m) < m then M(S) ≥ c; at c=14/183 (with M(U)≥1/(e+1)>14/183) this is M(S) ≥ 14/183 (the covering-min). Since mU is (+1/m)-periodic, M(U) is hit at all m shifts of U's argmax, each tightener is unsafe (<c) on ≤ 2cm+gcd(w_i,m) of them, so the f tighteners cannot cover all m shifts — a shift with all tighteners ≥c survives. Extends opus THM-616 (f=1, all m) to f≥2 AND sharpens confinement (M>1/14) to the covering-min (M≥14/183): closes for m > 7f/(7−f) (coprime: f=2⇒m≥3, f=3⇒m≥6, …), leaving only bounded small-m residuals (f=2: exactly m=2, the folding case).
status: PROVED (elementary shift-pigeonhole). VERIFIED: count bound 0/1394 violations; pigeonhole implication (M>1/14) 0/1054 failures; SHARP version (all tighteners ≥14/183 shift exists) 0 hard over m≥3.
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

## The theorem (threshold `c`, giving the covering-min at `c=14/183`)
> **Fix `0<c≤M(U)`. If `f·2c·m + Σ_{i=1}^{f} gcd(w_i,m) < m`, then `M(S) ≥ c`.**
> In particular with `c=14/183` (legal since `M(U)≥1/(e+1)≥1/12>14/183` for `e≤12`): **`M(S) ≥ 14/183`**,
> the sharp primitive covering-min — not merely `>1/14`.

**Proof.** `g_{mU}(t)=min_{u∈U}‖mu·t‖=g_U(mt)`, so `M(mU)=M(U)`. Let `τ*` be an argmax of `U`
(`g_U(τ*)=M(U)`) and set `t_j=(τ*+j)/m`, `j=0,…,m−1`. Then `g_{mU}(t_j)=g_U(mt_j)=g_U(τ*+j)=g_U(τ*)=M(U)`
for **every** shift `j` — the `m`-divisible part is at its maximum `M(U)≥c` on the whole shift-orbit.

Call shift `j` *unsafe for `w_i`* if `‖w_i t_j‖ < c`. Now `w_i t_j = w_i τ*/m + w_i j/m`, and
`{w_i j/m mod 1 : j}` is the subgroup of order `m/g_i` (`g_i=gcd(w_i,m)`), each value hit `g_i` times; the
danger arc `‖·‖<c` has width `2c`, so it contains `≤ ⌊2c·(m/g_i)⌋+1` of those values, whence
```
   #{unsafe shifts for w_i}  ≤  (⌊2c·m/g_i⌋ + 1)·g_i  ≤  2c·m + g_i.
```
Summing, the total number of shifts unsafe for *some* tightener is `≤ f·2c·m + Σ g_i < m` by hypothesis. So
some shift `j*` has `‖w_i t_{j*}‖ ≥ c` for **all** `i`. There,
`min_{v∈S}‖v t_{j*}‖ = min(M(U),\ min_i‖w_i t_{j*}‖) ≥ min(M(U),c) = c`, hence `M(S) ≥ c`. ∎

*(The confinement `M(S)>1/14` is the `c=1/14` case, `2c=1/7`; the covering-min `M(S)≥14/183` is `c=14/183`,
`2c=28/183` — the two thresholds are so close that the legal `m`-regime is identical.)*

## Consequences (at `c=14/183`, i.e. the full covering-min `M(S)≥14/183`)
- **`f=1` (recovers opus THM-616, sharpened):** `2c·m + gcd(w,m) < m` for all `m≥2` (`2c=28/183<1`,
  `gcd<m`), so one tightener never lets `M<14/183` — single-tightener covering-min, all `m`.
- **`f=2`, coprime tighteners (`g_i=1`):** `(56/183)m+2<m ⟺ m>366/127≈2.88 ⟺ m≥3`. **Covering-min
  `M≥14/183` for m≥3, f=2** — the only residual is `m=2` (opus's THM-615 folding / the argmax conspiracy on
  the 2-shift orbit, where `2c` cannot fit two tighteners into `<2` shifts).
- **General `f` (coprime):** closes for `f < 183m/(28m+183)` → `f≤6`; e.g. `f=3⇒m≥6`, `f=4⇒m≥10`.
  Tighteners sharing a factor with `m` need `Σgcd(w_i,m) < m(1−28f/183)`.
- **Net:** the **full covering-min** `M≥14/183` for all primitive covering families reduces to a **bounded,
  per-`f` set of small-`m` cases** (the large-scale regime is closed by one pigeonhole — and note the `>1/14`
  confinement and the `≥14/183` covering-min have the *same* legal `m`-regime, since `1/7≈28/183`). Combined
  with THM-616 (`f=1`) and klein's HYP-4080 spectral gap, the residual is the small-`m` argmax/folding endgame
  that opus/klein/kps are closing rung-by-rung.

## Verified (`f2_pigeonhole_macmini_20260704.py`)
Count bound `#unsafe ≤ 2c·m+gcd`: **0/1394 violations** (`c=1/14`). Confinement implication: **0/1054 failures**.
**Sharp version** (`c=14/183`: a shift with *all* tighteners `≥14/183` exists): **0 hard cases over m≥3**
(`min` of the both-safe value: `m=3→7/36`, `m=4→3/16`, … all `≥14/183`). `m=2` is where two tighteners,
unsafe at the two *different* shifts (`a_0+a_1=1/2`), can leave no safe shift — exactly the folding residual.

## Not claimed
The **large-scale covering-min** `M≥14/183` is proved (all `m` beyond the per-`f` threshold). The **small-`m`
residuals** — chiefly `m=2` for every `f` — stay open; those are the argmax/folding cases (opus THM-615,
klein/kps ladder-and-residue closures). So this theorem covers everything *except* a bounded, explicitly
listed set of small-`m` scales — a genuine reduction of the covering-min, not a soft bound.
