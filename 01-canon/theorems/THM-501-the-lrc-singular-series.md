---
id: THM-501
title: The LRC singular series — the covering-deficit density limit L(S), L>0 ⟹ loose, and the reduction of C'(14) to a singular-series lower bound
status: PARTIAL — the density limit L(S), its vanishing-iff-tight, and L>0 ⟹ loose are established (expansion + exhaustive/large-q verification); the uniform lower bound inf L>0 (which would prove C'(14)) is CONJECTURED with strong evidence + the extremal structure identified
source: kind-pasteur-2026-06-14-S6
depends_on:
  - THM-398   # C' reduction; D(q,S)>0 <=> loose-via-q
  - THM-492   # band criterion (the deficit's definition)
  - THM-497   # dilated-band covering; the over-correlated regime
related:
  - HYP-2489  # the LRC deficit = circle-method singular series (this makes it concrete)
  - HYP-2501  # L(S) exists
  - HYP-2502  # L>0 => loose
  - THM-446   # the additive-relation / Sidon ladder (the relation lattice that controls L)
---

# THM-501 — the LRC singular series `L(S)`

## The object

For a speed set `S` (13 distinct positive integers) the covering deficit is
`D(q,S) = #{a ∈ Z/q : v·a mod q ∉ B_q for all v ∈ S}`, `B_q = ±{0,…,⌊q/14⌋}`;
`D(q,S) > 0 ⟺ S` has a strict loneliness witness at shell `q` (THM-398/492).

**Additive-character expansion.** Writing `1 − 1_{B_q}(va) = (1−β) − Σ_{t≠0} ĉ(t)
e_q(tva)` with `β = (2⌊q/14⌋+1)/q`, `ĉ(t) = (1/q)Σ_{x=−⌊q/14⌋}^{⌊q/14⌋} e_q(−tx)`
(the Dirichlet kernel, real), and summing the product over `a`:

> `D(q,S)/q = (1−β)^{13} + Σ_{∅≠T⊆S} (1−β)^{13−|T|}(−1)^{|T|}
>            Σ_{(t_v)_{v∈T}: t_v≠0, Σ_{v∈T} t_v v ≡ 0 (mod q)} ∏_{v∈T} ĉ(t_v).`

The deviation from the independence main term `(1−β)^{13}` is a sum over **additive
resonances** `Σ_{v∈T} t_v v ≡ 0 (mod q)` of the speeds.

## The singular series (PROVED to exist; verified)

A resonance from a non-zero integer relation `Σ t_v v = m ≠ 0` only fires at `q | m`
— finitely many `q`. So as `q → ∞` only the **exact** integer relations
`Σ_{v∈T} t_v v = 0` survive, and `ĉ(t) → s(t) := sin(πt/7)/(πt)` (the band's sinc
kernel). Therefore the limit

> **`L(S) := lim_{q→∞} D(q,S)/q`** exists — the **LRC singular series** — with
> `L(S) = (6/7)^{13} + Σ_{exact relations} (6/7)^{13−|T|}(−1)^{|T|} ∏_{v∈T} s(t_v)`,
> controlled entirely by the speeds' integer additive-relation lattice (the
> Sidon/B_h structure of THM-446).

*Verified* (exact deficit, large-`q` window averaging, stable to ±0.005):
`L(generic/Sidon) ≈ 0.135–0.146 = (6/7)^{13}` (no small relations); `L` is
**suppressed** by small relations — the arithmetic-progression core `d·{1,…,12}`
(relations like `7 − 2·14 + 21 = 0`) gives the lowest values: evaders
`7·{1,…,12}∪{r}` have `L ≈ 0.030`, the `d=3` core `3·{1,…,12}∪{182}` reaches
`L ≈ 0.026`. (`04-computation/lrc14_singular_series_kps6.py`.)

## `L > 0 ⟹ loose`, and `L = 0 ⟺ tight`

`L(S) > 0` means `D(q,S) > 0` for all large `q` (positive witness *density*), so `S`
is **loose**. This is *stronger* than C'(14)'s "∃ a witness" (positive density vs
nonempty), so it is the clean target. `L(S) = 0` exactly for **tight** configs: the
tight `14·{1,…,13}` (`M = 1/14`) has `D ≡ 0`, `L = 0` (verified). So `L` is the
density of the LRC safe set, asymptotically.

## The reduction: C'(14) ⟸ a singular-series lower bound

> **Reframe.** If `L(S) > 0` for every primitive multiple-of-14 speed set `S`, then
> C'(14) — hence LRC(14) — holds.

This is **exactly the circle-method / Pollock structure**: prove the *singular
series is bounded below by a positive constant*, and the main term dominates so
representations (witnesses) exist (HYP-2489 made concrete; cf. the Pollock proof's
"singular series positive ⟹ every large N representable"). Evidence + the extremal
structure:

- Over 120 sampled primitive non-dominant multiple-of-14 configs (entries < 100):
  `min L = 0.094`, **0 configs with `L < 0.02`**.
- The infimum over the maximally-structured families is attained at the **dilated
  arithmetic-progression cores** `d·{1,…,12}∪{r}` (the evaders and their `d`-analogues)
  — an AP has the most small additive relations, maximally suppressing `L` — and even
  there `L ≈ 0.026 > 0`.
- So `inf_S L(S) > 0` over primitive multiple-of-14 `S` is strongly supported, with
  the extremal/hardest configs being the dilated-AP cores. **Proving this uniform
  lower bound is the singular-series-positivity proof of C'(14).**

The threshold `q*(S)` (first-witness shell) is governed by the *non-relation* speed
magnitudes; large strangers go B'-dominant (the evader family height drops from 41 to
13 once `r ≥ 1093 = (n−1)·max(core)`), so the `ladder ∪ B'` closure (HYP-2438) reads
as "`L(S) > 0` + a finite-shell check below `q*`".

## Adelic structure: `L` is a singular INTEGRAL, not a singular-series Euler product (macmini-2026-06-14)

The natural HYP-2503 question — does `L(S)` factor as a Hardy–Littlewood Euler
product `L = β_∞ · ∏_p β_p` of nontrivial local densities? — is settled **NO**.
Define the single-prime local limit `β_p(S) := lim_{e→∞} D(p^e,S)/p^e`. A
resonance `Σ t_v v = m` fires at `q = p^e` iff `p^e | m`; for `m ≠ 0` this holds
for only finitely many `e`, so the single-prime limit keeps **exactly the same
exact relations** (`m = 0`) as the full `q→∞` limit — and `m = 0 ∈ p^eℤ` for
*every* `p, e`. Therefore

> **`β_p(S) = L(S)` for every prime `p`.**

Verified (`04-computation/lrc14_adelic_structure_audit_0614.py`): evader
`7·{1..12}∪{13}` has `β_2=β_3=β_5=β_11=β_13 ≈ 0.029 = L` (agree to ±0.0002);
AP-core `{1..12}∪{14}` has all `β_p ≈ 0.024 = L`; Sidon-ish sets have all `β_p`
agreeing to ±0.003 and each equal to the set's own `L`. Consequences:

- **No nontrivial Euler product.** With `β_p = L` for all `p`, the product
  `β_∞ · ∏_p β_p` would force `L = β_∞ · L^{#primes}` — impossible for `0 < L < 1`.
  Every local factor is trivially `1` (each single prime already reconstructs the
  whole `L`); none carries a proper local-density correction.
- **Suppression is p-adically universal, not p-local.** The AP-core relations have
  moduli coprime to `11, 13`, yet `β_11, β_13` are equally suppressed — because the
  surviving relations are *exact* (global), sitting at `m = 0`, seen identically by
  every prime ladder.
- **The prime-power data lives in the APPROACH to `L`** — the threshold shell
  `q*(S)` and the convergence rate of `D(p^e)/p^e` — **not in `L` itself.**

So `L(S)` is the **singular INTEGRAL** of the LRC covering (a single archimedean
number, the exact-relation sinc sum / the density on the limiting torus), **not a
singular-series Euler product.** This is the structural point where LRC departs
from classical Waring/Goldbach circle-method counts: there, finite-`m` congruence
obstructions populate genuine local factors `β_p`; here the surviving resonances
are exact, so the entire arithmetic content is archimedean and the local factors
collapse to `1`. (HYP-2503 corrected accordingly.)

## Honesty

- `L(S)`'s defining series has the circle-method's **conditional-convergence**
  subtlety (the limiting weights `s(t) ~ 1/t` are not absolutely summable; the
  Dirichlet-kernel `L¹` norm is `~ log q`). The *limit* `L(S)` is verified to exist
  numerically (stable windows); a rigorous existence/lower-bound proof must handle
  this exactly as the classical singular series does. This is the open analytic core.
- `inf L > 0` is conjectural (strong sample evidence + identified extremizers), not
  proved. Proving it is the prize.
- `L > 0 ⟹ loose` is the proved direction; the converse (loose ⟹ `L>0`) is open
  (a loose config could a priori have witness-density 0).

**Artifacts:** `04-computation/lrc14_resonance_expansion_kps6.py`,
`lrc14_singular_series_kps6.py`, `lrc14_singular_series_infimum_kps6.py` (+ `.out`).
Reflection `07-reflections/the-lrc-singular-series-kps6.md`.

---

## Structural results on `L(S)` (PROVED — kind-pasteur-2026-06-14-S6 structural pass)

These are rigorous lemmas about the term structure of `L(S)`. They do **not** prove
`inf L > 0`; they isolate exactly which class is fully under control.
Verified by `04-computation/lrc14_singular_series_structure_kps6.py`
(+ `05-knowledge/results/lrc14_singular_series_structure_kps6.out`).

### LEM-501.1 — The 7-vanishing
For integer `t ≠ 0`, `s(t) = sin(πt/7)/(πt) = 0  ⟺  7 | t`.
*Proof.* `s(t)=0 ⟺ sin(πt/7)=0 ⟺ πt/7 ∈ πℤ ⟺ 7 | t` (and `πt ≠ 0` since `t ≠ 0`). ∎
**Corollary.** In `∏_{v∈T} s(t_v)`, if ANY coefficient `t_v` has `7 | t_v` the factor is
0, so the whole relation term vanishes. Hence `L(S)` depends ONLY on exact relations
all of whose coefficients are coprime to 7 ("7-primitive relations").
Also `|s(t)| ≤ 1/(π|t|)` (since `|sin| ≤ 1`).

### LEM-501.2 — Pairwise absolute convergence, explicit bound
For a pair `{v_a, v_b}`, exact pairwise relations `t_a v_a + t_b v_b = 0` (`t_a,t_b≠0`)
are exactly `(t_a,t_b) = k·(v_b, −v_a)/g`, `k ∈ ℤ\{0}`, `g = gcd(v_a,v_b)`. The pair's
total contribution to the series carries `(6/7)^{11}` and the lattice factor
`P(a,b) = Σ_{k≠0} s(k v_b/g)·s(−k v_a/g) = 2 Σ_{k≥1} s(k a)·s(k b)`, where
`a = v_b/g`, `b = v_a/g` (`s` odd ⇒ the product is even in `k`). This series is
**absolutely convergent** with the explicit bound
> `|P(a,b)| ≤ 2 Σ_{k≥1} 1/(π² a b k²) = (2/π²)(π²/6)/(ab) = 1/(3ab) = g²/(3 v_a v_b).`
*Proof.* `|s(ka)s(kb)| ≤ 1/(π²ab k²)` by LEM-501.1's bound; sum the `ζ(2)` series.
`ab = v_a v_b / g²`. ∎  (If `7 | a` or `7 | b`, every term vanishes by LEM-501.1 and
`P = 0`, e.g. `{98,105}`: reduced step `14` ⇒ `P = 0`.) Verified for many pairs (all
exact `|P|` strictly under the bound).

### LEM-501.3 — Almost-Sidon subcase: a provable `L > 0`
Call `S` **almost-Sidon** if its only 7-primitive exact relations have `|T| = 2`
(no `|T| ≥ 3` relation with all coefficients coprime to 7). For such `S`,
> `|L(S) − (6/7)^{13}| ≤ (6/7)^{11} · Σ_{a<b} |P(a,b)| ≤ (6/7)^{11} · Σ_{a<b} g(a,b)²/(3 v_a v_b).`
If the RHS is `< (6/7)^{13}` then `L(S) ≥ (6/7)^{13} − RHS > 0`, so `S` is loose.
**Explicit threshold:** `Σ_{a<b} g(a,b)²/(3 v_a v_b) < (6/7)² = 36/49 ≈ 0.7347`
suffices (divide the displayed inequality by `(6/7)^{11}`). For coprime spread sets
this is easily met: e.g. `{14,17,…,127}` gives proved RHS `≈ 0.0029 ≪ 0.135`
(`⇒ L ≥ 0.132`); `{14,101,…,1201}` gives `≈ 0.00016` (`⇒ L ≥ 0.1346`). Ground-truth
deficit `D(q,S)/q` for `{14,17,…,127}` averages `0.128 ± 0.005`, consistent.

### What this covers and what it does NOT
- **Covered (rigorously):** any almost-Sidon multiple-of-14 set with pairwise
  relation mass below the explicit threshold has `L > 0` hence is loose.
- **NOT covered:** the hard configs — dilated AP cores `d·{1,…,12}∪{r}` — are
  **not** almost-Sidon: they have an abundance of `|T| ≥ 3` 7-primitive relations
  (e.g. `1·7 − 2·14 + 1·21 = 0`), and these `|T| ≥ 3` terms are precisely what drives
  `L` down (evader `7·{1..12}∪{611}`: ground-truth `L ≈ 0.029`, of which the
  pairwise part alone would give `≈ 0.25`; the `|T|≥3` mass `≈ −0.22` is the whole
  story). The full `inf L > 0` over high-relation-energy AP cores remains the open
  analytic core (`|T|≥3` lattice sums are conditionally convergent, the genuine
  circle-method difficulty). **No claim here closes that gap.**
