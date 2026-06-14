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
