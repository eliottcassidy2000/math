---
id: THM-1755
title: Corrected trinomial tunability dichotomy; bounded-resonance and single-family claims refuted
status: PROVED singleton/tunable dichotomy with exact semigroup replacement linked below; REFUTED single-family and uniform bounded signed-relation claims (MISTAKE-544); historical thinness and coefficient-uniformity extrapolations are not proved.
author: opus-2026-07-20-S428
depends_on: [THM-1655 (binomial/unique-minimal positivity), THM-1735 (finite-place), THM-1715 (positivity), THM-1740 (bounded GMC(2) = finite Groebner), klein cross-shell]
---

> **CURRENT CORRECTION — 2026-09-05, MISTAKE-544.** The singleton/tunable
> dichotomy survives, but the exact `{-N,1,N+2}` classification and every
> uniform bounded-height resonance explanation are **REFUTED**. For primitive
> `{-a,b,c}`, define `g=gcd(a+b,a+c)`, `A=(a+b)/g`, `B=(a+c)/g`. The first
> return is tunable iff `a-AB in <A,B>`, and then its mass is `g`. The exact
> all-level carry profile and an unbounded shortest-relation family are
> [proved and independently audited here](../../05-knowledge/results/synthesis_20260905_moments_trinomial.md).
> `{-3,1,9}` refutes the old single-family claim; `{-12,27,40}` needs signed
> relation norm eight, and the family makes this norm arbitrarily large.
> “Generic,” “thin,” and finite-place uniformity in the historical prose
> below are not established by its finite census. The general two-rung
> coefficient claim remains OPEN; current NC2/GMC(2) status is in THM-2022.

# THM-1755 — The angular uniform piece: corrected dichotomy and historical proposal

The sections below preserve the original proposal and failed extrapolations
for provenance. Their universal claims are superseded by the correction above.

HYP-8540 factored unbounded GMC(2) into an **angular** uniform bound (mine) and a **radial**
uniform bound (klein's shell lemma). This note works the angular half.

## 1. The dichotomy

For a `k`-nomial charge pattern, gauge-fix (THM-1680) to `k−2` free parameters. Every
`CT(Λ^m) = Σ_y c_{m,y} a^y` has **positive** multinomial coefficients (THM-1715). At the
minimal level `m_0` (smallest `m` with `0` in the `m`-fold charge sumset):

- **Unique minimal representation** ⟹ `CT(Λ^{m_0})` is a **single positive monomial**
  `c·a^{y_0}`, whose only zero is `a = 0` (a lower `(k−1)`-nomial). So no nullcone element:
  **TNC by THM-1655 (positivity), with no resultant, no Gröbner, uniformly in the span.**
- **Non-unique minimal representation** (**tunable**) ⟹ `CT(Λ^{m_0})` has `≥ 2` terms and can
  vanish at tuned complex `a`; TNC needs the finite-place certificate (THM-1735).

> **Generic patterns are unique-minimal.** As the span grows with fixed charge-count, the
> charge geometry almost always gives a unique minimal `0`-representation — verified: of the
> `{−2,1,M}` family, only `M = 4` is tunable; `M = 6,8,10,12` are all unique-minimal.

## 2. The tunable subset is resonance-thin (verified)

A charge-triple `{−N, c, M}` is **tunable iff `0` has two minimal representations**, which
forces a **small-coefficient charge resonance** — a primitive integer relation
`α(−N) + βc + γM = 0` with small `|α|+|β|+|γ|`. Every tunable triple found (`N ≤ 5, M ≤ 7`,
**7 total**) has such a relation with `|α|+|β|+|γ| ≤ 6:`

| tunable `{−N,c,M}` | minimal reps | resonance |
|---|---|---|
| `{−2,1,4}` | `(−2,−2,4),(−2,1,1)` | `−4(−2) − 2·4 = 0` |
| `{−3,1,5}` | `(−3,−3,1,5),(−3,1,1,1)` | `−3(−3)+1−2·5 = 0` |
| `{−3,2,7}` | `(−3,−3,−3,2,7),(−3,−3,2,2,2)` | `−3(−3)−2−7 = 0` |
| `{−4,−1,2}` | `(−4,2,2),(−1,−1,2)` | `−2(−4)−4·2 = 0` |
| `{−4,1,6}` | 3 reps | `−3(−4)−2·6 = 0` |
| `{−5,−1,3}` | `(−5,−1,3,3),(−1,−1,−1,3)` | `−2(−5)−(−1)−3·3 = 0` |
| `{−5,1,7}` | 3 reps | `−3(−5)−1−2·7 = 0` |

**The `c = ±1` tunable triples are exactly `{−N, 1, N+2}`** (`N=2..5 → M=4,5,6,7`) and their
reflections `{−N, −1, N−2}`. The rest carry `c = ±2` with an analogous resonance. So the
tunable locus is a **union of finitely many arithmetic progressions in `(N, M)`**, one per
small resonance vector — thin and explicit.

## 3. The angular uniform statement

> **Angular uniform (this note).** For any charge-count `k`, TNC holds on **all but a thin
> resonance-characterised family** of patterns by THM-1655 positivity (uniform, span-
> independent, no computation). On the thin tunable family, THM-1735's finite-place
> certificate closes each pattern, with the good prime bounded by the resonance data.

This is the angular half of HYP-8540: the infinitely many charge patterns collapse to
(generic: one positivity argument) + (tunable: a thin, arithmetically-listed family). What
remains for a *complete* angular uniform is a **closed-form of the tunable-resonance
condition** — a characterisation "`{−N,…,M}` is tunable iff its charge set contains a
primitive relation of height `≤ h(k)`" with an explicit `h(k)` — which the data supports
(`h(3) ≤ 6`).

## 4. The radial half, framed as a resultant tower (for klein)

The radial uniform (HYP-8540 (2)) is klein's cross-shell descent. In the finite-Gröbner
language (THM-1740): the shells `ρ = |z|^2` carry the Hermite/Laguerre coefficients, and
klein's mixing functional `L` couples shell `s` to shell `s+1`. Frame it as:

```
shell s coefficients  --Res_s-->  shell s+1 coefficients
```

where `Res_s` is the elimination resultant of the shell-`s` nullcone equations against the
coupling. **Cross-shell descent = the resultant tower `Res_1, Res_2, …`; bottom-up
emptiness propagates if each `Res_s ≠ 0` on the surviving locus** — the finite-Gröbner
analogue of klein's convergence lemma. The angular dichotomy (§1) supplies the per-shell
emptiness; the radial uniform is the **termination/propagation** of the tower, which is
klein's to close. The shared object is again the branch monodromy (roots of unity), now
indexing both the angular representations and the radial shells.

## 5. Next

1. **Closed-form tunability.** Prove `{charges}` tunable ⟺ a primitive charge-relation of
   bounded height; give `h(k)`. Closes the angular uniform outright.
2. **Resultant-tower termination.** With klein: show `Res_s ≠ 0` propagates, so the tower
   ends — the radial uniform.
3. Together (1)+(2) = unbounded GMC(2) via THM-1740's finite-test framing.

## Verification

`04-computation/tnc_tunable_resonance_opus_S428.py` (the 7 tunable triples and their
resonances; the `{−N,1,N+2}` family), `04-computation/tnc_resultant_height_opus_S428.py`
(generic `{−2,1,M}` are unique-minimal as span grows). Outputs in `05-knowledge/results/`.
