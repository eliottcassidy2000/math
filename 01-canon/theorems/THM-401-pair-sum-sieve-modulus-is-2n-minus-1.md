---
id: THM-401
title: The pair-sum sieve's modulus is C = 2n−1 — two independent identities (Farey, odd-square)
status: PROVED (the modulus identity, parts i–iii); the transversal characterization is VERIFIED-with-exceptions
source: opus-2026-06-03-S590
depends_on:
  - HYP-2059   # pinch lemma (S557): optimal witness at a pair-sum modulus
  - HYP-2128   # 8·C(n,2)+1 = (2n-1)^2
related:
  - HYP-2135   # labelled sumset-support calculus after the modulus identity
  - HYP-2138   # composite C=2n-1 nonunit shell sporadic branch
  - HYP-2136   # chi=2 AP/interval regular orbit
  - HYP-2141   # tournament-level interval/additive beat, not Paley/QR
  - HYP-2137   # n=14 multiple-of-3 shell-swap branch
  - HYP-2084   # antipodal transversals mod 2n-1
  - HYP-2088   # clock-blocker ledger / second-gap exceptions (S573)
---

# THM-401 — the pair-sum sieve's modulus is `C = 2n−1`

## Goal

Make rigorous, **purely arithmetically (sidestepping the resonance-energy / measure
route)**, that the natural modulus of the pair-sum (pinch) sieve at the LRC floor
`δ = 1/n` is literally `C = 2n−1` — the additive "summand-shell" modulus. This turns
"the multi-sieve is the additive face" from a slogan into an identity.

## Setup

`M(S) = max_t min_i ‖v_i t‖`; floor `δ = 1/n`. By the **pinch lemma** (HYP-2059 / S557),
`M(S)` is attained at a *pair-sum* time `t = m/(v_a+v_b)`, and equals `r/s` with
`s = (v_a+v_b)/gcd(v_a,v_b)` (a reduced pair-sum). So the sieve's witnesses live at
**pair-sum denominators**, and its resolution at the floor is set by the smallest
pair-sum value just above `δ`.

## Theorem (the modulus identity)

Let `C = 2n − 1`.

**(i) Farey-neighbor identity.** `1/n` and `2/(2n−1)` are adjacent in the Farey sequence
`F_{2n−1}`:
```
det( 1  2 ; n  2n−1 ) = 1·(2n−1) − 2·n = −1,   so |ad − bc| = 1.
```
*Proof:* direct. ∎ Consequently `2/(2n−1)` is the **immediate Farey successor of `1/n`**
at denominator level `2n−1`; the open interval `(1/n, 2/(2n−1))` contains no fraction of
denominator `< 3n−1` (its first interior fraction is the mediant `3/(3n−1)`). So the
**pinch value one step above the floor is `2/(2n−1)`, of denominator `C = 2n−1`.**

**(ii) Summand-shell structure.** `C = 2n−1` is odd, and modulo `C` the nonzero residues
partition into `n−1` antipodal **summand shells**
```
P_a = { a, C − a } = { a, −a }   (a = 1, …, n−1),
```
with no fixed midpoint (no `a ≡ −a`, since `C` odd). Addition defines the shells:
`a + (C − a) = C ≡ 0`. ∎

**(iii) Odd-square identity (HYP-2128).** `C = 2n − 1 = √( 8·C(n,2) + 1 )`, since
`8·C(n,2) + 1 = 4n(n−1) + 1 = (2n−1)²`. So `C` is the **odd-square companion of the
triangular pair-count** `C(n,2) = T_{n−1}`. ∎

**(iv) The identity.** The sieve's floor-companion modulus (from (i), the Farey
successor `2/(2n−1)`) and the additive summand-shell modulus (from (ii), the shells
`{a, C−a}`) and the triangular-pair modulus (from (iii), `√(8·C(n,2)+1)`) are **the same
number `C = 2n−1`** — established three independent ways, none using measure. Hence:

> **The pair-sum (pinch) sieve's modulus is literally `C = 2n − 1`, the additive face.**
> The witness family `t = k/(2n−1)` is the summand graph at the odd node `C`, on whose
> shells `P_a = {a, C−a}` "addition" acts (the antipodal pairing) and the multiplicative
> units act by inverse clocks `k = a^{−1}` (`gcd(a,C)=1`).

**(v) Odd part / the 2-adic split.** `C = 2n−1` is odd; the `2`-adic (doubling) part of
the sieve factors off separately (the `⟨×2⟩` Frobenius/fragmentation, HYP-2126), so the
**odd part of the pair-sum sieve's modulus is exactly `2n−1`** — the additive shell
modulus — while the `2`-part is the multiplicative/dynamical face. This is the
add/mult split of S586–S589 at the level of the sieve modulus.

## Honest scope (what is NOT claimed)

- **(i)–(v) are proved** (elementary arithmetic + the pinch lemma + HYP-2128).
- The **transversal characterization** — "every config with `M(S) < 2/(2n−1)` is
  floor-tight and a perfect transversal of the shells `P_a`" — is **verified in bounded
  boxes (S570/S572) but is NOT a theorem**: S573 (HYP-2088) exhibits second-gap-lifted
  rows, e.g. `(1,5,6,11,16,17)` at `n=7` with `M = 5/33 ∈ (1/7, 2/13)`. So the *modulus*
  is pinned exactly; the *complete classification at that modulus* still needs the
  three-clock blocker ledger (HYP-2088), not a global second-gap-empty claim.

## Why this matters (sidestepping resonance energy)

The resonance-energy / measure bound (S550) is *Vitali-blind* at the floor (S551o) — it
cannot reach the measure-zero worry-set. THM-401 replaces it with **pure arithmetic**:
the floor `1/n` has an exact Farey companion `2/(2n−1)`, whose denominator `2n−1` is the
additive summand-shell modulus and the odd-square root of the triangular pair-count. The
pair-sum sieve therefore *is* the additive face — an identity, not an analogy — and the
remaining work (the transversal classification, the second-gap exceptions) is a finite
arithmetic ledger at the fixed modulus `2n−1`, not a measure estimate.

**Artifacts:** `04-computation/lrc_pairsum_modulus_2nm1_s590.py` (+`.out`). Builds on
HYP-2059 (pinch), HYP-2128 (`8C(n,2)+1=(2n-1)²`), HYP-2084 (transversals), HYP-2088
(second-gap ledger), S572o (summand graph at `C=2n-1`). See reflection
`07-reflections/lrc-the-pair-sum-sieve-modulus-is-2n-1-an-identity-s590.md`.
HYP-2135 adds the next labelled support calculus layer: speed shells, pair
shells, denominator shields, unit visibility, and lift denominators. HYP-2141
adds the tournament-level reading: the additive interval circulant is the beat,
while multiplicative units are witness symmetries rather than Paley/QR beats.
HYP-2138 isolates the composite-`C` nonunit-shell sporadic branch, refined at
n=14 by the HYP-2137 multiple-of-3 shell swap.
