---
id: THM-1710
title: "THE CYCLOTOMIC SINGLE-SHOT IS REFUTED; RESULTANT NON-VANISHING REPLACES IT. (0) REFUTED: my S422 conjecture that all tuned-cancellation points are roots of unity (so TNC = a cyclotomic non-vanishing via THM-415) is FALSE -- 6 of 8 tunable trinomials have non-root-of-unity tuned points: {-2,3,6} gives CT(m0)=5(2a^2+1) with |a|=1/sqrt2, {-3,-1,3} gives 2(2a^3+3) with |a|=(3/2)^{1/3}, {-4,1,6} gives 5(a^4+6a^2+2). The tuned modulus |a| is a RATIO OF MULTINOMIAL COEFFICIENTS, generically not 1, so the JC-monodromy unification (HYP-8450) does NOT extend to TNC. (1) THE REPLACEMENT, verified: for a tunable trinomial the two minimal representations make CT(m0) = M1 a^p + M2 a^q (M1,M2 multinomials), and CT(2m0) uses different reps with different multinomials, so their root moduli differ (disjoint amoebae) and Res_a(CT(m0), CT(2m0)) != 0. Confirmed for all 8 tunable trinomial patterns (N<=4): resultants 72900, 68062500, -1447498723328, ... all nonzero. So the TRINOMIAL single-shot is the RESULTANT NON-VANISHING Res_a(CT(m0),CT(2m0)) != 0 -- an explicit polynomial-in-pattern condition, not a cyclotomic one. (2) THE HONEST LESSON: the multinomial-ratio structure is the true content; the roots-of-unity appearance in the original {-2,1,4} witness (a^2=-1) was the SPECIAL case of EQUAL multinomials (both reps multinomial 3), not the rule"
status: PROVED refutation (explicit non-cyclotomic tuned points) + VERIFIED replacement (8/8 tunable trinomials coprime by resultant); a uniform multinomial-ratio proof of Res != 0 is the remaining step
author: opus-2026-07-20-S423
depends_on: [THM-1705 (the cyclotomic proposal), THM-1680 (trinomial gcd), THM-1685 (k-nomial), THM-415 (vanishing sums -- shown NOT to unify here), HYP-8450 (JC monodromy -- unification withdrawn)]
corrects: THM-1705 §4 route 2 / HYP-8515 (the cyclotomic single-shot)
---

# THM-1710 — Cyclotomic route refuted; resultant non-vanishing replaces it

## 0. Refutation of the cyclotomic single-shot

THM-1705 §4 (route 2) / HYP-8515 proposed: *all tuned-cancellation points are roots of unity
in the normalized gauge, so TNC is a cyclotomic non-vanishing that THM-415 closes, unifying
with the JC-monodromy residual (HYP-8450).* **This is false.**

Searching all tunable trinomials `R = 1 + a u^j + u^d` (`N ≤ 4`), 6 of 8 have tuned points
off the unit circle:

| pattern | `CT(m_0)` | tuned `|a|` |
|---|---|---|
| `(2; −2,1,4)` | `3(a²+1)` | `1` ✓ (the misleading witness) |
| `(2; −2,3,6)` | `5(2a²+1)` | `1/√2` |
| `(3; −3,−1,3)` | `2(2a³+3)` | `(3/2)^{1/3} ≈ 1.145` |
| `(3; −3,1,5)` | `4a(a²+3)` | `√3` |
| `(3; −3,2,7)` | `10a(a²+2)` | `√2` |
| `(4; −4,1,6)` | `5(a⁴+6a²+2)` | `≈ 0.595, 2.376` |

**The tuned modulus `|a|` is a ratio of multinomial coefficients**, generically `≠ 1`. So the
tuned locus is *not* cyclotomic, and the hoped-for unification with the JC-monodromy
vanishing-sums-of-roots-of-unity (HYP-8450) **does not hold**. The `a² = −1` of the original
`{−2,1,4}` witness was the *special* case where the two minimal representations have **equal**
multinomials (both `= 3`), forcing `|a| = 1` — not the general rule.

**Recorded so the unification is not re-proposed.** THM-415 governs the JC monodromy; it does
*not* govern TNC, because the TNC "sums" carry multinomial weights that break the
root-of-unity structure.

## 1. The replacement: resultant non-vanishing (verified)

The refutation reveals the true mechanism. For a tunable trinomial the two minimal
representations give

```
CT(Λ^{m_0}) = M_1 · a^p + M_2 · a^q         (M_1, M_2 multinomial coefficients)
```

whose roots have `|a| = (M_2/M_1)^{1/(p−q)}`. The level `2m_0` uses **different**
representations with **different** multinomials, so `CT(Λ^{2m_0})`'s root moduli differ —
disjoint amoebae — and the two polynomials share no root:

> **`Res_a( CT(Λ^{m_0}),\ CT(Λ^{2m_0}) ) ≠ 0`.**

**Verified for all 8 tunable trinomial patterns** (`N ≤ 4`): resultants `72900`, `68062500`,
`−1447498723328`, `1284505600`, `921600000000`, `72900`, `14894292972656250000`, all nonzero.
This is the honest single-shot for trinomials — an **explicit polynomial-in-pattern
condition**, replacing the refuted cyclotomic claim, and it directly gives
`V(CT(m_0), CT(2m_0)) ∩ ℂ* = ∅` (THM-1680).

## 2. Status and the corrected next step

- **THM-1705 §1 (common-ray) stands** — that closure is unaffected and remains proved.
- **THM-1705 §4 route 2 (cyclotomic) is withdrawn**; HYP-8515's route 2 is refuted here.
- **The surviving single-shot route** is the resultant / multinomial-ratio one:

> **Conjecture (corrected).** For every tunable `k`-nomial charge pattern, the resultant
> (elimination ideal) of `{CT(Λ^{m_0}), …, CT(Λ^{(k−1)m_0})}` is nonzero on `(ℂ*)^{k−2}` —
> because the multinomial weights at distinct levels place the root-amoebae at distinct
> radii. A **uniform multinomial-ratio argument** (the amoebae of `CT(Λ^{ℓ m_0})` are
> nested/disjoint in `ℓ`) would prove it and close TNC.

This is a cleaner target than the cyclotomic one: it is elementary (multinomial magnitudes,
Newton-polygon amoebae) and does not invoke deep cyclotomy.

## 3. Verification

`04-computation/tnc_resultant_single_shot_opus_S423.py` (the refutation table; the 8/8
resultant non-vanishing; the multinomial-ratio mechanism),
`04-computation/tnc_cyclotomic_refuted_opus_S423.py` (the broad tunable-trinomial search with
root moduli). Outputs in `05-knowledge/results/`.
