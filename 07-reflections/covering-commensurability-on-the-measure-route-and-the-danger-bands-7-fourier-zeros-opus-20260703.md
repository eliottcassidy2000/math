# Covering-commensurability, on the measure route: μ = (6/7)¹³ + gcd-controlled resonances, and the danger band's 7-Fourier-zeros

*opus-2026-07-03-S57. The owner asked me to pursue the covering-commensurability angle further. On the census
it was a red herring (mac-mini-S29: compressed families are loose, the rational witness is irrelevant). On
the MEASURE route it is exact and real — it gives the safe measure its resonance form, a clean heptagon
Fourier structure, and a verified correlation — but it does not close the crux. Careful and honest, given
MISTAKE-097 (two prior overclaims on this exact problem).*

## The right object: the safe measure, not the census witness

mac-mini-S29 pivoted correctly: the compressed lcm-blockers are **loose** (`M(covmin) ≈ 0.25–0.32 ≈ 4× the
danger radius`), so their large rational-witness denominator `q* ~ log M` is irrelevant — they are lonely on
a fat positive-measure set. The crux is not "find a small rational witness" but **`μ(safe) > 0`** for every
covering family, equivalently `M(covmin) ≥ 14/183 > 1/14`. So the covering-commensurability idea belongs on
the **measure**, where it becomes exact:

```
   μ(safe) = Leb{ t : ‖vᵢt‖ ≥ 1/14 ∀i } = (6/7)¹³ + Σ_{integer resonances Σ mᵢvᵢ = 0, ≠0} ∏ᵢ c(mᵢ),
   c(0) = 6/7,   c(m) = −D̂(m),   D̂(m) = ∫_{−1/14}^{1/14} e(−mx)dx = sin(πm/7)/(πm).
```

This is the *continuous* twin of the S56 census count — but the resonances `Σ mᵢvᵢ = 0` are now **exact
integer relations**, and this is the same `μ = (6/7)^k(1+R)` decomposition as my S48 (R2) resonance.

## Structural fact 1 — the danger band's 7-Fourier-zeros (the heptagon)

`D̂(m) = sin(πm/7)/(πm)` **vanishes exactly when `7 | m`** (verified: `D̂(7)=D̂(14)=0`). The band `1/14 =
1/(2·7)` carries the heptagon `7` as the *gap in its Fourier support*: any resonance supported on multiples
of 7 contributes **nothing** to the measure. This is a clean, new face of the ubiquitous `7` (the union-bound
wall `13·(2/14) = 13/7 > 1`, the covering threshold, THM-503's seven-vanishing): here it is where the
danger-band Fourier dies. It says the resonances that *can* move the measure avoid `7ℤ`.

## Structural fact 2 — the pair resonance is gcd-controlled (the commensurability)

The smallest pair resonance `mᵢvᵢ + mⱼvⱼ = 0` is `(mᵢ,mⱼ) = k·(vⱼ/g, −vᵢ/g)`, `g = gcd(vᵢ,vⱼ)`, and its
term scales like **`g²/(vᵢvⱼ)`**. So **commensurate pairs (large `g`) dominate the resonance**, and a
**covering family shares small factors `2,3,…`**, giving larger gcds — a structurally *stronger* resonance
than a generic 13-tuple has. Verified correlation (`sum gcd²` vs `μ`):

| family | covering | `μ(safe)` | `R = μ/(6/7)¹³ − 1` | `Σ gcd²` |
|---|---|---|---|---|
| loose lcm-block | ✓ | 0.108 | **−0.20** | 8313 |
| tight cov `{…126}` | ✓ | 0.105 | −0.22 | 4799 |
| `{2..14}` | ✓ | 0.061 | −0.55 | 321 |
| deep well `{1..12,182}` | ✓ | 0.024 | −0.82 | 321 |
| random cov `{…168}` | ✓ | 0.009 | **−0.94** | 351 |

Higher commensurability ⟹ larger `μ` (less negative `R`, safer). **The angle is real.**

## Why it does not close the crux (honest)

The families with `R` near `−1` — the ones where `μ` is barely positive, the actual crux — are exactly the
ones with **low** commensurability (`Σ gcd² ≈ 350`, deep well and random-cov). Commensurability *helps*, but
the crux is precisely where the family has *little* of it, and there the pair-resonance signs are delicate
(`D̂(m)` changes sign at `m > 7`). So the measured correlation does not give a *lower bound* `R > −1` for all
covering families — and `R > −1 ⟺ μ > 0 ⟺ M ≥ 1/14 ⟺` **LRC(14)**. The commensurability angle supplies the
resonance mechanism and the heptagon Fourier structure; it does not supply the sign control that the tight,
low-commensurability covering families need. That control is the measure floor (klein) / the covering-min
`≥ 14/183` — the hard core.

## What this contributes

- **The exact measure-resonance form** `μ = (6/7)¹³ + Σ_{Σmᵢvᵢ=0} ∏c(mᵢ)` — the object the measure route
  (klein, mac-mini-S29) must lower-bound, made explicit and connected to my S48 `μ=(6/7)^k(1+R)`.
- **The danger band's 7-Fourier-zeros** — a clean structural constraint (resonances avoid `7ℤ`), a new face
  of the heptagon that the measure floor can use.
- **The gcd-control of the pair resonance + the verified correlation** — commensurability is the
  measure-side content of the covering hypothesis.
- **Honest boundary**: not a closure; the crux is the low-commensurability tight covering families, where the
  angle is weakest.

## Status

- **Verified (correct):** `μ = (6/7)¹³ + Σ resonances`; `D̂(m)=0 iff 7|m`; the pair term `~ g²/(vᵢvⱼ)`;
  commensurability–measure correlation. (`measure_commensurability_route_opus_20260703_S57.py`.)
- **Aligned:** mac-mini-S29 (census red herring, measure route), klein (measure floor), my S48 (R2 resonance),
  the heptagon `7` (THM-503, the union-bound wall). Not a proof.
- **Open (= LRC(14)):** `R > −1` for the tight, low-commensurability covering families.

Related: mac-mini-S29 (the pivot), HYP-4055/kps (`q*≤13 ln M`, now seen as the red-herring side), HYP-4057/S56
(the census circle-method twin), HYP-4013/S48 (R2 = `μ=(6/7)^k(1+R)`), THM-503 (seven-vanishing), MISTAKE-097.
Script: `04-computation/measure_commensurability_route_opus_20260703_S57.py`.
