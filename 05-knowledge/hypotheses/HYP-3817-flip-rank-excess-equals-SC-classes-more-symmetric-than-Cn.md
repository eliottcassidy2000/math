---
id: HYP-3817
title: THE FLIP-RANK EXCESS = #{SELF-COMPLEMENTARY CLASSES MORE SYMMETRIC THAN THE ROTATION C_n} -- the fixed-point-sensitive instrument (|Aut| = commutant of U, refined by SC = reflection fixed points) predicts the covering excess the SPECTRUM could not see. Answers the S82 lead. INSTRUMENT: the automorphism group |Aut(T)| is the COMMUTANT of the Cayley U (perms P with PU=UP) -- a NON-spectral invariant (S82: the spectrum is reflection-blind) -- and its complement-extension Aut*(T) (auto- AND anti-automorphisms) detects SC exactly ([Aut*:Aut]=2 <=> SC = reflection FIXED point). VERIFIED (exact, n=3..7): the flip-rank covering excess = rho(n) - ceil(log2|G_n|) = 0,0,0,1,3 EQUALS #{iso classes C : C self-complementary AND |Aut(C)| > n} = 0,0,0,1,3. Mechanism unifies TWO prior findings the spectrum missed: (a) opus-S15 -- the covering OBSTRUCTION is the max-|Aut| class (n=7 = Paley heptagon |Aut|=21, argmax); high-|Aut| = few labeled reps (n!/|Aut|) = covering-hard; (b) HYP-3810 -- the SC (self-complementary) classes CARRY the excess (T-join parity). The RAW symmetry census #{|Aut|>n} = 0,0,0,1,5 OVERSHOOTS at n=7 (5 super-rotational classes: four |Aut|=9 + one |Aut|=21); filtering to SELF-COMPLEMENTARY (the reflection fixed points, HYP-3810's excess-carrier) drops the 2 non-SC |Aut|=9 classes, giving exactly 3 = excess(7). So the excess is counted by classes that are BOTH super-symmetric (|Aut|>n = more symmetric than the cyclic rotation) AND reflection-fixed (SC) -- precisely the two features invisible to the skew/Cayley spectrum (S82). RESOLUTION: (spectrum,|Aut|) resolves ~2x the spectrum alone (2,2,4,12 vs 1,2,2,6 at n=3..6) but still << |G_n|. Also confirmed: |Aut| all ODD (Moon); fiber f(C)=H(C)/|Aut(C)|, Sum_C f = 2^{C(n-1,2)} (the tiling-cube size). HONEST: VERIFIED on 5 points (0,0,0,1,3) incl 3 zeros + the raw-census-overshoot correction; but only 2 nonzero informative values, so it is a well-motivated CONJECTURE; the falsifiable prediction is excess(8) = #{SC on 8 vtx with |Aut|>8} (needs rho(8), currently infeasible, to confirm)
status: VERIFIED n=3..7 (exact) as a data match; CONJECTURE as a general law. EXACT: excess = 0,0,0,1,3 (rho=1,2,4,7,12 from S71+opus-S15; LB=ceil(log2 A000568)=1,2,4,6,9). #{SC & |Aut|>n}: n=3,4,5,6 computed by full iso-class enumeration (|Aut| dist n=6={1:41,3:12,5:2,9:1}, matches opus); n=7 by targeted fix(sigma) enumeration (sigma=(012)(345)(6): captures EVERY |Aut|>7 class since |Aut|>7 on 7 pts forces |Aut| in {9,21}, both with a 3+3+1 order-3 element) -> |Aut| dist {3:7,9:4,21:1}, SC among them {3:3,9:2,21:1} => #{SC & |Aut|>7}=3=excess. Raw #{|Aut|>n}=0,0,0,1,5 (overshoots n=7). HONEST: 2 nonzero data points; a mechanistically-grounded conjecture (opus symmetry + HYP-3810 SC-carries-excess), NOT proven; predicts excess(8)=#{SC,|Aut|>8 on 8 vtx}. Answers 'does |Aut|-resolution predict the excess': YES with the SC refinement, NO as a raw |Aut| census.
source: klein-2026-07-01-S83
depends_on:
  - HYP-3816   # S82: the excess is invisible to the spectrum (reflection-blind) -> need a non-spectral instrument
  - HYP-3810   # S77: SC classes CARRY the flip-rank excess (T-join parity) -- the SC filter
  - HYP-3814   # S81: complement = reflection; SC = fixed points; Aut = commutant of U
related:
  - HYP-3803   # S71: flip-rank rho(n)=1,2,4,7; excess first at n=6
  - HYP-3805   # opus-S15: obstruction = argmax|Aut| = Paley heptagon (n=7); max|Aut|=1,3,3,5,9,21
  - HYP-3804   # S72: skew-spectrum weakness (now: |Aut| is the missing non-spectral refinement)
results:
  - 04-computation/aut_instrument_excess_prediction_klein.py
  - 05-knowledge/results/aut_instrument_excess_prediction_klein.out
---

# HYP-3817 — the excess counts the self-complementary, super-symmetric classes

## The result
The flip-rank covering excess equals the number of tournament iso classes that are **both** self-complementary
**and** more symmetric than the cyclic rotation `C_n`:
> **`excess(n) = rho(n) - ceil(log2|G_n|) = #{ C : C self-complementary AND |Aut(C)| > n } = 0,0,0,1,3`** (n=3..7).

## The instrument (the S82 lead answered)
The spectrum is reflection-blind (S82), so the excess — carried by the SC *fixed points* — is invisible to it.
The **fixed-point-sensitive** instrument is the automorphism group:
- `|Aut(T)|` is the **commutant of the Cayley `U`** (perms `P` with `PU = UP`) — a genuinely **non-spectral**
  invariant (cospectral tournaments can differ in `|Aut|`).
- Its **complement-extension** `Aut*(T)` (auto- and anti-automorphisms) satisfies `[Aut* : Aut] = 2 <=> SC` —
  the reflection **fixed-point detector**.
- `(spectrum, |Aut|)` resolves ~2x the spectrum alone (`2,2,4,12` vs `1,2,2,6` at n=3..6), still `<< |G_n|`.

## Why this exact combination (two mechanisms unified)
- **Symmetry (opus-S15):** the covering obstruction is the **max-`|Aut|`** class — high symmetry means few
  labeled reps (`n!/|Aut|`), hence hardest to fit in a small subcube. `max|Aut| = 1,3,3,5,9,21` (n=2..7); the
  n=7 obstruction is the **Paley heptagon** (`|Aut|=21`).
- **Fixed points (HYP-3810):** the **self-complementary** classes *carry* the excess (T-join parity).
- The **raw census** `#{|Aut|>n} = 0,0,0,1,5` **overshoots** at n=7 (five super-rotational classes: four with
  `|Aut|=9`, one with `|Aut|=21`). Filtering to **SC** removes the two non-SC `|Aut|=9` classes, leaving
  `2 + 1 = 3` = excess(7). So the excess counts exactly the classes that are **both** super-symmetric **and**
  reflection-fixed — the two features the skew/Cayley spectrum cannot see (S82).

## Data (exact)
| `n` | `|G_n|` | `rho` | LB | excess | `#{|Aut|>n}` | `#{SC & |Aut|>n}` |
|---|---|---|---|---|---|---|
| 3 | 2   | 1  | 1 | 0 | 0 | 0 |
| 4 | 4   | 2  | 2 | 0 | 0 | 0 |
| 5 | 12  | 4  | 4 | 0 | 0 | 0 |
| 6 | 56  | 7  | 6 | **1** | 1 | **1** |
| 7 | 456 | 12 | 9 | **3** | 5 | **3** |

n=7 via targeted `fix(sigma)`, `sigma=(012)(345)(6)` (captures every `|Aut|>7` class): `|Aut|` dist
`{3:7, 9:4, 21:1}`, SC among them `{3:3, 9:2, 21:1}` → `#{SC & |Aut|>7} = 2+1 = 3`.

## Honest scope
Verified exactly on n=3..7 (five points, including three zeros and the raw-census-overshoot correction), but
only **two nonzero** informative values, so this is a **mechanistically-grounded conjecture**, not a theorem.
It combines two established findings (opus symmetry-obstruction + HYP-3810 SC-carries-excess). Falsifiable
**prediction**: `excess(8) = #{SC classes on 8 vertices with |Aut|>8}` — testable once `rho(8)` is computed
(currently infeasible: `|G_8|=6880`, cube `Q_21`).

## Net
Answers the S82 question. **YES**, the `|Aut|`-resolution predicts the excess — but only with the **SC
refinement**: `excess = #{SC & |Aut|>n}`, the count of self-complementary classes more symmetric than the
rotation. The raw `|Aut|` census overshoots. The right instrument is `|Aut|` (symmetry, commutant of `U`)
*refined by* SC (the reflection fixed points) — exactly the two things the spectrum is blind to, now doing the
counting the spectrum could not.
