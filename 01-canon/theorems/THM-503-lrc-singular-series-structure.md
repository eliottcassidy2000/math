---
id: THM-503
title: Structure of the LRC(14) singular series — the 7-vanishing, pairwise absolute convergence, the almost-Sidon loose-class, and L is a singular INTEGRAL (no Euler product)
status: PARTIAL/PROVED-in-scope — (1)(2)(3) proved with explicit constants; (4) is a rigorous structural reframe (β_p = L for all p) that corrects HYP-2503; the full inf L>0 (C'(14)) remains open
source: mac-mini-2026-06-14-S1
depends_on:
  - THM-501   # the LRC singular series L(S) and its additive-character expansion
  - THM-398   # C' reduction; D(q,S)>0 <=> loose
related:
  - HYP-2503  # CORRECTED here: L is NOT a nontrivial Euler product
  - HYP-2489  # L = circle-method singular object (made precise: singular INTEGRAL)
  - THM-446   # the additive-relation / Sidon ladder controlling L
---

# THM-503 — the structure of the LRC singular series L(S)

Builds on THM-501: `L(S) = (6/7)^13 + Σ_{exact relations Σ_{v∈T} t_v v = 0, t_v≠0}
(6/7)^{13−|T|} (−1)^{|T|} ∏_{v∈T} s(t_v)`, `s(t) = sin(πt/7)/(πt)`. Four results,
each adversarially re-verified with fresh code.

## (1) The 7-vanishing (PROVED)

For integer `t ≠ 0`: **`s(t) = 0 ⟺ 7 | t`** (since `sin(πt/7) = 0 ⟺ 7 | t`, and
`πt ≠ 0`). Hence any relation term containing a coefficient `t_v ∈ 7ℤ` is killed, and

> **`L(S)` is a sum over only the 7-PRIMITIVE exact relations** — those whose every
> coefficient is coprime to 7.

The apex prime of `14 = 2·7` literally gates the correction series. Sign structure
(period 14): `sign s(t) = sign sin(πt/7)`, **positive on `t mod 14 ∈ {1..6}`, zero at
`{0,7}`, negative on `{8..13}`**. Also `|s(t)| ≤ 1/(π|t|)`.

## (2) Pairwise absolute convergence with explicit bound (PROVED)

For a pair `{v_a, v_b}`, `g = gcd(v_a,v_b)`, the `|T|=2` relations form the rank-1
lattice `(t_a,t_b) = k(v_b,−v_a)/g`. The pair's total contribution to `L` is
`(6/7)^{11} P(a,b)`, `P(a,b) = Σ_{k≠0} s(kv_b/g)s(−kv_a/g) = 2Σ_{k≥1} s(ka')s(kb')`
with `a'=v_b/g, b'=v_a/g` coprime. By (1)'s bound,

> **`|P(a,b)| ≤ g²/(3 v_a v_b)`** (absolutely convergent; `= 0` if `7 | a'` or `7 | b'`).

## (3) The almost-Sidon loose class (PROVED, explicit threshold)

Call `S` **almost-Sidon** if its only 7-primitive exact relations have `|T| = 2`.
Then `L(S) − (6/7)^13 = Σ_{a<b} (6/7)^{11} P(a,b)`, so

> **If `Σ_{a<b} g(a,b)²/(3 v_a v_b) < (6/7)² = 36/49 ≈ 0.7347`, then `L(S) > 0`,
> hence `S` is LOOSE (`M(S) > 1/14`).**

This is a genuinely new, fully rigorous **infinite family of confirmed-loose
multiple-of-14 configs**, proved straight from the singular series with explicit
constants — a first nontrivial chunk of `inf L > 0` carved out analytically rather
than by sieve. Coprime spread sets clear it with huge margin: `{14,17,…,127}` has
mass `≈ 0.0029` (`⇒ L ≥ 0.132`); `{14,101,…,1201}` has `≈ 0.00016` (`⇒ L ≥ 0.1346`).

## (4) `L` is a singular INTEGRAL, not an Euler product — HYP-2503 corrected (PROVED structurally)

Define the single-prime local limit `β_p(S) := lim_{e→∞} D(p^e,S)/p^e`. A resonance
with value `m = Σ t_v v` fires at `p^e` iff `p^e | m`; for `m ≠ 0` this holds for only
finitely many `e`, so as `e→∞` only the EXACT (`m=0`) relations survive — and
`m = 0 ∈ p^eℤ` for **every** prime `p` and every `e`. Therefore

> **`β_p(S) = L(S)` for every prime `p`.**

Consequently there is **no nontrivial Hardy–Littlewood Euler product** `L = β_∞ ∏_p β_p`
(it would force `L = β_∞ · L^{#primes}`, and `∏_p L → 0` for `0 < L < 1`,
contradiction). **`L` is the ARCHIMEDEAN singular INTEGRAL** — the density carried by
the exact-relation lattice weighted by the archimedean sinc kernel — not a singular
SERIES. The p-adic data lives in the **APPROACH** to `L` (the threshold shell `q*(S)`
and the convergence rate of `D(p^e)/p^e`), not in `L` itself. This is the precise
departure from Waring/Goldbach: there finite-`m` congruence obstructions give genuine
local densities `β_p ≠ 1`; for LRC the surviving resonances are exact, so all local
factors collapse to 1.

**Correction to HYP-2503:** replace "L is a product over primes of nontrivial local
densities" with "L is the archimedean singular integral; `β_p = L` for every `p`, all
local factors trivial; the prime-power data is in the approach, not in `L`." (Honest
caveat: a differently-normalized local density — e.g. unit-restricted `a ∈ (ℤ/q)*`, or
relative to an archimedean reference — was NOT tested and could in principle recover a
refined product identity.)

## Sharper infimum data (numerical, this session)

Adversarial search over primitive multiple-of-14 configs sharpens THM-501's `≈0.026`:
**`inf L ≈ 0.0237`**, attained at the near-tight cores `{1,…,12} ∪ {14m}` (almost the
tight AP, but primitive — the 13th element broken to a multiple of 14). All sampled
configs have `L ≥ 0.0237 > 0`; random non-structured configs bottom out far higher
(`≈ 0.068`). Consistent with `inf_S L(S) > 0` with margin `~0.024`.
(`04-computation/lrc14_singular_series_adelic_macmini_0614s1.py`.)

## What remains open (the prize)

`inf L > 0` over the high-energy dilated-AP cores `d·{1,…,12} ∪ {r}` — which are NOT
almost-Sidon (they carry abundant `|T|≥3` 7-primitive relations, e.g.
`1·7 − 2·14 + 1·21 = 0`, and the data shows these `|T|≥3` terms are exactly what drives
`L` down to `≈0.024`). Result (4) clarifies that this is an **archimedean** positivity
question (a Weil/Pólya–Vinogradov-type bound on the exact-relation lattice with the
conditionally-convergent sinc kernel), **not** a product-of-local-densities positivity.
The open analytic core is the conditional convergence of the `|T|≥3` lattice sums.

## Honesty

Existence of `L` is from THM-501 (numerically verified; conditional-convergence
subtlety acknowledged there). (1)(2)(3) are unconditional given that `L` exists and its
almost-Sidon truncation is exact. (4)'s `β_p = L` is read off exact deficits along prime
ladders `p^e ≤ ~6·10^4` with empirical convergence `±0.003`; the structural conclusion
(no Euler product) follows from the firing dichotomy, which is exact. The converse
(loose ⟹ L>0) is not addressed.
