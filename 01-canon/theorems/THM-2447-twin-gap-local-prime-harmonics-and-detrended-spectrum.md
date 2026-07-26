---
id: THM-2447
title: "Twin-gap local prime harmonics and the detrended spectrum"
status: >
  PROVED (local factor lemma: for p >= 5 the admissible residues
  m mod p for the double twin pattern at gap g number exactly
  p - c_p(g), c_p(g) = 2 if p|g, 3 if p|9g^2-1, 4 otherwise; primes
  2, 3 impose no condition; exhaustively re-checked for all
  p in [5,199], g in [1,2p+2], 8,532 cases) + VERIFIED-EXACT
  (census to centers 10^8: after detrending by the fitted
  cumulative-availability model, the empirical gap counts recover
  every measured prime's harmonic amplitudes (p-2)/(p-4) and
  (p-3)/(p-4) within 2.4% across p in {5,7,11,13,17,19}, twelve
  measurements, with the wrong amplitude (p-1)/(p-4) rejected at
  every p) + VERIFIED (model layer: dyadic-window Cramer fits
  N_j(g)/w(g) ~ A_j exp(-lambda_j W(g)) at R^2 >= 0.988 with
  window-decaying lambda_j). The naive class-mean comparison is
  PROVEN-BIASED by demonstration: it fails by up to 23% at p = 19
  because availability decay is confounded with class placement.
  The theorem does not prove any all-n statement about twin gaps,
  the full Hardy-Littlewood asymptotic, or the convergence rate of
  the truncated weight product.
source: kind-pasteur-2026-07-26-S131b
depends_on:
  - THM-2443-twin-rank-parent-parity-margins-and-boundary-crossing-bridge
related:
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
  - HYP-9025-twin-gap-singular-series-partner-law
  - HYP-1994-twin-goldbach-necklace
script: 04-computation/twin_gap_prime_harmonics_thm2447.py
output: 05-knowledge/results/twin_gap_prime_harmonics_thm2447.out
script_sha256: 46a079453699f8fa38e53f4434a039385320cefb5a19217fa85a7f5e4488eea6
output_sha256: 9ccbed151df646f914e3f282ab4e34e43d9c1afdc3b0f17da2850a5a1d9ff4be
hash_basis: working-tree bytes (LF)
---

# THM-2447 -- each prime plays one line in the twin-gap spectrum

**PROVED + VERIFIED-EXACT + VERIFIED** as itemized in the status.

Typing per MISTAKE-268: `K = A002822` (twin ranks). This theorem
makes HYP-9025's "singular-series partner law" exact at the local
layer and measurable at the global layer: the distribution of
consecutive twin-rank gaps decomposes, empirically and to stated
precision, into a **continuous availability spectrum** times
**discrete prime harmonics** with exactly computable amplitudes.

## 1. The local factor lemma (PROVED)

Fix a prime `p >= 5`, a gap `g >= 1`, and let `a = 6^{-1} mod p`.
The pattern "both `m` and `m+g` are twin ranks" requires the four
numbers `6m -+ 1, 6(m+g) -+ 1` to be nonzero mod `p`, i.e. `m mod p`
must avoid the set

```text
F_p(g) = { a, -a, a - g, -a - g }  (mod p).                     (1)
```

`a != -a` since `p` is odd and `a != 0`. The only possible
coincidences are: `a = a - g <=> p | g` (two coincidences, pairing
both columns); `a = -a - g <=> p | 3g + 1`; `-a = a - g <=>
p | 3g - 1` (one coincidence each; jointly `p | 9g^2 - 1`); and
`p` cannot divide two of `{g, 3g-1, 3g+1}` at once. Hence

```text
|F_p(g)| = c_p(g) = 2 if p | g;  3 if p | 9g^2-1;  4 otherwise, (2)
```

and the admissible density is `(p - c_p(g))/p`. Primes `2, 3`
impose no condition (`6m -+ 1` is a unit mod 6). Relative to the
generic class, the **p-th harmonic amplitude** is

```text
h_p(g) = (p - c_p(g)) / (p - 4)
       = (p-2)/(p-4) on p | g,   (p-3)/(p-4) on p | 9g^2-1.     (3)
```

The companion re-proves (2) by exhaustive residue enumeration for
every `p in [5,199]` and every `g in [1, 2p+2]` (8,532 cases,
covering every class of every listed prime).

## 2. The measurement design, and why the naive one fails

Empirically the gap counts `N(g)` (census: centers `<= 10^8`,
440,310 gaps) carry two spectra at once: the smooth availability
decay in `g` (gaps compete for total mass -- the Cramer layer) and
the discrete arithmetic lines (3). The naive estimate of `h_p` --
ratio of class means of `N(g)/w_{p-omitted}(g)` -- is **biased
down, increasingly with `p`** (measured `0.87` vs exact `1.13` at
`p = 19`, a 23% error), because the class `p | g` contains only
gaps `g >= p`, which sit lower on the availability curve. The
companion keeps this failure visible as a negative control.

The correct design detrends first: fit
`log(N(g)/w_{p-omitted}(g)) = -lambda W_{p-omitted}(g) + const` on
the generic class only (`W` = cumulative weight, the availability
coordinate), then compare class means of the detrended residuals.

## 3. Results (VERIFIED-EXACT at stated tolerance)

```text
p   class      naive     detrended   exact     rel.err
5   p|g        2.8115    3.0244      3.0000    0.81%
5   p|9g^2-1   2.0049    2.0213      2.0000    1.06%
7   p|g        1.5390    1.6604      1.6667    0.38%
7   p|9g^2-1   1.3470    1.3356      1.3333    0.17%
11  p|g        1.2321    1.3090      1.2857    1.81%
11  p|9g^2-1   1.1511    1.1489      1.1429    0.53%
13  p|g        1.0530    1.2024      1.2222    1.62%
13  p|9g^2-1   1.1188    1.0982      1.1111    1.16%
17  p|g        0.9185    1.1576      1.1538    0.32%
17  p|9g^2-1   1.1167    1.0967      1.0769    1.83%
19  p|g        0.8693    1.1063      1.1333    2.39%
19  p|9g^2-1   1.0527    1.0735      1.0667    0.64%
```

Hostile control: the wrong amplitude `(p-1)/(p-4)` fits worse than
the exact `(p-2)/(p-4)` at every measured prime.

Window layer (HYP-9025 prediction 1, cumulative form): in each
dyadic window `[2^j, 2^{j+1})`, `j = 20..23`, the fit
`N_j(g)/w(g) ~ A_j exp(-lambda_j W(g))` achieves `R^2 >= 0.988`,
with `lambda_j` decreasing `0.0127 -> 0.0101` as the window grows --
the availability rate thins like the twin density itself, as the
Cramer model demands.

## 4. Scope and reading

The lemma (2)-(3) is exact and elementary; the global statements
are census-verified fits, not theorems for all `n`. Truncating the
weight at `p <= 199` leaves `1 + O(1/p)` tail factors absorbed into
the fitted constants. Nothing here proves the Hardy-Littlewood
asymptotic. The reading for the repo's harmonics programme: **the
"which partner combines" law of the owner's A014574 puzzle is a
spectral statement -- one exactly computable line per prime, visible
only after the continuous spectrum is removed.** The same grammar
(discrete character lines over a smooth carrier) is what the LRC(14)
91-stalk machinery manipulates exactly (mixed `F_7 x F_13`
characters, THM-2436; klein's nonzero `C_7 x C_13` tensor cells),
with the primes 7 and 13 in place of the sieve primes here.

## 5. Reproduction

```bash
python 04-computation/twin_gap_prime_harmonics_thm2447.py
```

Deterministic (sieve + exact counts + least squares); every check a
hard failure; final line `ALL CHECKS PASSED`.
