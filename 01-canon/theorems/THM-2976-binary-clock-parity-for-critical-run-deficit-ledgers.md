---
id: THM-2976
title: "Binary-clock parity for critical-run deficit ledgers"
status: >
  PROVED + VERIFIED-EXACT. In the THM-2966 spine normal form with any depth
  law d_M, the choice-independent forced parity of the doubled homogenized
  deficit ledger, beta_M(x) = (1+x)^{A_M} + (1+x^{M+1})(1+x)^{d_M} mod 2
  (A_M = M+d_M+1), obeys a binary clock: (T1) beta_M vanishes identically
  whenever M+1 is a power of two, for every d_M; (T2) otherwise its minimal
  positive odd position is exactly o*(M) = 2^{v_2(M+1)}; (T3) o*(M) sits on
  the isolated corner cell (o = d_M+1) iff d_M = 2^{v_2(M+1)}-1; (C4) the
  ladder rates gamma = 1/(2^k-1), D0=0 hit corner timing at every
  M = 2^r-2^{r-k}-1; (C5) for D0=0 a linear law floor(gamma*M) has
  infinitely many corner-timed levels iff gamma = 1/J with J odd (J=1 and
  the even-J unit fractions degenerate to the vacuous T1 clock), so
  gamma = 1/3 is the largest corner-clocked rate below the classical
  gamma = 1.
source: death-star-2026-07-30-coinC2
depends_on:
  - THM-2966
related:
  - THM-2160
  - THM-2225
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_binary_clock_parity_thm2976.py
output: 05-knowledge/results/amm12592_binary_clock_parity_thm2976.out
lean: >
  04-computation/lean/TournamentH7/TournamentH7/BinaryClockParity.lean —
  T1 (checkpoint_vanishing) and T2 (clock_coeff_lt / clock_coeff_top)
  kernel-checked in the TournamentH7 project via the expand-bracket
  factorization; axiom audit standard (propext, Classical.choice,
  Quot.sound); single-module `lake build TournamentH7.BinaryClockParity`.
---

# THM-2976 -- binary-clock parity for critical-run deficit ledgers

Frame: THM-2966 spine normal form for AMM 12592 fair extractors. Doubled
deficits `delta_{z,o}` (integers, `|delta| <= binom(d_m,k)`,
`delta == binom(d_m,k) mod 2` by Lucas) sit at monomials `p^z q^o`;
`D_M(p) = (1/2) sum_{m <= M} delta p^z q^o`; row `m` occupies the
anti-diagonal `z + o = A_m = m + d_m + 1`.

## Lemma 0 (forced-parity closed form; lane D, re-proved here)

Homogenize `2 D_M` to degree `A_M` by `(p+q)`-padding and mark the
`q`-exponent by `x`. The parity vector of the resulting integer coefficient
list is choice-independent and equals

```text
beta_M(x) = (1+x)^{A_M} + (1 + x^{M+1}) (1+x)^{d_M}   (mod 2).      (B)
```

*Proof.* Telescoping (THM-2966 (4)) gives
`2 D_M = 2 S_M - 1 + p^{M+1} + q^{M+1}`, where `2 S_M` is the doubled
scheme sum, an integer polynomial with all coefficients even (each cell
contributes `2 w p^z q^o`). Hence mod 2 the homogenized coefficients equal
those of `-(p+q)^{A_M} + p^{M+1}(p+q)^{d_M} + q^{M+1}(p+q)^{d_M}`
(`A_M - M - 1 = d_M`), which in the `x` convention is (B). Choices enter
only through `2 S_M`, i.e. not at all. QED.

## T1 (dyadic checkpoint vanishing)

**If `M + 1 = 2^r` then `beta_M == 0` identically — for every depth value
`d_M >= 0`.**

*Proof.* `A_M = 2^r + d_M`, and over `F_2` the Frobenius identity gives
`(1+x)^{2^r + d} = ((1+x)^{2^r})(1+x)^d = (1+x^{2^r})(1+x)^d
= (1+x^{M+1})(1+x)^{d_M}`, so the two summands of (B) coincide. QED.

At dyadic-timed levels the ledger carries **no forced parity obstruction at
all**: whatever half-integer debt the interior rows created, its
homogenized shadow at `M = 2^r - 1` is parity-free. This is the structural
reason the classical ratio-2 construction closes its books on dyadic
boundaries (THM-2966 sec 3.5's exact checksum identity at dyadic `M`).

## T2 (the clock)

**If `M + 1` is not a power of two, the minimal `o >= 1` with
`[x^o] beta_M = 1` is exactly `o*(M) = 2^{v}, v = v_2(M+1)`.**

*Proof.* Write `M + 1 = 2^v c`, `c >= 3` odd. For `o <= M` the shift term
in (B) is absent, so `[x^o]beta_M = binom(A_M, o) + binom(d_M, o) mod 2`.
Since `A_M = d_M + 2^v c`, we have `A_M == d_M (mod 2^v)`, so by Lucas
`binom(A_M, o) == binom(d_M, o) mod 2` for all `o < 2^v`: parity vanishes
below `2^v`. At `o = 2^v`, `binom(n, 2^v) mod 2` is bit `v` of `n`; adding
`2^v c` (whose bits below `v` vanish and whose bit `v` is 1) flips bit `v`
of `d_M` with no incoming carry, so `binom(A_M, 2^v) + binom(d_M, 2^v)
= 1 mod 2`. QED.

## T3 (corner timing)

**`o*(M)` lands on the isolated corner cell (`o = d_M + 1`, the unique
capacity-1 forced-odd cell of row `M`, lane D corner isolation) iff
`d_M = 2^{v_2(M+1)} - 1`.** Immediate from T2. In that case parity below
the corner is completely silent: the whole forced obstruction of the level
concentrates in the corner monomial.

## C4 (ladder membership)

For `gamma = 1/(2^k - 1)`, `D0 = 0`, the levels `M = 2^r - 2^{r-k} - 1`
(`r > k`) satisfy `d_M = floor(gamma M) = 2^{r-k} - 1` and
`v_2(M+1) = r - k`: corner timing at every scale. (`k = 1` degenerates:
`M + 1 = 2^{r-1}` is dyadic, T1's vacuous case — the classical `gamma = 1`
scheme lives entirely on parity-free checkpoints.)

*Proof.* `gamma M = (2^r - 2^{r-k} - 1)/(2^k - 1) = 2^{r-k} - 1/(2^k-1)`,
whose floor is `2^{r-k} - 1`; `M + 1 = 2^{r-k}(2^k - 1)` has `v_2 = r-k`.
QED.

## C5 (corner-clocked rates are the odd unit fractions)

**For `D0 = 0`, a linear law `d_M = floor(gamma M)` has infinitely many
corner-timed levels iff `gamma = 1/J` with `J >= 3 odd`.** All other
rational rates (and even-`J` unit fractions) have finitely many, all
transient.

*Proof.* A corner-timed level means `floor(gamma M) = 2^v - 1`,
`v = v_2(M+1)`, i.e. `M + 1 = 2^v j` with `j` odd and
`gamma(2^v j - 1) in [2^v - 1, 2^v)`. Dividing by `2^v`:
`gamma j in [1 - (1 - gamma)/2^v, 1 + gamma/2^v)`, so along any infinite
sequence of hits `gamma j -> 1` with `j` odd integers; hence eventually
`j = 1/gamma` exactly, forcing `gamma = 1/J`, `J` odd. Conversely for
`gamma = 1/J`, `J >= 3` odd, every `M = J 2^v - 1` works:
`v_2(M+1) = v` and `floor((J 2^v - 1)/J) = 2^v - 1`. For `J` even, `1/J
* (M+1) = 2^v` forces `M + 1 = J 2^v` with `v_2(M+1) = v + v_2(J) != v`.
QED.

So `gamma = 1/3` (`C = 4/3`) is the **largest corner-clocked rate below
the classical `gamma = 1`** — the first rung of the odd ladder. Referee
counts (window `M < 200000`): `1/3`: 17 hits up to `196607`; `1/5`: 18;
`1/2`, `1/4`, `2/5`, `3/8`, `2457/6592`: none past `M = 5`.

## Consequences for the minimal-C frontier (scope notes, not new claims)

1. **Extraction depth is clocked.** Any evaluation-type dual that converts
   forced parity into a 2-adic lower bound at a bias with
   `q`-numerator `2^s u` can extract valuation at most
   `s * o*(M) = s * 2^{v_2(M+1)}` from level `M`: deep extraction is
   possible only along dyadically-timed subsequences, and maximal (corner)
   extraction only at rates/levels classified by T3/C5.
2. **Certificate resonances (recorded, unproved):** the confirmed AMM
   12592 certificate (27) has `q_B/p_B ~ 0.33619` and rapidity ratio
   `r_A/r_B ~ 0.33076` straddling the top odd rung `1/3`, while its
   coefficient `2457/6592 ~ 0.37272` is *not* corner-clocked (transient
   hits only). Whether the true threshold `gamma*` relates to the odd
   ladder is open (HYP-9061).
3. **Constructions pay no parity toll at dyadic checkpoints** (T1); a
   `gamma < 1` cascade closing its books at `M = 2^r - 1` faces only
   archimedean/capacity obstructions there — the band capacity desert
   (lane D) remains the binding constraint.

All of T1-T3, C4, C5(i,ii) are refereed exactly (Lucas + independent
`F_2` convolution cross-check, exhaustive boxes plus randomized sweeps and
hostile controls) in the companion script. QED.
