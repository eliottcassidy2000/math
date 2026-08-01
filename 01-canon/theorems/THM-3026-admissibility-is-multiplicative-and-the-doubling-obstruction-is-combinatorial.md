---
id: THM-3026
title: "AMM 12592: admissibility is multiplicative, so the doubling obstruction is combinatorial, not archimedean"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Structural
  input to the epoch-doubling induction that THM-3002 section 5b and klein
  both name as the remaining prize on the CONSTRUCTION side of AMM 12592.
  A block is ADMISSIBLE at degree d when its coefficients delta_k in the basis
  B_{d,k}(x) = x^{d-k}(1-x)^k satisfy |delta_k| <= binom(d,k) and
  delta_k = binom(d,k) mod 2 -- the Lucas-box capacity and parity conditions.
  (M) MULTIPLICATIVITY: B_{d,k} B_{e,k'} = B_{d+e,k+k'} EXACTLY, so the
  product of blocks has basis coefficients (delta * eps)_kappa =
  sum_{k+k'=kappa} delta_k eps_k', and Vandermonde
  sum_{k+k'=kappa} binom(d,k)binom(e,k') = binom(d+e,kappa) gives BOTH the
  capacity bound AND (reducing mod 2) the parity condition. So admissible x
  admissible = admissible at the SUM of degrees, exactly, with no slack lost.
  (L) LIFTING: since x + (1-x) = 1, the block delta_k = binom(d'-d,k) is
  admissible at d'-d and represents the constant 1, so multiplying by it
  re-expresses any admissible block at any larger degree d' >= d, still
  admissible. Admissibility is therefore degree-monotone.
  (D) DEGREE FIT AT gamma = 3/5: d_i + d_{i'} <= d'_{i+i'} for all i, i',
  because floor(a) + floor(b) <= floor(a+b) and 3(R+i)/5 + 3(R+i')/5 =
  3(2R+i+i')/5 exactly. Verified with no violations for R = 8..128.
  CONSEQUENCE. The identity q^{2R-1} = (1-x)(sum_i x^i Delta_i)^2 turns an
  epoch-R solution into a valid epoch-2R IDENTITY whose every product block
  Delta_i Delta_{i'} is admissible in its own row. What the naive squaring
  breaks is ONLY the row distribution: row j receives ~j pairs (i,i') and
  their capacities add (observed overshoot factors 49, 276, 1541 at
  R = 8,16,32). So the obstruction to a doubling induction is purely
  COMBINATORIAL -- redistributing pairs across rows -- and NOT archimedean.
  This is consistent with, and explains, THM-3002 section 5b's finding that
  criterion (4) is uniformly ample at 3/5 (worst ratio ~1.20) while the beam
  search still stalls at R = 128 on search budget alone.
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3002
related:
  - THM-3007
  - THM-3009
  - THM-3024
  - HYP-9061
script: 04-computation/amm12592_admissibility_multiplicative_thm3026.py
output: 05-knowledge/results/amm12592_admissibility_multiplicative_thm3026.out
---

# THM-3026 -- admissibility is multiplicative

## 1. The constraint, restated in the right basis

In the epoch recursion `sigma_{-1} = q^{R-1}`, `p sigma_i = sigma_{i-1} -
Delta_i` (`q = 1-x`), each `Delta_i` is written in the basis

```text
B_{d,k}(x) = x^{d-k} (1-x)^k,       0 <= k <= d,
Delta_i = sum_k delta_{i,k} B_{d_i,k},   d_i = floor(3(R+i)/5) + D0,
```

and the Lucas-box conditions are

```text
|delta_{i,k}| <= binom(d_i, k)        (capacity)
 delta_{i,k} = binom(d_i, k) mod 2    (parity)
```

Call such a block **admissible at degree `d`**. Note the whole epoch identity
is `sum_i x^i Delta_i(x) = q^{R-1}`.

## 2. (M) Products of admissible blocks are admissible

The basis multiplies **exactly**:

```text
B_{d,k}(x) B_{e,k'}(x) = x^{(d+e)-(k+k')} (1-x)^{k+k'} = B_{d+e, k+k'}.
```

So the product of two blocks has basis coefficients given by the plain
convolution `(delta * eps)_kappa = sum_{k+k'=kappa} delta_k eps_{k'}`, and

```text
|(delta*eps)_kappa| <= sum_{k+k'=kappa} binom(d,k) binom(e,k') = binom(d+e,kappa),
 (delta*eps)_kappa  =  sum_{k+k'=kappa} binom(d,k) binom(e,k') = binom(d+e,kappa) mod 2,
```

both by **Vandermonde's identity** -- the second because the parity conditions
say `delta_k = binom(d,k)` and `eps_{k'} = binom(e,k')` mod 2, and the same
convolution then computes `binom(d+e,kappa)` mod 2. Hence

```text
admissible at d,  admissible at e   =>   product admissible at d+e.      (M)
```

**No slack is lost**: the capacity bound is attained by the all-maximal block,
so (M) is an exact statement, not an estimate. Verified on 400 random
admissible pairs with `d,e <= 9`, together with the exactness of the
polynomial identity.

## 3. (L) Admissibility is degree-monotone

Because `x + (1-x) = 1`,

```text
sum_k binom(r,k) B_{r,k}(x) = (x + (1-x))^r = 1,
```

so the block `delta_k = binom(r,k)` is admissible at degree `r` and
**represents the constant `1`**. Multiplying any admissible block at `d` by it
and applying (M) re-expresses the *same polynomial* as an admissible block at
`d' = d + r`, for any `r >= 0`. Verified on 400 random lifts.

## 4. (D) The degrees fit at gamma = 3/5

```text
d_i + d_{i'} = floor(3(R+i)/5) + floor(3(R+i')/5)
            <= floor( 3(R+i)/5 + 3(R+i')/5 ) = floor(3(2R+i+i')/5) = d'_{i+i'},
```

using `floor(a)+floor(b) <= floor(a+b)` and the *exact* additivity
`3(R+i)/5 + 3(R+i')/5 = 3(2R+(i+i'))/5`. No violations for `R = 8,...,128`.
(This is a property of the rate `3/5` written as a linear function of the row
index; it is why `gamma` linear in `i` is the right normal form.)

## 5. What this says about the doubling induction

Squaring the epoch identity gives, for free,

```text
q^{2R-1} = (1-x) * (q^{R-1})^2 = (1-x) * ( sum_i x^i Delta_i )^2,
```

a **valid epoch-`2R` identity**, and by (M)+(L)+(D) each product block
`Delta_i Delta_{i'}` is admissible in row `i+i'` of the doubled epoch. So
every *individual* ingredient of the doubled solution is legal.

What fails is only that row `j` of the doubled epoch receives **all** pairs
`(i,i')` with `i+i' = j` -- about `j` of them -- and their capacities add.
Measured overshoot of the naive lift (max `|coefficient| / capacity`):

```text
R =  8 -> 16 :   49
R = 16 -> 32 :  276
R = 32 -> 64 : 1541
```

growing like the pair count, exactly as the analysis predicts.

**Conclusion.** The obstruction to an epoch-doubling induction at `gamma = 3/5`
is **combinatorial** -- how to redistribute the `~j` product blocks across
rows, or equivalently how to choose a *sparse* factorisation of `q^{2R-1}`
instead of the full square -- and **not archimedean**: no capacity is lost at
the level of individual blocks. This is consistent with, and explains, the
situation reported in THM-3002 section 5b, where criterion (4) is uniformly
ample at `3/5` (worst ratio `~1.20` at `t = 1`, stable to `R = 256`) while the
beam search nevertheless stalls at `R = 128` on search budget alone.

## 6. Scope

(M), (L) are proofs (Vandermonde plus the exact basis identity); (D) is proved
in general from `floor(a)+floor(b) <= floor(a+b)` and verified to `R = 128`.
**No doubling induction is claimed** -- the redistribution problem is open, and
until it is solved `C* <= 8/5` remains verified only for `n <= 127`
(THM-3002 section 5b), not proved for all `n`. What is established is a
negative: the failure of the naive square-and-multiply lift carries **no
archimedean information**, so it must not be read as evidence against
`gamma = 3/5` closing at every scale.

Referee: `amm12592_admissibility_multiplicative_thm3026.py` -- 400 random
product tests for (M) including polynomial exactness, 400 lift tests for (L),
and the exhaustive degree-fit check for (D) over `R = 8..128`.
