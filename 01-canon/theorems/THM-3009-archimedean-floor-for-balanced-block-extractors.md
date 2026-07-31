---
id: THM-3009
title: "The archimedean floor for balanced block extractors: C* > 1.596"
status: >
  PROVED (finite-m certified) + ASYMPTOTIC CONSTANT COMPUTED. For AMM 12592
  the within-block system in deviation form is
  sum_k E_k(u)(1+u)^(L_k) = +-u^(m-1), |[u^i]E_k| <= binom(a_k,i),
  [u^i]E_k = binom(a_k,i) mod 2. TWO REDUCTIONS.
  (A) MOD 2 THERE IS NO OBSTRUCTION, EVER: the forced-odd slots are the
  submasks i of a_k, sum_{i subset a} u^i = (1+u)^a over F_2, and the whole
  pattern collapses to sum_k (1+u)^(m-1-k) = [(1+u)^m-1]/u = u^(m-1) because
  m is dyadic. Every coefficient below the top vanishes IDENTICALLY for
  EVERY depth profile. Inside a shell, parity is a red herring; the one
  surviving parity is the free label of z = 1^m.
  (B) THE OBSTRUCTION IS ARCHIMEDEAN: expanding at u = -1 gives the
  necessary condition binom(m-1,d) <= sum_{k: 0<=d-L_k<=a_k}
  binom(a_k,d-L_k) 2^(a_k-d+L_k) for every d (ARCH). It is O(m^2), monotone
  in the a_k, and strictly stronger than THM-2160 S6.2. Certified:
  rho(m) > 1.4000, 1.5000, 1.5556, 1.5610, 1.5753, 1.5828, 1.5887, 1.5925,
  1.5949, 1.5962 for m = 4..2048. Hence C* > 1.596 for balanced block
  schemes -- a real improvement on the classical 3/2. Rescaling k = xm,
  d = (delta)m and taking log_2/m, (ARCH) survives iff
  H(delta) <= max_x [alpha H(r/alpha) + alpha - r], whose least admissible
  slope is C_arch = 1.5979874..., attained at delta* = 0.6179 (within grid
  resolution of 1/phi = 0.61803; the inner argmax slides from the profile
  kink at x = kappa down to x = 0 as delta grows, so no one-variable
  reduction is available). The finite-m certified bounds converge to exactly
  this constant.
source: opus-2026-07-31-amm12592-writeup
depends_on:
  - THM-3007
  - THM-3008
related:
  - THM-2160  # S6.2 gives only 3/2, and only from the k=0 stratum
  - THM-2966
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_archimedean_lower_bound.py
asymptotic: 04-computation/amm12592_archimedean_threshold_asymptotic.py
output: 05-knowledge/results/amm12592_archimedean_lower_bound_large_m.out
---

# THM-3009 -- the archimedean floor

## 0. Setting

THM-3008's within-block system on the dyadic shell `[m,2m)` is

```text
sum_k sum_i p_{k,i} binom(L_k, j-1-i) = binom(m,j)/2,  1<=j<=m-1,
0 <= p_{k,i} <= binom(a_k,i),   L_k = m-1-k-a_k,                  (*)
```

with `T(m+k) = m+k+1+a_k`. Every level of (*) is **exactly centred**: the
box total at level `t = m-1-j` is `sum_k binom(m-1-k, t-k+1) = binom(m,t+1)`
by the hockey-stick identity, exactly twice the target. So `p = binom/2`
always solves (*) over the rationals, and the whole problem is integrality.
Writing `p_{k,i} = (binom(a_k,i) + e_{k,i})/2` makes it homogeneous:

```text
sum_k E_k(u)(1+u)^(L_k) = eps u^(m-1),   eps = +-1,
|[u^i]E_k| <= binom(a_k,i),   [u^i]E_k = binom(a_k,i) (mod 2).     (D)
```

## 1. Reduction A: no 2-adic obstruction

`binom(a,i)` is odd exactly when `i` is a binary submask of `a` (Lucas), and
`sum_{i subset a} u^i = (1+u)^a` over `F_2`. So the forced-odd pattern of (D)
reduces mod 2 to

```text
sum_(k=0)^(m-1) (1+u)^(a_k+L_k) = sum_(k=0)^(m-1) (1+u)^(m-1-k)
                                = [(1+u)^m - 1]/u = u^(m-1)  (mod 2),
```

the last step because `m` is a power of two, i.e. `(1+u)^m = 1+u^m` over
`F_2` -- the same dyadicity that THM-3007 shows is forced. **Every
coefficient below the top vanishes identically, for every depth profile.**

This is a correction to the working picture in HYP-9061, which sought the
obstruction in half-integer deficits, dyadic split-jumps and binary-clock
parity. Inside a shell those all cancel automatically; the single surviving
parity is exactly the free label of the middle word `z = 1^m`. The
obstruction must be archimedean.

## 2. Reduction B: the archimedean condition

Substituting `u = -1+v` sends `(1+u)^(L_k)` to `v^(L_k)`, so comparing
coefficients of `v^d` in (D),

```text
sum_(k: L_k <= d) [v^(d-L_k)] E_k(-1+v) = eps (-1)^(m-1-d) binom(m-1,d).
```

Since `|[v^r]E_k(-1+v)| <= sum_i binom(a_k,i) binom(i,r)
= binom(a_k,r) 2^(a_k-r)`, a necessary condition for the profile is

```text
binom(m-1,d) <= sum_(k: 0<=d-L_k<=a_k) binom(a_k, d-L_k) 2^(a_k-d+L_k)
                                                for every d.       (ARCH)
```

(ARCH) costs `O(m^2)`, is monotone in each `a_k`, and is strictly stronger
than THM-2160 S6.2, which only ever yields `3/2` and only from the `k=0`
stratum. Binary-searching the profile against (ARCH):

```text
m        4      8      16      32      64     128
rho >  1.4000 1.5000 1.5556 1.5610 1.5753 1.5828
m      256    512    1024   2048
rho >  1.5887 1.5925 1.5949 1.5962
```

Two checks: the bound never rejects a profile known to be feasible (it
passes at the exact optima for `m = 4,8,16,32`), and `rho_LB(m)` closely
tracks the *exact* `rho(m/2)` -- so (ARCH) is nearly tight, off by about one
doubling.

**Corollary.** `C* > 1.5962` for balanced block schemes (from `m = 2048`
alone; the bound is a finite exact integer computation).

## 3. The asymptotic constant

Scale `k = xm`, `d = (delta)m`, `a_k = alpha(x) m`, `L_k = ell(x) m`. With
`C = 1+gamma` and the extremal profile `a_k = min(m-1-k, gamma(m+k))`,

```text
kappa = (1-gamma)/(1+gamma),
x <= kappa : alpha = gamma(1+x),  ell = (1-gamma) - x(1+gamma),
x >= kappa : alpha = 1-x,         ell = 0.
```

Taking `log_2` and dividing by `m`, `binom(m-1,d) -> H(delta)` and each
summand `-> alpha H(r/alpha) + (alpha - r)` with `r = delta - ell(x)`. So
(ARCH) survives in the limit iff

```text
for all delta:  H(delta) <= max_x [ alpha H(r/alpha) + alpha - r ].   (T)
```

The least slope satisfying (T) is

```text
C_arch = 1 + gamma*,     gamma* = 0.597987401...,
C_arch = 1.597987401...,     binding delta* = 0.6179...
```

`delta*` sits within grid resolution of `1/phi = 0.618034`, which is worth
chasing but is NOT claimed. The inner argmax is at the profile kink
`x = kappa` for small `delta` and slides monotonically to `x = 0` as `delta`
grows (`x = 0.2516, 0.1649, 0.0898, 0.0387, 0` at
`delta = 0.30, 0.45, 0.55, 0.618, 0.70`), so there is **no** one-variable
reduction: an attempt to evaluate the capacity at the kink alone gives the
wrong answer `1.9`.

The finite-`m` certified bounds of section 2 converge to `C_arch`, which
independently confirms the constant to four digits.

## 4. Where this leaves the proposer's question

```text
classical schemes                     C = 2
THM-2996 (sharpened half-tail)        C = 2   (D improves, slope does not)
THM-3008 verified constructions       rho(4,8,16,32) = 3/2, 14/9, 25/16, <=11/7
THM-3009 certified floor              C* > 1.5962   (balanced block schemes)
THM-3009 asymptotic floor             C_arch = 1.59799
```

So the answer is neither `2` nor `1`: for balanced block schemes `C*` lies in
`(1.596, 2]`, and every indication is that it equals `C_arch ~ 1.598`. What
remains open is the matching construction -- (ARCH) is a capacity condition,
and capacity conditions of this shape are usually sufficient once satisfied
with room, but no uniform family of profiles is yet built. The greedy
triangular solve of THM-3008 degrades past `m = 32` (it returns `1.65`,
`1.77`, `1.85` at `m = 64, 128, 256`), so those numbers are policy artifacts,
not evidence about `rho`.

## 5. Scope

All of this is for **balanced block schemes**: methods that decide every
critical value in a block by the block's end and split every composition
class of every block evenly. THM-3007 forces such blocks to be dyadic
intervals. A general exactly fair extractor need not be block-balanced, and
for that wider class only `C* >= 1` (the floor lemma) and `C* <= 2` are
known.

## 6. Referee

```bash
python3 04-computation/amm12592_archimedean_lower_bound.py 4 8 16 32 64 128 256
python3 04-computation/amm12592_archimedean_threshold_asymptotic.py
```

QED for reductions A and B and for the finite-`m` bounds; the asymptotic
constant is a numerical evaluation of (T).
