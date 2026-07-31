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
  1.5949, 1.5962, 1.5970 for m = 4..4096. Hence C* > 1.597 for balanced block
  schemes -- a real improvement on the classical 3/2. Rescaling k = xm,
  d = (delta)m and taking log_2/m, (ARCH) survives iff
  H(delta) <= max_x [alpha H(r/alpha) + alpha - r], whose least admissible
  slope is C_arch = 1.5979874..., attained at delta* = 0.6179 (within grid
  1/phi EXACTLY). SOLVING THE STATIONARITY SYSTEM IN CLOSED FORM:
  delta* = 1/phi and p* = r/alpha = 1/sqrt5 exactly, the denominator of
  gamma* collapses to (1/2)log_2 5, and
      gamma* = log_2(phi) / ((1/2) log_2 5) = log_5(phi^2) = log_{sqrt5}(phi),
      C_arch = 1 + log_5(phi^2) = log_5(5 phi^2) = 1.5979874356654401497...
  The inner argmax slides from the profile kink x = kappa down to x = 0 as
  delta grows, so no one-variable reduction is available. The finite-m
  certified bounds converge to exactly this constant.
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
m      256    512    1024   2048   4096
rho >  1.5887 1.5925 1.5949 1.5962 1.5970
```

The gap to `C_arch = 1.5979874...` reads `0.00929, 0.00551, 0.00311, 0.00176,
0.00099`: it roughly halves per doubling, so the certified sequence converges
to the closed form of section 3.1, confirming it independently to four
digits.

Two checks: the bound never rejects a profile known to be feasible (it
passes at the exact optima for `m = 4,8,16,32`), and `rho_LB(m)` closely
tracks the *exact* `rho(m/2)` -- so (ARCH) is nearly tight, off by about one
doubling.

**Corollary.** `C* > 1.5970` for balanced block schemes (from `m = 4096`
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

### 3.1 Closed form

At the threshold three conditions hold at once: the max over `x` is interior
(`dg/dx = 0`), capacity touches requirement (`g = H(delta)`), and the touch
is tangential (`dg/ddelta = H'(delta)`). Tangency alone gives

```text
(1-p)/p = 2(1-delta)/delta,   i.e.   p = delta/(2-delta),
```

and the `x`-stationarity, using `log_2((1-p)/p) = 1 + H'(delta)`, gives

```text
gamma = -H'(delta) / [ H(p) + (1 + H'(delta))(1-p) ].
```

**The remaining equation is a quadratic** (this replaces the 40-digit
root-find of the first draft; the whole chain is now algebraic).

*Lemma.* The closure `alpha_1 = alpha_2` is equivalent to

```text
H(delta) = -H'(delta) (2 - delta).                                (CLO)
```

*Proof.* With `D = H(p) + (1+H')(1-p)` and `gamma = -H'/D`,

```text
1/gamma + 1 - p = [D + (1-p)(-H')]/(-H')
                = [H(p) + (1-p)(1 + H' - H')]/(-H')
                = [H(p) + 1 - p]/(-H'),
```

so `alpha_2 = (2-delta)(-H')/[H(p)+1-p]` while
`alpha_1 = H(delta)/[H(p)+1-p]`. The bracket cancels. QED

*Lemma.* (CLO) is equivalent to `delta^2 = 1 - delta`.

*Proof.* Put `L = log_2 delta`, `M = log_2(1-delta)`. Then
`H = -delta L - (1-delta) M` and `-H' = L - M`, so (CLO) reads
`-delta L - (1-delta)M = (2-delta)(L-M)`. Collecting,
`-2L = -M`, i.e. `2L = M`, i.e. `delta^2 = 1-delta`. QED

The unique root in `(1/2,1)` is `delta* = (sqrt5-1)/2 = 1/phi`. Then

```text
2 - delta* = (5 - sqrt5)/2 = sqrt5 (sqrt5-1)/2,
p*  = delta*/(2-delta*) = 1/sqrt5      (exactly),
-H'(delta*) = log_2(delta*/(1-delta*)) = log_2(1/delta*) = log_2 phi.
```

So `p* = 1/sqrt5` follows *algebraically* from `delta* = 1/phi`:
`delta/(2-delta) = ((sqrt5-1)/2)/((5-sqrt5)/2) = 1/sqrt5`. With
`delta* = 1/phi` one has `H'(1/phi) = -log_2 phi`, and the denominator
collapses exactly:

```text
H(1/sqrt5) + (1 - log_2 phi)(1 - 1/sqrt5)
   = [F + (1-s)(L-1)] + (1-L)(1-s) = F = (1/2) log_2 5,
   F = (1/2)log_2 5,  L = log_2 phi,  s = 1/sqrt5.
```

Hence

```text
gamma* = log_2(phi) / ((1/2)log_2 5) = 2 log_5 phi = log_5(phi^2)
       = log_{sqrt5}(phi),

C_arch = 1 + log_5(phi^2) = log_5(5 phi^2) = log_{sqrt5}(sqrt5 * phi)
       = 1.59798743566544014974502650205...
```

with `5 phi^2 = (15+5 sqrt5)/2` and `sqrt5 * phi = (5+sqrt5)/2 = 2+phi`.
Every step above is an identity; the numerics (51 digits) are only a check.

### 3.1a General alphabet: only q = 2 is golden

The `2` in (CLO) is the alphabet size (it enters through
`|Delta^j P(0)| <= 2^j`). Repeating the second lemma with `q` in place of `2`
gives `-qL = (1-q)M`, i.e.

```text
delta^q = (1-delta)^(q-1).
```

`q = 2` is the golden quadratic `delta^2 = 1-delta`; no other `q` gives a
metallic equation. Roots: `0.6180339887 (=1/phi), 0.5698402910,
0.5497004779, 0.5385972572` for `q = 2,3,4,5`, matching the independent
full-system numerics of section 9.3 to every printed digit (and resolving its
`q=3` "no root", which was a root-finder artifact). This is the *reason*
behind the metallic-ratio negative: the golden ratio is not the first member
of a metallic family here, it is the `q=2` member of the family
`delta^q = (1-delta)^(q-1)`.

**The extremal profile is therefore a Sturmian (Beatty) sequence of slope
`log_5(phi^2)`:** `a_k = min(m-1-k, floor(gamma*(m+k)))`.

### 3.2 Numerics

```text
C_arch = 1 + gamma*,     gamma* = 0.5979874356654401...,
C_arch = 1.5979874356654401...,     binding delta* = 1/phi.
```

The inner argmax is at the profile kink
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
`(1.596, 2]`, and every indication is that it equals
`C_arch = log_5(5 phi^2)`. What remains open is the matching construction -- (ARCH) is a capacity condition,
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

QED for reductions A and B, for the finite-`m` bounds, and -- since
section 3.1 -- for the closed form of the asymptotic constant.

## 10. Exactly what is still missing for `C* >= C_arch`

Note first that **no asymptotics are needed for a rigorous bound**: (ARCH) at
any single `m` is a finite exact-integer computation, and `m = 4096` already
certifies `C* > 1.597`. What the limit statement `C* >= C_arch` still needs is
only:

1. **The scaling limit.** `log_2 binom(m-1, delta m)/m -> H(delta)` and
   `log_2[binom(a,r) 2^(a-r)]/m -> alpha H(r/alpha) + (alpha - r)` uniformly
   on compacts -- routine Stirling, but it must be stated with the error term,
   because the sum over `k` has `O(m)` terms and the per-term polynomial
   factors must not accumulate.
2. **Interiority.** That the binding `delta` is interior (not at `0` or `1`)
   and the inner `max_x` is attained in the interior of the admissible range,
   so that the three stationarity conditions of section 3.1 are the correct
   optimality system rather than one branch of a boundary case. The inner
   argmax is known to slide from `x = kappa` to `x = 0` as `delta` grows
   (section 3.2), so this needs a genuine argument, not an appeal to
   smoothness.

Neither is deep, but neither is written. Everything else in the chain --
(ARCH) itself, the reduction to (CLO), and the evaluation
`(CLO) => delta^2 = 1-delta => C_arch = log_5(5 phi^2)` -- is now proved.


## 7. The construction: what fails, and why (recorded negative)

The natural candidate is a **division ladder**. Because `L_0 > L_1 > ...`
along the constrained strata and `L_k = 0` afterwards, one can set
`R := eps u^(m-1)` and repeatedly take `E_k :=` quotient of `R` by
`(1+u)^(L_k)`, `R :=` remainder. This terminates exactly and the quotient
degrees fit (`deg E_k <= a_k` automatically). Two things go wrong.

1. **Wrong parity class.** The raw ladder returns `E_0 = eps u^(m-1)` and
   `E_k = 0` otherwise, which violates `[u^i]E_k = binom(a_k,i) (mod 2)`.
   This is repairable: by reduction A the canonical parity representative
   `F_k(u) = sum_{i subset a_k} u^i` satisfies
   `sum_k F_k(u)(1+u)^(L_k) = u^(m-1) (mod 2)`, so
   `W := [eps u^(m-1) - sum_k F_k(1+u)^(L_k)]/2` has integer coefficients and
   one can run the ladder on `W` and set `e = f + 2g`.

2. **It violates the boxes at the corners, catastrophically.** The repaired
   ladder is exact and parity-correct but overshoots at `i = 0` (and `i = 1`),
   where the box `binom(a_k,0) = 1` is tightest. Measured overshoot factor
   at the first violated slot:

   ```text
   m         4      8       16          32              64            128
   C=8/5    3      21      3087        1.4e8           1.5e17        6.2e35
   C=2      3       7        15          31              63             127
   ```

   It succeeds only at `m = 4` (`C = 7/4, 15/8`).

The failure is structural and diagnostic: dividing from the top pushes mass
into the LOW-order coefficients, exactly where the binomial box is smallest,
whereas the boxes are enormous in the middle (`i ~ a_k/2`). Any successful
construction must keep the deviation mass near the middle of each stratum --
which is what the exact optima found for `m <= 32` do. So the outstanding
problem is a **middle-weighted** decomposition of `eps u^(m-1)` into
`sum_k E_k(u)(1+u)^(L_k)` with `E_k` dominated coefficientwise by
`(1+u)^(a_k)`; note `u^c (1+u)^(a-c)` is always legal, since
`binom(a-c, i-c) <= binom(a,i)`, so such blocks are the natural atoms.
Referee: `04-computation/amm12592_division_ladder_construction_attempt.py`.


## 8. The Catalan substitution, and what the true atoms are

### 8.1 The substitution kills the L_k

Put `w = u/(1+u)`, so `1+u = 1/(1-w)`, `u = w/(1-w)` -- the Catalan/Riordan
substitution. Write `E_k(u) = (1+u)^(a_k) Lam_k(w)`. Since
`a_k + L_k = m-1-k` identically,

```text
E_k(u)(1+u)^(L_k) = (1+u)^(m-1-k) Lam_k(w) = Lam_k(w)(1-w)^(-(m-1-k)),
u^(m-1) = w^(m-1) (1-w)^(-(m-1)),
```

so multiplying the system by `(1-w)^(m-1)` gives

```text
sum_(k=0)^(m-1) Lam_k(w) (1-w)^k = eps w^(m-1),   deg Lam_k <= a_k.   (W)
```

**The `L_k` disappear entirely.** (W) is a `(1-w)`-adic digit expansion of
`eps w^(m-1)` with per-digit degree bounds -- a much cleaner object than the
`u`-side system.

### 8.2 Why the golden ratio was inevitable here

The Catalan generating function `C(w) = (1 - sqrt(1-4w))/(2w)` satisfies

```text
C(-1)              = (1-sqrt5)/(-2) = (sqrt5-1)/2 = 1/phi = delta*,
sqrt(1-4w)|_(w=-1) = sqrt5                        = 1/p*.
```

Both threshold constants of section 3.1 are Catalan-GF data at `w = -1`,
which is the natural evaluation point of (W) (`1-w = 2` there). The user's
sequences `C(2n,n-1) = 1,4,15,56,210`, `C(2n+1,n-1) = 1,5,21,84,330`,
Catalan `1,2,5,14,42,132` and the central binomials are the Catalan-triangle
entries that (W) manufactures.

### 8.3 The exact invariant: Lam_k(1) = +-1

Evaluating (W) and its derivatives at `w = 1` gives `Lam_0(1) = eps` and the
carry recursion `R_(k+1)(1) = Lam'_k(1) - R'_k(1)`, with `Lam_k(1) = R_k(1)`
forced. Moreover the box constraint at the top coefficient `i = a_k` reads
exactly `|Lam_k(1)| <= 1`. Extracting `Lam_k` from the verified optima
confirms this with no exceptions:

```text
m = 8 :  Lam_k(1) = +1,-1,+1,-1,-1,-1,+1,+1
m = 16:  Lam_k(1) = -1,+1,... all +-1 for every k
```

The `i = a_k - 1` constraint then reproduces `a_0 >= (m-1)/2`, i.e.
THM-2160 S6.2, for a third independent time.

### 8.4 The unit-atom construction, and why it is not enough

If the positive coefficients of `Lam_k` sum to at most `1` and the negative
ones to at least `-1`, then `|[u^i]E_k| <= binom(a_k,i)` automatically
(because `0 <= binom(a_k-c, i-c) <= binom(a_k,i)`). Over the integers those
are exactly the **middle-weighted unit atoms**

```text
Lam_k  in  { 0,  +- w^c,  w^c - w^(c') },    c, c' <= a_k,
```

and (W) becomes a carry ladder whose digit is the exponent `c`
(`Lam'_k(1) = sigma c`, or `c - c'`). This is implemented and VERIFIED, but
the one-sided digit range (when `Lam_k(1) = +-1` only `Lam'_k(1)` of one sign
is reachable) makes the greedy myopic: it attains `C = 7/4` at `m = 4`,
`15/8` at `m = 8` and only `C = 2` from `m = 16` on.

### 8.5 What the true atoms actually are

Extracting `Lam_k` from the exact optima shows the unit-atom family is far
too narrow -- the real solutions use massive cancellation (`m = 16, k = 4`
has positive part `5126` and negative part `5127`). Their shape is instead

```text
Lam_4 = -1 + 20 w(1-w)^9 + (small corrections),
```

and since `(1+u)^a * w(1-w)^(a-1) = u` identically, this is

```text
E_k = sigma_k (1+u)^(a_k) + (low-degree monomial corrections),
```

i.e. **the full binomial -- the maximally middle-weighted object -- plus
adjustments in the low-order coefficients, where the box `binom(a_k,i)` still
has room.** The legality condition is just
`|sigma_k binom(a_k,i) + d_{k,i}| <= binom(a_k,i)`.

So the correct atom family for a matching construction is
`sigma(1+u)^a + D(u)` with `D` low-degree, NOT the unit atoms
`u^c(1+u)^(a-c)`; the latter are legal but cannot carry the required mass.
Building the ladder over this larger family is the remaining step.

Referees: `04-computation/amm12592_catalan_w_construction.py` (the reduction
(W), the unit-atom ladder, and its verification) and
`04-computation/amm12592_extract_lambda_atoms.py` (the `Lam_k` of the exact
optima).


## 9. The full atom family, and a metallic-ratio negative

### 9.1 Domination is a modulation; alternation saturates

`E_k` is dominated by `(1+u)^(a_k)` exactly when

```text
[u^i] E_k = binom(a_k,i) * P_k(a_k - i),      |P_k| <= 1,
```

so the legal atoms are the binomial profile modulated by any function
bounded by `1`. The two extremes are `P_k = +-1`, giving
`E_k = +-(1+u)^(a_k)` (`Lam_k = +-1`), and `P_k(y) = +-(-1)^y`, giving
`E_k = +-(u-1)^(a_k)` (`Lam_k = +-(2w-1)^(a_k)`) -- the **alternating** atom.
Since `Lam_k(1) = P_k(0)` and `Lam'_k(1) = -a_k(P_k(1)-P_k(0))`, one has
`|Lam'_k(1)| <= 2 a_k` with equality exactly at the alternating atom: twice
the reach of the unit atoms `u^c(1+u)^(a-c)`. It is `|Delta^j P(0)| <= 2^j`,
saturated by alternation, that produces (ARCH). So alternation is the
extremal atom, and the unit atoms of section 8.4 are the wrong family.

### 9.2 What the ladder over the full family achieves

Taking `P_k` integer-valued in `{-1,0,1}` makes every coefficient integral
automatically and turns each level into a signed subset-sum with coins
`binom(a_k, d-k)` over `{-1,0,+1}`. This is implemented and verified and
improves the small cases from `7/4, 15/8` to

```text
m = 4 : C = 8/5     m = 8 : C = 8/5     m >= 16 : only C = 2.
```

Allowing general integer `e_{k,i}` returns exactly the triangular solver of
THM-3008 (fresh variables always enter with coefficient `1`), now run with an
**alternating-extreme** splitting policy in place of the centred one. That
recovers the exact optima `3/2, 14/9, 25/16` for `m <= 16` and gives
`52/33` at `m = 32`, `109/66` at `m = 64` -- still degrading.

**The uniform construction at `C < 2` remains OPEN.** What is now known is
where it must live: at the alternating extreme of the box, not near its
centre, and the obstruction to the greedy is that the level residual must be
split among fresh unit-coefficient variables without destroying later levels.

### 9.3 Metallic ratios: refuted along the alphabet axis

The `2` in the tangency relation `p = delta/(2-delta)` traces to
`|Delta^j P(0)| <= 2^j`, i.e. to the binary alphabet. Replacing it by `q` and
re-solving the stationarity system gives

```text
q = 2:  delta* = 0.6180339887 = 1/phi,     1/p* = 2.2360679 = sqrt5
q = 3:  no root
q = 4:  delta* = 0.5497004779,             1/p* = 6.2766900
q = 5:  delta* = 0.5385972572,             1/p* = 8.2833744
```

None of the `q > 2` values is `1/x` for a metallic ratio `x^2 = nx+1`, and
none of the `1/p*` is `sqrt(n^2+4)`. **The golden ratio here is not the
`n = 1` member of a metallic family in this parameter**; it is specific to
the binary alphabet, which is exactly what the Catalan reading of section 8.2
predicts (`sqrt(1-4w)` at `w = -1`). Referee:
`04-computation/amm12592_metallic_generalization_test.py`.
