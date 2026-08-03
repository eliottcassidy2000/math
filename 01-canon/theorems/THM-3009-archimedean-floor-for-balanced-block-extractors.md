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
  AUDIT NOTE 2026-08-03 (boxeph,
  05-knowledge/results/amm12592-golden-floor-audit-boxeph.md): independent
  exact re-referee PASSED — full (ARCH) ladder re-derived from scratch for
  m <= 1024 (all values match), m = 2048/4096 refuted candidates re-verified,
  ladder monotonicity swept exhaustively at m = 8,16,32, and every refuted
  rate certified below gamma* by Fibonacci/Lucas integer comparison
  (certified bracket 115939/193882 < gamma* < 105183/175895). Residual debt:
  the sec 10.3 Stirling transfer lemma is the one unwritten analytic step of
  the asymptotic floor; after the audit's concavity/gamma-monotonicity
  upgrades the delta-scan is load-bearing only for threshold
  exactness-from-above. Scope: balanced-block class; the general-class
  promotion attempt is demoted (THM-3024, MISTAKE-361).
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
THM-3032 (sharpened half-tail)        C = 2   (D improves, slope does not)
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

## 10. The two remaining gaps, closed

### 10.1 Reparametrisation: eliminate x in favour of ell

Both gaps become tractable after replacing the stratum coordinate `x` by
`ell = L_k/m`. On the constrained branch `ell = (1-gamma) - x(1+gamma)`, so

```text
alpha = gamma(2-ell)/(1+gamma),   ell in [0, 1-gamma],
```

and the free branch is `ell = 0`, `alpha in [0, 2gamma/(1+gamma)]`; the two
meet at the junction `ell = 0`, `alpha = 2gamma/(1+gamma)`. With
`r = delta - ell` and `p = r/alpha`, `g = alpha H(p) + alpha - r`.

### 10.2 Interiority

**(I1) The free branch is dominated.** `dg/dalpha = H(p) - p H'(p) + 1
= 1 - log_2(1-p) > 0` for `p in [0,1)`. So on the free branch `g` increases
with `alpha` and is maximal at the junction. The max over `x` is therefore a
max over `ell in [0, min(delta, 1-gamma)]`.

**(I2) The `r = 0` endpoint is never optimal.** On the constrained branch

```text
dg/dell = -(gamma/(1+gamma))(1 - log_2(1-p)) - H'(p) + 1,
```

and as `ell -> delta` we have `p -> 0`, `H'(p) -> +infinity`, so
`dg/dell -> -infinity`.

**(I3) Hence the maximiser is the junction or interior; at the threshold it
is interior.** With `alpha* = H(delta*)/(H(p*)+1-p*)` and
`ell* = delta* - p* alpha*`,

```text
ell*        = 0.34027368653552164201,
1 - gamma*  = 0.40201256433455985025,
```

so `0 < ell* < 1-gamma*` with margins `0.3403` and `0.0617` -- strict
inequalities between explicit constants, with room to spare.

**Global minimality of `delta*`.** Scanning the deficiency
`Phi(delta) = max_x g - H(delta)` at `gamma = gamma*` over
`delta in [0.30, 0.95]` gives `Phi >= 0` throughout with minimum
`4.6e-5` at the grid point nearest `1/phi`: the tangency is at `1/phi` and
`Phi` is positive elsewhere. (A fully rigorous statement needs a Lipschitz
certificate on this scan; the margins away from `delta*` are of order
`10^-2`, so the certificate is routine but is NOT yet written.)

### 10.3 The Stirling estimate, and the convergence rate it predicts

From `n H(j/n) - log_2(n+1) <= log_2 binom(n,j) <= n H(j/n)` and the fact
that (ARCH) has at most `m` summands,

```text
log_2 LHS >= (m-1) H(d/(m-1)) - log_2 m,
log_2 RHS <= log_2 m + max_k [ a_k H(r_k/a_k) + a_k - r_k ].
```

So (ARCH) **fails** as soon as

```text
H(delta) - max_x g(x,delta)  >  (2 log_2 m + O(log m))/m,
```

the `O(log m)` absorbing the floors in `a_k, L_k` (each shifts `a H(r/a)` by
`O(log m)`). Hence for every `gamma < gamma*` the slope-`(1+gamma)` profile is
infeasible for all large `m`, i.e.

```text
C* >= 1 + gamma* = C_arch = log_5(5 phi^2).
```

The same estimate **predicts the convergence rate** of the certified
finite-`m` bounds, `gamma* - gamma_m ~ c log_2(m)/m`. Measured:

```text
m        256      512     1024     2048     4096
gap    0.00934  0.00551  0.00311  0.00176  0.00099
c      0.2987   0.3135   0.3186   0.3281   0.3367
```

`c` is near-constant, confirming the error term is of the right order. The
slow upward drift is the expected quadratic-tangency correction: near
`gamma*` the deficiency behaves like `-c_1(gamma*-gamma) + c_2(delta-delta*)^2`
(a fold), so the negative window has width `~sqrt(gamma*-gamma)` and depth
`~(gamma*-gamma)`.

## 11. The three objects, worked

### 11.1 The deficiency fold DERIVES the convergence constant

Near the threshold `Phi_gamma(delta) = max_x g - H(delta)` is a fold,
`Phi ~ -A(gamma*-gamma) + (B/2)(delta-delta*)^2`, so (ARCH) fails exactly when
the depth `A(gamma*-gamma)` exceeds the Stirling error `E`, giving
`gamma* - gamma_m = E/A`. By the envelope theorem `A = dg/dgamma` at the
maximiser, and at fixed `ell` only `alpha = gamma(2-ell)/(1+gamma)` depends on
`gamma`, so with `dg/dalpha = 1 - log_2(1-p)`,

```text
A = (2 - ell*)/(1+gamma*)^2 * (1 - log_2(1-p*)) = 1.20581835346289864.
```

Sharpening the Stirling constant matters here. The central term is
`2^{nH}/sqrt(2 pi n p(1-p))`, and the `O(m)` summands of (ARCH) form a
Gaussian of width `~sqrt m`, so the sum costs `+1/2 log_2 m` against the
left side's `-1/2 log_2 m`: `E = (1/2)log_2(m)/m`, not the crude
`2 log_2(m)/m` of section 10.3. Hence the prediction

```text
c := (gamma* - gamma_m)/(log_2(m)/m) = 1/(2A) = 0.414656153279.
```

The crude bound would give `2/A = 1.6586`, four times too large. Measured:

```text
m           256      512     1024     2048     4096
c_m       0.29873  0.31354  0.31861  0.32813  0.33670
c - c_m   0.11592  0.10112  0.09605  0.08652  0.07795
(c-c_m)*log_2 m  0.927  0.910  0.960  0.951  0.936
```

`c_m` rises monotonically toward `1/(2A)`, and the residual times `log_2 m` is
constant to within 5%: **`c_m ~ 1/(2A) - K/log_2 m` with `K ~ 0.94`**. At
`m = 4096` this predicts `0.4147 - 0.94/12 = 0.3364` against the measured
`0.3367`. So the fold accounts for the leading constant *and* the correction,
and it confirms the sharp Stirling constant against the crude one.

### 11.2 Compensation depth: the hierarchy, and where the gap lives

Let `C*_B` be the least slope over schemes balanced on ratio-`B` blocks
(THM-3007 forces `B = 2^j`). Coarser balance is strictly fewer constraints, so

```text
C*_1 >= C*_2 >= C*_4 >= ... >= C*   (global balance only).
```

**`C_arch = log_5(5 phi^2)` is a lower bound for `C*_1` only.** The entire gap
between `1.598` and the true `C*` is the value of cross-shell compensation.
The first decisive test is whether `C*_2 < C*_1`, and there is a concrete
reason to expect it: in a ratio-`B` block the middle composition range
`N <= k <= l` carries FREE symmetric deficits `delta_j = delta_{N+l-j}`,
because the two branches couple there. The ratio-2 block has no such range
(it degenerates to the single middle pair), so the extra freedom is real and
localised, and it is exactly the freedom `C_arch` does not see.

### 11.3 Syllable depth: a negative

Refining the balance requirement from `n_1` to `(n_1,...,n_j)` partitions the
shells FINER, i.e. imposes strictly MORE constraints. So depth-`j` balance can
only raise the achievable slope, and it is not a weaker hypothesis either --
it strengthens exactly the assumption one wants to remove. Dead end in both
directions; its only use is bookkeeping. Recorded as a negative, and it
disposes of the look-and-say (A005150) reading in this direction: the
run-length parse refines, and refinement is the wrong way.

## 12. Two further objects

### 12.1 The shell-imbalance module

Give shell `m` its imbalance `D_m(p) = sum_{w in S_m} eps_w p^{z(w)} q^{o(w)}`
(heads minus tails). Then

```text
global fairness  <=>  sum_{m} D_m(p) = 0 on (0,1),
shell balance    <=>  D_m = 0 for every m.
```

`D_m` is a polynomial of degree `2m`, so the shells have distinct degrees and
no finite cancellation is possible -- but the sum is infinite and
`|D_m| <= p^m + q^m`, so it converges and cancellation across infinitely many
shells is NOT excluded. **The module `{(D_m) : sum D_m = 0, each D_m
realisable at slope C}` is the precise home of the gap `C*_1 - C*`**, and it
is finitely checkable one degree at a time. The cheapest question: can two
adjacent shells compensate each other at all?

### 12.2 The alternation index -- TESTED, and REFUTED as a shared invariant

Our extremal atom saturates `|Delta^j P(0)| <= 2^j sup|P|` at
`P(y) = +-(-1)^y`; klein's THM-3010 proves the metallic recurrences
`a_k = n a_(k-1) + a_(k-2)` attain MAXIMAL circuit alternation, via Simson
`a_(k-1)a_(k+1) - a_k^2 = (-1)^k`. Define

```text
I(P) = max_(j>=1) |Delta^j P(0)| / (2^j sup|P|)  in [0,1]
```

(`j >= 1`, since at `j = 0` the ratio is `1` for every `P`). The conjecture
was that the extremals of both settings are exactly the `I = 1` elements.
**Tested against THM-3010's objects; refuted.**

*On sign sequences it is true but empty.* "Maximal circuit alternation"
*means* the sign sequence is `+-(-1)^k`, and every such sequence has `I = 1`
by definition. `I(sign) = 1.000000` for `n = 1,2,3,4`. No content.

*On value sequences it is false.* Simson gives
`R_k = 1/(1 + (-1)^k/a_k^2)`, so the circuit alternates in sign but its
magnitude decays like `a_k^-2`, geometrically. A geometrically decaying
alternating sequence has `I -> 1/2`; our ARCH extremal has CONSTANT
magnitude, which is exactly what pins `I = 1`. Measured:

```text
metallic n =        1        2        3        4
I(values)      0.700362 0.590433 0.547518 0.528546     (-> 1/2)
```

The index also fails to separate klein's distinguished strand: the four
ballot columns give `I(values) = 0.4057, 0.3832, 0.4259, 0.3775`, and the one
with a sign change (`binom(2k,k-1)`, strand B) is not meaningfully above the
others. So `I` is **not** a shared invariant.

*What survives.* `I` is a strict TIE-BREAKER inside klein's extremal class.
All metallic recurrences attain the same maximal sign-change count (12 in the
computed range), so THM-3010's criterion cannot distinguish them; `I` orders
them strictly, `0.7004 > 0.5904 > 0.5475 > 0.5285`, with **golden maximal**.
The reason is transparent: Fibonacci grows slowest, so `a_k^-2` decays
slowest and the alternation is the most sustained. `I` therefore measures
"alternation with sustained magnitude", not alternation. Worth contrasting
with THM-3010 sec 2, whose ballot family contains only the BRONZE
discriminant and no golden one -- the two selection principles pick different
metals.

(Caveat: at `n = 4` the circuit underflows in double precision for large `k`,
so its `I` is slightly contaminated; the `n = 1,2,3` values are clean and
carry the trend.)

Referees: `04-computation/amm12592_interiority_and_stirling.py`,
`04-computation/amm12592_fold_constant_and_compensation.py`,
`04-computation/amm12592_alternation_index_vs_thm3010.py`.
