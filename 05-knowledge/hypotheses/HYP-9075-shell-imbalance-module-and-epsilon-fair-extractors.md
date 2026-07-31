---
id: HYP-9075
title: "The shell-imbalance module for AMM 12592, and epsilon-fair extractors"
status: >
  PARTIAL. The module is now exactly characterised and the first question is
  settled negatively: NO two ADJACENT shells can compensate, verified
  exhaustively for m = 1, 2, 4. Structure forced: c_0 = c_2m = 0 and every
  c_k even (dyadic m), and C^(m)(1) = 0, hence
  C^(m)(u) = 2 u (1-u) G(u) with G integer of degree <= 2m-3. Adjacent
  compensation is the exact convolution c^(2m) = -c^(m) * (1+u)^(2m) landing
  in its own box, and there is a clean evaluation family
  |C^(m)(w)| (1+w)^(2m) <= N_(2m)(w) for all w > 0. A guessed obstruction
  (the two-word middle class) is REFUTED: at m = 2 the binding index is the
  low-order i = 2, and the middle constraint is satisfiable. General m and
  non-adjacent pairs are OPEN. Separately proposes epsilon-fair extractors as
  a new route, prompted by arXiv 2607.27088 (exact TV is hard, approximate TV
  is linear time).
source: opus-2026-07-31-amm12592-writeup
related:
  - THM-3007  # balanced blocks are exactly the dyadic intervals
  - THM-3008  # within-block deadline profiles
  - THM-3009  # the archimedean floor C_arch = log_5(5 phi^2)
  - THM-2966  # spine normal form
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
  - "K. Anand, A. Benford, H. Guo, Linear time approximation of the TV distance between product distributions, arXiv:2607.27088."
script: 04-computation/amm12592_shell_imbalance_module.py
output: 05-knowledge/results/amm12592_shell_imbalance_module.out
---

# HYP-9075 -- the shell-imbalance module

## 1. The module

Give shell `S_m = {b^m y : y != b^m}` a colouring `eps: S_m -> {+-1}` and set

```text
D_m(p) = sum_(w in S_m) eps_w p^(z(w)) q^(o(w)) = sum_k c_k p^(2m-k) q^k,
c_k = sum_(o(w)=k) eps_w.
```

Because `{p^(2m-k) q^k}` is the BERNSTEIN basis, `D_m = 0` iff `c = 0`, and
the realisable `c` are exactly

```text
|c_k| <= N_k,   c_k = N_k (mod 2),
sum_k N_k u^k = (1+u^m)(1+u)^m - 1 - u^(2m)      (THM-3007's class GF).
```

Global fairness is `sum_m D_m = 0`; shell balance is every `D_m = 0`. **The
gap `C*_1 - C*` of THM-3009 sec 11.2 lives exactly in the nonzero solutions.**

## 2. Forced structure

`N_0 = N_(2m) = 0`, so `c_0 = c_(2m) = 0`. For dyadic `m` every `N_k` with
`1 <= k <= 2m-1` is even (Lucas), so every `c_k` is even. And summing the
box constraint against `u = 1` in the convolution of section 3 forces
`C^(m)(1) = 0`. Hence

```text
C^(m)(u) = 2 u (1-u) G(u),    G in Z[u],  deg G <= 2m-3.
```

For `m = 2` that is only two free integers; the module is small at the bottom
and this is why the low cases are decidable by hand.

## 3. Adjacent compensation is an exact convolution

Re-expanding the degree-`2m` Bernstein basis in the degree-`4m` one
(multiplying by `(p+q)^(2m) = 1`) turns `D_m + D_(2m) = 0` into

```text
c^(2m)(u) = - C^(m)(u) (1+u)^(2m),
```

so `c^(m)` DETERMINES `c^(2m)` and the only question is whether it lands in
its own box. Since coefficientwise domination implies the weighted-`l^1`
bound, for every `w > 0`

```text
|C^(m)(w)| (1+w)^(2m) <= N_(2m)(w) = (1+w^(2m))(1+w)^(2m) - 1 - w^(4m),
```

a clean necessary family. At `w = 1` it gives `|C^(m)(1)| <= 2 - 2/4^m < 2`,
which with evenness forces `C^(m)(1) = 0` as used above.

## 4. Result: no adjacent compensation at the bottom

Exhaustive over the admissible box:

```text
m = 1:      3 vectors   -> none
m = 2:     27 vectors   -> none
m = 4:  91875 vectors   -> none
```

**A guessed obstruction is refuted.** We expected the two-word middle class
(`N_(2m) = 2`, THM-3007's cancelling corner pair) to be the bottleneck. It is
not: the minimum of `|sum_j c_j binom(2m,j)|` over admissible nonzero `c` is
`0` at both `m = 2` and `m = 4`, so the middle constraint is satisfiable. The
actual binding index at `m = 2` is `i = 2`, at the LOW-order edge: for
`c = (0,-2,0,2,0)` the convolution gives `|c^(4)_2| = 8` against
`N^(4)_2 = binom(4,2) = 6`. The obstruction is the edge of the row, where the
box `binom(2m,i)` grows more slowly than convolution with `(1+u)^(2m)`
amplifies -- not the middle, where the box is smallest but the weighted sum
can be cancelled.

## 5. Open

1. General `m`: is adjacent compensation impossible for every dyadic `m`?
   The evaluation family of section 3 is the natural tool; it has not been
   pushed past `w = 1`.
2. NON-adjacent pairs `(m, 2^j m)`, and finite sets of shells. The convolution
   becomes `(1+u)^(2^j m - ...)`, with far more room.
3. Infinite compensation. `deg D_m = 2m` are distinct, so no FINITE
   cancellation exists, but `|D_m| <= P(S_m) = p^m - p^(2m) + q^m - q^(2m)`
   makes the infinite sum converge: cross-shell cancellation over infinitely
   many shells is not excluded by any argument here.

## 6. A new route: epsilon-fair extractors

Prompted by arXiv 2607.27088 (Anand-Benford-Guo, linear-time approximation of
the TV distance between product distributions). Our problem is the EXACT
vanishing of a signed sum against a product measure, uniformly in the
parameter; that paper's theme is that the exact version of such a quantity is
hard while the approximate version is cheap. So relax:

```text
epsilon-fair:   |P_p(H) - 1/2| <= epsilon   for all p in (0,1),
C*(epsilon) := the least slope of an epsilon-fair extractor.
```

`C*(0)` is the original problem. The ARCH obstruction of THM-3009 came from
needing exact cancellation of FORCED half-integer corner deficits, so a slack
`epsilon` should relax it directly, and the object of interest is the rate

```text
C*(epsilon) -> C*(0)   as   epsilon -> 0.
```

Two reasons this is worth doing: it is the natural home of the
exact-vs-approximate dichotomy that the cited paper exemplifies, and a fast
decay would say the exact problem's rigidity is an integrality artifact
rather than a real information barrier. NOT yet computed; the (ARCH)
derivation should carry over with `binom(m-1,d)` replaced by
`binom(m-1,d) - epsilon 2^m`-ish, which is a one-line change to the existing
referee.
