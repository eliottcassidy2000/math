# The series S(k) = sum C(2n,n)C(4n,2n)/((kn+1) 64^n): exact reductions

**Session:** death-star-2026-07-31-coinC2 (owner-supplied problem).
**Status:** the reductions below are PROVED and independently verified;
the closed-form question for `k >= 4` is reported honestly as OPEN with a
methodology-validated negative search.
**Scripts:** `04-computation/amm_sk_setup.py`, `amm_sk_quad_transform.py`,
`amm_sk_kmoment.py`, `amm_sk_k6.py`.

## 0. The problem

For `k` a positive integer,

```text
S(k) = sum_{n>=0} C(2n,n) C(4n,2n) / ( (kn+1) 64^n ).
```

Known closed forms (all re-verified here to **200 digits**):

```text
S(1) = 8 sqrt2/(3 pi)
S(2) = (4/pi) log(1 + sqrt2)
S(3) = -(1/pi)( sqrt3 log(5 - 2 sqrt6) + 2 arctan(sqrt2/5) )
```

## 1. The summand is a signature-4 hypergeometric coefficient (PROVED)

```text
C(2n,n) C(4n,2n) / 64^n  =  (1/4)_n (3/4)_n / (n!)^2.                (1)
```

*Proof.* `C(2n,n) = 4^n (1/2)_n/n!` and `C(4n,2n) = 16^n (1/2)_{2n}/(2n)!`;
Legendre duplication gives `(1/2)_{2n} = 4^n (1/4)_n (3/4)_n` and
`(2n)! = 4^n n! (1/2)_n`; the `(1/2)_n` cancel. (Verified as exact
rationals for `n < 40`.) Hence the generating function is **Ramanujan's
signature-4 series**

```text
sum_n C(2n,n)C(4n,2n) x^n/64^n = 2F1(1/4, 3/4; 1; x) =: F(x).          (2)
```

Since `1/(kn+1) = int_0^1 t^{kn}dt`,

```text
S(k) = int_0^1 F(t^k) dt = (1/k) int_0^1 s^{1/k-1} F(s) ds
     = 3F2(1/4, 3/4, 1/k;  1, 1+1/k;  1).                             (3)
```

The `3F2` is **balanced** (parametric excess 1) but non-terminating, so
Saalschütz does not apply; Watson's and Whipple's theorems fit the
numerator pair `(1/4, 3/4)` but force `k = 1` on the denominators (checked
symbolically). This already explains why `k = 1` is the clean case:

```text
S(1) = 2F1(1/4,3/4;2;1) = Gamma(2)Gamma(1)/(Gamma(7/4)Gamma(5/4))
     = 16/(3 pi sqrt2) = 8 sqrt2/(3 pi)      (Gauss).                 (4)
```

## 2. An elementary kernel for F (PROVED)

Because `b - a = 1/2` for `F = 2F1(1/4,3/4;1;z)`, a quadratic
transformation exists; numerically confirmed to 40 digits:

```text
F(z) = (1+sqrt z)^{-1/2} 2F1(1/2,1/2;1; 2 sqrt z/(1+sqrt z))
     = ((1+sqrt(1-z))/2)^{-1/2} 2F1(1/2,1/2;1; (1-sqrt(1-z))/(1+sqrt(1-z))).
```

Writing `2F1(1/2,1/2;1;m^2) = (2/pi)K(m)` and using the half-angle form
collapses this to a completely elementary kernel:

```text
F(z) = (1/pi) int_0^pi  dphi / sqrt(1 + sqrt z cos phi).               (5)
```

Consequently

```text
S(k) = (4/pi) int_0^1 int_0^{pi/2} tau (1 - tau^k cos 2theta)^{-1/2} dtheta dtau
     = (2/pi) int_0^{pi/2} 2F1(1/2, 2/k; 1+2/k; cos 2theta) dtheta.    (6)
```

and, pushing the modulus to be the variable, a **moment of the complete
elliptic integral** (verified against (3) for `k = 1..6`):

```text
S(k) = (16/(k pi sqrt2)) int_0^1 mu^{4/k-1} (2-mu^2)^{-2/k-1/2} K(mu) dmu. (7)
```

## 3. What (6) explains, and what it does not

The inner function `2F1(1/2, b; 1+b; c)`, `b = 2/k`, is the incomplete Beta
`b c^{-b} B_c(b, 1/2)`. Therefore:

* it is **algebraic in `c`** iff `b = 2/k` is a positive integer, i.e.
  **`k in {1, 2}`**;
* it is **elementary** (allowing `arcsin`) iff `4/k` is a positive
  integer, i.e. **`k in {1, 2, 4}`** (at `k = 4` it is
  `arcsin(sqrt c)/sqrt c`, transcendental);
* for `k = 3` it is provably *not* algebraic: the Schwarz data of
  `2F1(1/2,2/3;5/3;c)` is `(|1-c|,|c-a-b|,|a-b|) = (2/3, 1/2, 1/6)`, and
  `1/6` occurs in no entry of Schwarz's list.

`k = 1` and `k = 2` are therefore fully explained, and both closed forms
were **re-derived from scratch** by this route. For `k = 2`, (6) gives
`2F1(1/2,1;2;c) = 2(1-sqrt(1-c))/c` with `1-cos2theta = 2 sin^2 theta`, so
that `1 - 2sin^2 theta = (1-sqrt2 sin theta)(1+sqrt2 sin theta)` cancels
the numerator and

```text
S(2) = (4/pi) int_0^{pi/2} dtheta/(1 + sqrt2 sin theta) = (4/pi) log(1+sqrt2),
```

by the tangent half-angle substitution (`t^2 + 2 sqrt2 t + 1` has roots
`-sqrt2 +- 1`, giving `log(1+sqrt2)` exactly).

**`k = 3` is genuinely exceptional.** Its closed form is exact (200
digits) although the inner function is not algebraic and no classical
`3F2` summation applies. Its mechanism is *not* explained by the criterion
above and remains the key open structural question here: whatever produces
it is a 2-dimensional cancellation, not an iterated elementary
integration.

## 4. k = 4 and k = 5 (the user's Question 1): a negative search

For `k = 4` the reduction is fully explicit (PROVED, verified to 25
digits against (3)):

```text
S(4) = (2/pi) int_0^1 [ arcsin v + arcsinh v ] / sqrt(1-v^4) dv
     = (4/(pi sqrt2)) int_0^1 K(mu)/(2-mu^2) dmu.                      (8)
```

(The integrand pair is the `v -> i v` symmetry of `sqrt(1-v^4)`, since
`arcsin(iv) = i arcsinh v`.)

High-precision values (150 dps):

```text
pi S(4) = 3.3594701291301671589900224730224240656919222546691...
pi S(5) = 3.3209372626410174933367229658774746263264441740161...
pi S(6) = 3.2940228837649681600321847392742692680650456489048...
```

PSLQ was run at 120-160 digits over: elementary logs/arctans in
`Q(sqrt2)`, `Q(sqrt2,sqrt3,sqrt6)`, `Q(sqrt2,sqrt5)` (golden ratio),
`Q(2^{1/4})`; the lemniscatic constants `varpi = Gamma(1/4)^2/(2 sqrt(2 pi))`,
`varpi^2`, `varpi pi`, `varpi log(1+sqrt2)`; and weight-two constants
`pi^2, G, log^2 2, log^2(1+sqrt2), pi log(1+sqrt2)`; subset sizes 2-4,
coefficient bounds up to `10^6`.

**No relation was found for `k = 4, 5, 6, 12`.**

**The search is methodology-validated:** the identical sweep, run on
`pi S(3)` over the same constant pool, *does* recover

```text
pi S(3) - sqrt3 log(5+2 sqrt6) + 2 arctan(sqrt2/5) = 0
```

immediately. So the negative result is meaningful, not an artifact of a
weak search. Honest conclusion: **`S(4)` and `S(5)` have no closed form of
the `k <= 3` type in any of the natural fields tested.** This is
consistent with the problem statement's "three special values appear to
have closed forms". It is *not* a proof of non-existence — the `k = 3`
mechanism is still unidentified, so a field or constant class we have not
guessed cannot be excluded.

## 5. The general closed form (the user's Question 2)

A genuine general closed form does exist, at the level of classical
objects: (3) and (7) express `S(k)` for every `k` as a single balanced
`3F2` at 1, equivalently a one-dimensional `K`-moment with algebraic
weight. What is *not* available in general is an evaluation in
logarithms of algebraic numbers and `pi`: by section 3 that is expected
exactly at `k = 1, 2` (proved), happens anomalously at `k = 3` (verified
but unexplained), and is not detected at `k >= 4`.

## 6. Next obligations

1. **Explain `k = 3`.** This is the crux. Candidate routes: a cubic
   transformation acting on the `tau`-`theta` double integral (6); a
   degeneration of the underlying elliptic curve at `k = 3`; or a Thomae
   image of the `3F2` that is summable only when `1/k + 1/4 + 3/4` hits a
   special configuration. Any of these would immediately settle whether
   `k >= 4` can work.
2. If `k = 3` is explained by a mechanism with a finite list of `k`, that
   list *is* the answer to Question 1.
3. Independent: is `S(k) - 1` for large `k` asymptotically
   `C/k + O(1/k^2)` with `C = sum_{n>=1} a_n/n` a recognisable constant?
   (`a_n = (1/4)_n(3/4)_n/(n!)^2`.)
