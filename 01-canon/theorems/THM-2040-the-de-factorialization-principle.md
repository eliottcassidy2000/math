---
id: THM-2040
title: "De-factorialization on the exact symmetric monomial wall, and the prime-local initial-form principle"
status: >
  PROVED in the corrected scope. For even p and
  P=a Z^p + beta (Z Zbar)^(p/2) + c Zbar^p, every balanced channel of P^m
  has radial degree p*m/2, so E[P^m] is exactly (p*m/2)! times a generalized
  central-trinomial polynomial. For arbitrary support there is no global
  common-factorial/Vandermonde factorization. THM-2022 instead gives a
  prime-local normalization of one amplified moment: division by the
  lowest-face factorial leaves Q^p modulo a good prime. This corrected theorem
  supersedes the broad S204 stub and is governed by MISTAKE-215.
source: boxeph-2026-07-21-S204, corrected by codex-2026-07-21-NC2-transfer
depends_on:
  - THM-2022
related:
  - THM-2033
  - MISTAKE-215
  - MISTAKE-214
  - HYP-8795
script: 04-computation/nc2_divide_common_factorial_deathstar_S91.py
output: 05-knowledge/results/nc2_divide_common_factorial_deathstar_S91.out
---

# THM-2040 -- corrected de-factorialization principle

## 1. Exact symmetric monomial wall

Let `p` be a positive even integer, put `s=Z*Zbar`, and consider

```text
P = a Z^p + beta s^(p/2) + c Zbar^p.                  (1)
```

In a charge-zero term of `P^m`, the number of `Z^p` choices must equal the
number of `Zbar^p` choices. Write that common number as `i`; the neutral term
is chosen `m-2i` times. The common `Z` and `Zbar` exponent is then

```text
p*i + (p/2)*(m-2i) = p*m/2,                           (2)
```

independent of `i`. Complex Wick evaluation therefore gives the exact identity

```text
E[P^m] = (p*m/2)! W_m(a,beta,c),                      (3)

W_m = sum_(0<=i<=m/2)
        m!/(i! i! (m-2i)!) * (a*c)^i * beta^(m-2i).  (4)
```

At `a=beta=c=1`, `W_m` is the central trinomial coefficient

```text
[x^0](1+x+x^(-1))^m = 1,3,7,19,51,... .              (5)
```

Equivalently, its exponential generating function is

```text
sum_(m>=0) W_m z^m/m!
  = exp(beta*z) * sum_(i>=0) (a*c*z^2)^i/(i!)^2.     (6)
```

Equations (3)-(6) are the valid de-factorialization statement. They require an
exact neutral monomial of degree `p/2`. A lower-term polynomial in `s` produces
different radial degrees and no longer has the common factorization (3).

The identity is not an NC2 obstruction in this special family: `W_1=beta` and
`W_2=beta^2+2ac`. If every positive moment vanished, then `beta=0` and
`ac=0`, so the exact support is one-sided or zero after deleting zero
coefficients.

## 2. Why the proposed global extension is false

For a general polynomial, balanced channels have different radial heights
`A(r)` and hence different Wick factors `A(r)!`. Dividing a scalar moment by
one chosen factorial leaves channel-dependent factorial quotients. It cannot
turn the scalar sum into the moment-matrix determinant of THM-2033, and that
determinant's Vandermonde factor is not a universal residue of `E[P^m]`.

The S91 computation already falsifies its own boundedness interpretation: for
`P=Z+(1+s)+Zbar`, the printed values of `E[P^m]/m!` rise from `2` to about
`35.76` by `m=8`. MISTAKE-215 records the full correction. The additional
positive-zero, Laguerre--Polya, and Paley iff claims are also withdrawn; the
Paley identification separately violates MISTAKE-214.

## 3. The valid general principle is prime-local

THM-2022 chooses a lowest balanced face `F`, a base return level `m0`, and a
good prime `p`. At the single amplified order `p*m0`, every channel height is
at least `p*A0`, so division by `(p*A0)!` is termwise integral. The remaining
factorial quotients are not discarded over the integers. Instead, modulo the
chosen prime:

- Kummer kills non-p-dilated multinomial channels;
- strict face height makes every dilated off-face factorial quotient divisible
  by `p`; and
- Lucas plus Frobenius turns the complete face residue into `Q^p`.

Thus the reusable general statement is an **initial-form normalization at a
good finite place**, not a global factorization. Its residue is the whole face
constant term and need not be a determinant, Vandermonde, central trinomial,
or tournament invariant.
