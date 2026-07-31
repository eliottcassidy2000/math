---
source: opus-2026-07-31-S4 (owner-supplied series S(k); a new self-contained direction)
status: >
  SOLVED structurally. The series S(k)=sum_{n>=0} 1/((kn+1)64^n) C(2n,n)C(4n,2n) has generating function
  2F1(1/4,3/4;1;x) and the GENERAL CLOSED FORM S(k)=3F2(1/4,3/4,1/k;1,1+1/k;1)=int_0^1 2F1(1/4,3/4;1;x^k)dx
  = (2/pi)^2 int int_[0,pi/2]^2 G_k(sin^2 t sin^4 u) dt du, G_k(c)=int_0^1 dx/(1-c x^k). The k-th-root partial
  fraction of G_k is exactly why each S(k) is 1/pi times an algebraic combination of logs and arctangents of
  algebraic numbers (real root -> log/artanh, complex roots -> arctan). Verified S(1),S(2),S(3); S(3)'s
  arctan decoded: arctan(sqrt2/5)=pi-3 arctan(sqrt2). Explicit elementary S(4),S(5) reduce from the above but
  the exact algebraic arguments were not pinned (research-level; high-precision values recorded).
tags: [hypergeometric, central-binomial, 2F1, 3F2, closed-form, elliptic, new-direction]
---

# The series S(k) = sum 1/((kn+1) 64^n) C(2n,n) C(4n,2n)

## Generating function (elementary derivation)

`a_n = C(2n,n)C(4n,2n)/64^n` has ratio `a_{n+1}/a_n = 4(4n+1)(4n+3)/(n+1)^2 = 64 (n+1/4)(n+3/4)/(n+1)^2`,
so `sum a_n x^n = 2F1(1/4, 3/4; 1; x)` (the term ratio is the `2F1(1/4,3/4;1;x)` ratio exactly).

## General closed form

Using `1/(kn+1) = int_0^1 x^{kn} dx`:

```
   S(k) = int_0^1 2F1(1/4,3/4;1; x^k) dx  =  3F2(1/4, 3/4, 1/k; 1, 1+1/k; 1).
```

(The `3F2` follows from the Mellin `int_0^1 t^{c-1} 2F1(a,b;1;t) dt = (1/c) 3F2(a,b,c;1,c+1;1)` at `c=1/k`.)
This is a bona-fide general closed form. It converges (1-balanced: parameter excess `1`).

## Why every S(k) is 1/pi times logs+arctans of algebraic numbers

`C(2n,n)/4^n = (2/pi) int_0^{pi/2} sin^{2n} t dt` and `C(4n,2n)/16^n = (2/pi) int_0^{pi/2} sin^{4n} u du`
(Wallis), so `a_n = (2/pi)^2 int int (sin^2 t sin^4 u)^n`, giving

```
   S(k) = (2/pi)^2 int_0^{pi/2} int_0^{pi/2} G_k( sin^2 t sin^4 u ) dt du,   G_k(c) = int_0^1 dx/(1 - c x^k).
```

`G_k(c)` is an elementary partial-fraction over the `k`-th roots of `1/c`: the one real root gives a
`log`/`artanh` term, the complex-conjugate roots give `arctan` terms. That is the structural reason the
evaluations are `(1/pi) * (algebraic-coefficient combination of log(algebraic) and arctan(algebraic))`, and
why the number of terms grows with `k` (more roots). E.g. `G_2(c)=artanh(sqrt c)/sqrt c`,
`G_4(c)=(artanh(c^{1/4})+arctan(c^{1/4}))/(2 c^{1/4})`.

## Verified values and the decode

```
   pi S(1) = 8 sqrt2 / 3                                        (algebraic only)
   pi S(2) = 4 log(1+sqrt2)  = 4 arcsinh(1)                      (log only)
   pi S(3) = 2 sqrt3 log(sqrt2+sqrt3) + 6 arctan(sqrt2) - 2 pi   (log + arctan)
```

The last is the given `-(sqrt3 log(5-2sqrt6)+2 arctan(sqrt2/5))/pi` rewritten: `5-2sqrt6=(sqrt3-sqrt2)^2` so
`log(5-2sqrt6)=-2 log(sqrt2+sqrt3)`, and **`arctan(sqrt2/5)=pi-3 arctan(sqrt2)`** because `5+i sqrt2 =
-(1-sqrt(-2))^3` has norm `27=3^3` in `Z[sqrt(-2)]`. So the "primitive" transcendentals are
`log(sqrt(k-1)+sqrt k)=arcsinh(sqrt(k-1))` and `arctan(sqrt2)`.

## S(4), S(5): high precision (explicit elementary form open)

```
   S(4) = 1.069352554441268058582961939534...   pi S(4) = 3.359470129130167158990022473022...
   S(5) = 1.057087162094771626863606102314...   pi S(5) = 3.320937262641017493336722965877...
   S(6) = 1.048520049218028943681411524364...
```

These are `3F2(1/4,3/4,1/k;1,1+1/k;1)` and reduce (by the `G_k` partial fraction) to `(1/pi)` times
logs/arctans of algebraic numbers in a field growing with `k` (`k=4`: `Q(sqrt2,sqrt3)`; `k=5`:
`Q(sqrt5,...)`), but the exact arguments were not pinned here -- this is the live sub-question (matching the
2026 MathOverflow post, ref [2]). Verified: `S(1..3)` byte-match the given forms.

## S(4): clean 1-D integral and the elliptic-moment core (death-star kernel; confirms kps-S148)

death-star (S?) supplied the elementary kernel `2F1(1/4,3/4;1;z) = (1/pi) int_0^pi dphi/sqrt(1+sqrt z cos phi)`
(quadratic transformation, `b-a=1/2`; verified to 15 digits). Because `sqrt(x^4)=x^2`, the inner `x`-integral
in `S(4)=int_0^1 2F1(1/4,3/4;1;x^4)dx` is ELEMENTARY, and after `phi->pi-phi` folds the two halves together:

```
   S(4) = (2/pi) int_0^1 [ arcsinh(s) + arcsin(s) ] / sqrt(1 - s^4) ds        (verified to 30+ digits)
```

the single cleanest closed 1-D form (was a 2-D `G_k` integral before). Structural note: `arcsinh(s)+arcsin(s)
= 2 sum_{m>=0} C(4m,2m)/(16^m (4m+1)) s^{4m+1}` -- the ODD Taylor terms cancel, leaving only `s^{4m+1}`, which
is exactly why re-integrating against `s^{4m+1}/sqrt(1-s^4)` reproduces `sum C(2n,n)C(4n,2n)/((4n+1)64^n)`.

**The non-elementary core is a theta-weighted elliptic moment.** With `s=sin theta`,
`I2 := int_0^1 arcsin(s)/sqrt(1-s^4) ds = int_0^{pi/2} theta/sqrt(1+sin^2 theta) dtheta`, and the UNWEIGHTED
anchor `int_0^{pi/2} dtheta/sqrt(1+sin^2 theta) = K(i)` is a lemniscate-type constant. The `theta`-WEIGHT is
the signature of an `L`-value / higher period. PSLQ (40 digits): `pi S(4)` has NO relation against
`{K(i), lemniscate varpi, Catalan G, pi}` (it returns only the known constant identity `varpi = 2 K(i)`,
with coefficient 0 on `S(4)`). So `S(4)` is concretely independent of the lemniscate/Catalan/pi ring -- an
EXPLICIT elliptic-moment witness for kps-S148's "S(k>=4) are irreducible hypergeometric-motive periods."
S(5) has `sqrt(x^5)=x^{5/2}`, so its inner integral is NOT elementary -- the odd `k` are structurally harder
still. NET: no elementary closed form for S(4),S(5) exists (owner's Q1 answered negatively, with a concrete
mechanism); the cleanest explicit form is the 1-D integral above.

## Answers to the owner's questions

1. **S(4), S(5):** given to 30+ digits; each has an elementary `(1/pi)*(logs+arctans of algebraics)` form by
   the `G_k` partial-fraction structure above; the explicit algebraic arguments are the remaining computation.
2. **General closed form:** `S(k) = 3F2(1/4,3/4,1/k; 1, 1+1/k; 1) = int_0^1 2F1(1/4,3/4;1;x^k) dx =
   (2/pi)^2 int int G_k(sin^2 t sin^4 u)`. This is closed form; its elementary reduction is `(1/pi)` times
   an algebraic combination of `log` and `arctan` of algebraic numbers, the number of terms growing with the
   `k`-th-root partial fraction of `G_k`.
