---
id: THM-3260
title: "Ternary-Cantor Bessel--Cartier residual exclusion"
status: >
  PROVED + VERIFIED-EXACT.  Let D>=3 be odd and divisible by 3.  If every
  ternary digit of D-1 is 0 or 2, then the resonant exact-quadratic
  factorial-moment pair at D is coprime.  Both full moment polynomials have
  their terminal Bessel lower faces; their only shared Newton slope is 1/2.
  The residual edges form a two-sector Cartier module over F_3, and exact
  Frobenius recurrences prove the sectors coprime for every such ternary
  word.  This gives an infinite mixed-prime closure, not a finite
  extrapolation.
source: root/factorial-composite-newton-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3204-odd-prime-power-resonance-laguerre-bessel-residual-exclusion
related:
  - THM-3200-even-resonance-bessel-residual-newton-exclusion
  - THM-3131-prime-resonance-newton-slope-separation
script: 04-computation/factorial_ternary_cantor_bessel_residual_thm3260.py
output: 05-knowledge/results/factorial_ternary_cantor_bessel_residual_thm3260.out
script_sha256: 62059117695a038748099f8a8da173adb4472a85c6386e74791d114b59e67d4c
output_sha256: cd63c5859ff2cd3e44c5a2eac7651dbb7e01e5c36d323888dbb65af33fcef715
hash_basis: LF-normalized bytes
---

# THM-3260 -- ternary-Cantor Bessel--Cartier residual exclusion

**PROVED + VERIFIED-EXACT.**

Let

```text
L(t^k)=k!,                    q(t)=a_0+b t+c t^2,             (1)
```

with `a_0 b c!=0`.  Suppose a three-moment window begins at `r>=1`, put

```text
D=r+2,                         N=D-1,                         (2)
```

and assume

```text
D is odd,  3|D,  and every base-3 digit of N is 0 or 2.      (3)
```

Then

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                              (4)
```

cannot all vanish.

Condition `(3)` has three equivalent and useful descriptions:

1. `N=D-1` is a ternary Cantor integer and its final digit is `2`;
2. if `D=3^a U`, then the base-3 word of `U` has exactly one digit `1`,
   and that `1` is its final digit;
3. the base-3 word of `D` is obtained from a word over `{0,2}` by replacing
   its terminal run of one or more `2` digits by `10...0`.

Thus this is an infinite chamber.  Its first genuinely mixed resonances
past THM-3124's finite range include

```text
D=219=22010_3,  237=22210_3,  489=200010_3,  507=200210_3.  (5)
```

The theorem closes all word lengths at once; `(5)` is not the universe of
the argument.

## 1. The resonant pair and its terminal expansion

THM-3124 shows that a bad window forces `b/a_0=-1/D`.  Normalize by `a_0`,
put `v=Dc/a_0`, and define

```text
A_n(v)=L((D-t+v t^2)^n)=sum_(j=0)^n a_(n,j)v^j.              (6)
```

It is enough to exclude a common root of `A_(D-2)` and `A_(D-1)`.  Put

```text
C(n,j)=(n+j)!/[j!(n-j)!].                                   (7)
```

Changing from the number `ell` of selected `-t` terms to
`k=n-j-ell`, the coefficient formula becomes the exact terminal expansion

```text
a_(n,j)=n! sum_(k=0)^(n-j)
  (-1)^(n-j-k) D^k/k! C(n-k,j).                              (8)
```

The `k=0` term is exactly the scaled terminal Bessel coefficient used in
THM-3204.  Unlike the prime-power argument, individual off-face correction
terms need not be integral relative to that terminal coefficient.  The
right invariant is their absolute height above the target face.

## 2. The Cantor-window valuation lemma

Write `v` for `v_3`, `s` for the ternary digit sum, and let `a=v(D)`.  The
following lemma is the bridge from the full Laguerre sum `(8)` to the
Bessel boundary.

### Lemma 2.1 (thickened Cantor window)

Let `N` have only `0,2` as ternary digits and suppose `N` ends in `2`.  Put
`D=N+1` and `a=v(D)>=1`.  In the indicated ranges,

```text
v(D^k/k! C(N-k,j))       >= j/2,                             (9)
v(D^k/k! C(N-1-k,j))     >= (j-1)/2       (j>=1),           (10)
v(D^k/k! C(N-1-k,0))     >= 0.                              (11)
```

Here the left sides are integers, so a half-integral right side means its
ceiling.  If `k>0`, every inequality is strict as an inequality of rational
numbers.  At `k=0`, equality in `(9)` can occur only for even `j`, and
equality in `(10)` only for odd `j`.

### Proof

Legendre's formula gives

```text
2v(C(m,j))=j+s(j)+s(m-j)-s(m+j).                             (12)
```

For completeness, here is the digit mechanism behind the lemma.  If
`B_N(x)` is the total number of borrows when subtracting `x` from `N`, and
`K_N(x)` the total number of carries when adding `x` to `N`, columnwise
subtraction and addition give the exact identities

```text
s(N-x)=s(N)-s(x)+2B_N(x),
s(N+x)=s(N)+s(x)-2K_N(x).                                   (13)
```

The point of the alphabet `{0,2}` is that a nonzero input digit must be
paid for on one of the two sides: over a center digit `0` subtraction
borrows, while over a center digit `2` addition carries.  Incoming
borrow/carry bits only strengthen this assertion.  Scanning the digits of
`k` and `j` simultaneously, and charging a borrow chain at its first digit,
gives

```text
(2a-1)k+s(k)+s(j)+s(N-k-j)-s(N-k+j) >= 0,                   (14)

(2a-1)k+s(k)+s(j)+s(N-1-k-j)-s(N-1-k+j)+1 >= 0.             (15)
```

In `(15)` use only `j>=1`; `j=0` is immediate from
`ak-v(k!)>0` for `k>0`.  A zero in `(14)` or `(15)` forces no digit of `k`
to spend the positive credit supplied by the trailing `a`-digit borrow
chain, hence forces `k=0`.  At `k=0`, the final digit parity gives the stated
even/odd equality restrictions.

Here is an equivalent induction audit of the charging step.  Write the
last digits of `(N,k,j)` as `(epsilon,u,w)`, with
`epsilon in {0,2}` and `u,w in {0,1,2}`; pass the subtraction borrow and
addition carry to the next column.  Substitution of

```text
s(3x+r)=s(x)+r,               r=0,1,2,                       (16)
```

in `(14)`--`(15)` reduces each of the `18` digit triples, for each incoming
bit, to the same inequality for the shortened words plus a nonnegative
integer.  The added integer can vanish with `u!=0` only while the initial
run `epsilon=2` is being scanned; the term `(2a-1)k` pays exactly that run,
and the outgoing column is then positive.  This proves `(14)`--`(15)` by
induction on the word length and also proves the equality statement.

Finally, substituting `m=N-k` and `m=N-1-k` in `(12)`, and using

```text
2v(D^k/k!)=(2a-1)k+s(k),                                    (17)
```

turns `(14)` and `(15)` into `(9)` and `(10)`.  This proves the lemma.

### Consequence for the full moment pair

Let `h=v((D-1)!)=v((D-2)!)`.  Lemma 2.1 and `(8)` show that every
`k>0` summand lies strictly above the following lower-face lines, while the
`k=0` Bessel term supplies their residual coefficients:

```text
A_(D-1):       y=h+j/2,
A_(D-2):       y=h                 on 0<=j<=1,
               y=h+(j-1)/2        on 1<=j<=D-2.             (18)
```

Thus no cancellation or denominator argument is hidden here: the full
moment faces are literally the terminal Bessel faces.

## 3. The two residual sectors

We now calculate those Bessel faces.  For every even ternary Cantor integer
`M`, define polynomials over `F_3` by

```text
P_M(X)=sum_(k=0)^(M/2)
  overline(3^(-k) C(M,2k)) X^k,

Q_M(X)=sum_(k=0)^((M-2)/2)
  overline(3^(-k) C(M-1,2k+1)) X^k.                         (19)
```

A coefficient in `(19)` is declared `0` when its valuation is strictly
larger than the displayed power.  The `k=0` case of the digit argument in
Lemma 2.1, which does not require the final digit to be `2`, proves that
all displayed divisions are integral.

For `M=N`, up to separate nonzero row scalars, `P_N` is the residual edge
of `A_(D-1)` and `Q_N` the slope-`1/2` residual edge of `A_(D-2)`.

## 4. Exact Cartier recurrences

The pair `(P_M,Q_M)` contracts one ternary digit at a time.  If `M` is a
ternary Cantor integer, then

```text
P_(3M)   = P_M^3-X Q_M^3,             Q_(3M)   =-X Q_M^3,   (20)

P_(3M+2) = P_M^3+X R_M^3,             Q_(3M+2) =-P_M^3,     (21)
```

where

```text
R_M=P_M-Q_M       if M==0 (mod 3),
R_M=-(P_M+Q_M)    if M==2 (mod 3).                           (22)
```

These are Frobenius cubes in `F_3[X]`, not suggestive fits.

To prove them, let `U(m)` be the unit part of `m!` modulo `3`.  Complete
blocks of three give

```text
v((3m+r)!)=m+v(m!),
U(3m+r)=(-1)^m U(m) eta_r,
eta_0=eta_1=1, eta_2=-1.                                   (23)
```

Apply `(23)` to the three factorials in `C(n,j)`.  Sorting the residual
exponent `k` modulo `3` gives the complete contraction table

```text
index k       3ell             3ell+1                         3ell+2
-----------------------------------------------------------------------
P_(3M)        p_(M,ell)       -q_(M,ell)                      0
Q_(3M)        0               -q_(M,ell)                      0
P_(3M+2)      p_(M,ell)        r_(M,ell)                      0
Q_(3M+2)     -p_(M,ell)        0                              0,        (24)
```

with `r_(M,ell)` the coefficient prescribed by `(22)`.  Terms whose
normalized valuation is positive give the zeros in `(24)`.  Reading the
three residue classes as Cartier components is exactly `(20)`--`(21)`.

The base pair is

```text
P_0=1,                         Q_0=0.                         (25)
```

## 5. Coprimality is an invariant of the digit tree

We prove by induction that

```text
gcd(P_M,Q_M)=1                                                     (26)
```

for every ternary Cantor `M`.

If `M=3m`, equations `(20)` say that a common nonzero root must make both
`Q_m` and `P_m` vanish.  The root `X=0` is impossible because
`P_(3m)(0)=P_m(0)^3=1`.

If `M=3m+2`, equation `(21)` first forces `P_m=0`.  At such a root,
`R_m` is a nonzero scalar multiple of `Q_m`, so the first equation in
`(21)` forces `Q_m=0`; again induction excludes it.  Thus `(26)` follows
from `(25)`.

When `N` ends in `2`, `(21)` also shows that `Q_N` has the full endpoint
degree `(N-2)/2`, while `P_N` has degree `N/2`.  Consequently the exact
lower hulls in `(18)` are

```text
A_(D-1): (0,h)----------------------------(D-1,h+(D-1)/2),

A_(D-2): (0,h)---(1,h)--------------------(D-2,h+(D-3)/2).  (27)
```

Their only shared slope is `1/2`.  Its denominator is prime to `3`, so the
ordinary residual-polynomial Newton theorem applies.  Equation `(26)` says
that the two residual edges are coprime.  The two full polynomials therefore
have no common algebraic root, completing the proof of the theorem.

## 6. Why the tempting broader chamber fails

The earlier empirical phrase “after stripping powers of `3`, the leading
ternary digit is `2` and the trailing digit is `1`” is false.  The first
hostile inside that broader language is

```text
D=201=21110_3,                    D/3=2111_3.                (28)
```

At the common slope `1/2`, the exact full residual polynomials are

```text
Q=1+X^81=(1+X)^81,
P=1+X+X^81+X^82=(1+X)Q.                                   (29)
```

Thus their gcd has degree exactly `81`.  The obstruction is not a random
large factor: the extra ternary `1` glues the two Cartier sectors at the
single root `X=-1`, and Frobenius magnifies it to multiplicity `3^4`.

Even forbidding adjacent `11` is insufficient: `D=1731`, whose stripped
unit is `210101_3`, already fails the simple Bessel chamber.  The precise
clean condition is that `D-1`, rather than merely the two ends of `D/3^a`,
lies in the ternary Cantor set.

## 7. Verification and nonclaims

The companion independently checks:

1. the factorial-unit contraction `(23)`;
2. every Bessel face and Cartier recurrence for all `256` Cantor words of
   length at most eight;
3. direct residual gcds for all `64` chamber words of length at most seven;
4. `6,412,320` individual correction terms in Lemma 2.1 for every chamber
   word of length at most six;
5. equality of the **full** moment and terminal Bessel faces at
   `D=21,57,75,165,183,219,237`; and
6. the exact factorization `(29)` for `D=201`.

Reproduce with

```bash
python3 04-computation/factorial_ternary_cantor_bessel_residual_thm3260.py
python3 -O 04-computation/factorial_ternary_cantor_bessel_residual_thm3260.py
```

The computation is an audit of a symbolic proof, not evidence substituted
for unbounded quantifiers.

No claim is made here for:

- odd composite `D` for which `D-1` contains a ternary digit `1`;
- the non-simple but occasionally coprime polygons at such values (for
  example `D=117,129`);
- exact quadratics outside the resonance forced by THM-3124; or
- the full multivariate factorial conjecture.
