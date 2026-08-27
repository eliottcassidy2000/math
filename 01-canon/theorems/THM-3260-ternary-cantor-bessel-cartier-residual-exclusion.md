---
id: THM-3260
title: "Ternary-Cantor Bessel--Cartier residual exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let D>=3 be odd and divisible by 3.  If every
  ternary digit of D-1 is 0 or 2, then the resonant exact-quadratic
  factorial-moment pair at D is coprime.  Both full moment polynomials have
  their terminal Bessel lower faces; their only shared Newton slope is 1/2.
  The residual edges form a two-sector Cartier module over F_3, and exact
  Frobenius recurrences prove the sectors coprime for every such ternary
  word.  This gives an infinite mixed-prime closure, not a finite
  extrapolation.
audit: >
  An independent audit rederived the 12-state carry potential, all local
  induction inequalities and strictness, the factorial-unit contraction,
  both Cartier digit recurrences, residual gcd induction, Newton-edge use,
  and the D=201 Frobenius hostile.  Fresh normal and optimized runs byte-match
  the stored transcript and declared hashes after the THM-3260 repair now
  recorded as MISTAKE-523.
source: root/factorial-composite-newton-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3204-odd-prime-power-resonance-laguerre-bessel-residual-exclusion
related:
  - THM-3200-even-resonance-bessel-residual-newton-exclusion
  - THM-3131-prime-resonance-newton-slope-separation
script: 04-computation/factorial_ternary_cantor_bessel_residual_thm3260.py
output: 05-knowledge/results/factorial_ternary_cantor_bessel_residual_thm3260.out
script_sha256: 94daa134712b930a752e273365b6b93ce4c21590c0ac2833640884ba1cae3a10
output_sha256: 2cfd78dbcd49f989637090d501bf357f1f7a00cc15a630321b8db53ca7e6fa8f
hash_basis: LF-normalized bytes
---

# THM-3260 -- ternary-Cantor Bessel--Cartier residual exclusion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

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

We give the complete finite-state certificate for the digit step.  For
nonnegative integers `A,B,C`, a Cantor integer `M`, and

```text
A+B+C+rho=M,                   rho in {0,1,2},               (12)
```

run two base-`3` additions in parallel:

1. double `B`, with incoming carry `d in {0,1}`;
2. add `A` to the output digits of that doubling, with incoming carry
   `g in {0,1}`.

Let `D_d(B)` and `G_(d,g)(A,B)` be the total numbers of outgoing carries in
the first and second additions.  Define the twelve-state potential

```text
H(rho,d,g)=0                     if rho=0,
H(rho,d,g)=rho-d-g               if rho=1 or 2.              (13)
```

The exact carry inequality is

```text
s(B)-C-2(D_d(B)+G_(d,g)(A,B)) <= H(rho,d,g).                (14)
```

Here is a fully explicit induction proof.  Write the low digits of
`(A,B,C,M)` as `(alpha,beta,gamma,epsilon)`, where every digit is in
`{0,1,2}` and `epsilon in {0,2}`.  The next state is determined by

```text
alpha+beta+gamma+rho = epsilon+3rho',
2beta+d                  = f+3d',
alpha+f+g                = h+3g'.                            (15)
```

After deleting the low column, the left side of `(14)` is

```text
beta-gamma-2(d'+g')
 + [s(B')-C'-2(D_(d')(B')+G_(d',g')(A',B'))] -2C'.          (16)
```

Thus induction follows from the local inequality

```text
beta-gamma-2(d'+g')+H(rho',d',g') <= H(rho,d,g).            (17)
```

There are `12` incoming states.  For each state, exactly `18` quadruples
`(epsilon,alpha,beta,gamma)` satisfy the first equation of `(15)`.  The
following table gives the maximum of the left side of `(17)` over all
those transitions.  Equality with the displayed potential in every row is
the complete transition certificate.

```text
(rho,d,g)   H    transition maximum       (rho,d,g)   H    transition maximum
--------------------------------------------------------------------------
(0,0,0)     0           0                 (0,0,1)     0           0
(0,1,0)     0           0                 (0,1,1)     0           0
(1,0,0)     1           1                 (1,0,1)     0           0
(1,1,0)     0           0                 (1,1,1)    -1          -1
(2,0,0)     2           2                 (2,0,1)     1           1
(2,1,0)     1           1                 (2,1,1)     0           0.        (18)
```

The base case has `A=B=C=M=0`, hence `rho=0`; `(14)` is immediate for all
four incoming pairs `(d,g)`.  Equations `(15)`--`(18)` therefore prove
`(14)` for arbitrary word length.  The companion enumerates all `216`
valid transitions and checks `(17)` exactly; this is a finite certificate
of the induction, not a bounded integer census.

We now apply the carry inequality.  Kummer's theorem and

```text
C(m,j)=((2j)!/j!) binom(m+j,2j)                              (19)
```

give, with `A=m-j` and `B=j`,

```text
2v(C(m,j))
 =j-s(j)+2(D_0(j)+G_(0,0)(m-j,j)).                          (20)
```

First take `(A,B,C,rho)=(N-k-j,j,k,0)` in `(14)`.  Equations
`(14)` and `(20)` give the stronger shifted-Bessel bound

```text
2v(C(N-k,j)) >= j-k.                                        (21)
```

Next take `(A,B,C,rho)=(N-1-k-j,j,k,1)`.  Since
`H(1,0,0)=1`, the same argument gives

```text
2v(C(N-1-k,j)) >= j-k-1.                                   (22)
```

Finally Legendre's formula gives

```text
2v(D^k/k!)=(2a-1)k+s(k).                                    (23)
```

For `k>0`, the difference between the right side of `(23)` and `k` is

```text
2(a-1)k+s(k)>0.                                             (24)
```

Adding `(21)` or `(22)` to `(23)` proves `(9)` and `(10)`, including
strictness for every `k>0`.  Equation `(11)` is immediate.  At `k=0`, an
integer valuation cannot equal `j/2` for odd `j`, or `(j-1)/2` for even
`j`; this gives the equality parity assertion and completes the proof.

### Consequence for the full moment pair

Let `h=v((D-1)!)=v((D-2)!)`.  Lemma 2.1 and `(8)` show that every
`k>0` summand lies strictly above the following lower-face lines, while the
`k=0` Bessel term supplies their residual coefficients:

```text
A_(D-1):       y=h+j/2,
A_(D-2):       y=h                 on 0<=j<=1,
               y=h+(j-1)/2        on 1<=j<=D-2.             (25)
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
  overline(3^(-k) C(M-1,2k+1)) X^k.                         (26)
```

A coefficient in `(26)` is declared `0` when its valuation is strictly
larger than the displayed power.  The `k=0` case of the digit argument in
Lemma 2.1, which does not require the final digit to be `2`, proves that
all displayed divisions are integral.

For `M=N`, up to separate nonzero row scalars, `P_N` is the residual edge
of `A_(D-1)` and `Q_N` the slope-`1/2` residual edge of `A_(D-2)`.

## 4. Exact Cartier recurrences

We first give the coefficient statement from which the polynomial
recurrences follow.  If an integer `z` satisfies `v(z)>=e`, write

```text
[z]_e = overline(3^(-e)z) in F_3.                            (27)
```

Thus `[z]_e=0` when `v(z)>e`.  Put

```text
p_(M,l)=[C(M,2l)]_l,
q_(M,l)=[C(M-1,2l+1)]_l,                                   (28)
```

and set either symbol to zero outside its natural range.

### Lemma 4.1 (coefficient Cartier contraction)

For every even ternary Cantor integer `M`, including both endings
`M==0,2 (mod 3)`, the following twelve identities hold:

```text
residual index k     3ell                 3ell+1               3ell+2
----------------------------------------------------------------------------
p_(3M,k)             p_(M,ell)           -q_(M,ell)            0
q_(3M,k)             0                    -q_(M,ell)            0
p_(3M+2,k)           p_(M,ell)            r_(M,ell)            0
q_(3M+2,k)          -p_(M,ell)            0                     0,      (29)
```

where

```text
r_(M,ell)=p_(M,ell)-q_(M,ell)       if M==0 (mod 3),
r_(M,ell)=-p_(M,ell)-q_(M,ell)      if M==2 (mod 3).          (30)
```

### Proof

Let `F(m)=v(m!)` and let `U(m)` be the unit part of `m!` modulo `3`.
Complete blocks of three give the exact one-digit contraction

```text
F(3m+r)=m+F(m),
U(3m+r)=(-1)^m U(m) eta_r,
eta_0=eta_1=1, eta_2=-1.                                   (31)
```

Apply `(31)` to the numerator factorial and the two denominator factorials
in

```text
C(n,j)=(n+j)!/[j!(n-j)!].                                  (32)
```

In particular, the bracket is computed without any choice of lift by

```text
[C(n,j)]_e
 = 1_(F(n+j)-F(j)-F(n-j)=e)
   U(n+j) U(j)^(-1) U(n-j)^(-1),                            (32a)
```

because the shifted-Bessel bounds exclude a valuation below `e` in every
entry used below.

For auditability, the low remainders of those three factorial arguments
are displayed below.  Each entry is `(n+j,j,n-j) mod 3`; the column is the
residue class of the residual index `k`, so the Bessel index is `j=2k` in
a `p` row and `j=2k+1` in a `q` row.

```text
target row       k=3ell             k=3ell+1           k=3ell+2
-----------------------------------------------------------------------
p_(3M,k)          (0,0,0)            (2,2,1)            (1,1,2)
q_(3M,k)          (0,1,1)            (2,0,2)            (1,2,0)
p_(3M+2,k)        (2,0,2)            (1,2,0)            (0,1,1)
q_(3M+2,k)        (2,1,0)            (1,0,1)            (0,2,2).       (33)
```

The shifted-Bessel bounds `(21)`--`(22)`, with `k=0`, show before division
that every bracket in `(28)` and `(29)` is integral.  Substitution of
`(31)` using the twelve remainder triples `(33)` gives, respectively,

```text
(p,-q,0), (0,-q,0), (p,r,0), (-p,0,0).                     (34)
```

Here is the only entry in `(34)` not obtained by immediate cancellation of
the three `eta` factors.  For the middle component of `p_(3M+2,k)`, the
remainders are `(1,2,0)`.  After the first contraction, reduce the final
digit of `M`; `(31)` gives

```text
[C(3M+2,6ell+2)]_(3ell+1)
 = p_(M,ell)-q_(M,ell)              if M==0 (mod 3),
 =-p_(M,ell)-q_(M,ell)              if M==2 (mod 3).         (35)
```

For the six zero entries, the valuation part of `(31)` is at least one
larger than the displayed normalizing exponent; if the neighboring
quotient lies on its face, the low remainder in `(33)` supplies that extra
unit, and otherwise the shifted-Bessel bound already supplies it.  This is
also read directly from `(34)`, where a bracket of positive valuation
reduces to zero.  Thus `(34)` is precisely the twelve identities `(29)`,
for both possible final digits of `M`.  This proves the lemma.

The pair `(P_M,Q_M)` therefore contracts one ternary digit at a time.  In
polynomial form Lemma 4.1 says

```text
P_(3M)   = P_M^3-X Q_M^3,             Q_(3M)   =-X Q_M^3,   (36)

P_(3M+2) = P_M^3+X R_M^3,             Q_(3M+2) =-P_M^3,     (37)
```

where

```text
R_M=P_M-Q_M       if M==0 (mod 3),
R_M=-(P_M+Q_M)    if M==2 (mod 3).                           (38)
```

These are Frobenius cubes in `F_3[X]`, not suggestive fits.

The base pair is

```text
P_0=1,                         Q_0=0.                         (39)
```

## 5. Coprimality is an invariant of the digit tree

We prove by induction that

```text
gcd(P_M,Q_M)=1                                                     (40)
```

for every ternary Cantor `M`.

If `M=3m`, equations `(36)` say that a common nonzero root must make both
`Q_m` and `P_m` vanish.  The root `X=0` is impossible because
`P_(3m)(0)=P_m(0)^3=1`.

If `M=3m+2`, equation `(37)` first forces `P_m=0`.  At such a root,
`R_m` is a nonzero scalar multiple of `Q_m`, so the first equation in
`(37)` forces `Q_m=0`; again induction excludes it.  Thus `(40)` follows
from `(39)`.

When `N` ends in `2`, `(37)` also shows that `Q_N` has the full endpoint
degree `(N-2)/2`, while `P_N` has degree `N/2`.  Consequently the exact
lower hulls in `(25)` are

```text
A_(D-1): (0,h)----------------------------(D-1,h+(D-1)/2),

A_(D-2): (0,h)---(1,h)--------------------(D-2,h+(D-3)/2).  (41)
```

Their only shared slope is `1/2`.  Its denominator is prime to `3`, so the
ordinary residual-polynomial Newton theorem applies.  Equation `(40)` says
that the two residual edges are coprime.  The two full polynomials therefore
have no common algebraic root, completing the proof of the theorem.

## 6. Why the tempting broader chamber fails

The earlier empirical phrase “after stripping powers of `3`, the leading
ternary digit is `2` and the trailing digit is `1`” is false.  The first
hostile inside that broader language is

```text
D=201=21110_3,                    D/3=2111_3.                (42)
```

At the common slope `1/2`, the exact full residual polynomials are

```text
Q=1+X^81=(1+X)^81,
P=1+X+X^81+X^82=(1+X)Q.                                   (43)
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

1. the `216` transitions in the unbounded carry certificate `(17)`--`(18)`;
2. the factorial-unit contraction `(31)`;
3. every Bessel face and Cartier recurrence for all `256` Cantor words of
   length at most eight;
4. direct residual gcds for all `64` chamber words of length at most seven;
5. `6,412,320` individual correction terms in Lemma 2.1 for every chamber
   word of length at most six;
6. equality of the **full** moment and terminal Bessel faces at
   `D=21,57,75,165,183,219,237`; and
7. the exact factorization `(43)` for `D=201`.

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
