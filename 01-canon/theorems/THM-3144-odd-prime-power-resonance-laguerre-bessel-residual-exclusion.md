---
id: THM-3144
title: "Odd-prime-power resonance Laguerre--Bessel residual exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every odd prime
  p and every a>=2, the resonant quadratic factorial-moment pair at d=p^a
  has exactly one shared p-adic Newton slope, namely 1/(p-1).  A negative-
  parameter Laguerre coefficient identity and a Kummer carry bound show that
  the full coefficients have the terminal Bessel lower faces.  The two
  residual edge polynomials obey a scalar-prefix-plus-top-monomial identity
  and are coprime.  Thus every odd prime-power resonance is excluded
  uniformly.  The proof is symbolic and is not a finite extrapolation.
audit: >
  An independent audit rederived the Laguerre identity, Kummer equality
  case, both digit-sum polygons, the p=3 endpoint, residual scaling and gcd,
  and the Newton application.  It extended hostile probes through 3^6,
  5^4, 7^3, and 23^2, confirmed both exceptional lifts are exactly one,
  and replayed normal and optimized transcripts and hashes byte-for-byte.
source: root/factorial-composite-newton-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3140-even-resonance-bessel-residual-newton-exclusion
related:
  - THM-3131-prime-resonance-newton-slope-separation
script: 04-computation/factorial_odd_prime_power_laguerre_residual_thm3144.py
output: 05-knowledge/results/factorial_odd_prime_power_laguerre_residual_thm3144.out
script_sha256: 6bcb9584197639cbf6e0a3d76ecec057a4b4ef7533da6230ba1aa1c9b5fb36fe
output_sha256: 5b26e4f8999189c128dc85bd015338fcab057123e0f6228845c7bb5f9531b789
hash_basis: LF-normalized bytes
---

# THM-3144 -- odd-prime-power Laguerre--Bessel residual exclusion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Let

```text
L(t^k)=k!,                    q(t)=a_0+b t+c t^2,             (1)
```

with `a_0 b c!=0`.  Suppose a three-moment window begins at `r>=1` and its
THM-3124 resonance

```text
d=r+2                                                        (2)
```

is an odd prime power `p^a`, with `p` odd and `a>=2`.  Then

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                               (3)
```

cannot all vanish.

Together with THM-3131 and THM-3140, this leaves only odd composite
resonances having at least two distinct prime divisors.  Together also with
THM-3124's finite-exact census, any still-open exact-quadratic bad window
must satisfy

```text
r>=201,   d=r+2 odd,   omega(d)>=2,   and   b/a_0=-1/d,       (4)
```

where `omega(d)` is the number of distinct prime divisors.

## 1. Resonant pair and a Laguerre coefficient identity

Put `D=p^a`.  As in THM-3124, normalize by `a_0`, put `u=c/a_0`, make the
invertible change `v=Du`, and define

```text
A_n(v)=L((D-t+v t^2)^n)=sum_(j=0)^n a_(n,j)v^j.              (5)
```

It is enough to prove that

```text
A_(D-2) and A_(D-1) have no common root.                     (6)
```

Direct expansion gives

```text
a_(n,j)=binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)D^(n-j-ell)(-1)^ell(2j+ell)!.                (7)
```

There is a useful exact negative-parameter Laguerre form.  With

```text
L_m^(alpha)(x)=sum_(k=0)^m binom(m+alpha,m-k)(-x)^k/k!,      (8)
```

putting `m=n-j` and using

```text
binom(-2j-1,m-k)=(-1)^(m-k)binom(2j+m-k,m-k)                 (9)
```

in `(7)` yields

```text
a_(n,j)=n!(2j)!/j! * L_(n-j)^(-n-j-1)(-D).                 (10)
```

Thus the cancellation problem is a concrete congruence for generalized
Laguerre values outside the usual positive-parameter PF cone.

For the valuation proof it is better to isolate the `k=0` summand.  Define

```text
b_(n,j)=(n+j)!/(2^j j!(n-j)!)
       =binom(n+j,2j)(2j-1)!!,                               (11)
T_(n,j)=(-1)^(n-j)n!2^j b_(n,j).                            (12)
```

Then `T_(n,j)` is the terminal Bessel contribution to `a_(n,j)`, and `(7)`
can be written

```text
a_(n,j)=T_(n,j)(1+sum_(k=1)^(n-j) R_(n,j,k)),               (13)

R_(n,j,k)=((-D)^k/k!)
  *binom(n+j-k,2j)/binom(n+j,2j)
 =((-D)^k/k!)*binom(n-j,k)/binom(n+j,k).                    (14)
```

The ratios in `(14)` are rational, but the next lemma proves they are
`p`-integral in exactly the range needed here.

## 2. The Kummer ratio lemma and its only exceptions

Let

```text
n=D-delta,       delta in {1,2},
N=n+j,           m=n-j,           1<=k<=m.                  (15)
```

Since `N<2D<p^(a+1)`, Kummer's theorem shows

```text
v_p(binom(N,k))<=a:                                          (16)
```

there can be at most one base-`p` carry at each of the `a` lower digit
positions, and no carry out of the top digit.  Also

```text
ak-v_p(k!)>=a,                                               (17)
```

Indeed, `v_p(k!)<=k-1` proves `(17)`.  If `k>=2`, oddness of `p` gives the
stronger bound `v_p(k!)<k-1<=a(k-1)`, so `(17)` is strict.
Equations `(14)`--`(17)` give

```text
v_p(R_(n,j,k))
 =ak-v_p(k!)+v_p(binom(m,k))-v_p(binom(N,k)) >=0.            (18)
```

Equality in `(18)` forces `k=1`.  It then forces

```text
v_p(N)=a,       v_p(m)=0,       hence N=D and j=delta.       (19)
```

Conversely `(19)` does give equality.  At this single coefficient,

```text
1+R_(n,delta,1) ==2delta+1                    (mod p),       (20)
```

while every `k>=2` term vanishes modulo `p`.  Therefore:

- for every `j!=delta`, `a_(n,j)` has the same valuation and normalized
  residue as `T_(n,j)`;
- at `j=delta`, its valuation can only increase;
- the increase occurs modulo `p` precisely at the hostile pairs
  `(delta,p)=(1,3)` and `(2,5)`.

Both hostile coefficients are strictly above the lower Newton faces found
below.  They cannot change a face or a residual edge polynomial.  This is
the precise reason the termwise-Bessel argument, which is false before
normalization, nevertheless survives on the tropical boundary.

## 3. The two exact Bessel Newton polygons

Because `p` is odd, the power `2^j` in `(11)` is a unit.  Legendre's formula
gives, with `s_p` denoting base-`p` digit sum,

```text
v_p(b_(n,j))
 =[j+s_p(j)+s_p(n-j)-s_p(n+j)]/(p-1).                       (21)
```

### 3.1. The `n=D-1` polygon

For `1<=j<=D-1`, digit complementation and the absence of a carry across
`D+(j-1)` give

```text
s_p(D-1-j)=a(p-1)-s_p(j),
s_p(D-1+j)=1+s_p(j-1).                                      (22)
```

Hence

```text
v_p(b_(D-1,j))
 =j/(p-1)+[a(p-1)-1-s_p(j-1)]/(p-1)
 >=j/(p-1).                                                  (23)
```

Equality holds exactly when `j-1` has digit sum `a(p-1)-1`, namely

```text
j=D-p^rho,                   0<=rho<a,                       (24)
```

as well as the constant index `j=0`.  Therefore the Bessel lower hull is the
single edge

```text
(0,0)----------------(D-1,(D-1)/(p-1)),   slope 1/(p-1).    (25)
```

### 3.2. The `n=D-2` polygon

The coefficients at `j=0,1` are units.  For `2<=j<=D-2`, use

```text
s_p(D-2-j)=a(p-1)-s_p(j+1),
s_p(D-2+j)=1+s_p(j-2),
s_p(j+1)=s_p(j)+1-(p-1)v_p(j+1).                            (26)
```

Subtracting the line through `(1,0)` of slope `1/(p-1)` from `(21)` gives

```text
v_p(b_(D-2,j))-(j-1)/(p-1)
 =a+v_p(j+1)-[1+s_p(j-2)]/(p-1) >=0.                        (27)
```

Equality in `(27)` holds exactly at

```text
j=D-p^rho+1,                1<=rho<a,                       (28)
```

and `j=1` is also on the line.  The last equality point is `j=D-p+1`.  If
`p>3`, equation `(21)` gives the leading valuation `(D-1)/(p-1)` at
`j=D-2`; the integral valuations and `(27)` put each of the intervening
points above the chord to that leading coefficient.  If `p=3`, the last
equality point is already `j=D-2`, with valuation `(D-3)/2`.  Thus for
`p>3` the lower hull is

```text
(0,0)--(1,0)--(D-p+1,(D-p)/(p-1))
                   --(D-2,(D-1)/(p-1)),                     (29)
```

where the final edge has slope `1/(p-3)`.  For `p=3`, omit the final vertex
and edge from `(29)`; the third vertex is then the leading vertex.

By Section 2, the full coefficients `a_(n,j)` lie at least as high as these
Bessel points and agree at every vertex and every face index.  If

```text
h=v_p((D-2)!)=v_p((D-1)!),                                  (30)
```

the actual Newton polygons of the two polynomials in `(6)` are exactly
`(25)` and `(29)`, translated vertically by `h`.  Their only common slope is

```text
lambda=1/(p-1).                                              (31)
```

Its denominator is prime to `p`, so this is a tame residual face.

## 4. Scalar-prefix residual separation

Put

```text
e=p-1,                         M=(D-1)/e.                    (32)
```

For the common face, write the standard residual polynomials as

```text
R_1(Z)=sum_(k=0)^M u_(1,k)Z^k          for A_(D-1),
R_0(Z)=sum_(k=0)^(M-1) u_(0,k)Z^k      for A_(D-2),          (33)
```

where `R_1` uses coefficient indices `j=ke` and `R_0` uses `j=1+ke`.
Equations `(24)` and `(28)` show that their nonzero supports align through
degree `M-1`; `R_1` has one additional nonzero top coefficient in degree
`M`.

At an aligned nonzero index `j=0` or `j=D-p^rho` with `1<=rho<a`, Section 2
allows us to use the terminal units.  Directly from `(11)`--`(12)`,

```text
u_(0,k)/u_(1,k)
 ==[(D-1-j)(D-2-j)]/[(D-1)(j+1)]
 ==-2                                                   (mod p). (34)
```

Thus, for a nonzero `gamma in F_p`,

```text
R_1(Z)=(-1/2)R_0(Z)+gamma Z^M.                              (35)
```

The constant coefficient of `R_0` is nonzero.  Any common divisor in `(35)`
would divide `Z^M`, but `Z` does not divide `R_0`.  Therefore

```text
gcd(R_0,R_1)=1 in F_p[Z].                                   (36)
```

The tame residual-Newton lemma of THM-3140 now excludes a common root of the
only shared valuation `(31)`.  The raw slope sets exclude every other
valuation.  Hence the two polynomials in `(6)` are coprime over `Q`, proving
the theorem.  QED.

## 5. Hostile boundary and scope

The two exceptional equal-valuation perturbations are not suppressed:

```text
n=D-1, j=1, p=3:     1+R_(n,j,1) ==0 (mod 3),
n=D-2, j=2, p=5:     1+R_(n,j,1) ==0 (mod 5).               (37)
```

They demonstrate why a claim of coefficientwise equality with the Bessel
polynomial would be false.  What survives is exactly the lower-face equality
needed by the proof.  At odd composites with distinct prime divisors this
prime-power digit complement breaks.  The smallest stubborn control `D=15`
has nonconstant residual gcds at both `p=3` and `p=5`; it is a certificate
failure, not a common-root witness, and its window is already finite-exact
empty by THM-3124.

No claim is made for odd composites with `omega(D)>=2`, translated supports,
higher support size, arbitrary `SFC(1)`, ambient `SFC(3)`, or full `FC(3)`.

## 6. Exact verification

The companion independently checks the Laguerre identity, the ratio-valuation
lemma and its unique equality positions, the two hostile cancellations, the
predicted hulls, the scalar-prefix residual identity, and residual gcds.  Its
finite universe contains every odd prime power through `400`, including
`3^5`, `5^3`, and `7^3`; small cases also use a repeated-product expansion.
These checks audit the implementation.  The theorem is proved uniformly by
`(10)`--`(36)`, not by the cutoff.  Run

```text
python3 04-computation/factorial_odd_prime_power_laguerre_residual_thm3144.py
python3 -O 04-computation/factorial_odd_prime_power_laguerre_residual_thm3144.py
```

and compare both byte-for-byte with the declared output.

**End of proof.**
