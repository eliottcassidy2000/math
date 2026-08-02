---
id: THM-3140
title: "Even-resonance Bessel residual Newton exclusion for quadratic factorial windows"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the resonant
  quadratic factorial-moment pair, every even resonance d=r+2 has one
  2-adic Newton edge of slope one.  After the exact scaling v=X/2, its
  residual polynomial is the Bessel polynomial B_n(X) modulo 2.  The
  Bessel three-term recurrence makes consecutive residuals coprime, so
  every even resonance is excluded uniformly.  Consequently a still-open
  bad resonance must be odd composite.  The proof is symbolic and is not a
  finite extrapolation.
audit: >
  An independent audit rederived the coefficient normalization, endpoint
  valuations, Newton-root sign, unit-residue scaling, Bessel recurrence,
  and exact-quadratic reduction.  It replayed the companion in normal and
  optimized modes, checked the normalization through d=200 and at d=256,510,
  and caught and repaired the general residual lemma's zero-root caveat.
source: root/factorial-composite-newton-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3131-prime-resonance-newton-slope-separation
script: 04-computation/factorial_even_resonance_bessel_residual_thm3140.py
output: 05-knowledge/results/factorial_even_resonance_bessel_residual_thm3140.out
script_sha256: fd1d7f7c1c87d51ee4fe69355571847547e23545de2d928325471771da0a28cf
output_sha256: 45431aa35b7a192026fabede50fa1a380d9e6b1202276d45bd661d6965882c21
hash_basis: LF-normalized bytes
---

# THM-3140 -- even-resonance Bessel residual Newton exclusion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Let

```text
L(t^k)=k!,                    q(t)=a+bt+ct^2,                 (1)
```

with `abc!=0`.  If `r>=1`, `d=r+2` is even, and

```text
L(q^r)=L(q^(r+1))=L(q^(r+2))=0,                              (2)
```

then a contradiction results.  Thus no exact quadratic has a bad
three-moment window at an even resonance.  Combining this with THM-3131,
every still-open bad resonance is odd composite.

The mechanism repairs the precise information loss left open by THM-3131.
At a composite resonance the raw Newton slopes can coincide.  Here they do:
at `p=2` both polynomials have the single slope `1`.  Their residual initial
polynomials, however, are consecutive Bessel polynomials modulo `2`, and
those are coprime.

## 1. A residual-Newton exclusion lemma

We first record the general local step, including the scaling convention.
Let `K` be a complete discretely valued field with uniformizer `pi`, residue
field `k` of characteristic `p`, and normalized valuation `v(pi)=1`.  Suppose
`f,g in K[V]` have a common lower-Newton slope

```text
lambda=a/e,                    gcd(a,e)=1,                    (3)
```

with `e>0` and `p` not dividing `e`.  In the tame totally ramified extension
obtained by adjoining `varpi` with `varpi^e=pi`, normalize the extended
valuation by `w(varpi)=1`.  For `f=sum_i c_i V^i`, put

```text
mu_f=min_i(e v(c_i)-a i),
in_lambda(f)(Y)
  = reduction of varpi^(-mu_f) f(varpi^(-a)Y).                (4)
```

The coefficients in `(4)` are integral.  Only the coefficient indices on
the `lambda`-face survive reduction.  These indices are congruent modulo
`e`, so, for the left endpoint `i_f`,

```text
in_lambda(f)(Y)=Y^(i_f) R_(f,lambda)(Y^e),                   (5)
```

where `R_(f,lambda)` is the usual residual edge polynomial.  Define the
analogous objects for `g`.

If `alpha` is a common root of valuation `-a/e`, then

```text
y=varpi^a alpha
```

is a unit.  Reducing the two scaled root equations gives a **nonzero**
common residue root `bar(y)` of the two initial forms, and hence
`bar(y)^e` is a common root of the two residual edge polynomials.  Therefore

```text
gcd(R_(f,lambda),R_(g,lambda))=1                             (6)
```

excludes a common root of that valuation.  Applying this to every shared
finite slope excludes every common **nonzero** root.  If additionally
`f(0)g(0)!=0`, it excludes every common factor; without that side condition a
common factor `V` would be invisible to finite slopes.  The nonzero-residue
condition is essential: it follows here from the exact valuation of `alpha`,
not merely from integrality.  In the application below `lambda=1`, `e=1`,
and `(13)`--`(14)` give nonzero constant coefficients, so the extension is
trivial, the zero-root caveat is absent, and no wild-ramification issue occurs.

## 2. The exact coefficient and its Bessel residue

THM-3124 reduces `(2)` to the resonance

```text
b/a=-1/d.                                                     (7)
```

After dividing `q` by `a`, putting `u=c/a`, and making the invertible change
`v=du`, define

```text
A_n(v)=d^n L((1-t/d+(v/d)t^2)^n)
      =L((d-t+v t^2)^n) in Z[v].                             (8)
```

The two equations still to be excluded are

```text
A_(d-2)(v)=A_(d-1)(v)=0.                                    (9)
```

Write `A_n(v)=sum_(j=0)^n a_(n,j)v^j`.  Choosing `j` quadratic factors and
then reindexing by `k`, the number of constant `d` factors, gives

```text
a_(n,j)
 =(-1)^(n-j) binom(n,j) sum_(k=0)^(n-j)
     binom(n-j,k)(-d)^k(n+j-k)!

 =(-1)^(n-j)n!2^j sum_(k=0)^(n-j)
     ((-d)^k/k!) binom(n+j-k,2j)(2j-1)!!.                    (10)
```

Here `(-1)!!=1`.  The last binomial-times-double-factorial expression is an
integer.  If `d` is even and `k>=1`, Legendre's formula gives

```text
v_2(d^k/k!)
  =k v_2(d)-(k-s_2(k))
  =(v_2(d)-1)k+s_2(k) >=1,                                  (11)
```

where `s_2(k)` is the binary digit sum.  Hence every `k>=1` term in `(10)`
vanishes after normalized reduction modulo `2`; only `k=0` survives.

Put

```text
h_n=v_2(n!),
b_(n,j)=(n+j)!/(2^j j!(n-j)!)
       =binom(n+j,2j)(2j-1)!! in Z.                          (12)
```

Because `n!/2^h_n` is odd, `(10)`--`(12)` prove, for every `0<=j<=n`,

```text
v_2(a_(n,j)) >= h_n+j,
2^(-(h_n+j))a_(n,j) == b_(n,j)                 (mod 2).      (13)
```

This is the exact coefficient lemma.  It is stronger than merely observing
the lower hull in a finite range.

Both endpoint residues are nonzero:

```text
b_(n,0)=1,                    b_(n,n)=(2n-1)!! ==1 (mod 2). (14)
```

Therefore the entire `2`-adic Newton polygon of `A_n` is the single edge

```text
(0,h_n)----------------(n,h_n+n),          slope 1.          (15)
```

After the slope-one substitution, the integral scaled polynomial

```text
F_n(X)=2^(-h_n)A_n(X/2) in Z_2[X]                            (16)
```

has residual initial polynomial

```text
beta_n(X)=sum_(j=0)^n b_(n,j)X^j in F_2[X].                  (17)
```

The integer polynomial with the coefficients in `(12)` is the classical
Bessel polynomial `B_n(X)` in the normalization

```text
B_0=1,                         B_1=1+X,
B_(n+1)=(2n+1)X B_n+B_(n-1).                                 (18)
```

Reducing `(18)` modulo `2` gives

```text
beta_(n+1)=X beta_n+beta_(n-1).                              (19)
```

Consequently

```text
gcd(beta_n,beta_(n+1))
 =gcd(beta_n,beta_(n-1))
 =...=gcd(beta_1,beta_0)=1                                  (20)
```

for every `n>=0`.  This is the uniform residual separation.

## 3. Exclusion of every even resonance

Take `n=d-2` in `(15)`--`(20)`.  If the two polynomials in `(9)` had a
common complex root, they would have a nonconstant common factor over `Q`.
After passing to `Qbar_2`, choose a common root `alpha`.  The single Newton
edge `(15)` forces

```text
v_2(alpha)=-1.                                                (21)
```

Thus `x=2 alpha` is a unit.  The two equations in `(9)`, scaled as in `(16)`
and reduced in the finite residue extension containing `x`, would make the
nonzero residue `bar(x)` a common root of `beta_(d-2)` and `beta_(d-1)`.
This contradicts `(20)`.  Hence `(9)`, and therefore `(2)`, is impossible.
QED.

Combining THM-3124, THM-3131, and this theorem, a still-open exact-quadratic
bad window must now satisfy

```text
r>=201,              d=r+2 odd composite,              b/a=-1/d. (22)
```

## 4. Hostile boundary and scope

Odd composite resonance is a real boundary of this proof, not an omitted
parity case.  The place `p=2` no longer divides `d`, so `(11)` does not erase
the perturbation terms.  At the smallest stubborn odd controls, the naive
edge-residual observer also retains common factors at every prime divisor:

```text
d=15, p=3: shared slopes 1/3,4/9; residual gcd degrees 1,1;
d=15, p=5: shared slope  1/5;   residual gcd degree  2;
d=33 and d=35: every prime-divisor observer likewise has a
               nonconstant residual gcd on at least one shared face. (23)
```

The denominators in these displayed odd controls are divisible by the
residue characteristic, so they are also wild faces; the tame lemma in
Section 1 is deliberately not overextended to them.  These are failures of
this certificate, not common-root witnesses.  THM-3124 already excludes the
displayed small windows by exact gcd.

No claim is made for odd composite resonances, translated supports, higher
support size, arbitrary `SFC(1)`, ambient `SFC(3)`, or full `FC(3)`.

## 5. Exact verification

The companion independently expands small bivariate powers, verifies the
reindexed coefficient identity, checks `(13)` and the exact hull `(15)` for
every even `4<=d<=160`, verifies the integer Bessel recurrence, computes the
consecutive residual gcds over `F_2`, and freezes the odd hostile controls in
`(23)`.  These finite checks audit the implementation; the theorem is proved
by `(10)`--`(21)`, not by the cutoff.  Run

```text
python3 04-computation/factorial_even_resonance_bessel_residual_thm3140.py
python3 -O 04-computation/factorial_even_resonance_bessel_residual_thm3140.py
```

and compare both byte-for-byte with the declared output.

**End of proof.**
