---
id: THM-3050
title: "Rational product-Gamma radial nullcone and critical Borel order"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  If the squared radius of one circular complex variable is a positive
  rational scale times any finite product of independent
  positive-rational-shape Gamma variables, its polynomial moment nullcone is
  exactly the two strict one-sided charge loci and its polynomial Mathieu
  implication holds.  The scalar moment sequence has sharp Borel order equal
  to the number of Gamma factors, with an explicit critical generalized
  hypergeometric transform.  THM-3047 supplies a new formal-corner family of
  these laws for rational t.  This is one scalar radial grade, not a
  multi-coordinate Wick product or a physical moving-lower resultant.
source: kind-pasteur-2026-08-01-product-gamma-radial-nc2
audit: >
  Two independent hostile audits ACCEPTED the exact-support algebraic
  descent, scalar product-Pochhammer normalization, Kummer/Lucas/Frobenius
  residue, rational-parameter boundary, THM-3047 specialization, Borel-order
  asymptotic and hypergeometric normalization, exponential-integrability
  hostile, and scalar-versus-vector-grade scope.  One independently swept
  288,834 integrality and 44,892 strict-divisibility cells.  Both replayed
  normal, optimized, and stored output after LF normalization, matched the
  declared hashes, and passed the documentation checker.
depends_on:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
related:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3031-the-exponential-period-bridge-to-FC2-requires-value-one-not-zero
script: 04-computation/gmc_rational_product_gamma_radial_nc2_borel_thm3050.py
output: 05-knowledge/results/gmc_rational_product_gamma_radial_nc2_borel_thm3050.out
script_sha256: 4061524cde4485208ae2e0c6a65c579c6f2b05f92eac8d09a3d0c1b33b23aef5
output_sha256: f539139a2284cb12f68141e59579d41a2da479df28834769d90665c92846b1fc
hash_basis: LF-normalized bytes
---

# THM-3050 -- rational product-Gamma radial nullcones

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2022 extends the two-dimensional Gaussian nullcone proof to one rational
Gamma radial factor.  THM-3047 unexpectedly supplies products of many Gamma
factors as the exact positive-real model of the formal resultant width flag.
The two mechanisms fit because all factors see the same scalar radial grade.
Their product strengthens, rather than obstructs, the prime-block separation
at the lowest balanced face.

The same moments also explain why the ordinary exponential-period operation
is mismatched here.  A product of `J` Gamma factors has Borel order `J`, not
order one.

## 1. Product-Gamma circular laws

Fix an integer `J>=1` and positive rational numbers

```text
c, alpha_1,...,alpha_J in Q_(>0).                        (1)
```

Let `X_j` be independent `Gamma(alpha_j,1)` variables, let `U` be independent
and uniform on the unit circle, and define

```text
T=c product_(j=1)^J X_j,             Z=sqrt(T) U.        (2)
```

For `A,B>=0`, angular orthogonality and independence give

```text
E[Z^A Zbar^B]=0                         if A!=B,
E[Z^A Zbar^A]=w(A),

w(n)=c^n product_(j=1)^J (alpha_j)_n.                   (3)
```

Here `(alpha)_n` is the rising factorial.

> **Product-Gamma radial nullcone theorem.**  For the moment functional in
> `(3)` and every polynomial `P in C[Z,Zbar]`,
>
> ```text
> E[P^m]=0 for every m>=1
> ```
>
> holds if and only if every charge `A-B` in the exact support of `P` is
> strictly positive, or every charge is strictly negative.  Consequently,
> for every fixed polynomial `Q`,
>
> ```text
> E[Q P^m]=0
> ```
>
> for all sufficiently large `m` whenever all positive moments of `P`
> vanish.

The converse and the Mathieu conclusion are immediate from charge addition.
The content is exclusion of every support whose charge convex hull contains
zero.

## 2. Prime-block proof

The moment polynomials for `(3)` have rational coefficients.  Therefore the
algebraic specialization and exact-support localization in THM-2022 Section 2
apply unchanged: a complex torus null point on a fixed support would give one
over a number field, with every support coefficient still nonzero.

Apply THM-2022's lowest balanced face construction.  It supplies a face `F`,
an integer return length `m_0`, a nonzero face constant term `Q`, and one
common scalar radial grade `A_0` for every balanced channel contributing to
`Q`.  At length `p m_0`, every balanced channel has grade at least `pA_0`.

Write each shape as

```text
alpha_j=h_j/k_j,
```

and choose a sufficiently large rational prime `p>m_0` which is good for the
specialized coefficients and `Q`, is coprime to all `k_j`, and is a unit for
the numerator and denominator of `c`.  For every integer `n>=pA_0`,

```text
w(n)/w(pA_0)
 =c^(n-pA_0)
  product_(j=1)^J product_(ell=pA_0)^(n-1)
       (h_j+k_j ell)/k_j                               (4)
```

is `p`-integral.  More strongly, if `A'>A_0`, each interval

```text
pA_0 <= ell < pA'                                      (5)
```

is a union of `A'-A_0` complete blocks of length `p`.  Since `k_j` is a
unit modulo `p`, each block contains exactly one solution of
`h_j+k_j ell=0 mod p`.  Hence

```text
v_p(w(pA')/w(pA_0)) >= J(A'-A_0)>0.                    (6)
```

Normalize the moment of order `pm_0` by `w(pA_0)`.  Formula `(4)` ensures
that every channel term is integral.  Kummer's multinomial carry theorem
kills every channel whose allocation vector is not componentwise divisible
by `p`.  Among the dilated channels, `(6)` kills every strict off-face grade.
The face channels have ratio one.  Lucas and Frobenius therefore leave
exactly the nonzero residue `Q^p`, exactly as in THM-2022 equations
`(13)--(15)`.  This contradicts nullity and proves the theorem.

The factor count `J` only strengthens `(6)`.  THM-2022 Section 9's hostile
concerns a vector of independent factorial **grades**.  Here the internal
Gamma factors all depend on the same scalar grade `A`; their product is the
single scalar weight `w(A)`.  No coordinate can move down while another moves
up.

## 3. THM-3047 becomes a radial NC2 family

For the `k`-slot character of THM-3040/3047, with `k>=2`, put

```text
A_k=k!(H_k-1),
B_k=k!(k+1-2H_k),
I_k=A_k+B_k=k!(k-H_k).                                  (7)
```

Fix rational `t>0`.  Specialize `(1)--(2)` by taking

```text
J=I_k,
c=t^I_k,
alpha=1/t       with multiplicity A_k,
alpha=1/t+1     with multiplicity B_k.                  (8)
```

Then `(3)` is exactly THM-3047's universal width factor

```text
w(M)=F_M^(k)(t).                                        (9)
```

Thus every rational positive-real formal width flag gives one explicit
circular product-Gamma law with the exact one-sided nullcone and Mathieu
property.  At `k=2`, `(A,B,I)=(1,0,1)`, recovering the single-Gamma family in
THM-2022 up to rational scale.  For `k>=3`, already `I_3=7`; these are genuine
multi-Gamma radial laws not stated in THM-2022.

Rationality of `t` is load-bearing for this elementary finite-place proof.
For transcendental `t`, fixed-parameter algebraic specialization need not
produce number-field coefficients.  For algebraic irrational `t`, a number
field is available, but the rational arithmetic-progressions in `(4)--(6)`
must be replaced by a separate prime-ideal/splitting argument.  No such
extension is claimed here.

## 4. Sharp Borel-order spectrum

The analytic statement is valid for arbitrary positive real `c,alpha_j`,
without the rationality used above.  For real `sigma`, define

```text
B_sigma(z)=sum_(n>=0) w(n) z^n/(n!)^sigma.              (10)
```

Since

```text
(alpha)_n = Gamma(n+alpha)/Gamma(alpha)
          =Theta(n^(alpha-1)n!),                        (11)
```

the `n`th root of the coefficient in `(10)` is asymptotic to

```text
c (n/e)^(J-sigma).                                      (12)
```

Therefore the radius of convergence is exactly

```text
0             if sigma<J,
1/c           if sigma=J,
infinity      if sigma>J.                               (13)
```

At the critical order one has the explicit generalized hypergeometric
function

```text
B_J(z)
 = _J F_(J-1)(alpha_1,...,alpha_J;
              1,...,1; c z),                           (14)
```

with `J-1` denominator parameters equal to one.  All coefficients after the
constant term are positive, so

```text
B_J(z)>1                 for 0<z<1/c.                   (15)
```

For the THM-3047 family, the critical order is `I_k`, the critical radius is
`t^(-I_k)`, and `(14)` has `A_k` copies of `1/t` and `B_k` copies of
`1/t+1`.  At `k=2`, this reduces to

```text
B_1(z)=(1-tz)^(-1/t).                                   (16)
```

For `k>=3`, the ordinary exponential/Borel-one series has radius zero.

## 5. Why an ordinary exponential-period argument does not replace Frobenius

THM-3031's value-one bridge for the factorial conjecture uses a polynomial on
a compact interval, so its exponential integral is entire and its moment
series may be integrated termwise.  That analytic operation is unavailable
for the present radial nullcone.

The obstruction already occurs for the algebraic one-sided polynomial

```text
P=Z^d,                  d>=3.                            (17)
```

All positive polynomial moments of `P` vanish by angular charge.  Yet for
every `s!=0`,

```text
E[|exp(sP)|]=infinity.                                  (18)
```

Choose an angular arc on which `Re(sU^d)` has a fixed positive lower bound,
and restrict `X_2,...,X_J` to compact positive intervals.  The remaining
`X_1` integral dominates

```text
integral exp(C x^(d/2)-x) x^(alpha_1-1) dx,
```

which diverges because `d/2>1`.  Thus the formal ordinary exponential series
equals one coefficientwise but does not define an absolutely integrable
exponential observable.  The matched transform is `(14)`, not the Borel-one
period, and the nullcone proof genuinely uses prime-local Frobenius.

## 6. Scope

This theorem gives new circular **one-complex-variable** radial laws.  It does
not prove a higher-dimensional Gaussian conjecture, identify the internal
Gamma factors with independent Wick coordinates, or bypass the explicit
GMC(3) counterexamples.  It also does not turn THM-3047's formal corner into a
physical resultant width or control the extra transport from moving lower
offsets.

The nullcone theorem uses rational parameters.  The Borel spectrum uses only
positive real parameters.  These quantifiers must not be interchanged.

## 7. Exact companion

The dependency-free rational companion checks:

- `2272` exact prime-block integrality cells and `192` strict off-face
  divisibility cells for one through four Gamma factors, including
  nonintegral shapes and scales;
- `52` critical hypergeometric coefficient recurrences, and records the
  proof-level three-way radius law;
- `84` exact THM-3047 specialization cells for `k=2..5` and three rational
  positive values of `t`;
- balanced positive and strict one-sided null controls, and records the
  proof-level ordinary exponential-integrability hostile.

Reproduce with

```text
python 04-computation/gmc_rational_product_gamma_radial_nc2_borel_thm3050.py
python -O 04-computation/gmc_rational_product_gamma_radial_nc2_borel_thm3050.py
```

Both runs equal the stored six-line transcript after LF normalization.

**QED.**
