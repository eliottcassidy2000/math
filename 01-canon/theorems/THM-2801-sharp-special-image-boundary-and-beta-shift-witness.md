---
id: THM-2801
title: "Sharp Special Image Conjecture boundary and beta-shift witness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The Special
  Image Conjecture is true in dimension one and false in every dimension
  n>=2.  An explicit bidegree-(4,4) polynomial on the rank-one 2x2 Segre
  cone has every positive power in ker(E_2), while multiplication by the
  coordinate z_2 leaves the kernel at every power.  Its complete marked
  contraction is a beta integral, and every coordinate power z_2^s is a
  persistent observer once m>=s.  Independently, a four-term
  bidegree-(2,2) polynomial already refutes SIC(3), inside an all-r family
  with a closed mixed-output formula.  These witnesses refute the full
  Mathieu-space route but do not lie in Zhao's xi-linear cubic or separable
  Hessian sectors and prove nothing about JC(2).
source: root-sic-jc-lrc-2026-07-28
attribution: >
  The two displayed r=3 formulas and the two-pair formula were supplied by
  the user without a source.  The cited arXiv:2607.23450 is unrelated.
  This file proves and sharpens the formulas but makes no discovery-priority
  claim.  The beta-shift ladder, coordinate-only multiplier, all-r family,
  and r=1 four-term reduction are project-internal deductions.
depends_on:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
related:
  - THM-1435-zhao-vc-witness-transport-machinery-and-the-closed-shortcut
  - THM-1490-the-gaussian-moment-counterexample-verified-proved-shortened-and-obstructed
  - THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate
script: 04-computation/sic_sharp_boundary_beta_shift_thm2801.py
output: 05-knowledge/results/sic_sharp_boundary_beta_shift_thm2801.out
script_sha256: bdaf2a070aa3495bb710916820cc6be873e660aa0aa64f2eaca067226194e035
output_sha256: faf02e5c10cca6d9dbb7654b82e1adf47fe6c123e01a082d64d9d6aeba345b52
independent_script: 04-computation/sic_sharp_boundary_beta_shift_independent_audit_thm2801.py
independent_output: 05-knowledge/results/sic_sharp_boundary_beta_shift_independent_audit_thm2801.out
independent_script_sha256: cd7f76f266eafabde173431bcff9424311e7a8728eb9ad06bdd4dead00ef1aa3
independent_output_sha256: b639844cec5976e1ed5ef321133e0677f0ec6db1b42edd6ba24ae9069aafb430
hash_basis: LF-normalized bytes
---

# THM-2801 -- sharp Special Image Conjecture boundary and beta-shift witness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The new mechanism is not merely that a scalar moment vanishes.  A whole
power tower lies in the operator kernel, but multiplying by one coordinate
shifts a coefficient extraction by one place and exposes a nonzero endpoint
jet.  This gives the sharp dimensional status

```text
SIC(1) is true;
SIC(n) is false for every n>=2.                         (1)
```

The dimension-one assertion uses THM-2022.  The counterexample and every
all-\(m\) identity below are proved directly.

## 1. Operator convention and the scalar-output gate

Put

```text
A_n=C[xi_1,...,xi_n,z_1,...,z_n]
```

and define

```text
E_n(P(xi)Q(z))=P(partial_z)Q(z).                       (2)
```

The Special Image Conjecture `SIC(n)` says that

```text
M_n=ker(E_n)
```

is a Mathieu--Zhao subspace: if `a^m` belongs to `M_n` for every
`m>=1`, then `g a^m` belongs to `M_n` for every fixed `g` and all
sufficiently large `m`.

Zhao's kernel theorem (`arXiv:0902.0210`, Theorem 3.1) identifies the same
space with the special differential image:

```text
M_n=ker(E_n)
   =sum_i (xi_i-partial_(z_i))A_n.                       (3)
```

Changing every displayed generator's sign gives the equivalent
`partial_(z_i)-xi_i` convention.

For a monomial,

```text
E_n(xi^alpha z^beta)
 = beta!/(beta-alpha)! z^(beta-alpha)  if alpha<=beta,
 = 0                                  otherwise.        (4)
```

Consequently, if a polynomial is bihomogeneous of bidegree `(d,d)`, its
`E_n` image has degree zero.  On that balanced slice `E_n` equals the
scalar contraction

```text
F_n(xi^alpha z^beta)=alpha! delta_(alpha,beta).          (5)
```

This gate is load-bearing.  A generic Gaussian counterexample controls only
`F_n`; the witnesses below control `E_n` because their bidegrees balance.

## 2. The rank-one two-pair witness

In `A_2`, define

```text
R=xi_1 z_1+xi_2 z_2,      Z=xi_1 z_2,
W=2 xi_2 z_1,             T=xi_1 z_1-xi_2 z_2,          (6a)

F=(R+Z)(R^2 W-(1/2)(2R+Z)T^2).                         (6b)
```

These are coordinates on the rank-one `2 x 2` matrix
`(xi_i z_j)`.  Its determinant relation is

```text
R^2=T^2+2ZW.                                            (7)
```

It gives the useful exact reframe

```text
F=W(R+Z)^3-(1/2)R^2(R+Z)(2R+Z).                        (8)
```

The polynomial `F` has bidegree `(4,4)`, while `Z` has bidegree `(1,1)`.
For every `m>=1`,

```text
E_2(F^m)=0,                                             (9)

E_2(ZF^m)
 =(4m+2)! m!/(2m+1)!! !=0.                             (10)
```

### 2.1 Hopf/Segre reduction

Let `g_1,g_2` be independent standard circular complex Gaussians and set
`xi_i=conj(g_i)`, `z_i=g_i`.  Equation `(5)` says that `E_2` on a
balanced polynomial is its Gaussian expectation.  Write

```text
rho=|g_1|^2+|g_2|^2,
t=(|g_1|^2-|g_2|^2)/rho,
u=Z/rho.                                                (11)
```

The radial variable has the `Gamma(2,1)` law,

```text
E(rho^d)=(d+1)!,                                        (12)
```

and is independent of the uniform Hopf direction.  On that direction,
`t` is uniform on `[-1,1]`, the phase of `u` is uniform, and

```text
W/rho=(1-t^2)/(2u).                                     (13)
```

Thus, with `h=F/rho^4`,

```text
h=(1+u)/(2u) [1-t^2(1+u)^2].                           (14)
```

The apparent division by `u` is only a Laurent-coordinate expression away
from the measure-zero Hopf poles; it came from the polynomial identity
`(7)`.

For

```text
A_m(a)=integral_0^a (1-y^2)^m dy,                       (15)
```

phase constant-term extraction followed by the uniform `t` integral gives

```text
E_Hopf(h^m)
 =2^(-m)[u^m](1+u)^(m-1) A_m(1+u),                    (16)

E_Hopf(u h^m)
 =2^(-m)[u^(m-1)](1+u)^(m-1) A_m(1+u).                (17)
```

Indeed,

```text
(1/2) integral_(-1)^1 [1-t^2(1+u)^2]^m dt
 =A_m(1+u)/(1+u).                                      (18)
```

But

```text
A_m'(1+u)=(-u(2+u))^m,
```

so

```text
A_m(1+u)=A_m(1)+O(u^(m+1)).                            (19)
```

The polynomial `(1+u)^(m-1)` has degree `m-1`.  Hence the coefficient
in `(16)` is zero, while the coefficient in `(17)` is exactly `A_m(1)`.
Finally,

```text
2^(-m)A_m(1)
 =2^(-m) integral_0^1(1-y^2)^m dy
 =m!/(2m+1)!!.                                         (20)
```

Equations `(12)`, `(16)`, and `(20)` prove `(9)--(10)`.

An independent audit expands `(16)` as an alternating finite difference:
the pure coefficient is the `m`th difference of a polynomial of degree
`m-1`; after shifting one coefficient, the surviving remainder is the beta
sum in `(20)`.  This proof does not use the flat-jet argument `(19)`.

### 2.2 A complete observer ladder

The shift works at every finite depth.  For fixed `s>=1`,

```text
E_2(Z^s F^m)=0                                      if m<s,

E_2(Z^s F^m)
 =(4m+s+1)! m!/(2m+1)!! binom(m-1,s-1)             if m>=s. (21)
```

Multiplication by `u^s` changes the extraction in `(16)` from `u^m` to
`u^(m-s)`.  Equation `(19)` then leaves

```text
[u^(m-s)](1+u)^(m-1)=binom(m-1,s-1),
```

and the radial degree becomes `4m+s`, proving `(21)`.

More importantly, no `xi` variable is needed in the MZ multiplier.  The
operator identity

```text
E_2(xi_1^s H)=partial_(z_1)^s E_2(H)                  (22)
```

applied to `H=z_2^sF^m` shows that the coefficient of `z_1^s` in
`E_2(z_2^sF^m)` is

```text
E_2(Z^sF^m)/s!.                                        (23)
```

It is nonzero for every `m>=s`.  At `s=1` there are no other channels:

```text
E_2(z_2F^m)
 =[(4m+2)!m!/(2m+1)!!] z_1.                            (24)
```

To see this, the `z_1` coefficient is `(10)` by `(22)`.  The `z_2`
coefficient is the contraction with `xi_2z_2=rho(1-t)/2`; the angular
integrand `h^m` is even in `t`, its odd part integrates to zero, and its
constant part is a radial multiple of the zero angular moment `(16)`.

Thus `(9)` and `(24)` alone are a direct Mathieu-space counterexample:

```text
F^m in M_2 for every m>=1,
z_2F^m notin M_2 for every m>=1.                        (25)
```

Equation `(21)` says more: every coordinate power `z_2^s` is a persistent
observer after the sharp delay `m=s`.

## 3. A four-term three-pair witness and its all-r family

Pair the derivative and coordinate variables as

```text
(tau,t),             (w,z),             (v,y).
```

For every integer `r>=1`, put

```text
f_r=tau^r(t+z)[w t^r-v y(t+y)^(r-1)].                   (26)
```

It has bidegree `(r+1,r+1)`.  For every `m>=1`,

```text
E_3(f_r^m)=0,                                           (27)

E_3(zf_r^m)
 =(rm+1)!m! t
  +1_(m=1) r! z
  -1_(m=2) 4(r-1)(2r)! y.                              (28)
```

In particular, `(28)` is never zero.  The user's displayed polynomial is
the case `r=3`.  The smallest member,

```text
f_1=tau(t+z)(wt-vy),                                    (29)
```

has only four terms and bidegree `(2,2)`, and satisfies

```text
E_3(zf_1^m)=(m+1)!m! t+1_(m=1)z.                       (30)
```

### Proof

Write `k` for the number of `w` factors and `ell=m-k`.  Expansion gives

```text
f_r^m
 =tau^(rm) sum_(k=0)^m binom(m,k)(-1)^ell w^k v^ell
    t^(rk)y^ell(t+y)^((r-1)ell)(t+z)^m.                 (31)
```

For the scalar output in `(27)`, balance forces `z^k` from `(t+z)^m`
and the zeroth `y` term from `(t+y)^((r-1)ell)`.  The `k`th contribution
is

```text
(rm)!m!(-1)^(m-k)binom(m,k).                            (32)
```

Their sum is zero.

For the coefficient of the output monomial `t` in `(28)`, balance instead
selects `z^(k-1)`, for `k>=1`.  After the derivative factorials cancel,
the sum is

```text
(rm+1)!m! sum_(k=1)^m
  (-1)^(m-k)binom(m,k-1)
 =(rm+1)!m!.                                           (33)
```

The output-`z` and output-`y` channels reduce respectively to alternating
differences of a linear and a quadratic polynomial in `k`.  They vanish
above orders one and two and give exactly the two exceptional terms in
`(28)`.  This proves the full formula, not only its nonvanishing.

## 4. The sharp bidegree-one boundary

No balanced bilinear polynomial can be the base of an SIC counterexample,
in any dimension.  Let

```text
f=xi^T A z.                                               (34)
```

The formal complex-Gaussian/MacMahon identity is

```text
sum_(m>=0) E_n(f^m)t^m/m!
 =1/det(I-tA).                                           (35)
```

It follows either by expanding the exponential and using `(4)`, or by
formal Gaussian elimination.  If `E_n(f^m)=0` for every `m>=1`, then the
left side of `(35)` is one, so

```text
det(I-tA)=1.
```

Thus `A` is nilpotent.  Under the pairing-preserving gauge

```text
z=S z',                  xi=S^(-T)xi',
```

the matrix changes by conjugation.  Put it in strictly upper-triangular
form and assign weights

```text
wt(z_j)=j,               wt(xi_i)=-i.                  (36)
```

Every monomial `xi_i z_j` occurring in `f` then has weight at least one.
Fix a monomial `g=xi^gamma z^delta`.  Every monomial of `g f^m` has weight
at least

```text
m+wt(g).                                                 (37)
```

If such a monomial survives `E_n`, its leftover exponent

```text
e=(z exponent)-(xi exponent)
```

is componentwise nonnegative and has the fixed total

```text
D=|delta|-|gamma|.
```

If `D<0`, survival is impossible.  If `D>=0`, its weight is

```text
sum_i i e_i <=nD.                                       (38)
```

Equations `(37)--(38)` are incompatible once

```text
m>nD-wt(g).
```

Taking the maximum over the finitely many monomials of an arbitrary fixed
polynomial `g` proves

```text
E_n(gf^m)=0 for all sufficiently large m.               (39)
```

Therefore every bidegree-`(1,1)` base in the radical of `ker(E_n)` already
has the Mathieu property.  The minimal balanced base degree of an SIC
counterexample is at least two.  The four-term `f_1` in `(29)` attains
degree two for `n=3`.  For `n=2`, the witness `(6b)` has degree four, leaving
the precise minimization problem:

```text
Can bidegree (2,2) or (3,3) fail in two pairs,
or is four minimal?                                      (40)
```

## 5. Sharp dimension boundary

Equation `(25)` proves `not SIC(2)`.  Adjoining unused pairs of variables
preserves every displayed `E_n` value, so it proves

```text
not SIC(n) for every n>=2.                              (41)
```

Dimension one goes the other way.  THM-2022 proves `GMC(2)`.  Propositions
3.2 and 3.1 of Derksen--van den Essen--Zhao give, in this fixed dimension,

```text
GMC(2) => ker(F_1) is MZ => ker(E_1) is MZ.             (42)
```

Independently, van den Essen--Wright--Zhao prove the one-variable Image
Conjecture when the coefficient ideal is radical; taking
`A=C[xi_1]` and the prime ideal `(xi_1)` gives `SIC(1)` directly.  Thus the
phrase "false in all dimensions" is correct only if it means all previously
open dimensions; literally, dimension one is the sharp exception.

## 6. What this changes for the Jacobian program

The full-SIC route is now closed, but no reverse implication gives a planar
Keller counterexample.  Zhao's exact JC-bearing restrictions are much
smaller:

```text
xi-linear sector:       f=sum_i xi_i H_i(z),
                        H_i cubic homogeneous,
                        multiplier g=z_i;

separable Hessian sector:
                        f=(sum_i xi_i^2)P(z),
                        P quartic homogeneous.           (43)
```

The first witness has bidegree `(4,4)`, `xi`-degree four, and grading drift
zero.  The two sectors in `(43)` have bidegrees `(1,3)` and `(2,4)`,
respectively, and grading drift two.  Simultaneous linear changes of the
`z` variables and contragredient changes of `xi` preserve the two separate
degrees.  Hence this `F` cannot be linearly conjugated into either
JC-bearing sector.

The coordinate-only multiplier `(24)` removes one superficial difference:
the observer now has exactly Zhao's `g=z_i` type.  The remaining obstruction
is entirely in the base polynomial.  The repaired open question is:

> Is the image kernel MZ on the `xi`-linear cubic nilpotent-Jacobian locus,
> or on the separable Hessian-nilpotent locus, even though the full kernel
> is not?

THM-1435's explicit homogeneous Hessian-nilpotent witness remains open.
The present theorem neither constructs that witness nor proves or refutes
`JC(2)`.

## 7. Structural interpretation and transfer boundary

The first failure already lives on the rank-one Segre cone.  After removing
the radial coordinate it is a `CP^1` coefficient problem:

```text
power tower invisible       <=> endpoint Taylor coefficient u^m is absent;
marked descendant visible   <=> multiplication shifts to u^(m-1);
observer ladder             <=> shift by s exposes binom(m-1,s-1). (44)
```

This is the precise reusable lesson.  A collapsed scalar or fibrewise
spectrum can miss an adjacent marked descendant forever.  It supports
retaining owner, endpoint, root, address, or determinant-fibre sidecars in
LRC carrier calculations.  It does not itself supply an LRC map:

```text
SIC target: linear averaged kernel membership;
LRC target: pointwise existential avoidance;
Hopf quotient loses: clock, owner, endpoint chronology, ancestry. (45)
```

Likewise, a KLS/Poincare inequality can bound variance only after a lawful
positive state graph and measure have been built.  The complex endpoint
currents and discontinuous address walls are not log-concave measures, so
the beta/needle picture licenses a four-edge marked test, not a direct KLS
transfer.

## 8. Exact audits and source boundary

Run

```bash
python 04-computation/sic_sharp_boundary_beta_shift_thm2801.py
python -O 04-computation/sic_sharp_boundary_beta_shift_thm2801.py
python 04-computation/sic_sharp_boundary_beta_shift_independent_audit_thm2801.py
python -O 04-computation/sic_sharp_boundary_beta_shift_independent_audit_thm2801.py
```

The main companion:

1. applies `E_2` literally through `m=8`;
2. checks the Taylor/Beta proof through `m=30` and `174` observer-ladder
   identities;
3. applies `E_3` to the user's `r=3` witness through `m=8`;
4. checks the full all-`r` formula for `1<=r<=4`, `1<=m<=4`;
5. checks the finite-difference identities through `m=100`; and
6. includes radial and low-order hidden-channel positive controls.

The independent companion uses a separate sparse-dictionary contraction
and the alternating-difference proof.  It retains the exceptional `z` and
`y` channels in `(28)` and the dimension-one gate.

The user cited `arXiv:2607.23450`; as checked on 2026-07-28, that identifier
is Marco Tomamichel's quantum-information preprint *A strong converse for
stabilizer codes over Pauli channels via the blowing-up lemma*, not a paper
about SIC.  The nearest public structural source located is Long's
`arXiv:2607.18186`, whose three-real-Gaussian counterexample has the
dehomogenized shape visible from `(6b)` on `R=1`; no provenance inference is
made from that similarity.  The operator definition and JC restrictions
used here are from Zhao `arXiv:0902.0210` and
Derksen--van den Essen--Zhao `arXiv:1506.05192`.
