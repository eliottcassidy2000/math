---
id: THM-2245
title: "Degree-fourteen spectral quartic discriminant reduction"
status: >
  PROVED + VERIFIED-EXACT. In the nonsplit terminal quartic exact-square
  branch, every reduced degree-fourteen Keller mate forces an explicit
  depressed spectral quartic H(y) to have a multiple root. The first two
  Faber fluxes complete to
  (235298 T^2+B(y))^2=(10976 y)^2 H(y). The exceptional y=0 branch already
  makes T constant and contradicts the nonsplit deck; if H were squarefree,
  its smooth genus-one double cover would make every rational trajectory
  constant and give the same contradiction. Every survivor is therefore on
  the explicit quartic discriminant and admits a double-root/conic
  normalization. The exact third flux and Keller one-form are retained for
  that singular branch. This is a strict degree-fourteen reduction, not a
  closure of that branch or a proof of planar JC.
source: codex-2026-07-25-degree14-spectral-quartic
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-2217-square-prefix-pole-alternative-and-odd-leading-degree-terminal-wall
  - THM-2240-dc2-grade-response-gauge-is-not-a-continuation-state
script: 04-computation/jc2_degree14_spectral_quartic_discriminant_thm2245.py
output: 05-knowledge/results/jc2_degree14_spectral_quartic_discriminant_thm2245.out
script_sha256: fdad44cf1923cf99d20c45823784f3741ad8c1a29307e639899ae66a5423792d
output_sha256: 401e56f6e069b23e4aec77eea767fd784b1d703b5bd11f3886b4c800cb1b0678
hash_basis: working-tree bytes (LF)
---

# THM-2245 -- degree fourteen lands on a singular spectral quartic

THM-2214 closes nonsplit terminal quartic mates through reduced fibre degree
ten. THM-2230 shows that the reduced degree and the nonresonant Faber
coefficients are independent of the chosen mate. The first open degree is
therefore fourteen. At that degree, the two constant Faber fluxes still
produce a curve of positive genus generically; the new point is that the
curve has a compact quartic model and that every survivor must lie on its
discriminant.

## 1. Inherited nonsplit coordinates

Use the terminal exact-square-prefix coordinates of THM-2214 over
`K=C(x)`. In the genuine quadratic extension obtained by adjoining a square
root of the nonsquare leading coefficient, write

```text
P=w^4+2d w^2+q w+(d^2-s),             q^2=T.          (1)
```

The deck involution sends

```text
(w,q) -> (-w,-q)
```

and fixes `d,s,T`. The constant field of the extension is `C`, and `q` is
nonzero. Normalize the degree-fourteen coefficient. Before removing the
harmless target shear, write

```text
Q=J_0(P)+E_14+alpha E_10+beta E_6+gamma E_2,
alpha,beta,gamma in C.                               (2)
```

The summand `J_0(P)` contributes no Jacobian or flux. THM-2230's full Faber
gauge sets it to zero and makes `alpha,beta,gamma` intrinsic in the fixed
fibre coordinate.

THM-2129's Hamiltonian identity and nonsplit parity give

```text
Phi_Q=0,              Psi_Q=Psi0 in C,
R_Q'=kappa/U,          kappa!=0.                     (3)
```

Here a prime differentiates the coefficient field with `w` fixed. As in
THM-2214, `q=A/U` and `T=A^2/V` for the polynomial coefficients of the
linear remainder and the quartic leading term.

## 2. The degree-fourteen Laurent bank

The new degree-fourteen row of the exact Laurent recurrence is

```text
Phi_14
 =-(7/64)q(-40dTs+T^2+40s^3),

Psi_14
 =(35/64)(2d^2T^2-12dTs^2+T^2s+2s^4),

R_14
 =(35/128)q(-8d^2Ts+dT^2+8ds^3-4Ts^2).             (4)
```

The companion obtains (4) both from the fibre-derivative recurrence and
from an independently coded finite multinomial expansion. It similarly
rechecks the degree `2,6,10` bank consumed from THM-2214.

Put

```text
y=7s-2alpha.                                         (5)
```

On the open set `y!=0`, the first equation in (3) is linear in `dT` and
gives

```text
dT=
 [343T^2-640alpha^3-480alpha^2y
  +2688alpha beta+1344beta y-6272gamma+40y^3]
 /(1960y).                                           (6)
```

The excluded denominator is not a lost branch. If `y=0`, the same first
flux is exactly

```text
343T^2-640alpha^3+2688alpha beta-6272gamma=0.         (7)
```

Thus `T^2` is constant. Since `T` lies in the algebraic function field
whose constant field is `C`, it is constant. Then `q^2=T` makes `q`
constant as well. This contradicts the fact that the nonzero `q` is
anti-invariant under the nonsplit deck. Hence

```text
y!=0                                                  (8)
```

in every degree-fourteen survivor.

## 3. Completed-square spectral equation

Define

```text
B(y)
 =219520y^3-439040alpha^3
   +1843968alpha beta-4302592gamma,                  (9)

H(y)
 =425y^4
  +(840beta-300alpha^2)y^2
  +(-1600alpha^3+6720alpha beta-15680gamma)y
  +13720Psi0+1200alpha^4-6720alpha^2beta
  +7840alpha gamma+7056beta^2.                      (10)
```

Substitute (6) into the second equation of (3). Exact completion of the
resulting quadratic in `T^2` gives

```text
(235298T^2+B(y))^2=(10976y)^2 H(y).                 (11)
```

The constants are structural rather than decimal approximations:

```text
117649=7^6,             10976=32*7^3.
```

Since `y!=0`, put

```text
W=(235298T^2+B(y))/(10976y).                         (12)
```

Then the coefficient trajectory lies on

```text
W^2=H(y).                                            (13)
```

## 4. The quartic must be singular

Suppose `H` were squarefree. It has degree four and leading coefficient
`425`, so the smooth projective completion of (13) has genus one.
The rational functions `(y,W) in C(x)^2` define a rational map from
`P^1` to that completion, and properness extends it across their poles.
Riemann--Hurwitz forbids a nonconstant map:

```text
-2=deg(f)(2*1-2)+Ramification>=0.
```

Thus `y` and `W` are constant. Equation (12) then makes `T^2` constant,
and the constant-field/deck argument following (7) again gives the
contradiction that the nonzero anti-invariant `q` is constant.

Therefore every degree-fourteen survivor satisfies

```text
gcd(H,H')!=1,

equivalently                  Disc_y(H)=0.           (14)
```

This is a genuine restriction on the normalized Faber coefficients and
the constant second flux; it is not a choice of mate or target shear.

## 5. Exact normalization of the remaining discriminant

The quartic (10) has no cubic term. Let `e` be any multiple root in the
remaining branch. Polynomial division and the missing cubic coefficient
force the exact factorization

```text
H(y)=425(y-e)^2((y+e)^2+D)                           (15)
```

for some `D in C`. If `h_2,h_1,h_0` denote the coefficients of
`y^2,y,1` in (10), then (15) is equivalently

```text
h_2=425D-850e^2,
h_1=-850De,
h_0=425e^2(D+e^2).                                  (16)
```

After removing the double factor, the first normalization is the conic

```text
v^2=(y+e)^2+D.                                       (17)
```

The original degree-fourteen curve is still the additional double cover
obtained by recovering `T` from (11). Thus singularity of `H` does not by
itself rationalize the full Keller trajectory. Equations (15)--(17)
identify the next exact task: compute the branch divisor of that remaining
`T`-cover and retain the Keller one-form below.

## 6. The third flux is the continuation sidecar

Eliminating `dT` from the third observable in (3) gives

```text
R_Q/q
 =-T[-343T^2+640alpha^3-320alpha^2y
      -2688alpha beta+896beta y+6272gamma+160y^3]
   /(8960y)
 =:F(y,T).                                           (18)
```

Using `q^2=T`, `q=A/U`, and `R_Q'=kappa/U`, equation (18) becomes the
exact Keller one-form

```text
A(2T F'+F T')=2kappa T.                              (19)
```

The spectral discriminant without (19) is only a response state, not a
continuation state. This is the planar analogue of the kernel-memory
warning made explicit for the Ore problem in THM-2240: the next proof must
carry both the normalized singular curve and its differential.

## 7. Scope and reproduction

The nonsplit deck, exact polynomial square prefix, terminal Faber normal
form, and reduced degree fourteen are load-bearing. On a split deck an
anti-invariant flux may be a nonzero constant, while without the exact
square prefix the coordinates (1) and the bank (4) are unavailable.

The theorem proves only the singular-discriminant reduction (14). It does
not show that the factorizations (15) are impossible, close higher reduced
degrees, or prove the planar Jacobian conjecture.

Run

```bash
python3 04-computation/jc2_degree14_spectral_quartic_discriminant_thm2245.py
python3 -O 04-computation/jc2_degree14_spectral_quartic_discriminant_thm2245.py
```

The exact companion checks the degree-fourteen bank by two independent
coefficient engines, the first-flux elimination and `y=0` branch, the
completed-square identity, its quadratic discriminant factorization, and
the third observable. The genus and constant-field arguments are the
mathematical proof of the global implication. Both execution modes produce
the stored transcript byte for byte. QED.
