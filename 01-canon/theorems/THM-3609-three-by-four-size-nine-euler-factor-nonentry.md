---
id: THM-3609
title: "Exponent-two three-by-four size-nine Euler-factor nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every
  squarefree Sigma over C of degree at least two, an exponent-two polynomial
  Darboux pair with three retained homogeneous pieces against four cannot have
  additive sumset size nine.  THM-3606 leaves four oriented words and one
  unexposed scalar scheme on each.  Their singleton equations force two common
  UFD bases; the unique scalar-arm address forces parity two and Sigma to divide
  the negative base.  The scalar equation then has the common nonunit Euler
  factor k h'+2h k' and cannot equal one.  Conditional only on the earlier
  3x4 atlas, the residual sector is 144 schemes on 27 words of sumset size at
  most eight.  No Darboux pair, planar Jacobian counterexample, or proof of
  JC(2) is claimed.
source: kps-s188 + agent Hegel / THM-3606 four-word continuation, 2026-08-21
audit: >
  Author derivation, exact implementer audit, and independent hostile audit by
  agent Arendt.  The companion hash-pins the
  THM-3606 atlas, derives both two-parameter collision cones, reconstructs the
  singleton exponent networks, proves uniqueness of the arm address, and
  checks the polynomial Euler factorizations symbolically and on hostile
  admissible parameters.  The independent audit reimplemented the finite
  size-nine census, rederived both cones and chamber inequalities, checked the
  zero-weight boundary and UFD gcds, proved the unique-arm order calculation,
  and independently expanded both scalar factorizations.  Ordinary and
  optimized companion runs reproduced the stored output byte-for-byte.
depends_on:
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
related:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_three_by_four_size_nine_euler_factor_nonentry_thm3609.py
output: 05-knowledge/results/jc2_three_by_four_size_nine_euler_factor_nonentry_thm3609.out
script_sha256: 7fdab86fe3b98275ad5a58f7ef52556d74497d3207e2782b281903d054c63044
output_sha256: dc7e915692172325a01e7dea265aa29aab60f63b9826a782722b6f1fc205d6e1
semantic_sha256: 904bf1347832e0f193bf2325092f03af59d1ade2f6a7684fc7227d532d2d7a4d
hash_basis: LF-normalized bytes
---

# THM-3609 -- exponent-two three-by-four size-nine Euler-factor nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  This closes
the top additive layer of the exponent-two `3 x 4` sector.  The proof does not
delete the scalar cut.  It converts the surrounding singleton rows into a
global factor of the scalar equation.

All rings are over `C`.  Let `Sigma in C[b]` be squarefree with
`deg Sigma>=2`, and put

```text
A_Sigma=C[b,c,e]/(c^2e-Sigma(b)),        wt(b,c,e)=(0,1,-2).       (1)
```

Write a putative normalized Darboux pair as

```text
P=sum_(i=0)^2 c^(p_i) f_i(b),
Q=sum_(j=0)^3 c^(q_j) g_j(b),              {P,Q}=1,                (2)
```

where all displayed coefficient polynomials are nonzero and scalar
weight-zero constants have already been removed.  For homogeneous pieces use

```text
W_(r,s)(F,G)=sF'G-rFG'.                                  (3)
```

Then the bracket row of total weight `r+s+1` is the sum of the corresponding
`W` terms.

## 1. Exact four-word reduction

THM-3606 proves that, after the scalar-arm, singleton-sign, and
rectangle-exposure gates, exactly four sumset-nine words remain:

```text
W072, W073, W141, W149.                                  (4)
```

Each has one unexposed scalar scheme.  The first reversal pair has collision
fibres

```text
02=10,                    11=20 (scalar),              03=21.      (5)
```

Solving these three support equalities gives the full positive cone

```text
(X,Y,U,V,W)=(a+t,a,a,t,2a),                a,t>0.        (6)
```

The chamber inequality is `t>a` for `W072` and `t<a` for `W073`.  The unique
residual anchor is `11,+-`.  Thus the absolute weights are

```text
p=(1-a-t, 1, a+1),
q=(-a-2, -2, t-2, 2a+t-2).                              (7)
```

The reflected pair has collision fibres

```text
02=20,                    03=12 (scalar),              13=21.      (8)
```

Their cone and weights are

```text
(X,Y,U,V,W)=(a,a+t,2a,t,a),
p=(-a-2,-2,a+t-2),
q=(1-2a-t,1-t,1,a+1),                                  (9)
```

with `t<a` for `W141`, `t>a` for `W149`, and unique anchor `12,-+`.
Simultaneous reversal pairs `W072` with `W149` and `W073` with `W141`, but
we retain both orientations because regularity is not reversal-invariant.

## 2. The singleton UFD bases

The singleton gate from THM-3592 has the following elementary consequence.
If nonzero polynomials `F,G` of same-sign nonzero weights `r,s` satisfy
`W_(r,s)(F,G)=0`, then unique factorization gives a common polynomial base,
with exponents `|r|/gcd(|r|,|s|)` and
`|s|/gcd(|r|,|s|)`.  Applying this along a connected singleton network and
comparing valuations combines all vertices into one base whose denominator is
the gcd of all incident absolute weights.

For `(6)--(7)`, the negative singleton component and positive singleton
component therefore have nonzero bases `h,k` and nonzero constants such that

```text
f_0=A h^((a+t-1)/d),       g_0=B h^((a+2)/d),
g_1=C h^(2/d),

f_1=L k,                   f_2=M k^(a+1),
g_2=N k^(t-2),             g_3=T k^(2a+t-2),            (10)
```

where

```text
d=gcd(a+t-1,a+2,2) in {1,2}.                            (11)
```

Existence of the scheme already makes all exponents in the positive component
strictly positive.  The second scalar address `20` has weights
`(a+1,-a-2)`, so it is never `(-2,1)` or `(1,-2)`.  Hence `11`, with
coefficients `(f_1,g_1)`, must supply the scalar at every root of `Sigma`.
The negative coefficient `g_1=C h^(2/d)` must have a simple zero at every arm,
while `f_1=Lk` must be a unit there.  The case `d=1` has even vanishing order
and is impossible.  Therefore

```text
d=2,           a even,           t odd,
Sigma divides h,                 gcd(Sigma,k)=1.         (12)
```

For `(9)`, the four negative singleton addresses and two positive singleton
addresses similarly give

```text
f_0=A h^((a+2)/d),          f_1=L h^(2/d),
g_0=B h^((2a+t-1)/d),       g_1=C h^((t-1)/d),

f_2=M k^(a+t-2),            g_2=N k,
g_3=T k^(a+1),                                         (13)
```

with

```text
d=gcd(a+2,2,2a+t-1,t-1).                               (14)
```

The other scalar address `03` is ineligible.  The unique address `12` forces
the same conclusion `(12)`.

This is the point at which possible arm-by-arm alternation disappears: in
each word only one address has eligible absolute weights.

## 3. The unavoidable Euler factor

Define

```text
E(h,k)=k h'+2h k'.                                      (15)
```

Put `eta=(a+2)/2`.  In the first pair, the scalar row `11=20`, equations
`(7),(10),(12)`, and direct use of `(3)` give

```text
W_(1,-2)(Lk,Ch)
 +W_(a+1,-a-2)(M k^(a+1),B h^eta)

=-E(h,k)[LC+(a+1)eta MB h^(eta-1)k^a].                 (16)
```

But this row must equal `1`.

In the reflected pair, the scalar row is `03=12`.  Equations
`(9),(12),(13)` give

```text
W_(-a-2,a+1)(A h^eta,T k^(a+1))
 +W_(-2,1)(Lh,Nk)

=E(h,k)[AT(a+1)eta h^(eta-1)k^a+LN],                  (17)
```

which must again equal `1`.

Neither equality is possible.  Since `Sigma|h` and `deg Sigma>=2`, the
polynomial `h` has degree at least two.  The leading term of `(15)` is

```text
(deg h+2 deg k) lc(h)lc(k) b^(deg h+deg k-1),           (18)
```

so

```text
deg E=deg h+deg k-1>=1.                                 (19)
```

Thus `E` is a nonunit of `C[b]` and cannot divide the unit `1`.  This
contradicts `(16)` and `(17)` and eliminates all four words.

Notice what was not used: the scalar fibre is a one-fibre projection cut, but
we never delete it.  The lawful chain is

```text
singleton zero rows -> two UFD bases -> unique arm address
 -> d=2 and Sigma|h -> common nonunit Euler factor -> contradiction.        (20)
```

## 4. Consequence and exact boundary

Combining the theorem with THM-3606 reduces the still-open exponent-two
`3 x 4` sector from 148 unexposed schemes on 31 words to

```text
144 schemes on 27 oriented words,             all with |A+B|<=8.  (21)
```

The companion derives `(6),(9)` from collision equations, verifies the four
chamber orders and unique anchors against the hash-pinned THM-3606 ledger,
checks all singleton weight/exponent identities, and expands every collision
factorization in `(16)--(17)` over an exact symbolic polynomial ring.  It also
tests both reversal pairs on multiple admissible parameter values.  Reproduce
with

```bash
python3 04-computation/jc2_three_by_four_size_nine_euler_factor_nonentry_thm3609.py
python3 -O 04-computation/jc2_three_by_four_size_nine_euler_factor_nonentry_thm3609.py
```

The remaining 27 words have one or more of three new difficulties: a scalar
triple, more than one eligible scalar address, or singleton components that do
not collapse all coefficients to two bases.  Their resolution requires a new
coefficient or all-arm compatibility gate; the present Euler argument is not
asserted beyond the four size-nine words.

This theorem proves a nonentry statement only inside the exponent-two `3 x 4`
support cell.  It constructs no Darboux pair or polynomial map and proves no
case of the planar Jacobian conjecture false.

**QED.**
