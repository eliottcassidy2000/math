---
id: THM-3047
title: "Formal-corner width product-Gamma moment and strict Hankel positivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The universal positive-real width flag isolated by THM-3040 is exactly the
  moment sequence of a scaled product of independent Gamma variables.  It is
  therefore a strict Stieltjes moment sequence: every generalized minor of
  its infinite Hankel matrix is positive.  Its adjacent curvature and every
  higher finite log difference have exact alternating signs.  This is a
  coefficientwise formal-corner theorem, not positivity of a physical width,
  raw chart, wall-stripped core, or a corner with moving lower offsets.
source: kind-pasteur-2026-08-01-product-gamma-width
audit: >
  An independent immutable-file hostile audit ACCEPTED the all-k character
  signs, exact product-Gamma/Mellin representation, continuous full-support
  product law, Andreief generalized-Hankel strict total positivity, adjacent
  curvature, all alternating finite-log signs, and every stated scope
  boundary.  It replayed normal, optimized, and stored output byte-for-byte,
  matched both LF hashes, and passed the documentation checker.
depends_on:
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-3030-eighth-resultant-jet-and-corner-constant-closed-forms
script: 04-computation/gmc_formal_corner_product_gamma_hankel_thm3047.py
output: 05-knowledge/results/gmc_formal_corner_product_gamma_hankel_thm3047.out
script_sha256: 1310c2e95bb9539af09c7a5fc97589254784f794456db76986f859cdf5f69e15
output_sha256: de53c0b8052ecb47f0dea5f4cecbe3e8dfe9ba44e81c19219d784dd06ce52125
hash_basis: LF-normalized bytes
---

# THM-3047 -- the formal-corner width flag is a product-Gamma moment sequence

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3040 turns the factorial-resultant corner into a one-dimensional width
quotient.  The quotient is stronger than a convenient formula: on the
positive real `t`-axis it is the moment recurrence of one explicit positive
random variable.  Thus the entire universal width direction, at every slot
number and every order, lies inside the strict Stieltjes moment cone.

This positivity does not cross the formal-corner boundary by itself.  The
point is instead to identify precisely which part of the GMC resultant already
has a positive measure model, and therefore cannot be responsible for a
Hankel-sign obstruction.

## 1. Universal width flag

Fix a slot number `k>=2` and put

```text
A=A_k=k!(H_k-1),
B=B_k=k!(k+1-2H_k),
I=I_k=A+B=k!(k-H_k).                                  (1)
```

Here `H_k=sum_(j=1)^k 1/j`.  THM-3040 proves that the normalized formal-corner
width quotient is

```text
R_(M+1)^C(t)/R_M^C(t)
 =(1+Mt)^A(1+(M+1)t)^B.                              (2)
```

The universal flag singled out by this recurrence, normalized by `F_0=1`, is

```text
F_M^(k)(t)
 = product_(s=1)^(M-1)(1+st)^I (1+Mt)^B,             (3)
```

with the empty-product convention.  In the standard normalization of
THM-3040, `(3)` is the width-dependent factor left after division by the fixed
lower resultant.  For the four-slot corner,

```text
(A,B,I)=(26,20,46).                                  (4)
```

The signs in `(1)` are intrinsic.  Since `H_k>1`, `A>0`.  Also
`k+1-2H_k=0` at `k=2`, while its increment from `k` to `k+1` is
`1-2/(k+1)>0`; hence `B>=0`, with equality only at `k=2`.

## 2. Exact product-Gamma representation

Assume `t>0` and write `a=1/t`.  Let

```text
X_1,...,X_A  be independent Gamma(a,1) variables,
Y_1,...,Y_B  be independent Gamma(a+1,1) variables,   (5)
```

and let all `X` and `Y` variables be mutually independent.  If `B=0`, the
second family is empty.  Define

```text
Z=t^I product_(i=1)^A X_i product_(j=1)^B Y_j.         (6)
```

Then for every integer `M>=0`,

```text
F_M^(k)(t)=E[Z^M].                                    (7)
```

Indeed, using the rising factorial `(x)_M`, independence gives

```text
E[Z^M]=t^(IM)(a)_M^A(a+1)_M^B.                        (8)
```

The two elementary cancellations are

```text
t^M(a)_M     =product_(s=1)^(M-1)(1+st),
t^M(a+1)_M   =product_(s=1)^M(1+st).                  (9)
```

Raising `(9)` to the powers `A,B` proves `(7)`.  More generally, for
`Re(z)>-a`, the Mellin transform is

```text
E[Z^z]
 =t^(Iz) [Gamma(a+z)/Gamma(a)]^A
          [Gamma(a+1+z)/Gamma(a+1)]^B.                (10)
```

Thus THM-3040's quotient is exactly the Mellin shift recurrence of `(10)`.

## 3. Strict Stieltjes and generalized Hankel positivity

The law `mu_(k,t)` of `Z` is supported on `(0,infinity)`.  It is non-atomic
and has positive mass in every nonempty open subinterval: take logarithms in
`(6)` and convolve the continuous full-support log-Gamma densities.  Equation
`(7)` therefore makes `(F_M)_(M>=0)` a Stieltjes moment sequence with infinite
support.

In fact its whole infinite Hankel matrix is **strictly totally positive**.
For every `n>=1` and strictly increasing nonnegative integer sequences

```text
r_1<...<r_n,             c_1<...<c_n,                 (11)
```

one has

```text
det [F_(r_i+c_j)^(k)(t)]_(i,j=1)^n >0.                (12)
```

Andreief's identity rewrites the determinant as

```text
1/n! integral det[x_j^(r_i)] det[x_j^(c_i)]
             product_j dmu_(k,t)(x_j).                (13)
```

On the chamber `0<x_1<...<x_n`, both generalized Vandermonde determinants in
`(13)` are strictly positive.  Their product is symmetric under permutation
of the `x_j`, and the complement of the distinct-point chambers has measure
zero.  Since `mu_(k,t)` charges every open interval, `(13)` is strictly
positive.  This proves `(12)`, not merely positivity of the contiguous
principal minors.

## 4. Exact curvature and the complete alternating log hierarchy

Put

```text
Q_M=F_(M+1)/F_M=(1+Mt)^A(1+(M+1)t)^B.                 (14)
```

For every `M>=1`,

```text
F_(M-1)F_(M+1)/F_M^2
 =((1+Mt)/(1+(M-1)t))^A
  ((1+(M+1)t)/(1+Mt))^B
 >1.                                                   (15)
```

So adjacent log-convexity is strict with an exact multiplicative curvature.
There is also an all-order statement.  For every `r>=1` and `M>=0`,

```text
(-1)^(r-1) Delta^(r+1) log F_M >0.                    (16)
```

To prove it, let `f(x)=log(1+tx)`.  From `(14)`,

```text
Delta^(r+1)log F_M=A Delta^r f(M)+B Delta^r f(M+1).   (17)
```

Repeated use of the fundamental theorem of calculus gives

```text
Delta^r f(x)
 =integral_[0,1]^r f^(r)(x+u_1+...+u_r)du,

f^(r)(x)=(-1)^(r-1)(r-1)!t^r/(1+tx)^r.               (18)
```

Because `A>0`, `B>=0`, and `t>0`, equations `(17)--(18)` prove the strict
sign in `(16)`.

## 5. Sharp boundaries

### `t=0`

At `t=0`, `F_M=1` for every `M`.  The representing measure degenerates to
the point mass at `1`, and every Hankel minor of order at least two vanishes.
Thus strictness in `(12)` genuinely needs `t>0`.

### Negative `t`

No universal nonnegative moment statement survives for `t<0`.  Already for
`k=2`, where `(A,B,I)=(1,0,1)`, one has

```text
F_2^(2)(-2)=1-2=-1.                                    (19)
```

### Moving lower offsets

For fixed lower offsets, THM-3040 separates one `M`-independent nonzero lower
resultant from `(3)`.  If a lower offset moves with `M`, the full corner has
an additional lower-resultant transport factor.  Nothing here controls its
Hankel minors or finite log differences.  The product-Gamma theorem applies
to the universal factor `(3)`, not automatically to that full moving-lower
corner.

### Physical and raw objects

The substitutions `2^M=3^M=0` are coefficientwise formal corner operations.
Consequently `(7)--(16)` do not assert positivity for an integer physical
width, the raw selected Macaulay minor, the wall-stripped norm core, or a
global jet polynomial `P_j`.  Those objects carry factors and transport data
which `(3)` deliberately forgets.

## 6. Exact companion

The dependency-free companion uses rational arithmetic and explicit
determinants.  It checks:

- the integral character ledger `(A,B,I)` for every `k=2..12`;
- `187` all-`k` formal factor-profile identities and `272` materialized
  rational Gamma-moment identities;
- `240` exact adjacent-curvature cells and `1056` alternating finite-log
  difference cells;
- `13,720` positive generalized Hankel minors of sizes `1..3`, with
  independently chosen row and column index sets;
- the rank-one `t=0` boundary and the negative `k=2,t=-2` hostile.

Reproduce with

```text
python 04-computation/gmc_formal_corner_product_gamma_hankel_thm3047.py
python -O 04-computation/gmc_formal_corner_product_gamma_hankel_thm3047.py
```

Both runs equal the stored seven-line transcript byte-for-byte.

The theorem identifies a large strictly positive universal width sector, but
it does not prove the Gaussian Moment Conjecture or any physical-width
nonvanishing statement.

**QED.**
