---
id: THM-3835
title: "A polynomial marked-root ratio cannot support the nonlinear cubic plane atlas"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every dominant plane atlas of the THM-3811 surface, the marked
  cubic-root ratio z=h/k is genuinely rational and nonpolynomial.  If z were
  polynomial, the SL2 determinant would make k scalar; the genus-three
  square sidecar then makes h scalar, contradicting dominance.  The result
  is all-degree and does not exclude a nonconstant denominator k.
source: jc_quartic_c3_construct / triangular root-ratio polynomialization lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The companion verifies the degree-seven
  sidecar, the determinant factorization, the constant-row lift, and an exact
  Groebner certificate [1] excluding H(T,c)=(T-a)S(T)^2 for every c!=0.
  It retains c=0 as a hostile positive control.  Normal and optimized runs
  byte-match the frozen transcript; independent hostile audit remains.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
related:
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
script: 04-computation/jc2_polynomial_marked_root_ratio_nonentry_thm3835.py
output: 05-knowledge/results/jc2_polynomial_marked_root_ratio_nonentry_thm3835.out
script_sha256: 130e542df0a86ea57008d0c0fca909852449d2872e0d14cae2441f1778be94d9
output_sha256: d3aaf8d4e8f87af6dccb3bc69d92acc7ae499ad48cf9a9b9fc7792f64f106a90
semantic_sha256: 02fdbd5deb402cfb98842cf84318b6088f0a7d100d66cbf10004d8a8362f3b3b
hash_basis: raw LF bytes
---

# THM-3835 -- the marked-root ratio must keep its denominator

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**  Work over an algebraically closed field `K` of characteristic zero.
Let

```text
psi:A2_(x,y) -> U                                                (1)
```

be a dominant morphism to the THM-3811 nonlinear cubic surface.  Write
`R=K[x,y]` and denote the pulled-back intrinsic generators by

```text
h,k,m,C,D in R,                 Ck-mh=1.                         (2)
```

Every dominant plane atlas, and hence every prospective planar Keller
factorization through `U`, satisfies these hypotheses.  In the function
field put

```text
z=h/k.                                                            (3)
```

Then

```text
z is not in R.                                                    (4)
```

Equivalently, the denominator `k` in the THM-3832 triangular root-ratio
chart is structural: it cannot be removed by a polynomial source
parametrization.  No bounded-degree or support assumption is present.

## 1. A polynomial ratio collapses the unimodular row

Suppose, toward a contradiction, that `z` belongs to `R`.  Then `h=zk`, so
the determinant law `(2)` becomes

```text
1=Ck-mh=k(C-mz).                                                  (5)
```

Thus `k` is a unit of `R`, hence

```text
k=c in K*.                                                        (6)
```

This is stronger than merely saying that `h/k` was written in lowest terms:
the actual `SL_2` completion makes a polynomial quotient impossible unless
its denominator is scalar.

## 2. The constant-row sidecar has at least two odd roots

THM-3822 proves that every actual pullback `(1)` has a polynomial `w` with

```text
w^2=H(h,k),                                                       (7)
```

where

```text
H(T,c)=
 84T^7 +(36c^2+196c)T^6 +(84c^3+36c^2)T^5
       +(49c^4+112c^3)T^4 -12c^5T^3
       +(-14c^6+12c^5)T^2 +c^8.                                 (8)
```

The following exact one-variable lemma is the key new gate:

```text
For every c in K*, H(T,c) is not (T-a)S(T)^2
for any a in K and S in K[T].                                    (9)
```

To prove it uniformly, write

```text
S=lT^3+uT^2+vT+w_0.                                              (10)
```

Equate all eight coefficients in

```text
(T-a)S^2-H(T,c)=0                                                (11)
```

and adjoin `c c_inv-1=0`.  Over `Q`, the exact grevlex Groebner basis of
these coefficient equations in

```text
(l,u,v,w_0,a,c_inv,c)                                            (12)
```

is `[1]`.  This proves `(9)` after every characteristic-zero scalar
extension.  The invertibility condition is load-bearing: at `c=0`,

```text
H(T,0)=84T^7=T(sqrt(84)T^3)^2.                                  (13)
```

The companion freezes both the universal certificate and this hostile
boundary.

Because `(8)` has odd degree seven, it has at least one root of odd
multiplicity.  If it had only one *distinct* such root `a`, all remaining
multiplicities would be even and algebraic closedness would give exactly the
forbidden form `(9)`.  Consequently, for every `c!=0`, there are two distinct
roots

```text
a_1!=a_2                                                         (14)
```

whose multiplicities in `H(T,c)` are odd.

## 3. UFD parity makes h constant

Factor `H(T,c)` over `K` and substitute `T=h`.  For distinct roots `a_i`,
the polynomials `h-a_i` are pairwise comaximal in the UFD `R`, since their
difference is a nonzero scalar.  Equation `(7)` therefore says that every
irreducible valuation occurring in `h-a_i`, for an odd-multiplicity root
`a_i`, has even order.  Hence

```text
h-a_i=gamma_i u_i^2,          gamma_i in K*, u_i in R.           (15)
```

Absorb `gamma_i` into `u_i^2`, using algebraic closedness.  Applying `(15)`
to the two roots in `(14)` gives

```text
u_1^2-u_2^2=a_2-a_1 in K*.                                      (16)
```

Thus `(u_1-u_2)(u_1+u_2)` is a unit of `R`.  Both factors are units, hence
scalars; characteristic zero then makes `u_1,u_2`, and finally `h`, scalar.
Together with `(6)`, the whole first row `(h,k)` is constant.

## 4. Constant first row contradicts dominance

The first intrinsic lift law is

```text
D(7h^2+3k^2)=1+2Ck,                    A=hD.                     (17)
```

If `7h^2+3k^2=0`, equation `(17)` forces `C=-1/(2k)` to be scalar.  If the
coefficient is nonzero, `(17)` expresses `D`, and hence `A`, as an affine
polynomial in `C`.  In either case `A,C` are algebraically dependent.  This
contradicts dominance of `(1)`, because the structural map
`U -> A2_(A,C)` is dominant and `psi` was assumed dominant.

Therefore `(4)` holds.  The theorem closes the polynomial-ratio boundary but
deliberately leaves the genuine construction lane

```text
z=h/k in K(x,y) but not in R,      k nonconstant.                 (18)
```

No plane atlas and no Jacobian counterexample is constructed.  **QED,
pending independent hostile audit.**
