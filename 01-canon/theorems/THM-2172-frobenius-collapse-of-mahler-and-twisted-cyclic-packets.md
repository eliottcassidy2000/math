---
id: THM-2172
title: "Frobenius collapse of Mahler recurrences and twisted cyclic packets"
status: >
  PROVED. Over F_p[[w]], every linear p-power Mahler recurrence is already
  an additive algebraic equation because G(w^(p^j))=G(w)^(p^j). The
  fibbinary and Moser--de Bruijn series from THM-2005 therefore satisfy the
  unique cubic equations G_F=1+wG_F^3 and (1+w)G_M^3=1 over F_2[[w]].
  In the twisted cyclic quotient F_p[y]/(y^n-a), Frobenius is an explicit
  weighted permutation of coefficient positions when p does not divide n
  and a is nonzero; this gives an exact support-packet law for Shunia's
  quotient polynomials at times t and pt. The hypotheses are sharp for
  general support preservation. Reduction still destroys coefficient
  magnitude, order, and natural division, so the packet law does not replace
  THM-2159/2165's Archimedean spectral gap and arithmetic quantum.
source: codex-2026-07-24-cross-frontier-face-transfer
depends_on:
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-2159-shunia-finite-root-spectral-extraction
related:
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2159-shunia-finite-root-spectral-extraction
  - THM-2165-shunia-linear-exponent-spectral-extraction
---

# THM-2172 -- Frobenius collapse of Mahler and twisted cyclic packets

This theorem records two exact transfers of the whole-layer Frobenius
mechanism. It also identifies the precise point at which neither transfer
can carry an ordered or Archimedean conclusion.

## 1. Prime-power Mahler recurrences become algebraic

Let `p` be prime and let

```text
G(w) in F_p[[w]].
```

Coefficientwise Frobenius gives, for every `j>=0`,

```text
G(w^(p^j))=G(w)^(p^j).                                (1)
```

Consequently, any linear `p`-power Mahler equation

```text
G(w)=sum_(j=0)^r a_j(w)G(w^(p^j)),
a_j(w) in F_p[w],                                    (2)
```

is exactly the additive polynomial equation

```text
G(w)=sum_(j=0)^r a_j(w)G(w)^(p^j).                   (3)
```

Unless the right side of (2) is tautologically just `G(w)`, equation (3) is
a nonzero polynomial equation over `F_p(w)`. Thus `G` is algebraic over
`F_p(w)`.

The proof is the formal freshman's dream. If

```text
G(w)=sum_(m>=0)c_m w^m,
```

then in characteristic `p`

```text
G(w)^(p^j)=sum_(m>=0)c_m^(p^j)w^(m p^j)
          =sum_(m>=0)c_m w^(m p^j)
          =G(w^(p^j)).
```

No convergence or coefficientwise isolation is involved: Frobenius
preserves the complete formal series.

## 2. Exact fibbinary and Moser cubics

Use the indicator generating series from THM-2005, including their zero
terms:

```text
G_F(w)=sum_(m in F)w^m,
G_M(w)=sum_(m in M)w^m,                               (4)
```

where `F` is the set of nonnegative binary integers with no adjacent ones
and `M` is the set of nonnegative integers whose base-four digits lie in
`{0,1}`. THM-2005 proves over `Z[[w]]` that

```text
G_F(w)=G_F(w^2)+wG_F(w^4),
G_M(w)=(1+w)G_M(w^4).                                (5)
```

Their coefficients are zero or one, so coefficientwise reduction modulo
two retains their exact supports. By (1), their reductions in `F_2[[w]]`
satisfy

```text
G_F=G_F^2+wG_F^4,
G_M=(1+w)G_M^4.                                      (6)
```

Both series have constant term one and are therefore units. Cancelling one
factor gives the sharper cubic equations

```text
G_F=1+wG_F^3,
(1+w)G_M^3=1.                                       (7)
```

These equations determine the two series uniquely in `F_2[[w]]`.

For the first equation, every solution has constant term one. If `X,Y` are
solutions and `X-Y` is divisible by `w^k`, then

```text
X-Y=w(X^3-Y^3)
```

is divisible by `w^(k+1)`. Starting at `k=1` proves `X=Y`.

For the second equation, every solution again has constant term one. If
`X,Y` are solutions, the unit `H=X/Y` obeys `H^3=1`. But

```text
H^3-1=(H-1)(H^2+H+1),
```

and the second factor has constant term one in characteristic two, so it is
a unit. Hence `H=1` and `X=Y`.

Thus (7) is an exact algebraic encoding of the two automatic supports, not
merely an asymptotic or sampled shadow.

## 3. Frobenius in a twisted cyclic quotient

Let

```text
R_(p,n,a)=F_p[y]/(y^n-a),
```

and write every element uniquely as

```text
f=sum_(r=0)^(n-1)c_r y^r.                            (8)
```

The Frobenius endomorphism has the exact coordinate formula

```text
f^p
 =sum_(r=0)^(n-1)c_r
    a^(floor(pr/n)) y^(pr mod n).                    (9)
```

Indeed, freshman's dream first gives

```text
f^p=sum_r c_r^p y^(pr)=sum_r c_r y^(pr),
```

and repeated use of `y^n=a` gives (9).

Assume now that

```text
p does not divide n,             a!=0 in F_p.        (10)
```

Then `r -> pr mod n` is a permutation and every weight in (9) is nonzero.
Therefore Frobenius is a weighted permutation of the coefficient positions:

```text
supp_coeff(f^p)
 ={pr mod n:r in supp_coeff(f)}.                     (11)
```

In particular it preserves the number of nonzero coefficient positions and
cannot create cancellation between distinct positions. Equivalently,
`R_(p,n,a)` is finite etale over `F_p` and Frobenius is an automorphism; (9)
gives the stronger labelled coordinate action needed here.

Both hypotheses in (10) are necessary for general support preservation.

- If `a=0`, then in `F_3[y]/(y^2)` the nonzero element `y` has `y^3=0`,
  even though `3` does not divide `2`.
- If `p|n`, then in `F_p[y]/(y^p-1)` the nonzero element `y-1` has
  `(y-1)^p=y^p-1=0`, even though `a=1`.

Thus the good-prime boundary is sharp, not a proof artifact.

## 4. Exact Shunia packet law

For an integer `a` and `t>=0`, let

```text
F_t(y)=(1+y)^t mod (y^n-a)
      =sum_(r=0)^(n-1)c_(t,r)y^r,                    (12)
```

as in THM-2159 and THM-2165. Reduce (12) modulo a prime `p`. In the quotient
ring of Section 3,

```text
F_(pt)=F_t^p.                                        (13)
```

Hence, whenever `p` does not divide `na`,

```text
supp{c_(pt,r) mod p}
 =p supp{c_(t,r) mod p} mod n,                       (14)
```

with the individual nonzero entries rescaled by the explicit powers of `a`
in (9). This is a genuine whole-packet Frobenius transfer. It propagates the
complete coefficient support without isolating a preferred Fourier mode or
coefficient.

It does not perform Shunia's integer-root extraction. That extraction needs
the ordered statement that one ratio lies strictly between two consecutive
integers. Reduction modulo `p` has no such order and cannot even determine
natural division in general: the positive pairs

```text
(U,V)=(1,1),             (U',V')=(p+1,1)             (15)
```

have identical residue vectors modulo `p`, while

```text
floor(U/V)=1,            floor(U'/V')=p+1.            (16)
```

Therefore THM-2159/2165's dominant-mode contraction and arithmetic distance
to the adjacent integers remain load-bearing. Equation (14) is an exact
finite-characteristic sidecar, not a replacement for their Archimedean
argument.

## 5. Transfer ledger and loss boundary

The maps proved above are:

| Source | Target map | Preserved predicate | Destroyed information |
|---|---|---|---|
| `p`-power Mahler recursion | `G(w^(p^j)) -> G(w)^(p^j)` | the complete formal equation; for zero-one coefficients, exact support | integer multiplicity outside the zero-one locus; order, norm, and analytic estimates |
| twisted cyclic quotient | `r -> pr mod n`, weighted by `a^(floor(pr/n))` | labelled nonzero coefficient support at primes satisfying (10) | coefficient magnitude and natural-number quotient; at bad primes even support |
| Shunia time packet | `F_t -> F_(pt)=F_t^p` | the complete residue-support packet | the leading-coordinate order and distance to an integer boundary |

The Archimedean loss is literal. The integer series

```text
A(w)=1,
B(w)=1+p sum_(m>=1)m! w^m                            (17)
```

have the same reduction in `F_p[[w]]`, while `A` is entire and `B` has
radius of convergence zero. Thus no theorem about coefficient growth,
spectral contraction, Abel tails, or ordered extraction follows from the
finite-characteristic image alone.

For the binary fibbinary and Moser series, reduction happens to retain every
support bit, so (7) is stronger than a generic parity shadow. Even there, the
block-tail and harmonic-mass conclusions of THM-2005 still require their
positive block-count estimates; algebraicity in characteristic two supplies
no inequality-preserving functor back to those conclusions.

This is exactly the safe portion of the THM-2022 analogy:

```text
retain a whole packet -> apply Frobenius -> preserve its nonzero address.
```

The NC2 exit is algebraic nonvanishing, so that packet is sufficient there.
The Shunia exit is an ordered integer quotient, and the THM-2005 tail is an
Archimedean summability statement, so each requires an additional sidecar.
