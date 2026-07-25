---
id: THM-2175
title: "Good-place certificates for factor-edge nonpowers"
status: >
  PROVED. Let R be a nonzero polynomial over a number field and M>=2.
  At every finite place where R is integral and its leading coefficient is
  a unit, a geometrically non-M-th-power reduction certifies that R is not
  an M-th power over the algebraic closure. Outside an explicit finite set
  determined by the squarefree multiplicity packet, discriminants, and
  pairwise resultants, reduction preserves the entire root-multiplicity
  packet and hence preserves the M-th-power predicate in both directions.
  Applied to an algebraic THM-2136 scalar factor edge, one good-place
  nonpower rules out THM-2134's coarsened-power branch and therefore forces
  exactly its terminally-short inequality. It does not exclude that
  inequality, glue different factors, or prove planar JC.
source: codex-2026-07-24-cross-frontier-face-transfer
depends_on:
  - THM-2134-preterminal-factor-edge-power-dichotomy
  - THM-2136-toric-scalarization-of-factor-edge-powers
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
---

# THM-2175 -- good-place certificates for factor-edge nonpowers

The relevant finite-place invariant is not coefficient support and not one
distinguished root. It is the complete packet of root multiplicities.

## 1. Geometric powers and the multiplicity packet

Let `K` be a number field, let

```text
0!=R(S) in K[S],                    M>=2,             (1)
```

and call `R` a **geometric `M`-th power** if

```text
R=Q^M                 for some Q in Kbar[S].         (2)
```

Factor `R` in characteristic zero and group its distinct irreducible factors
by multiplicity. There is a unique decomposition

```text
R=c product_(e in E) A_e(S)^e,                       (3)
```

where

```text
c in K^*,
each A_e is monic, nonconstant, and squarefree,
gcd(A_e,A_f)=1 for e!=f.                              (4)
```

The finite set `E` is the root-multiplicity packet. Over `Kbar`, every root
of `A_e` occurs in `R` with multiplicity exactly `e`. Since every nonzero
constant has an `M`-th root in an algebraically closed field,

```text
R is a geometric M-th power
  iff
M divides e for every e in E.                        (5)
```

This criterion remains valid over an algebraically closed field of arbitrary
characteristic. In particular, the residue characteristic is allowed to
divide `M`.

## 2. A one-place nonpower certificate

Let `pfrak` be a nonzero prime ideal of `O_K`. Assume

```text
R in O_(K,pfrak)[S],
lc(R) in O_(K,pfrak)^*.                              (6)
```

Thus coefficientwise reduction is defined and preserves the degree. Write
`Rbar` for the reduction in `k(pfrak)[S]`.

> **One-place certificate.** If `Rbar` is not an `M`-th power in
> `k(pfrak)bar[S]`, then `R` is not an `M`-th power in `Kbar[S]`.

Equivalently, geometric `M`-th powers remain geometric `M`-th powers at
every degree-preserving integral place.

### Proof

Suppose `R=Q^M` over a finite extension `L/K` containing the coefficients
and roots of `Q`, and choose a prime `Pfrak` of `L` above `pfrak`. Put

```text
a=lc(Q),                         lc(R)=a^M.           (7)
```

The second condition in (6) makes `a` a `Pfrak`-adic unit. The monic
polynomial

```text
Q_0=Q/a
```

has roots among the roots of the monic integral polynomial

```text
R/lc(R)=Q_0^M.                                        (8)
```

Those roots are integral over `O_(K,pfrak)`. Their elementary symmetric
functions, the coefficients of `Q_0`, are therefore integral over that
local ring. They lie in `L`, so they belong to the integrally closed
valuation ring `O_(L,Pfrak)`. Hence `Q` has `Pfrak`-integral coefficients
and unit leading coefficient. Reducing (2) gives

```text
Rbar=Qbar^M             in k(Pfrak)[S]
                       subset k(pfrak)bar[S].         (9)
```

The contrapositive proves the certificate.

The algebraic closure in the residue test is load-bearing. For example,
`2S^2` is not a square in `F_3[S]` but is a square in `F_3bar[S]`; over
`C` it is of course also a square. A native-field constant obstruction is
not a factor-edge obstruction.

## 3. Exact preservation outside a finite bad set

The converse to Section 2 can fail at a collision prime, but only at
finitely many primes.

Choose a finite set `Sigma` containing every prime of `K` at which:

1. `c` is not a unit;
2. some coefficient of an `A_e` is not integral;
3. `disc(A_e)` is not a unit for some `e`; or
4. `Res(A_e,A_f)` is not a unit for some `e!=f`.

All displayed discriminants and resultants are nonzero by (4), so `Sigma`
is finite. For every `pfrak` outside `Sigma`, reduction of (3) gives

```text
Rbar=cbar product_(e in E) Abar_e(S)^e,              (10)
```

where every `Abar_e` is squarefree and distinct `Abar_e,Abar_f` are
coprime. Thus (10) has exactly the same root-multiplicity packet `E` over
`k(pfrak)bar` as (3) has over `Kbar`. By (5),

```text
R is an M-th power in Kbar[S]
  iff
Rbar is an M-th power in k(pfrak)bar[S]              (11)
```

for every `pfrak` outside `Sigma`.

This proves more than the existence of one detecting prime: every
sufficiently good finite place gives the same geometric power verdict.
Nevertheless, only the nonpower direction from Section 2 is needed for the
Keller application.

## 4. Exact factor-edge consequence

Assume all hypotheses and notation of THM-2134 for a reduced proper-power
planar Keller pair, a fixed factor `pi` of the primitive top root `h`, and
its globally exposed nonradial factor edge:

```text
tau*d>e,
M=m/gcd(m,n).                                         (12)
```

Assume additionally that the coefficients involved are algebraic. Then the
THM-2136 scalar edge polynomial

```text
R_pi(S) in K[S]                                      (13)
```

lies in some number field `K`. THM-2136 proves the exact equivalence

```text
P_pi(T) is an M-th power in C(t)[T]
  iff
R_pi(S) is an M-th power in C[S].                    (14)
```

Let `pfrak` be any finite place of `K` satisfying (6) for `R_pi`. If

```text
Rbar_pi is not an M-th power in k(pfrak)bar[S],      (15)
```

Section 2 says that `R_pi` is not a geometric `M`-th power. Root
multiplicities do not change after extending `Kbar` to `C`, so the right
side of (14) fails. Therefore the power branch of THM-2134 fails.

The dichotomy in THM-2134 then forces exactly

```text
D <= L+r_0 floor(e*n/(tau*r_0)).                     (16)
```

This is the complete licensed consequence. The Keller pair itself is not
being reduced modulo `pfrak`; only the scalar edge polynomial is. Therefore
no finite-characteristic Jacobian or Poisson assertion is hidden in the
argument.

If the original coefficients and the selected factor are algebraic, the
scalarization data are algebraic as well, after enlarging `K` finitely.
Alternatively, (13) may simply be taken as the explicit arithmetic
hypothesis of the corollary.

## 5. Hostile controls and logical direction

Coefficient support does not detect the power branch. The two full-support
quadratics

```text
R_square=1+2S+S^2=(1+S)^2,
R_simple=1+S+S^2                                    (17)
```

have the same monomial support and all coefficients nonzero, but the first
is a square and the second has two simple roots and is not a square.

Nor does one power specialization prove a characteristic-zero power. The
second polynomial in (17) becomes

```text
1+S+S^2=(S-1)^2             over F_3,                (18)
```

because its discriminant `-3` vanishes there. This is precisely a root
collision excluded by `Sigma`. By contrast, one geometrically **nonpower**
degree-preserving reduction always certifies a characteristic-zero
nonpower by Section 2.

A scalar discriminant is also insufficient for general `M`: for example

```text
(S-1)^2(S-2)
```

has discriminant zero but is not a square because its multiplicity packet
is `{2,1}`. The complete packet modulo `M`, rather than squarefreeness alone,
is the invariant.

Finally, (16) is an alternative, not a contradiction. This theorem does
not:

1. exclude the terminally-short inequality;
2. make roots for different factors compatible;
3. lift local scalar roots to one polynomial approximate root;
4. treat transcendental factor-edge coefficients by specialization; or
5. prove `JC(2)` or `DC(2)`.

The exact transfer ledger is:

```text
source:       algebraic scalar edge R_pi;
map:          degree-preserving reduction at pfrak;
preserved:    geometric M-th-power failure, and outside Sigma the full
              root-multiplicity packet;
lost:         root labels, coefficient magnitude, and cross-factor gluing;
sidecar:      leading-coefficient unit plus discriminant/resultant exclusions;
exit:         only THM-2134's terminally-short inequality (16).
```

QED.
