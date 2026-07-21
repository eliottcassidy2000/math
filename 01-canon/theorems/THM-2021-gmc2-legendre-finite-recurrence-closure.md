---
id: THM-2021
title: "GMC(2) PROPORTIONAL CENTRAL RESONANCE: exact toral-times-radial factorization for every primitive charge pair, unconditional NC2 off an explicit countable exceptional set, and the symmetric exceptional set is exactly finite recurrence of Legendre zeros"
status: >
  PARTS A--C PROVED. For P=Z^p a(s)+b(s)+Zbar^q c(s), h=s^(pq/g)a^(q/g)c^(p/g),
  r=(p+q)/g, and h=kappa*b^r, every moment factors exactly as
  E[P^m]=A_m^(p0,q0)(kappa)L(b^m). Hence NC2 holds whenever A_m(kappa) is eventually
  nonzero; in particular it holds off an explicit countable exceptional set for every
  primitive charge pair. For p0=q0=1, A_m is a Legendre transform and every possible
  exceptional kappa is negative real. PART D is CONDITIONAL ON AN ANNOUNCED EXTERNAL
  RESULT: the February 2026 IAS announcement by Mangoubi--Kadets--Weller Weiser states
  that a fixed nonzero point is a zero of only finitely many Legendre polynomials; this
  would close the entire symmetric proportional hypersurface. No public proof/preprint
  of that new theorem was located, so the full symmetric corollary is not counted as a
  self-contained canon proof. The higher-charge finite-recurrence statement is HYP-8771.
source: codex-2026-07-21-gmc2-proportional-central-resonance
depends_on: [THM-2017, THM-1510, THM-1695]
related: [THM-2018, THM-2019, THM-2020, THM-1660, THM-1670, HYP-8766, HYP-8769, HYP-8771]
script: 04-computation/gmc2_proportional_legendre_finite_recurrence_thm2021.py
output: 05-knowledge/results/gmc2_proportional_legendre_finite_recurrence_thm2021.out
external:
  - "Dan Mangoubi, On Common Roots of Legendre Polynomials, IAS Special Year Learning Seminar, February 18, 2026, joint work with Borys Kadets and Adi Weller Weiser: https://www.ias.edu/math/events/special-year-learning-seminar-41"
  - "NIST Digital Library of Mathematical Functions, Chapter 18 (Legendre generating functions and zeros): https://dlmf.nist.gov/18.12 and https://dlmf.nist.gov/18.16"
---

# THM-2021 -- proportional central resonance and Legendre finite recurrence

This theorem fills a central part of HYP-8766 by combining three ideas that had
not previously been put together in the repository:

1. retain the **whole channel profile** before applying the scalar radial
   functional (the support-Dirichlet lesson of THM-2005);
2. recognize when that profile has rank one in the radial variable; and
3. ask for **finite recurrence of zeros**, not merely the consecutive-root
   descent used for Hermite/Legendre base cases.

The third distinction is load-bearing. A sequence with no two consecutive zeros
can still have infinitely many zeros, so the elementary three-term recurrence
alone does not give the eventual radial vanishing needed by EMP.

## A. Exact factorization for every primitive charge pair

Let

```text
P = Z^p a(s) + b(s) + Zbar^q c(s),       s=Z Zbar,
g = gcd(p,q),       p0=p/g, q0=q/g,      r=p0+q0,
h = s^(pq/g) a(s)^q0 c(s)^p0.
```

THM-2017 gives the exact primitive-return expansion

```text
E[P^m] = sum_(0<=k<=m/r)
  m!/((q0 k)!(p0 k)!(m-rk)!) L(h^k b^(m-rk)).       (1)
```

Assume the **proportional central-resonance identity**

```text
h = kappa b^r.                                       (2)
```

Then every radial word in the `k`th channel is the same power:

```text
h^k b^(m-rk) = kappa^k b^m.
```

Consequently

```text
E[P^m] = A_m^(p0,q0)(kappa) L(b^m),                  (3)

A_m^(p0,q0)(kappa)
 = sum_(0<=k<=m/r)
   m!/((q0 k)!(p0 k)!(m-rk)!) kappa^k.              (4)
```

This is exact, not asymptotic. It includes arbitrary radial degree and genuine
mixed-sign/complex coefficients. Concurrent THM-2018 found (3) independently
for `p0=q0=1`; (3)--(4) record the all-charge version and the correct zero-set
criterion.

## B. The eventual-nonzero criterion closes NC2

Suppose that, for the fixed `kappa`, `A_m(kappa)` is nonzero for every sufficiently
large `m`. If all Gaussian moments of `P` vanish, (3) gives

```text
L(b^m)=0                                               (m sufficiently large). (5)
```

Choose `N>=1` beyond the exceptional levels and set `B=b^N`. Then
`L(B^j)=L(b^(Nj))=0` for every `j>=1`. The exponential moment problem (EMP,
THM-1510) forces `B=0`, hence `b=0`. Equation (2) gives `h=0`; because `C[s]`
is a domain, `a=0` or `c=0`. Thus `P` is charge-one-sided.

Define

```text
E_(p0,q0) = union_(m>=1) {kappa in C : A_m^(p0,q0)(kappa)=0}.  (6)
```

Every `A_m` has constant term one, so it is not the zero polynomial. Therefore
`E_(p0,q0)` is countable. For `kappa` outside (6), all the factors in (3) are
nonzero, and the preceding argument proves NC2. This gives an unconditional,
arbitrary-degree, all-primitive-charge closure on the complement of an explicit
countable set inside the exact zero-offset resonance band.

More generally, even at a point of (6), NC2 follows as soon as the fixed point is
a zero of only finitely many members of the family. That finite-recurrence
property is precisely HYP-8771.

## C. Symmetric primitive charges are Legendre

For `p0=q0=1`, write `S_m=A_m^(1,1)`. Then

```text
S_m(kappa)=sum_(0<=k<=m/2) m!/(k!^2(m-2k)!) kappa^k, (7)

sum_(m>=0) S_m(kappa)t^m
 = 1/sqrt(1-2t+(1-4kappa)t^2).                       (8)
```

The Legendre generating function gives, when `1-4kappa != 0`,

```text
S_m(kappa)
 = (1-4kappa)^(m/2)
   P_m((1-4kappa)^(-1/2)),                            (9)
```

with the limiting value at `kappa=1/4` read from (8), and it is nonzero.
Equivalently,

```text
(m+1)S_(m+1)=(2m+1)S_m-m(1-4kappa)S_(m-1).          (10)
```

All zeros of `P_m` are real and lie in `(-1,1)`. A finite `kappa` for which
`S_m(kappa)=0` must therefore have

```text
kappa = (1-x^(-2))/4 < 0                             (11)
```

for a nonzero Legendre zero `x`. Hence:

> **Unconditional symmetric closure.** On `h=kappa b^2`, NC2 holds for every
> `kappa` not in the explicit countable set (11). In particular it holds for
> every `kappa` outside the negative real axis.

This phase localization is much sharper than coefficient positivity: arbitrary
complex radial `b` is allowed. Only an exactly opposite-phase proportionality
can even reach the residual.

## D. The announced 2026 Legendre theorem would close the entire symmetric hypersurface

The IAS announcement for joint work of Mangoubi--Kadets--Weller Weiser states:
for every fixed nonzero `x`, only finitely many Legendre polynomials `P_m` vanish
at `x`. Applied to (9), this says that `S_m(kappa)` is eventually nonzero for
every finite `kappa`. Section B would then prove:

> **Conditional corollary (on the announced Legendre finiteness theorem).** NC2
> holds on the full arbitrary-degree hypersurface
> `s a(s)c(s)=kappa b(s)^2`.

The official February 2026 seminar page states the result but identifies the
work as joint work rather than linking a public proof. We therefore record the
deduction exactly and do not relabel the external input as proved in-repository.
Stieltjes' stronger 1890 conjecture says that at most one `P_m` can vanish at a
fixed nonzero `x`; THM-2021 needs only finiteness, not that stronger conjecture.

## E. Why the elementary recurrence is insufficient

From (10), consecutive `S_m,S_(m+1)` have no common zero. That proves only that
the zero set in the index `m` has gaps. It does **not** imply that the zero set is
finite, and hence does not imply (5). This corrects the proposed closing sentence
in the reservation version of concurrent THM-2018. The distinction mirrors the
support-profile lesson: an endpoint or a local adjacency relation can erase the
tail property the proof actually needs.

## F. Verification and Tournament Analysis

The exact script verifies:

- direct Gaussian expansion = channel sum = factorized expression through
  `m=7` on realizable `(+1,-1)` and `(+1,-2)` polynomial slices;
- (8) and (10) exactly through levels 20 and 30;
- negative-real root localization numerically through `m=30`; and
- exact pairwise coprimality of the return polynomials on 1,379 level pairs:
  `(+1,-1)` through 36, `(+1,-2)` through 30, `(+1,-3)` through 26, and
  `(+2,-3)` through 22.

For Tournament Analysis the vertices are **moment levels**, not channels or
individual roots. The pairwise observable is
`deg gcd(A_m,A_n)`. The chronological edge `m->n` is flipped exactly when a
shared-root event occurs; chronological order is the tie Hamiltonian path. All
four computed tournaments have zero flips, zero directed triangles, singleton
SCCs, transitive score histogram, and one Hamiltonian path. This quotient
preserves exact pairwise root recurrence but destroys root identity,
multiplicity, and location. Channels were rejected as vertices because their
sets change with `m`; proof obligations were rejected because they do not retain
the recurrence predicate.

## G. Scope and next target

This does **not** prove full NC2. It closes a new infinite-dimensional rank-one
subvariety of the zero-offset three-weight resonance, unconditionally outside a
countable set, and shows that the symmetric remainder is already the subject of
a very recent external finiteness theorem. For higher primitive charges the
right target is HYP-8771. Concurrent THM-2020's finite-place channel separation
is especially relevant there: every exceptional `kappa` is algebraic, so the
archimedean recurrence question can be attacked through factorial valuations at
a finite place rather than by another magnitude estimate.
