---
id: THM-2021
title: "GMC(2) PROPORTIONAL CENTRAL RESONANCE REFINEMENT: exact toral-times-radial factorization, Legendre zero geometry, and finite-zero recurrence"
status: >
  REFINEMENT; PARTS A--C PROVED. The exact all-charge factorization and the
  symmetric Legendre transform are proved here. THM-2018 now gives the stronger
  NC2 conclusion on the entire proportional hypersurface h=kappa*b^r for every
  primitive charge pair: its root-of-unity EGF argument proves A_m(kappa) is
  nonzero arbitrarily far out, which is enough because EMP makes L(b^m)
  nonzero eventually. Thus neither the announced Legendre finite-zero theorem
  nor higher-charge HYP-8771 is an NC2 dependency. In the symmetric case the
  elementary no-consecutive-zero recurrence itself already gives a
  self-contained full closure. The external theorem remains a sharper result
  about the zero profile. Every possible symmetric zero parameter is negative
  real. See MISTAKE-213. THM-2022 subsequently proves full NC2 and GMC(2) for
  arbitrary support.
source: codex-2026-07-21-gmc2-proportional-central-resonance
depends_on: [THM-2018, THM-2017, THM-1510, THM-1695]
related: [THM-2018, THM-2019, THM-2020, THM-1660, THM-1670, HYP-8766, HYP-8769, HYP-8771, HYP-8772, MISTAKE-212, MISTAKE-213]
script: 04-computation/gmc2_proportional_legendre_finite_recurrence_thm2021.py
output: 05-knowledge/results/gmc2_proportional_legendre_finite_recurrence_thm2021.out
external:
  - "Dan Mangoubi, On Common Roots of Legendre Polynomials, IAS Special Year Learning Seminar, February 18, 2026, joint work with Borys Kadets and Adi Weller Weiser: https://www.ias.edu/math/events/special-year-learning-seminar-41"
  - "NIST Digital Library of Mathematical Functions, Chapter 18 (Legendre generating functions and zeros): https://dlmf.nist.gov/18.12 and https://dlmf.nist.gov/18.16"
---

# THM-2021 -- proportional central resonance and Legendre finite recurrence

> **Post-incoming-work synthesis.** The completed version of THM-2018 subsumes
> this theorem's NC2 application. Its exact EGF
> `exp(t) Phi_(p0,q0)(kappa*t^r)` and a nontrivial `r`th-root rotation show that
> the toral factor is nonzero at arbitrarily large levels. Since EMP makes the
> radial factor nonzero at every sufficiently large level, the two sequences
> intersect. THM-2021 is retained for its independent all-charge factorization,
> countable exceptional-set description, and sharp symmetric Legendre geometry.
> Moreover its no-consecutive-zero recurrence already closes the symmetric
> proportional slice directly; the contrary first reading is MISTAKE-213.
> THM-2022 now closes arbitrary support. The Legendre discussion remains a
> valid stronger zero-set refinement and its external attribution is not used
> in the Gaussian proof.

This theorem fills a central part of HYP-8766 by combining three ideas that had
not previously been put together in the repository:

1. retain the **whole channel profile** before applying the scalar radial
   functional (the support-Dirichlet lesson of THM-2005);
2. recognize when that profile has rank one in the radial variable; and
3. distinguish **finite recurrence of zeros** from both consecutive-root
   descent and the still weaker non-eventual-zero property sufficient for NC2.

The distinction is load-bearing. A sequence with no two consecutive zeros can
still have infinitely many zeros, so the elementary recurrence does not prove
finite recurrence. The completed THM-2018 also challenges the assumption that
finite recurrence is needed for NC2: unbounded toral nonzeros already intersect
EMP's cofinite radial nonzero tail.

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
mixed-sign/complex coefficients. Concurrent THM-2018 first reserved the
symmetric case and subsequently proved the all-charge version together with
the stronger EGF non-eventual-zero argument. Equations (3)--(4) retain an
independent derivation and the finer zero-set criterion.

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
property is precisely HYP-8771. It is stronger than NC2 needs: THM-2018 proves
that `A_m(kappa)` is nonzero for arbitrarily large `m`, and those levels must
eventually meet the cofinite nonzero tail of `L(b^m)`.

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

If `kappa!=1/4`, two consecutive zeros in (10) propagate backward to
`S_0=0`, contradicting `S_0=1`. At `kappa=1/4`, (8) reduces to
`(1-2t)^(-1/2)` and all coefficients are nonzero. Hence `S_m(kappa)` is
nonzero at arbitrarily large levels for every finite `kappa`. If all moments
vanished with `b!=0`, EMP would make `L(b^m)` nonzero on a cofinite tail, so
(3) would force `S_m(kappa)=0` on that entire tail, a contradiction. Therefore:

> **Self-contained full symmetric closure.** NC2 holds on every
> arbitrary-degree hypersurface `s*a(s)*c(s)=kappa*b(s)^2`.

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

## D. The announced 2026 Legendre theorem sharpens the symmetric zero profile

The IAS announcement for joint work of Mangoubi--Kadets--Weller Weiser states:
for every fixed nonzero `x`, only finitely many Legendre polynomials `P_m` vanish
at `x`. Applied to (9), this says that `S_m(kappa)` is eventually nonzero for
every finite `kappa`. Section B then gives another proof of:

> **NC2 corollary (already unconditional by Section C and THM-2018).** NC2
> holds on the full arbitrary-degree hypersurface
> `s a(s)c(s)=kappa b(s)^2`.

The official February 2026 seminar page states the result but identifies the
work as joint work rather than linking a public proof. We therefore record the
finite-zero conclusion exactly and do not relabel the external input as proved
in-repository. The NC2 corollary itself does not depend on it: THM-2018 only
needs infinitely many nonzero levels, proved by root-of-unity symmetry.
Stieltjes' stronger 1890 conjecture says that at most one `P_m` can vanish at a
fixed nonzero `x`; THM-2021 needs only finiteness, not that stronger conjecture.

## E. What the elementary recurrence does and does not prove

From (10), consecutive `S_m,S_(m+1)` have no common zero. This does **not** imply
that the zero set in the index `m` is finite, so it does not prove the stronger
Legendre statement in Section D or HYP-8771. It **does**, however, imply nonzero
levels arbitrarily far out, exactly what NC2 needs after EMP. The first version
of this theorem combined the supports backward and called consecutive-root
descent insufficient; that error is recorded as MISTAKE-213. The completed
THM-2018 proves the same sufficient non-eventual-zero property for all charges
by root-of-unity symmetry.

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

This theorem alone does **not** prove full NC2. Together with THM-2018 it closes the complete
infinite-dimensional rank-one subvariety `h=kappa*b^r` of the zero-offset
three-weight resonance for every charge pair. THM-2021 additionally identifies
the symmetric zero geometry and gives a stronger countable/off-axis statement
level by level. HYP-8771 remains a literal finite-zero-recurrence problem, not
an NC2 blocker. THM-2020's finite-place channel separation is still relevant to
that sharper sequence question: every exceptional `kappa` is algebraic, so its
recurrence can be attacked through factorial valuations at a finite place.
THM-2022 subsequently carries that finite-place philosophy to the whole
support: its Frobenius lowest-face theorem proves full NC2/GMC(2).
