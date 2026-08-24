---
id: THM-3966
title: "Euler rigidity forces split-A1-fibre affine-plane boundaries to have disjoint line normalizations"
status: >
  PROVISIONAL PROOF CANDIDATE / AWAITING INDEPENDENT TEXT AUDIT. Let X be a
  normal affine complex surface with Cl(X)=Z^n and compactly supported Euler
  characteristic 1+n. If X contains a dense open U isomorphic to A2, then
  X minus U is pure divisorial, has exactly n prime components, every
  component has normalization A1, and distinct components are disjoint. A
  split-A1 fibration with n exceptional fibres, each supported on two
  disjoint A1s, automatically has Euler characteristic 1+n. Consequently,
  in such a finite normalization every ramification prime must have A1
  normalization and ramification primes cannot meet. This is a necessary
  boundary invoice, not a sufficient construction of an affine-plane open.
source: jc-cohn3709 / post-THM-3964 class-Euler synthesis, 2026-08-24
depends_on:
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
  - THM-3964-polynomial-graph-hidden-double-root-normalization
---

# THM-3966 -- the class rank and Euler tariff rigidify the whole boundary

**PROVISIONAL PROOF CANDIDATE / AWAITING INDEPENDENT TEXT AUDIT.** Work first
over `C`; Section 4 records the characteristic-zero scope. Write `chi` for
compactly supported topological Euler characteristic.

## 1. Boundary Euler-rigidity theorem

Let `X` be a normal integral affine complex surface satisfying

```text
Cl(X) isomorphic to Z^n,                 chi(X)=1+n.       (1)
```

Suppose there is a dense open immersion

```text
U isomorphic to A2  ->  X.                                  (2)
```

Then the closed boundary `Z=X minus U` is pure of dimension one. It has
exactly `n` prime components

```text
Z=D_1 union ... union D_n,                                  (3)
```

the normalization of every `D_i` is `A1`, and

```text
D_i intersection D_j is empty whenever i!=j.                (4)
```

The conclusion concerns the normalizations: a component itself may still be
an intrinsically unibranch singular affine curve.

### 1.1 Affineness removes isolated boundary components

Let `D` be the union of the divisorial components of `Z`. Suppose a closed
point `p` were an isolated component of `Z` outside `D`. Since the residual
closed set `Z minus D` is zero-dimensional and hence finite, choose a normal
affine neighbourhood

```text
p in W subset X minus D,
W intersection (Z minus D)={p}.                             (5)
```

Both `W` and `U` are affine opens of the separated scheme `X`; therefore
their intersection is affine. But

```text
W intersection U=W minus {p}.                              (6)
```

Normality gives `S2`, and the missing point has codimension two, so Hartogs
extension gives

```text
Gamma(W,O_W) isomorphic to Gamma(W minus {p},O).             (7)
```

The open immersion in `(6)` is a morphism between affine schemes induced by
the isomorphism `(7)`, hence is itself an isomorphism. This contradicts the
missing point. Thus `Z=D` set-theoretically: no isolated point can be added
to, or hidden beside, the divisorial boundary.

### 1.2 The class group counts the curves exactly

THM-3922 says that the prime boundary classes form a `Z`-basis of `Cl(X)`.
Equation `(1)` therefore forces exactly `n` prime components, proving `(3)`.
This is an equality, not merely a lower bound: a hypothetical source cannot
delete an additional etale divisor once the `n` class slots are occupied.

### 1.3 Euler equality has no slack

THM-3920 makes each boundary prime rational and unibranch at every point of
`X`. Since `D_i` is affine, its normalization has the form

```text
D_i^nu isomorphic to P1 minus S_i,       s_i=#S_i>=1.       (8)
```

Unibranchness and algebraic closure make the finite normalization map
bijective on points. A finite bijective complex map is a topological
homeomorphism, so

```text
chi(D_i)=chi(D_i^nu)=2-s_i<=1.             (9)
```

For a point `p` lying on boundary components, let `k_p` be their number.
Only finitely many points have `k_p>=2`. The normalization of the reduced
union is the disjoint union of the component normalizations, and its fibre
over `p` has `k_p` points. Constructible additivity therefore gives the exact
incidence formula

```text
chi(Z)=sum_i (2-s_i) - sum_p (k_p-1).                     (10)
```

On the other hand, closed-open additivity, `(1)`, and `chi(A2)=1` give

```text
chi(Z)=chi(X)-chi(U)=n.                                   (11)
```

Combining `(10)` and `(11)` yields the nonnegative zero invoice

```text
sum_i (s_i-1) + sum_p (k_p-1)=0.                          (12)
```

Every summand vanishes. Thus `s_i=1` for every `i`, so every normalization
is `A1`, and no point lies on two distinct boundary primes. This proves
`(4)` and the theorem.

## 2. The split-A1-fibre calculator

Let `tau:X->A1_t` be a morphism from a complex affine surface. Suppose there
is a set `S` of `n` points such that

```text
X over (A1 minus S) isomorphic to (A1 minus S) times A1,   (13)
```

and the reduced fibre above every `a in S` is the disjoint union of two
copies of `A1`. Multiplicities in the scheme fibre are irrelevant to the
reduced complex topology. Stratifying by `S` gives

```text
chi(X)=(1-n)chi(A1)+n(2chi(A1))=1+n.                       (14)
```

Consequently, if `X` is normal and `Cl(X)=Z^n`, Section 1 applies to every
affine-plane open in `X`.

This packages a recurring normalization pattern. Inverting the degeneration
polynomial often exposes the product in `(13)`; each distinct root then
splits one affine-line fibre into two disjoint affine lines. The same integer
`n` appears both as the Nagata class rank and as the Euler excess. Equality
between those two ledgers is what removes all boundary slack.

## 3. Finite-map ramification corollary

Let `pi:X->A2` be finite with `X` as in Sections 1--2. If a same-function-
field Keller source existed, normalization-form Zariski Main would give an
open immersion

```text
A2 isomorphic to U -> X                                    (15)
```

on which `pi` is etale. Every ramification prime of `pi` must therefore be
a component of `X minus U`. Sections 1--2 impose both conditions

```text
each ramification-prime normalization is A1;
distinct ramification primes are disjoint.                 (16)
```

Thus either a second projective puncture on one ramification normalization,
or one finite incidence between two ramification primes, excludes the Keller
atlas. Notice that `(16)` is stronger than the class-primitivity test: a
primitive ramification class may still carry too many punctures.

## 4. Scope, descent, and hostile controls

The proof used ordinary complex Euler characteristic only for convenience.
A finite-type characteristic-zero Keller packet and its proposed open
immersion descend to a finitely generated subfield, which embeds in `C`.
For applications over another algebraically closed characteristic-zero
field, the class and fibre hypotheses must be checked geometrically after
that base change. Equivalently, one may run the proof with compactly
supported l-adic Euler characteristic.

Each input is load-bearing:

- Class rank alone counts components but does not control their punctures or
  incidences.
- Euler characteristic alone does not bound the number of deleted divisors.
- Unibranchness is needed to identify a component's Euler characteristic with
  that of its normalization; a multi-branch singularity carries an extra
  identification correction.
- Affineness of `U` is needed in Section 1.1. Deleting a point from a normal
  surface can produce a nonaffine open with the same regular functions.
- The exceptional-fibre calculation counts distinct roots, not their
  multiplicities.

Finally, the conclusion is necessary, not sufficient. A basis of disjoint
boundary primes with affine-line normalizations does not construct an open
`A2`, prove its affineness, or make a finite map etale. **QED, pending text
audit.**
