---
id: THM-2102
title: "Power-free weighted faces are tame and proper-power faces obey a first-defect descent"
status: >
  PROVED. If one component of a planar Keller pair has, for some positive
  weight, a leading face which is not a nontrivial perfect power, then that
  component is triangular and the pair is a polynomial automorphism. The
  proof repeatedly shears the mate by powers of the component until the
  leading weighted bracket is one, then applies THM-2113. Thus any hypothetical
  planar counterexample has a proper-power leading face for every positive
  weight. In the remaining common-root branch, the first weighted defect
  satisfies an exact Hamiltonian equation; a terminal first defect forces the
  primitive root to be weighted-linear and one top multiplicity to equal one.
  This sharply refines but does not terminate the all-proper-power branch or
  prove JC(2).
source: codex-2026-07-22-JC2-power-free-face-descent
related:
  - THM-2045
  - THM-2063
  - HYP-8950
  - HYP-8955
---

# THM-2102 -- power-free face descent and the first proper-power defect

Let `w=(w_1,w_2)` be positive integer weights and put

```text
W=w_1+w_2.                                              (1)
```

For a nonzero polynomial `A`, write `deg_w A` for its weighted degree and
`in_w(A)` for the sum of its terms of maximum weight. The Poisson bracket is

```text
{A,B}=A_x B_y-A_y B_x.                                 (2)
```

A polynomial is **power-free** here when it is not `lambda*h^m` for a scalar
`lambda!=0` and an integer `m>=2`. Equivalently, the gcd of the multiplicities
of its irreducible factors is one. This is weaker than squarefreeness; for
example `x^2 y^3` is power-free.

## 1. The homogeneous common-root lemma

We first isolate the algebraic engine.

> **Lemma.** Let `F,G` be nonzero positive-degree `w`-homogeneous polynomials
> of degrees `delta,epsilon`. If `{F,G}=0`, then after scalar normalization
> there is a power-free `w`-homogeneous polynomial `h` and positive integers
> `m,n` such that
> ```text
> F=h^m,                    G=c h^n.                    (3)
> ```
> If `F` itself is power-free, then `G=lambda F^k` for an integer `k>=1`.

Weighted Euler identities and the zero bracket give

```text
delta F G_x=epsilon G F_x,
delta F G_y=epsilon G F_y.                             (4)
```

Hence the differential of `G^delta/F^epsilon` vanishes in `C(x,y)`, so

```text
G^delta=c F^epsilon.                                   (5)
```

Unique factorization proves (3) by dividing all common irreducible
multiplicities by their gcd. If the multiplicity gcd of `F` is one, (5)
forces `delta|epsilon`, and `G=lambda F^(epsilon/delta)`. This proves the
lemma.

## 2. Finite descent of a power-free leading face

Let

```text
f,g in C[x,y],                   {f,g}=1,               (6)
F=in_w(f),                       delta=deg_w f,
G=in_w(g),                       epsilon=deg_w g,
```

and assume `F` is power-free. If

```text
delta+epsilon>W,                                        (7)
```

then the component of maximum weight in (6) is `{F,G}` and must vanish. The
lemma gives

```text
G=lambda F^k.                                          (8)
```

Replace

```text
g by g-lambda f^k.                                     (9)
```

This preserves `{f,g}=1` and strictly lowers `deg_w g`. Positive integral
weights make the descent finite. Repeating (9) reaches a mate `g_*` with
leading face `H` of degree `eta` and

```text
delta+eta<=W.                                          (10)
```

Strict inequality in (10) is impossible: every possible monomial in
`{f,g_*}` would have negative weighted degree. Therefore

```text
delta+eta=W,                     {F,H}=1.               (11)
```

All lower-face cross terms have weight below zero, so the second equality in
(11) is the entire weight-zero component, not just one summand.

THM-2113 applies to the quasi-homogeneous Keller pair `(F,H)`. After possibly
swapping `x,y`,

```text
F=a x+b y^q,                    a!=0.                   (12)
```

In particular `delta=w_1`. Any monomial of `f` involving `x` has weight at
least `w_1=delta`; the only such top monomial is the displayed `a x`, and no
such monomial can occur on a lower face. Consequently

```text
f=a x+P(y).                                             (13)
```

This already proves that `f` is a coordinate. It also proves that the original
pair is an automorphism, not merely that one component is a coordinate. In
the determinant-one coordinates

```text
u=f,                         v=y/a,                     (14)
```

equation (6) becomes `partial g/partial v=1`, so

```text
g=v+Q(u).                                               (15)
```

Thus `(f,g)` is triangular in polynomial coordinates and has a polynomial
inverse.

> **Power-free face theorem.** If `in_w(f)` is power-free for one positive
> weight and `f` has a polynomial Jacobian mate, then (13)--(15) hold after a
> variable swap. Hence every component of a hypothetical JC(2) counterexample
> has a nontrivial perfect-power leading face for **every** positive weight.

## 3. Brieskorn principal faces, with arbitrary lower terms

Take integers `p,q>=1`, the weight `w=(q,p)`, and suppose

```text
in_w(f)=a x^p+b y^q,                   ab!=0.           (16)
```

The binomial in (16) is squarefree: a common divisor of it and both partial
derivatives would divide both `x` and `y`. It is therefore power-free. The
theorem gives the exact classification

```text
f has a polynomial Jacobian mate   iff   min(p,q)=1.    (17)
```

Indeed, if `p,q>=2`, (12) contradicts (16). If `p=1`, every monomial involving
`x` has weight at least the top degree `q`, so

```text
f=a x+P(y)
```

and `y/a` is an explicit mate; `q=1` is symmetric. Thus `x^2+y^3`, with
**arbitrary lower `(3,2)`-weight terms**, has no Jacobian mate. No assertion
about one place at infinity or generic-fiber properness is used.

## 4. Exact first-defect law in the proper-power branch

The preceding descent stops being automatic precisely when the top face is a
proper power. The first failure nevertheless has a rigid form.

Suppose a Keller pair has normalized leading faces

```text
in_w(f)=h^m,                   in_w(g)=c h^n,            (18)
```

where `h` is power-free and homogeneous of degree `d`. Put

```text
D=(m+n)d-W>0,                 p=min(m,n).               (19)
```

Let `rho>0` be the smallest defect below either leading face. Write the full
faces at that defect as

```text
f=h^m+A+(faces of defect >rho),
g=c h^n+B+(faces of defect >rho),                       (20)
```

where `deg_w A=md-rho`, `deg_w B=nd-rho`, and a missing face is recorded as
zero. Define

```text
K=m h^(m-p) B-c n h^(n-p) A.                           (21)
```

The bracket contribution at the first defect is exactly

```text
h^(p-1){h,K},                                           (22)
```

of weight `D-rho`. Therefore

```text
rho<=D.                                                 (23)
```

If `rho<D`, the positive-weight component must vanish:

```text
{h,K}=0.                                                (24)
```

The common-root lemma and power-freeness of `h` then give

```text
K=0                                  if d does not divide rho,
K=lambda h^(m+n-p-rho/d)             if d divides rho, (25)
```

where the displayed exponent is required to be a nonnegative integer. For
example, when `n>=m`, every nonresonant first defect obeys the fractional-power
coefficient law

```text
B=c(n/m) h^(n-m) A.                                    (26)
```

If `rho=D`, the degree-zero component of the Keller equation is

```text
h^(p-1){h,K}=1.                                         (27)
```

Thus

```text
p=1,                         {h,K}=1.                  (28)
```

THM-2113 makes `h` weighted-linear. In particular, if the first component has
a genuine proper-power top face `m>=2`, terminality at the first defect forces

```text
n=1.                                                    (29)
```

Equations (23)--(29) are the first-defect descent. They do not assert that all
later resonant defects terminate.

## 5. Sharp controls and the surviving wall

The power-free hypothesis cannot simply be erased. With ordinary weights,

```text
f=(x+y)^2+x,                    g=x+y                    (30)
```

has `{f,g}=1`. Its top face is the proper square `h^2`, and its terminal first
defect has `m=2,n=1,rho=D=1,K=-x`, with `{h,K}=1` exactly as (27)--(29)
predict. The pair is tame after a linear change but is not forced into the
axis-triangular form by its original top face.

The theorem also recovers the `gcd(r+1,s)=1` slice of THM-2045: the monomial
`x^(r+1)q^s` is power-free exactly on that slice. THM-2045 remains genuinely
stronger when the gcd is greater than one, which is precisely this theorem's
proper-power wall.

The live planar-Jacobian target is now narrower: prove termination through
the resonant laws (25), or show that a hypothetical pair cannot maintain a
proper-power leading face under every positive-weight change. JC(2) and DC(2)
remain open. QED.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that one must first turn a face polynomial into
a statement about places at infinity. The descent uses only the graded
Hamiltonian bracket, UFD multiplicities, and polynomial shears. The quotient
from a face to its primitive root preserves the zero-bracket relation but
forgets multiplicities; `(m,n,rho)` and the defect face `K` are the necessary
sidecars.

Candidate tournament vertices were faces, irreducible factors, defects,
shears, primitive roots, and proof obligations. Orienting defects by decreasing
weight gives a tie Hamiltonian path and a useful scheduler, but its scores,
cycles, SCCs, edge flips, and path counts do not encode the resonance
`d|rho`. The faithful carrier is the ordered defect filtration

```text
(h;m,n; rho; A,B; K; divisibility d|rho),              (31)
```

which preserves the Keller predicate through (22) and records exactly why the
power-free descent terminates or stalls.
