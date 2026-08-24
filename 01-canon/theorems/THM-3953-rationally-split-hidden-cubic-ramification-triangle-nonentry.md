---
id: THM-3953
title: "Rationally split hidden cubic ramification forms a forbidden boundary triangle"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER MISTAKE-474
  QUANTIFIER REPAIR.
  Let the hidden ramification cubic of a normal depressed-cubic surface split
  into three distinct polynomial graph roots. The missing linear h-row
  forces the exact parametrization r0=c a(a+b), r1=c b(a+b), r2=-cab.
  If a/b is nonconstant, the three primitive pair-collision polynomials are
  nonconstant and pairwise coprime. They give three distinct source points
  joining the three ramification primes in a triangle. For each pair, its two
  formal primitive factors are respectively zero/smooth and nonzero
  `2:-1`/singular when they vanish; a factor may be a unit. The natural cubic
  is a domain with only finitely many singular points, hence already normal;
  normalization cannot separate the singular incidences. A same-function-field Keller open
  would delete all three ramification primes, contradicting THM-3951's
  boundary tree. Constant-ratio and globally repeated-root packets are routed
  separately and are not overclaimed.
source: jc-degree6-one-place / post-THM-3951 split-residual audit, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS AFTER MISTAKE-474 REPAIR
  (audit_triangle_3953, 2026-08-24). The companion verifies the root
  parametrization and converse, all pair differences and their pairwise
  Bezout-resultant controls, the six-row smooth/singular collision table,
  the monic-cubic domain obstruction, and explicit positive/hostile controls.
  The audit reconstructed the UFD parametrization, pairwise coprimality,
  monic-cubic irreducibility, finite-singular-locus R1+S2 normality,
  arbitrary-point common-resolution paths, Zariski Main and the ramification-
  boundary bridge. It found only the frontmatter's false 3/3 population
  overread; the repaired type-wise statement and MISTAKE-474 preserve the
  triangle proof. Normal and optimized runs byte-match the frozen 59-gate
  output, and all hashes agree.
depends_on:
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
related:
  - THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow
script: 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
output: 05-knowledge/results/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.out
script_sha256: f90897ea150c75ade3480e499eb686bdce6bccc31c07bb4cd4325f0b8245848a
output_sha256: 7bbd1fd9c2f4c6efeb68dca40c973caa523b1dc39262c96c003710b64a83b624
semantic_sha256: 507b6809acb59a52a8255221897706e8adfd7fa5359e98e96acfcf45c9975c05
hash_basis: raw LF bytes
---

# THM-3953 -- three polynomial ramification roots form a forbidden triangle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER MISTAKE-474
QUANTIFIER REPAIR.**
Work over an algebraically closed field `k` of characteristic zero. This is
the rationally split complement to THM-3950--3951's irreducible
equianharmonic residual.

## 1. Universe and exact root parametrization

Let `r0,r1,r2 in k[t]` be pairwise distinct polynomial functions satisfying

```text
r0 r1+r0 r2+r1 r2=0.                                      (1)
```

Assume that at least one root ratio is nonconstant. Define

```text
C=2(r0+r1+r2),                  E=2r0r1r2,
G(h,t)=E+C h^2-2h^3=-2(h-r0)(h-r1)(h-r2).                (2)
```

The identity `(1)` is exactly the missing linear `h`-row in `(2)`.

All solutions have the following UFD parametrization. Write

```text
r0=d a,                    r1=d b,                    gcd(a,b)=1.      (3)
```

Then `(1)` becomes

```text
d ab+(a+b)r2=0.                                             (4)
```

The relation `a+b=0` would force `ab=0`, contrary to the three roots being
distinct. Moreover `gcd(a+b,ab)=1`, so `(a+b)|d`. Consequently there is a
nonzero `c in k[t]` such that

```text
r0=c a(a+b),             r1=c b(a+b),             r2=-c ab.          (5)
```

Conversely `(5)` satisfies `(1)`. The pairwise differences are

```text
r0-r1=c D01,        D01=(a-b)(a+b),
r0-r2=c D02,        D02=a(a+2b),
r1-r2=c D12,        D12=b(2a+b).                         (6)
```

Pairwise distinctness says none of these polynomials is zero. Because
`a/b` is nonconstant, every `Dij` is nonconstant: if one of the displayed
products were a nonzero scalar, both of its factors, hence `a` and `b`,
would be scalar. Also

```text
gcd(D01,D02)=gcd(D01,D12)=gcd(D02,D12)=1.                (7)
```

Indeed a common irreducible factor in any pair, combined with the relevant
linear forms in `(6)`, would divide both `a` and `b`; the only coefficients
introduced are `2` and `3`, which are units in characteristic zero.

Thus one may choose roots

```text
t01 in V(D01),             t02 in V(D02),
t12 in V(D12),                                                 (8)
```

and `(7)` makes the three parameter values distinct. At `tij`, the graph
roots `ri` and `rj` coincide. A zero of the common factor `c` may make the
third root coincide as well, but it neither removes the chosen incidence nor
identifies the three distinct parameter values.

## 2. The natural cubic is an integral normal surface

Put

```text
F=T^3-3PT-(E+CP),
X0=Spec k[P,t,T]/(F),                 pi:X0 -> A2_(P,t).   (9)
```

First, `X0` is integral. If the monic cubic `F` were reducible, it would have
a root in `k[P,t]`: a root over the fraction field is integral and the UFD
`k[P,t]` is normal. Comparison of `P`-degrees makes that root independent of
`P`; the two coefficient rows then force

```text
root=-C/3,                         C^3+27E=0.             (10)
```

Substitution of `(5)` gives the sharper factorization

```text
C^3+27E=
 2c^3(a-b)^2(a+2b)^2(2a+b)^2.                           (11)
```

None of the three linear factors can vanish identically without making a
pair of roots in `(5)` identical. Thus `(11)` is nonzero and contradicts
`(10)`. Hence `F` is irreducible. Notice that integrality only needs three
distinct roots; the nonconstant-ratio hypothesis enters later, when each
pair-collision carrier must have a finite zero.

The singular locus of `X0` is finite. Indeed

```text
F_T=3(T^2-P).                                             (12)
```

At a ramification point put `h=-T`, so `P=h^2` and `h` is one of the roots
`ri`. Differentiating `(2)` in `h` gives

```text
G_h=2h(C-3h)=-2h F_P.                                    (13)
```

If the three roots are distinct at `t`, then their simple root `h` cannot be
zero: `(1)` would otherwise force another root to vanish. Thus `G_h!=0`
implies `F_P!=0`, and `X0` is smooth above that parameter. Every singular
point therefore lies above a zero of

```text
c D01 D02 D12,                                            (14)
```

a finite set; each fibre supplies only finitely many candidates from `(12)`
and `(2)`. The hypersurface domain `X0` is `S2`, and a finite singular locus
has codimension two. Serre's criterion proves

```text
X0 is normal.                                             (15)
```

This normality statement is essential: the singular pair collisions below
are not separated by passing to a maximal overorder.

## 3. Exact smooth-versus-singular collision table

The ramification divisor contains the three distinct graph primes

```text
Ei:       P=ri(t)^2,                 T=-ri(t),       i=0,1,2.         (16)
```

Their generic points are smooth on `X0` and ramified for `pi` by
`(12)-(13)`. At a primitive pair collision where `c!=0`, the six factors in
`(6)` split into the following exact table.

| pair | collision factor | repeated value | third value | surface |
|---|---|---:|---:|---|
| `E0,E1` | `a+b=0` | `0` | `c a^2` | smooth |
| `E0,E1` | `a-b=0` | `2c a^2` | `-c a^2` | singular |
| `E0,E2` | `a=0` | `0` | `c b^2` | smooth |
| `E0,E2` | `a+2b=0` | `2c b^2` | `-c b^2` | singular |
| `E1,E2` | `b=0` | `0` | `c a^2` | smooth |
| `E1,E2` | `2a+b=0` | `2c a^2` | `-c a^2` | singular |

This is a table of formal factor types, not an equal-population statement.
A displayed factor can be a unit and then contributes no collision point.

For two equal roots `x,x` and third root `y`, equation `(1)` reads

```text
x(x+2y)=0,                       F_P=-(x+2y).             (17)
```

This proves the last column. At a nonzero collision `x+2y=0`, the product
form `(2)` also gives `G_t=0`, hence `F_t=0`; together with `(12)` this is a
singular point. At a zero collision, `F_P=-2y!=0`. If `c=0`, all three roots
are zero and the collision is singular; this is the common-gcd override to
the table. All such singularities are among the finite set already covered
by `(14)-(15)`.

## 4. The three boundary primes make a cycle

Suppose the same function field admitted source coordinates `x,z` with

```text
k(x,z)=Frac(X0),              P,t in k[x,z],
Jac_(x,z)(P,t) in k*.                                      (18)
```

The Keller map `A2_(x,z) -> A2_(P,t)` is etale and quasi-finite. Since `X0`
is already the finite normalization of the target in its function field,
Zariski Main gives an open immersion

```text
A2_(x,z)=U  -->  X0.                                      (19)
```

Every `Ei` is generically a ramification prime for `pi`, whereas `(19)` is
etale. Hence all three primes lie in `X0 minus U`.

We use the arbitrary-point form of THM-3951's boundary-incidence forest.
On a common resolution of a normal completion of `X0` and `(P2,L_infinity)`,
the SNC boundary is a tree. For each point in `(8)`, connectedness of the
proper birational fibre joins the strict transforms of the corresponding
two boundary primes by an exceptional path; smoothness of the point is not
needed. The three fibres are disjoint because the three parameter values are
distinct. Thus the paths

```text
E0tilde -- E1tilde,          E0tilde -- E2tilde,
E1tilde -- E2tilde                                           (20)
```

are internally disjoint and form a subdivision of a triangle in the
boundary dual multigraph. This contradicts the tree property. Therefore no
same-function-field planar Keller chart `(18)` exists.

## 5. Equality and failure boundaries

The hypotheses above have three exact boundaries.

1. **Constant root ratios.** In `(5)`, `a,b in k*` gives
   `ri=lambda_i c(t)`. If `c` has at least two distinct zeros, any pair of
   graph primes already meets twice and the same connected-fibre forest
   argument applies; `(11)` and the finite-collision proof still make `X0`
   a normal domain. If `c` is a unit or has only one distinct zero, the
   present triangle construction is absent and these scalar packets require
   a separate atlas analysis.
2. **Duplicate roots.** If, say, `r0=r1=x`, then `(1)` gives
   `x(x+2r2)=0`. Thus either the repeated graph is identically zero or the
   third root is `-x/2`. This is a globally repeated ramification component,
   not three distinct boundary primes, and is outside the theorem. If all
   three roots coincide, `(1)` forces all of them to vanish and `(9)` is
   reducible.
3. **Rational rather than polynomial roots.** Clearing denominators can add
   vertical/infinity primes and can destroy the affine collision count.
   Nothing here treats that case.

Accordingly the theorem closes precisely the three-distinct-polynomial-root,
nonconstant-ratio split residual. Together with THM-3950--3951, it explains
both sides of the first hidden-cubic dichotomy:

```text
irreducible residual  -> fixed positive-genus shadow + two-prime cycle;
polynomially split    -> three graph primes + boundary triangle.       (21)
```

It does not construct or exclude arbitrary cubic fields, extra factor
distributions, rational-root denominators, or the scalar equality packets.
`JC(2)` remains open.

Run

```bash
python3 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
python3 -O 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
```

for the assertion-only exact companion.
