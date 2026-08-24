---
id: THM-3953
title: "Rationally split hidden cubic ramification forms a forbidden boundary triangle"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Let the hidden ramification cubic of a normal depressed-cubic surface split
  into three distinct polynomial graph roots. The missing linear h-row
  forces the exact parametrization r0=c a(a+b), r1=c b(a+b), r2=-cab.
  If a/b is nonconstant, the three primitive pair-collision polynomials are
  nonconstant and pairwise coprime. They give three distinct source points
  joining the three ramification primes in a triangle. Half of the primitive
  collisions are smooth and half singular, but the natural cubic is a domain
  with only finitely many singular points, hence already normal; normalization
  cannot separate the singular incidences. If a/b is constant, an invertible
  polynomial change identifies the surface with xy=c(t)^3. A unit c gives a
  nonconstant unit, while a nonunit c gives an exact Nagata class group with
  nonzero 3-torsion, contradicting THM-3922. Thus every
  three-distinct-polynomial-root split packet fails a same-field Keller
  atlas; only duplicate roots and rational denominators remain outside.
source: jc-degree6-one-place / post-THM-3951 split-residual audit, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-zero-debt-lift and jc-cohn3709,
  2026-08-24). Both audits reconstructed the UFD root parametrization,
  pair-carrier coprimality, monic-cubic domain and finite-singular-locus
  normality proofs, every collision row, and the arbitrary-point boundary
  triangle. They also checked the strengthened constant-ratio coordinate
  isomorphism to xy=c(t)^3, its unit case, all Nagata valuations and units,
  the exact class-group Smith quotient, and the THM-3922 contradiction.
  Normal and optimized runs byte-match the frozen 75-gate output, all hashes
  agree, documentation checks pass, and no repair was required.
depends_on:
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow
script: 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
output: 05-knowledge/results/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.out
script_sha256: 07c8fd4243f6696e04357849e4aebd09f7b45fe26f6c2e6b9d7c8cc392345ab3
output_sha256: 234cf1fb24482025a245b7737468a9de95db7edd75478e864478100d0915adf5
semantic_sha256: 593a316f92f4c3148a732de8d4c71c97eb09572c395905143db3334f42504559
hash_basis: raw LF bytes
---

# THM-3953 -- three polynomial ramification roots form a forbidden triangle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. This is
the rationally split complement to THM-3950--3951's irreducible
equianharmonic residual.

## 1. Universe and exact root parametrization

Let `r0,r1,r2 in k[t]` be pairwise distinct polynomial functions satisfying

```text
r0 r1+r0 r2+r1 r2=0.                                      (1)
```

Define

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

Pairwise distinctness says none of these polynomials is zero. If `a/b` is
nonconstant, every `Dij` is nonconstant: if one of the displayed products
were a nonzero scalar, both of its factors, hence `a` and `b`, would be
scalar. Also

```text
gcd(D01,D02)=gcd(D01,D12)=gcd(D02,D12)=1.                (7)
```

Indeed a common irreducible factor in any pair, combined with the relevant
linear forms in `(6)`, would divide both `a` and `b`; the only coefficients
introduced are `2` and `3`, which are units in characteristic zero.

In this nonconstant-ratio case one may therefore choose roots

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

## 4. Nonconstant ratios make a boundary triangle

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

## 5. Constant ratios reduce exactly to `xy=c(t)^3`

It remains to treat the case in which every root ratio is constant. Since
`gcd(a,b)=1`, this is exactly

```text
a,b in k*,
r0=lambda0 c(t),       r1=lambda1 c(t),       r2=lambda2 c(t).        (21)
```

The three constants are distinct. Put

```text
sigma1=a^2+ab+b^2,                 sigma3=-a^2b^2(a+b)^2,
K=8sigma1^3+54sigma3
 =2(a-b)^2(a+2b)^2(2a+b)^2 !=0.                         (22)
```

The following change of coordinates is polynomial and invertible:

```text
x=3T+2sigma1 c,
Y=27P-x^2+6sigma1 x c-12sigma1^2 c^2,
y=-Y/K.                                                   (23)
```

Its inverse is

```text
T=(x-2sigma1 c)/3,
P=(-K y+x^2-6sigma1 x c+12sigma1^2 c^2)/27.              (24)
```

Direct substitution gives

```text
F=(K/27)(xy-c(t)^3).                                     (25)
```

Thus `X0` is exactly the normal Danielewski surface

```text
Sc:       xy=c(t)^3.                                     (26)
```

If `c` is a nonzero constant, `(26)` is `Gm times A1` and `x` is a
nonconstant unit. No dense open `A2` can lie in `Sc`, because restriction
would send that unit to a unit of `k[A2]`, hence to a scalar.

Now suppose `c` is nonconstant and factor it as

```text
c(t)=kappa product_j (t-alpha_j)^(m_j),
r>=1,                         m_j>=1.                    (27)
```

The surface `(26)` is a hypersurface domain whose singular locus is the
finite set `(x,y,t-alpha_j)=0`, so it is normal. Localizing at `x` gives

```text
Sc[x^(-1)]=k[x^(+-1),t],                                 (28)
```

a UFD with units `k* x^Z`. The height-one primes removed by this localization
are exactly

```text
Q_j=(x,t-alpha_j).                                       (29)
```

At the generic point of `Q_j`, `y` and all other factors are units, hence

```text
v_Qj(t-alpha_j)=1,              v_Qj(x)=3m_j,
div(x)=sum_j 3m_j Q_j.                                  (30)
```

A unit of `Sc` restricts in `(28)` to `lambda x^n`; `(30)` forces `n=0`,
so `Sc*=k*`. The Nagata localization sequence is therefore exact as

```text
0 -> Z --(1 |-> (3m_j))--> Z^r -> Cl(Sc) -> 0,
Cl(Sc)=Z^(r-1) direct_sum Z/(3 gcd_j m_j).                (31)
```

In particular `Cl(Sc)` has nonzero `3`-torsion. If a same-function-field
Keller source existed, Zariski Main would make it a dense `A2` open in this
normal finite degree-three completion. THM-3922 says the class group of every
proper such completion is freely generated by its prime divisorial boundary.
Equation `(31)` is impossible. This closes every constant-ratio packet.

## 6. Equality and failure boundaries

The hypotheses above leave two exact boundaries.

1. **Duplicate roots.** If, say, `r0=r1=x`, then `(1)` gives
   `x(x+2r2)=0`. Thus either the repeated graph is identically zero or the
   third root is `-x/2`. This is a globally repeated ramification component,
   not three distinct boundary primes, and is outside the theorem. If all
   three roots coincide, `(1)` forces all of them to vanish and `(9)` is
   reducible.
2. **Rational rather than polynomial roots.** Clearing denominators can add
   vertical/infinity primes and can destroy the affine collision count.
   Nothing here treats that case.

Accordingly the theorem closes every three-distinct-polynomial-root split
residual. Together with THM-3950--3951, it explains both sides of the first
hidden-cubic dichotomy:

```text
irreducible residual  -> fixed positive-genus shadow + two-prime cycle;
polynomially split    -> boundary triangle or torsion-class obstruction. (32)
```

It does not construct or exclude arbitrary cubic fields, extra factor
distributions, rational-root denominators, or globally repeated
ramification. `JC(2)` remains open.

Run

```bash
python3 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
python3 -O 04-computation/jc2_rationally_split_hidden_cubic_boundary_triangle_thm3953.py
```

for the assertion-only exact companion.
