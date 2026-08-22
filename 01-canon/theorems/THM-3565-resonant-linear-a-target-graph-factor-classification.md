---
id: THM-3565
title: "Resonant linear-a target-graph factor classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For phi in
  C[a,b] with deg_a(phi)=1, the
  pullback of c+phi=0 under the fixed THM-1300 Keller map is reducible if
  and only if phi=-2h(b)^3a+b h(b)^2-2h(b) for a nonzero h in C[b].  Its
  generic core cubic has exactly one linear and one irreducible quadratic
  factor, and the source pullback has the explicit factor
  x-(1+xy)h(F2).  Consequently every genuinely quadratic target graph has
  irreducible pullback.  Among collision-compatible affine graphs, h=+-2
  recovers exactly THM-3559's two Kummer rows.
source: kps-s188
audit: >
  An independent assumption-free audit reproduced the rational-root
  dichotomy, saturated quadratic-denominator Groebner basis, Laurent hostile,
  polynomial descent, residual nonsquare gate, and generic-to-geometric
  component typing.  It also prompted the field-correct branch wording and
  removal of unnecessary nonzero assumptions on r,v,w in the companion.
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3564-target-graph-trisection-degree-resonance-irreducibility
related:
  - THM-3559-affine-target-coordinate-pullback-no-go
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
companion: 04-computation/jacobian_resonant_linear_a_factor_classification_kps_s188.py
output: 05-knowledge/results/jacobian_resonant_linear_a_factor_classification_kps_s188.out
script_sha256: f15ac70a58883fa27a91adbe166137a09ba00f4fe996a61ae5e6119a109e59ae
output_sha256: 2b679c8f2bbeb9c060a31de77ad64756c12214118b47d10e2590706bdedb643e
hash_basis: LF-normalized bytes
---

# THM-3565 -- resonant linear-a target-graph factor classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The first
resonance left by THM-3564 is not an unstructured exceptional cell.  It is
one explicit one-function family.

All fields and varieties are over `C`.  Target coordinates are `(a,b,c)`,
source coordinates are `(x,y,z)`, and `F=(F1,F2,F3)` is THM-1300's fixed
etale, quasi-finite degree-three map.

## 1. Statement

Let

```text
phi=A(b)a+B(b) in C[a,b],                 A!=0.             (1)
```

Then the following are equivalent:

1. the core cubic of `F` over the target graph `c+phi(a,b)=0` is reducible
   over `C(a,b)`;
2. the pullback hypersurface `V(F3+phi(F1,F2))` is reducible;
3. there is a nonzero polynomial `h in C[b]` such that

```text
phi=-2h^3 a+b h^2-2h.                                      (2)
```

For `(2)`, put

```text
D=3ah^2-bh+1,
C=12a^2h^2-4abh+16a-b^2.                                  (3)
```

The exact generic factorization is

```text
E_phi(x)
 =(Dx-h)
  [DCx^2+hCx+2(2ah^2-bh+2)],                              (4)
```

and the quadratic factor in `(4)` is irreducible.  On source space one
factor is especially simple:

```text
R_h=x-(1+xy)h(F2).                                         (5)
```

Thus the resonance produces exactly two irreducible pullback components,
not an arbitrary factor atlas.

## 2. Rational-root dichotomy

Put `K=C(b)` and write

```text
phi=A(a-r),                         A in K*, r in K.        (6)
```

THM-2473's graph-specialized core cubic is

```text
E_phi(x)=L_phi x^3+(4+3bphi)x+2phi,

L_phi=27a^2phi^2+18abphi+16a-b^3phi-b^2.                  (7)
```

It is primitive as a polynomial in `K[a][x]`, because

```text
gcd(phi,4+3bphi)=1.                                        (8)
```

If `(7)` is reducible over `K(a)`, it has a rational root `P/Q` in lowest
terms.  The polynomial rational-root theorem gives

```text
P divides 2phi,                    Q divides L_phi.         (9)
```

THM-3564's infinity valuation forces `v_infinity(P/Q)=1`.  Since `phi` is
linear and irreducible in `K[a]`, after scaling there are only two cases:

```text
I.   x=u/(a-v),

II.  x=u(a-r)/(a^2+va+w),          gcd(a-r,a^2+va+w)=1,    (10)
```

with `u in K*`.  This is exhaustive; the second denominator and its
coprimality are the branch that a leading-order argument alone would miss.

## 3. The constant-numerator branch

Substitute case I into `(7)` and clear `(a-v)^3`.  The leading coefficient,
then the next coefficient, give

```text
A=-2/(27u^3),                 v=bu/2+r/3.                  (11)
```

The coefficient of `a^2` then factors as the product

```text
R S=0,

R=-3bu+2r+18u^2,             S=3bu+2r-18u^2.              (12)
```

The `R=0` branch gives

```text
r=3bu/2-9u^2,                v=bu-3u^2.                   (13)
```

There is no second branch hidden in `S=0`.  On substituting `S=0`, the
remaining two coefficients, up to nonzero scalars, are

```text
(b-6u)^2(b+12u),             bu(b-6u)^2(b+6u).            (14)
```

In the field `K=C(b)`, the first equation in `(14)` forces `b=6u` or
`b=-12u`, while the second forces `b=6u`, `b=-6u`, or `b=0`.  Since `u!=0`
and `b` is nonzero in `K`, their only common possibility is `b=6u`; there
equation `R=0` also holds.  Thus `(13)` is the only solution.  It gives

```text
phi=-2a/(27u^3)+(b-6u)/(9u^2),

x=u/(a-bu+3u^2).                                           (15)
```

Direct substitution verifies that the displayed `x` is a root.

## 4. The coprime quadratic-denominator branch

Substitution of case II and successive leading-coefficient elimination give
again `A=-2/(27u^3)` and

```text
v=-bu/2-4r/3,

w=(9b^2u^2+18bru-108bu^3+8r^2+324u^4)/36.                (16)
```

Every one of the five remaining coefficients contains the factor `R` from
`(12)`.  If `R=0`, the denominator in case II factors exactly as

```text
a^2+va+w=(a-r)(a-bu+3u^2),                                (17)
```

contradicting the coprimality in `(10)` and merely reducing to case I.

Suppose `R!=0`.  Scale

```text
p=b/u,                         q=r/u^2.                    (18)
```

After dividing the five equations by `R`, their exact reduced Groebner basis
over `Q` is

```text
q-p^2/8+3p/4,

(p-12)(p-6)(p+6).                                         (19)
```

The three rows `(p,q,R/u^2)` are

```text
(12,9,0),                    (6,0,0),                    (-6,9,54).   (20)
```

The first two return to the cancelled branch `R=0`.  The only genuinely
coprime row has `b=-6u`, `r=9u^2`; hence

```text
phi=16a/b^3-4/b.                                           (21)
```

This is a rational target graph but not an element of `C[a,b]`.  Therefore
case II contributes no polynomial graph.  The exact Groebner basis `(19)`
is the finite algebraic certificate for the only elimination step in the
proof.

## 5. Polynomial descent gives one h-family

It remains to impose that both coefficients in `(15)` are polynomials in
`b`.  Write `u=p_0/q_0` in lowest terms in `C(b)`.  Since

```text
-2/(27u^3)=-2q_0^3/(27p_0^3) in C[b],                    (22)
```

coprimality forces `p_0` to be constant.  Put

```text
h=q_0/(3p_0) in C[b] minus {0}.                            (23)
```

Then `u=1/(3h)`, and `(15)` becomes exactly `(2)`.  Conversely, substituting
`(2)` into `(7)` gives the factorization `(4)`.  This proves the equivalence
of conditions 1 and 3.

The discriminant of the quadratic factor in `(4)` is

```text
-(6ah^2-3bh+4)^2 C.                                       (24)
```

The quadratic `C` in `a` has discriminant

```text
64(b^2h^2-2bh+4).                                         (25)
```

If `C` were a square in `C(b)[a]`, `(25)` would vanish.  Then `bh` would be
one of the two constant roots of `t^2-2t+4`, forcing `h=constant/b`, which
is not polynomial.  Hence `C` is nonsquare, `(24)` is nonsquare, and the
residual quadratic in `(4)` is irreducible.

Quasi-finiteness forces every pullback component to dominate the target
graph.  Therefore generic factorization is equivalent to geometric
factorization, proving conditions 1 iff 2 and showing that there are exactly
two components.

## 6. The source factor and collision boundary

There is an identity before choosing `h`.  If `H` is an indeterminate and
`u_src=1+xy`, direct expansion gives

```text
F3-2H^3F1+H^2F2-2H=(x-u_src H) Q_H                       (26)
```

for an explicit polynomial `Q_H in C[x,y,z,H]`.  Substitution
`H=h(F2)` proves `(5)` without denominators.

Let `q=(-1/4,0,0)` be the triple collision value.  The target graph contains
`q` precisely when

```text
phi(-1/4,0)=h(0)(h(0)-2)(h(0)+2)/2=0.                    (27)
```

Thus `h(0)` is `0`, `2`, or `-2`.  At the three collision sources, `(5)`
contains respectively only

```text
h(0)=0:  p0,                 h(0)=-2: p+,                h(0)=2: p-.  (28)
```

The linear component never carries a collision pair; the irreducible
quadratic component carries the other two.

If `n=deg h`, then

```text
deg_total(phi)=3n+1.                                         (29)
```

Consequently no genuinely quadratic polynomial `phi` belongs to `(2)`.
Combining with THM-3564 for `deg_a(phi)=0` or `2` proves:

> Every target graph with `deg_total(phi)=2` has irreducible pullback.

Among degree-at-most-two graphs through `q`, the only reducible rows are
`phi=0` and, from constant `h=+-2`,

```text
phi=-16a+4b-4,                  phi=16a+4b+4.              (30)
```

These are exactly THM-3559's pure-`F3` row and two affine Kummer exceptions.
Thus the earlier finite `F_7` observation upgrades to a characteristic-zero
factor theorem.  It does not yet rule out the possibility that an
irreducible quadratic complete pullback is itself a coordinate `A2`; that
requires the Euler/normalization gate, not factor classification.

## 7. Exact verification

Run

```bash
python3 04-computation/jacobian_resonant_linear_a_factor_classification_kps_s188.py
python3 -O 04-computation/jacobian_resonant_linear_a_factor_classification_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion keeps every
coefficient equation active without assuming `r`, `v`, or `w` nonzero,
recomputes `(19)` over `Q`, verifies the exceptional rational row `(21)`,
verifies `(4)`, `(24)--(27)`, and checks both polynomiality and the identity
of the source cofactor before substituting `H=h(F2)`.  An independent
assumption-free derivation reproduced every structural branch.

**QED.**
