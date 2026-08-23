---
id: THM-3784
title: "Rational Keller tower different, codifferent, norm, and trace duality"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the
  THM-3774 degree-(m+1) tower, the reciprocal boundary unit is exactly
  g=-f_m'(t)/m and the source coordinate x=1/g is a generator of the trace
  codifferent.  The norm of g is, up to the proved sign and scalar, the
  irreducible branch discriminant, cuspidal exactly when m>=2.  The
  polynomial codifferent ladder
  x,tx,...,t^m x has trace packet (0,...,0,-m) and a constant unit pairing
  determinant.  This unifies the axis pole, branch curve, and sheetwise
  trace cancellation.  The trace-zero repair lane remains open; no
  polynomial Keller pair or JC(2) counterexample is claimed.
source: jc_zero_debt_lift / different-codifferent reframe, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT by jc_sparse_direct_search.  The derivative and
  reciprocal-Jacobian signs, Lagrange trace packet, recurrence, full
  trace-codifferent module equality, anti-triangular determinant sign and
  unit claim, discriminant/norm signs, and the m=1 and m=2 boundary cases
  were rederived independently.  Resultant and multiplication-matrix norm
  calculations agree for m=1,...,10; direct source typing is checked for
  m=1,...,8.  Normal, optimized, and frozen runs agree byte for byte, and
  the recorded hashes match.  The audit repaired the scope wording: the
  irreducible branch curve is smooth for m=1 and cuspidal for m>=2.
depends_on:
  - THM-3774-three-component-rational-keller-cover-tower
  - THM-3780-rational-keller-tower-all-affine-plane-fillings-obstruction
related:
  - THM-3779-three-component-tower-maximal-danielewski-polynomial-observable
  - THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable
script: 04-computation/jc2_rational_keller_tower_different_codifferent_trace_thm3784.py
output: 05-knowledge/results/jc2_rational_keller_tower_different_codifferent_trace_thm3784.out
script_sha256: 8aaf279cbf32ac22dc7cf6c80ea46af99507b654bf4930cd3c548b3432e597b8
output_sha256: a67027da9c7d7b5a53ab195804ce34aea01323ca2b922af9829671b5ac33cca6
semantic_sha256: f3c7450cf9e2733919394572d1b18057ee976f6c2f6818203cdd8a6f73d61f7d
hash_basis: raw LF bytes
---

# THM-3784 -- the axis pole is the inverse different

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Three
features of the rational Keller tower that initially look unrelated are one
algebraic object:

```text
axis pole x,
branch discriminant (cuspidal for m>=2),
trace cancellation across the m+1 sheets.                              (1)
```

The common object is the different of the finite normalization.  Its
generator is the reciprocal Jacobian factor `g=1/x`; the source coordinate
`x` generates the codifferent.  The branch equation is its norm, and the
trace-zero polynomial ladder is the coefficient-dual power basis.

Let `k` be an algebraically closed field of characteristic zero, fix an
integer `m>=1`, and put

```text
B_0=k[U,P],
f_m(T)=T^(m+1)-mPT+mU,
S=B_0[t]=B_0[T]/(f_m),
L=Frac(S).                                               (2)
```

THM-3774 proves that `f_m` is irreducible and that `L=k(x,y)` for the tower.
THM-3780 identifies

```text
t=x^2(1+xy),
g=P-[(m+1)/m]t^m=1/x.                                  (3)
```

The new assertions are

```text
g=-f_m'(t)/m,
x=1/g=-m/f_m'(t),                                      (4)

different_(S/B_0)=gS,
codifferent_(S/B_0)=xS.                                (5)
```

Define the branch binomial

```text
D_m=(m+1)^(m+1)U^m-m^(m+1)P^(m+1).                    (6)
```

Then the exact field norms are

```text
N_(L/k(U,P))(g)=(-1)^(m+1) D_m/m,
N_(L/k(U,P))(x)=(-1)^(m+1) m/D_m.                      (7)
```

Finally, with `Tr=Tr_(L/k(U,P))`,

```text
Tr(t^r x)=0                 for 0<=r<m,
Tr(t^m x)=-m.                                             (8)
```

Equivalently, for every polynomial `q(T)` of degree at most `m`,

```text
Tr(x q(t))=-m lc(q).                                    (9)
```

The `B_0`-basis

```text
x, tx, ..., t^m x                                      (10)
```

of the codifferent has trace-pairing determinant

```text
(-1)^((m+1)(m+2)/2) m^(m+1) in k*.                    (11)
```

## 1. The derivative is the reciprocal boundary unit

Differentiating the monic cover equation in `(2)` gives

```text
f_m'(T)=(m+1)T^m-mP.                                   (12)
```

Substitution of `T=t` and comparison with `(3)` immediately yield `(4)`.
Thus the factorization from THM-3780 can be written

```text
J_(x,y)(t,P)=x=-m/f_m'(t),
J_(t,P)(U,P)=g=-f_m'(t)/m,
J_(x,y)(U,P)=xg=1.                                     (13)
```

The two nonconstant Jacobians are not merely reciprocal functions.  They
are the different and inverse-different factors of the same finite cover.

The ring `S` is finite free of rank `m+1` over `B_0`, with power basis

```text
1,t,...,t^m.                                           (14)
```

It is the normalization of `B_0` by THM-3780.  Since characteristic zero
makes the extension separable, the trace pairing on `L/k(U,P)` is
nondegenerate.

## 2. Lagrange interpolation gives the entire first trace packet

Work over a splitting field and write `t_1,...,t_(m+1)` for the distinct
roots of `f_m`.  The standard partial-fraction identity is

```text
1/f_m(Z)=sum_i 1/[f_m'(t_i)(Z-t_i)].                   (15)
```

Expanding at `Z=infinity`, the right side is

```text
sum_(r>=0) [sum_i t_i^r/f_m'(t_i)] Z^(-r-1).           (16)
```

Because `f_m` is monic of degree `m+1`, the left side begins with
`Z^(-m-1)`.  Comparing the first `m+1` coefficients gives

```text
sum_i t_i^r/f_m'(t_i)=0       for 0<=r<m,
sum_i t_i^m/f_m'(t_i)=1.                              (17)
```

Multiplication by `-m` and `(4)` prove `(8)`.  Linearity proves `(9)`.
Thus `-x/m` is the trace functional that extracts the leading coefficient
in the power basis `(14)`.

For later powers, multiplying the relation `f_m(t)=0` by
`t^(r-m-1)x` gives the exact recurrence

```text
Tr(t^r x)
 =mP Tr(t^(r-m)x)-mU Tr(t^(r-m-1)x),       r>=m+1.     (18)
```

In particular every trace pairing among the elements in `(10)` and `(14)`
belongs to `B_0`; there is no hidden denominator beyond the codifferent
generator itself.

## 3. The trace packet proves the codifferent equality

Let

```text
S^vee={z in L : Tr(zS) subset B_0}.                    (19)
```

Equation `(18)` shows `xS subset S^vee`.  Pair the basis `(10)` against
`(14)`.  Its matrix is

```text
M_(ij)=Tr(t^(i+j)x),                 0<=i,j<=m.         (20)
```

By `(8)`, every entry strictly before the anti-diagonal is zero and every
anti-diagonal entry is `-m`.  Entries after the anti-diagonal are governed
by `(18)` but do not affect the determinant.  Reversing `m+1` columns gives

```text
det M
 =(-1)^((m+1)m/2)(-m)^(m+1)
 =(-1)^((m+1)(m+2)/2)m^(m+1),                         (21)
```

which is a unit of `B_0`.  Therefore the trace pairing identifies `xS`
with the full dual module `Hom_(B_0)(S,B_0)`.  This proves

```text
S^vee=xS.                                              (22)
```

Because `x=-m/f_m'(t)` and `m` is a unit of `k`, equation `(22)` is the
standard monogenic formula

```text
S^vee=f_m'(t)^(-1)S.                                   (23)
```

Its inverse fractional ideal is the trace different

```text
(S^vee)^(-1)=f_m'(t)S=gS,                              (24)
```

proving `(5)`.  In this monogenic hypersurface, `(24)` is also the Kähler
different: `Omega_(S/B_0)` has presentation by the single relation
`f_m'(t)dt=0`, so its zeroth Fitting ideal is `(f_m'(t))=(g)`.

## 4. The norm of the different is the branch curve

Put `n=m+1`.  Since `f_m` is monic,

```text
N(f_m'(t))=Res_T(f_m,f_m'),
disc(f_m)=(-1)^(n(n-1)/2) Res_T(f_m,f_m').             (25)
```

Using `g=-f_m'(t)/m` gives

```text
N(g)=(-1)^n m^(-n)(-1)^(n(n-1)/2) disc(f_m).           (26)
```

THM-3774's discriminant formula is

```text
disc(f_m)=(-1)^(m(m+1)/2)m^m D_m.                     (27)
```

The total sign exponent in `(26),(27)` is

```text
n(n+1)/2+m(m+1)/2=(m+1)^2 == m+1 mod 2,               (28)
```

and the power of `m` is `m-(m+1)=-1`.  This proves the
first formula in `(7)`; the second follows from `x=g^(-1)`.

Thus the branch divisor is literally the norm divisor of the different.
At the quadratic-cover boundary `m=1`, it is the smooth parabola

```text
D_1=4U-P^2,
g=P-2t,
N(g)=4U-P^2,
(Tr(x),Tr(tx))=(0,-1).                                (29)
```

For the cubic centerpiece `m=2`, the packet reads

```text
f_2=t^3-2Pt+2U,
g=P-(3/2)t^2=-f_2'(t)/2,
N(g)=(8P^3-27U^2)/2,
(Tr(x),Tr(tx),Tr(t^2x))=(0,0,-2).                     (30)
```

The same `A_2` cusp that detects non-Galois branching therefore measures the
norm of the reciprocal factor whose inverse is the source axis coordinate.

## 5. What the trace-zero ladder preserves and forgets

Every member of

```text
x, tx, ..., t^(m-1)x                                  (30)
```

is an honest polynomial in the original source coordinates because both
`x` and `t=x^2(1+xy)` are polynomial.  Equations `(8)` say that all `m`
members have zero field trace, while the final rung `t^m x` is the unique
first coefficient detector.  This supplies a new sidecar for constructions
that leave the rational target field `k(U,P)`:

```text
source object     polynomial codifferent rung t^r x;
target shadow     zero trace for r<m;
preserved         sheetwise signed cancellation and source polynomiality;
forgotten         individual sheet, vertical component, and Poisson bracket;
needed sidecar    component valuations plus pairwise brackets;
cheap next test   mix trace-zero rungs with the THM-3779 observables.     (31)
```

Trace zero alone does not equalize the component spectrum and does not imply
a Darboux relation.  It is instead the first exact mechanism in this tower
that cancels across conjugate sheets while using polynomial functions outside
the base target field.  Whether such rungs can be combined with non-base
observables to repair the final axis address is **OPEN**.

The exact companion verifies the discriminant and norm independently by
resultants and multiplication matrices through `m=10`, the full trace packet
and recurrence, the unit determinant `(11)`, reciprocal norms, and direct
source typing.  No polynomial Keller pair and no planar Jacobian
counterexample is constructed.  **QED.**
