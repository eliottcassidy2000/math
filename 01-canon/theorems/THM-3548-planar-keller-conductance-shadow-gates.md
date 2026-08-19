---
id: THM-3548
title: "Planar Keller conductance-shadow gates"
status: >
  PROVED / HOSTILE-AUDITED.  Four exact necessary conditions are collected
  under one typed conductance dictionary.  A normalized entrywise-intensity
  table of a large planar Keller differential is within |kappa|/T of rank
  one; coefficient channels close fibrewise polygons but must also satisfy
  global Segre rank-one equations; a polynomially removable Hermitian
  target-shear minimizer forces tameness; and a specified Puiseux escape
  branch obeys m <= r(d-2).  These are necessary filters, not a reduction of
  JC(2), and the Puiseux conclusion is conditional on the displayed branch
  expansion.
source: boxeph-2026-08-18-jacobian-dephasing
depends_on: []
related:
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
---

# THM-3548 -- planar Keller conductance-shadow gates

**PROVED / HOSTILE-AUDITED.**  Let

```text
F=(P,Q):C^2 -> C^2,              Jac(P,Q)=kappa in C*.  (1)
```

The squared magnitudes of derivative or coefficient channels behave like
conductances, but they forget the phases that impose the determinant and
coefficient cancellations.  The four statements below quantify what
survives and name the missing sidecars.

## 1. Differential plaquette: quantitative rank-one shadow

The following matrix statement applies pointwise to `A=DF(z)`.  Let

```text
A = [ a  b ],       C=(|A_ij|^2),       T=sum_ij C_ij,
    [ c  d ]

r=|ad|,             s=|bc|,             |det A|=|kappa|.
```

Then

```text
|r-s| <= |kappa| <= r+s,                 (2)
r^2+s^2 >= |kappa|^2/2,                  (3)
|det C| <= |kappa| T/2.                  (4)
```

If `M=C/T`, the Frobenius distance to matrices of rank at most one satisfies

```text
dist_F(M,{rank<=1}) = sigma_2(M) <= |kappa|/T.          (5)
```

### Proof

Write `alpha=ad`, `beta=bc`.  The triangle and reverse-triangle
inequalities for `kappa=alpha-beta` give `(2)`, and

```text
|kappa|^2 <= (r+s)^2 <= 2(r^2+s^2),
```

giving `(3)`.  Moreover

```text
|det C|=|r^2-s^2|=|r-s|(r+s)
       <= |kappa| (|a||d|+|b||c|)
       <= |kappa| T/2,
```

which is `(4)`.  The best Frobenius rank-one approximation error is
`sigma_2(M)`.  Since the nonnegative entries of `M` sum to one,

```text
(1,1) M (1,1)^T / 2 = 1/2,
```

so `sigma_1(M)>=1/2`.  Finally

```text
sigma_2(M)=|det M|/sigma_1(M)
          <= 2 |det C|/T^2 <= |kappa|/T.
```

This proves `(5)`.

If a sequence has `T->infinity`, pass to a subsequence.  There are two
exact regimes.

1. If `r+s->infinity`, then `r/s->1`.  With
   `Phi=arg(alpha overline(beta))`,

   ```text
   |kappa|^2=(r-s)^2+4rs sin^2(Phi/2),                 (6)
   ```

   hence `Phi->0` modulo `2 pi`.  Two large determinant matchings have
   nearly equal magnitudes and nearly equal phases, so their difference is
   the fixed unit: a **dark plaquette**.
2. If `r+s` remains bounded, some entry of `A` diverges while the entry
   opposite to it in the same determinant matching tends to zero.  This is
   the **channel-degenerate** or triangular-shadow regime.

Both are sharp.  The determinant-one matrices

```text
[ N  N ]                         [ N    0   ]
[ N  N+1/N ]                     [ 0   1/N  ]
```

realize the dark and channel-degenerate regimes respectively.

The phase loss is load-bearing.  The matrices

```text
A0 = [ 1  1 ],
     [ 1  1 ]

A1 = [       1                    1              ]
     [ exp(-2 pi i/3)  exp(-2 pi i/3)exp(pi i/3) ]
```

have the identical intensity table of four ones, whereas
`det A0=0` and `det A1=1`.  The missing datum is the four-edge Wilson phase
`arg(ad overline(bc))`.

## 2. Coefficient fibres: polygons plus a Segre sidecar

Write finite expansions

```text
P=sum_(a in A) p_a X^a,        Q=sum_(b in B) q_b X^b,
```

where `a,b in N^2`, and put

```text
z_ab=det(a,b) p_a q_b.                                 (7)
```

For `s in N^2`, the coefficient of `X^(s-(1,1))` in the
Jacobian is exactly

```text
sum_(a+b=s) z_ab.                                      (8)
```

Thus every nonconstant coefficient fibre is a closed complex polygon.  If
it has `m` nonzero channels, `c_e=|z_e|^2`, and `E=sum_e c_e`, then

```text
max_e sqrt(c_e) <= sum_(f!=e) sqrt(c_f),                (9)
c_max <= (m-1)E/m.                                     (10)
```

For two channels, the conductances are equal and the amplitudes are
opposite.  A singleton nonconstant fibre is impossible.  Equality in `(10)`
holds exactly when the other `m-1` amplitudes have equal magnitude, are
aligned, and oppose the maximal amplitude.  The constant fibre obeys the
same rule after adjoining the fixed amplitude `-kappa`.

Indeed, closure gives `|z_max|<=sum_rest |z|`; Cauchy--Schwarz gives

```text
c_max <= (m-1)(E-c_max),
```

which is `(10)` together with its equality statement.

Fibrewise closure is not sufficient.  Define the full coefficient matrix

```text
W_ab=p_a q_b,
T(W)_s=sum_(a+b=s) det(a,b)W_ab.                        (11)
```

Then

```text
Jac(P,Q)=kappa
iff T(W)=kappa e_(1,1) and W lies in the Segre rank-one locus.  (12)
```

In particular all `2 x 2` minors of `W` vanish; equivalently, every even
cycle in the coefficient bipartite graph obeys its alternating-product
binomial.  Squaring the channel amplitudes forgets both the polygon phases
and these global factorization constraints.  Equations `(9)`--`(12)` are an
exact finite support filter, not a finite reduction of JC(2).

## 3. Target-shear frustration

At a point of `C^2`, put

```text
u=grad P,     v=grad Q,     alpha=||u||^2,
beta=<v,u>/alpha,
```

where the Hermitian inner product is linear in its first entry.  For every
complex scalar `h`, the Gram determinant identity gives

```text
||v+h u||^2 = |kappa|^2/alpha + alpha |h+beta|^2.       (13)
```

The pointwise energy-minimizing shear is `h=-beta`.  A legal polynomial
target shear `Q -> Q+H(P)` can realize it only if

```text
beta=B(P) for some B in C[t].                          (14)
```

If `(14)` holds, choose `H'=-B`.  The two sheared gradients are Hermitian-
orthogonal, and their squared norms multiply to `|kappa|^2`.  These norm
squares are polynomials in the doubled variables `(x,y,bar x,bar y)` whose
product is a unit, hence both are constant.  Equivalently, each holomorphic
gradient is a bounded polynomial and therefore constant.  The sheared pair
is affine, so the original pair is a target shear of an affine automorphism.

Consequently every nonautomorphic planar Keller pair has genuine global
shear frustration:

```text
beta notin C[P].                                       (15)
```

THM-2230 identifies the same target shears as the entire fibre of mates with
fixed first component and fixed constant Jacobian response.  The theorem is
not needed for identity `(13)`; it explains why no omitted polynomial mate
can evade the frustration test.

## 4. A conditional Puiseux escape bound

Let `d=max(deg P,deg Q)`.  Suppose a convergent Laurent--Puiseux branch in a
sector has leading expansions

```text
x(t)=b t^(-r)+higher powers,       b!=0, r>0,
F(x(t))=a+q t^m+O(t^(m+delta)),    q!=0, m,delta>0.     (16)
```

Set `T(t)=||DF(x(t))||_F^2`.  Then

```text
T(t) >= const |t|^(-2(r+m)),
T(t) <= const |t|^(-2r(d-1)),                           (17)
m <= r(d-2).                                            (18)
```

To prove the lower bound, differentiate `(16)` and use

```text
x'=DF(x)^(-1) (F o x)',
||DF^(-1)||_F=||DF||_F/|kappa|                         (19)
```

for a `2 x 2` matrix of determinant `kappa`.  The derivative orders are
`|t|^(-r-1)` and `|t|^(m-1)`, giving the first inequality in `(17)`.
Polynomial degree gives the upper inequality.  Comparing exponents proves
`(18)`.  Combining the lower bound with `(5)` yields

```text
dist_F(C(t)/T(t),{rank<=1})=O(|t|^(2(r+m))).            (20)
```

In particular, such a branch cannot occur at degree two; at degree three it
must satisfy `m<=r`.

The existence of an expansion `(16)` for every hypothetical nonproper map
is a separate curve-selection-at-infinity import and is not asserted here.
No part of this theorem proves properness or JC(2).  It says that any branch
presented in this standard form has the displayed quantitative profile.

## 5. Counterexample search consequence

A conductance-based search should therefore retain, not dephase away:

1. the Wilson phase of each differential plaquette;
2. the complex closing phases of every coefficient polygon;
3. the global Segre cycle binomials coupling different fibres; and
4. the variation of the pointwise Schur minimizer across a whole `P`-fibre.

The cheapest viable cells have at least three coefficient channels, a
balanced dark-plaquette escape, and nontrivial cycle holonomy.  Trees,
singleton fibres, isolated two-channel cancellations, and globally removable
shears have already lost one of the required sidecars.
