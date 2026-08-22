# Arbitrary-ratio Dirichlet block reduction for disconnected-low LRC14

**Status.** PROVED reduction, exact-arithmetic audited.  This note does **not**
certify the residual affine rays and does not close the disconnected-low
branch.  Its output is a finite bank plus a finite raw-level head.

## 1. The arbitrary-ratio one-tooth function

Fix a body-safe cell, labels `1<=e,f<=14`, ruler `L>=168` divisible by 14,
and levels `p<q<8p`.  Put

```
z=Lp-e,  w=Lq-f.
```

Body safety makes all `p` first-clause teeth full: if `r=ej mod L`, then
`r>=L/14` and `r+e<=13L/14`, so the nominal indices `0,...,p-1` are exactly
the physical teeth and none clips at 0 or 1.

Map one first tooth by the affine coordinate of the second clause.  Its image
has length

```
lambda = w/(7z).
```

Let `C={x mod 1: ||x||<=1/14}`, so `|C|=1/7`.  If
`lambda=m+r`, `m=floor(lambda)`, `0<=r<1`, the overlap contributed by that
first tooth, as a function of its circular phase `x`, is

```
f(x)=(L/w)(m/7 + H_r(x)),
H_r(x)=int_0^r 1_C(x+u) du.                              (1)
```

This is the full-turn decomposition of the periodized trapezoid in THM-3352
(formerly numbered THM-3350 during the live session).  It is valid at every
ratio; the large full-turn baseline is constant and therefore irrelevant to
oscillation.

Equation (1) immediately gives

```
mean(f)=L/(49z),
osc(f)<=L/(7w).                                           (2)
```

Moreover, almost everywhere,

```
H_r'(x)=1_C(x+r)-1_C(x).
```

Each indicator has two unit jumps around the circle, so the total variation
of the piecewise-constant derivative is at most four.  Thus

```
Var(f')<=4L/w.                                            (3)
```

This is the repair of the false cap-two carryover.  The *absolute height* of
`f` can be large when `q/p` is large; only (2)'s oscillation is used below.

For `G_d(x)=sum_{s=0}^{d-1} f(x+s alpha)`, translation and subadditivity give

```
mean(G_d)=dL/(49z),  osc(G_d)<=dL/(7w),
Var(G_d')<=4dL/w.                                        (4)
```

## 2. Dirichlet block and the many-turn estimate

Write `h=q-p`.  Pigeonhole approximation to `h/p in (0,7)` supplies

```
1<=d<=8,  0<=a<=7d,
c=dh-ap,  |c|<=p/9.                                      (5)
```

Put `D=d+a`.  For `alpha=(w-z)/z`,

```
d alpha-a=A/z,   A=Lc+eD-df.                             (6)
```

Since `D<=8d`, the exact label perturbation bound is

```
|eD-df|<=888.                                             (7)
```

Consequently, with

```
R(p,L)=(Lp/9+888)/(Lp-14),
rho=|A|/z,
n=floor(p/d),
T=n rho,
```

one has `rho<=R`; for `p>=264,L>=168`, `R<1/2`, so `rho` is the actual
short signed rotation step.

Group the first `dn` nonnegative samples into `d`-blocks and discard fewer
than `d` remaining nonnegative samples.  The exact composite-trapezoid
identity along the lifted `G_d` path has three terms.

1. If `T>=5`, its complete turns supply at least
   `(4/5)(p-7)/(49p)`.
2. The endpoint half-difference is bounded by oscillation, not height:
   `osc(G_d)/2<=dL/(14w)<=4L/(7w)`.
3. Each derivative jump contributes at most `rho/8` times its jump size.
   Using (4), `ceil(T)<=n rho+1`, and `d<=8`, the Peano loss is at most
   `L(p rho^2+8rho)/(2w)`.

Since `w>=L(p+1)-14`, every `T>=5` channel is therefore bounded below by

```
B(p,L)=4(p-7)/(245p)
       -4L/[7(L(p+1)-14)]
       -L[pR(p,L)^2+8R(p,L)]/[2(L(p+1)-14)].              (8)
```

The numerators of `dB/dp` and `dB/dL`, after substituting
`p=x+264,L=y+168`, have respectively 36 and 20 monomials, all positive;
their smallest coefficients are 6373 and 8013.  Thus `B` increases in both
variables on the declared domain.  At its corner,

```
B(264,168)=149164364417/46870975182570,
B(264,168)-Dmax/5
 =85330033783953387991/7166476998347435667648750 >0.      (9)
```

For the stronger individual target,

```
B(273,168)-1/294
 =3539548829/134745104471250 >0.                          (10)
```

Hence raw `p>=264` many-turn edges supply the exact five-edge average needed
by a cross-component tree, while raw `p>=273` many-turn edges individually
exceed `1/294`.

## 3. Few turns and the affine bank

If `T<5`, then `n>= (p-7)/d`, so

```
|A| < 5d z/(p-7) < 5dLp/(p-7).
```

Using (7), `d<=8`, and `L>=168`,

```
|c| < 40p/(p-7)+37/7.
```

At `p=264` the right side is `83429/1799<47` and decreases thereafter.
Thus every nonzero resonance has

```
1<=|c|<=46.                                               (11)
```

For fixed `(d,a,c)` and a residue representative `1<=p0<=d` satisfying
`a p0+c=0 mod d`, put

```
q0=p0+(a p0+c)/d,
p=p0+dn,
q=q0+(d+a)n.                                              (12)
```

These are the affine resonance rays.  Eventual membership in `p<q<8p`
requires `c>0` when `a=0` and `c<0` when `a=7d`; interior `a` admits both
signs.  The exact parameter counts by `d=1,...,8` are

```
644, 1288, 1918, 2548, 3206, 3794, 4452, 5040,
```

for a total of **22,890**.  This is a finite cover, not a disjoint
decomposition: one physical point can have more than one Dirichlet witness.

On a nonzero-`c` ray,

```
dq-(d+a)p=c,
```

so `gcd(p,q)` divides `c`.  In particular every ray point with `p>=264` is
automatically a high primitive channel, since
`(p+q)/gcd(p,q)>(2*264)/46>8`.  Mixed gcd values must still be retained by an
all-raw compiler; they must not be silently identified with a fixed common
scale.

If `c=0`, reduction to coprime `(P,Q)` gives `P|d`, hence `P<=8`.  There are
145 high primitive channels with `P<=8,Q<8P`, including `(3,5)`.  For the
current `gcd(p,q)<=3` target these have raw `p<=24`, so no `c=0` point occurs
in the `p>=264` tail.  An all-raw strengthening must route these fixed-ratio
dilation families separately.

## 4. Exact residual compiler schema

The open compiler has two inputs.

* **Finite head:** every feasible physical context and every high raw pair
  with `p<264`, `q<8p`, restricted to the desired common scales (currently
  `g<=3`, or `g<=20` in the broader probe).  Evaluate with the exact
  THM-3352 Euclidean mass engine.
* **Affine tails:** for each of the 22,890 rows (12), start at the least `n`
  satisfying `p>=264` and the strict cone.  Use the arbitrary-ratio periodic
  PL function (1), the exact signed Peano identity, and an affine-ray
  continuum/error certificate `I(n)>=J0-K/p(n)`.  Alternatively generalize
  THM-3200's floor-branch stabilization method to these moving-ratio affine
  rays; THM-3200 itself applies only to a fixed primitive channel under common
  dilation and cannot be invoked verbatim here.  Retain the moving gcd only as
  a filter, never as a ray coordinate.

An all-raw affine certificate is preferable: because `c!=0` bounds the gcd,
it simultaneously covers every common scale occurring on that ray.  If such
a certificate fails, split by divisors of `c` (equivalently gcd strata) and
then by the finite floor-event residue modulus.  The `c=0` fixed-ratio bank is
separate.

Repository artifacts:

* `07-reflections/lrc14-generalized-dirichlet-resonance-reduction-20260812.md`
* `04-computation/lrc14_generalized_dirichlet_reduction_20260812.py`
* `05-knowledge/results/lrc14_generalized_dirichlet_reduction_20260812.out`
