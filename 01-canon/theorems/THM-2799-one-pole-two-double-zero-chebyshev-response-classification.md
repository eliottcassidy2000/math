---
id: THM-2799
title: "One-pole two-double-zero Chebyshev response classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  balanced-response chamber e=2,h=1 is classified in every
  degree N>=4.  After affine normalization its maps are indexed by the
  simple unit-circle roots of an explicit reciprocal polynomial R_N, modulo
  inversion, giving exactly floor((N-2)/2) affine/target classes.  R_N is
  the critical-point polynomial of the Chebyshev U_(N-2) carrier.  N=4 is
  split; every N>=5 class is genuinely nonsplit.  This is a response-layer
  theorem, not Keller-chart entry, JC(2), or DC(2).
source: root/jc-e2-chebyshev-response-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2768-modular-c2-c3-a4-s4-bass-serre-quotient
script: 04-computation/jc2_one_pole_e2_chebyshev_response_thm2799.py
output: 05-knowledge/results/jc2_one_pole_e2_chebyshev_response_thm2799.out
script_sha256: 36ee6c93b0bb878bdfda2b772e23a5ecf5960a3128d5ee2829cc36a00e79c266
output_sha256: faf0382f2f7c93db904ef02bd3899f86b16da90a9f51674d03e106e4469ed827
hash_basis: LF-normalized bytes
---

# THM-2799 -- one-pole two-double-zero Chebyshev response classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2796 classifies every balanced square-potential response with at most
one extra double zero.  The first open chamber is `e=2`.  When there is one
pole, its apparently nonlinear Stieltjes system collapses to a sparse
four-term numerator and then to the critical points of one Chebyshev
polynomial.  This gives the first complete all-degree chamber beyond the
`e<=1` wall.

## 1. Complete normal form

Work over an algebraically closed field of characteristic zero.  Consider a
balanced nonconstant response from THM-2796 with

```text
e=2,                  h=1,                  p=(N),
s=N-4,                r=N-2,               N>=4.     (1)
```

Translate the unique pole to `0`, scale one double zero to `1`, write the
other as `t`, and normalize `F(infinity)=1`.  Thus

```text
t!=0,1.                                               (2)
```

Define

```text
J_N(t)=sum_(j=0)^(N-2) t^j,

R_N(t)=sum_(j=0)^(N-3) (j+1)(N-2-j)t^j.              (3)
```

Then every response in `(1)` has

```text
R_N(t)=0,                                             (4)

a=-(N/2)J_N(t),
b=N(J_N(t)-1),
c=N-1-(N/2)J_N(t),

B_t(x)=x^N+a x^2+b x+c,                              (5)
E_t(x)=(x-1)(x-t),
S_t(x)=B_t(x)/E_t(x)^2,
q_t=N(N-2)J_N(t)/2,                                  (6)

F_t(x)=B_t(x)/x^N,
G_t(x)=q_t E_t(x)/(2x^(N+1)),
V_t(x)=4x^(N+2)S_t(x)/q_t^2,             kappa=1.   (7)
```

Conversely, every root `t` of `R_N` gives the displayed response.  The
polynomial `S_t` is monic, squarefree, has degree `N-4`, and is disjoint
from `0,1,t`.  Hence `(3)--(7)` are a complete classification, not merely
necessary equations.

### Proof of sparsity and the ratio equation

Write `F=B/x^N`, where `B=S E^2` is monic of degree `N`.  The divisor of
`F'` in THM-2796 has numerator exactly the monic quadratic `E`, up to a
nonzero scalar.  Therefore

```text
xB'-NB=qE,                       q!=0.                (8)
```

For every coefficient `b_k x^k` of `B`, the left side of `(8)` contains
`(k-N)b_k x^k`.  Since its degree is two, all coefficients in degrees
`3,...,N-1` vanish.  Thus

```text
B=x^N+a x^2+b x+c.                                  (9)
```

The four equations

```text
B(1)=B'(1)=B(t)=B'(t)=0                              (10)
```

first give the coefficients in `(5)`.  After the derivative equations and
`B(1)=0` are substituted,

```text
2 [B(t)-B(1)]/(t-1)
 =-[t-1]^2 R_N(t).                                  (11)
```

Since `t!=1`, the last root equation is exactly `(4)`.

For a root of `R_N`, `J_N(t)` cannot vanish: otherwise
`t^(N-1)=1` and the left side of `(11)` before factorization is
`2-2N!=0`.  Hence `q_t!=0`.  Equation `(8)` becomes

```text
xB_t'-N B_t=q_t(x-1)(x-t).                          (12)
```

It proves both divisibility by `E_t^2` and the absence of hidden
degeneracies.  Indeed, every multiple root of `B_t` must be a root of the
quadratic on the right.  A root of multiplicity at least three would make
the right side have multiplicity at least two, impossible because
`t!=1`.  Thus `1,t` are exact double roots and all remaining roots are
simple.  The constant term is nonzero because otherwise `(12)` would force
`0` to be a root of `E_t`.

Finally, `(12)` gives `F_t'=q_tE_t/x^(N+1)=2G_t`, and direct substitution
in `(7)` gives

```text
F_t=V_t G_t^2.                                       (13)
```

This proves the response equation and the converse classification.

## 2. Chebyshev root theorem

Put `n=N-2`.  A direct summation of `(3)` gives

```text
(1-t)^3 R_N(t)
 =n-(n+2)t+(n+2)t^(n+1)-n t^(n+2).                 (14)
```

For `t=exp(2ix)`, `0<x<pi`, the right side of `(14)` is a nonzero scalar
factor times

```text
n sin((n+2)x)-(n+2)sin(nx)
 =2 sin(x)^2 d/dx [sin((n+1)x)/sin(x)].             (15)
```

The bracketed quotient is

```text
U_n(cos x),                                          (16)
```

the Chebyshev polynomial of the second kind.  It has `n` simple roots in
`(-1,1)`.  Rolle interlacing gives exactly `n-1` distinct roots of
`U_n'` there.  Since `x -> exp(2ix)` is injective on `0<x<pi`, equations
`(14)--(16)` exhibit all

```text
deg R_N=n-1=N-3                                     (17)
```

roots of `R_N`: they are simple, nonzero, unequal to one, and lie on the
unit circle.

The polynomial is reciprocal:

```text
t^(N-3) R_N(1/t)=R_N(t).                            (18)
```

Swapping the two double zeros changes `t` to `t^(-1)` and is implemented by
the source scaling `x -> tx`.  Conversely, the unique pole and unordered
double-zero pair show that these are the only identifications.  When `N`
is even, `t=-1` is the unique inversion-fixed root; when `N` is odd there
is none.  Therefore the exact number of affine/target classes is

```text
#classes(e=2,h=1,N)=floor((N-2)/2).                  (19)
```

This chamber first becomes passport-multiple at `N=6`, where `(19)` gives
two maps.  THM-2796's separate `e=2,h=2,p=(4,1)` chamber already gives the
first passport multiplicity overall at `N=5`.

## 3. Split boundary and the quartic/V4 lane

At `N=4`,

```text
R_4(t)=2(t+1),              t=-1,
B=(x^2-1)^2,               S=1,
V=x^6/4.                                            (20)
```

Thus the only degree-four member of this chamber is split.  For every
`N>=5`, `S_t` is nonconstant and squarefree, so `V_t` is genuinely
nonsplit by the squareclass criterion of THM-2796.

This is a useful stopping rule for the proposed quartic `V4`/resolvent-cubic
program: the `e=2,h=1` balanced response cannot supply a nonsplit
degree-four object.  It does **not** classify graph quartics of Keller maps,
their cubic resolvents, or the degree-four point-cap branch.  Any such
quartic mechanism must enter through a different response chamber or
through data discarded before the response reduction.

## 4. Exact controls

The companion performs, with exact rational polynomial arithmetic:

1. the sparse numerator, double-root, quotient, and logarithmic-derivative
   identities through `4<=N<=14`;
2. reduction modulo `R_N` of every algebraic coefficient identity;
3. reciprocity and nonzero discriminant of every tested `R_N`;
4. exact real-root counts for `U_(N-2)'` on `(-1,1)`;
5. the inversion-orbit count `(19)`;
6. the split `N=4` and first nonsplit `N=5` controls; and
7. normal, optimized, and stored transcript identity.

The script contains no Python `assert` node.  Run:

```text
python 04-computation/jc2_one_pole_e2_chebyshev_response_thm2799.py
python -O 04-computation/jc2_one_pole_e2_chebyshev_response_thm2799.py
```

The finite controls support but do not replace the all-degree proofs.

## 5. Scope and next exact target

This theorem classifies an abstract balanced response layer.  It does not
prove that one of these responses enters a polynomial Keller chart, satisfy
the two inherited Faber flux equations, or arise from a Weyl-algebra
endomorphism.  Consequently it proves neither JC(2) nor DC(2).

The next response target is now sharply isolated: `e=2,h=2`, beginning with
the two `N=5,p=(4,1)` Nielsen classes detected in THM-2796.  There the pole
cross-ratio is a genuine retained coordinate; a passport-only
classification is impossible.  The cheapest decisive test is to combine
the degree-two Stieltjes polynomial `E` with that cross-ratio and then
intersect the resulting finite atlas with the Faber quotient and flux
equations.
