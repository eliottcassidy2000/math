# Balanced response `e=2,h=1` derivation

Status: exact proof draft, not canon.

Let the unique pole be translated to zero.  Then `T=x`, `D=x^N`,
`deg E=2`, `deg S=N-4`, and THM-2796 equation (2) is

```text
2xSE' + xES' - NES = C != 0.
```

Put `P=SE^2`.  Multiplication by `E` gives

```text
xP'-NP=CE.                                                (A)
```

Write `E=x^2+ux+v`.  Since `P` is monic of degree `N`, (A) forces

```text
P=x^N-C[x^2/(N-2)+ux/(N-1)+v/N].                         (B)
```

The roots of `E` are nonzero and distinct.  Scale one to `1` and write the
other as `rho`; then `u=-(1+rho)`, `v=rho`.  Requiring `P(1)=P(rho)=0`
gives

```text
C=N(N-1)(N-2)/(N-(N-2)rho),                              (C)

H_N(rho)
 =rho^(N-1)(N-(N-2)rho)-N rho+(N-2)=0.                   (D)
```

The two possible denominators in (C)/(B) cannot vanish at a root of (D).
Conversely (B)--(D) give `E|P`.  At a root `z` of `E`, (A) and `z!=0`
give `P'(z)=0`, so `E^2|P`.  Differentiating (A) gives

```text
zP''(z)=CE'(z)!=0,
```

so both roots are exactly double.  Any other repeated root of `P` would,
by (A), also be a root of `E`; hence `S=P/E^2` is squarefree and disjoint
from `E`.  Also `P(0)=-Cv/N!=0`, so the pole is disjoint.  This proves the
converse factor typing.

For `N>=4`, `rho=1` is a root of `H_N` of exact multiplicity three because

```text
H_N'''(1)=-N(N-1)(N-2)!=0.
```

Every other root is simple.  Indeed `H_N'=0` says

```text
rho^(N-2)[(N-1)-(N-2)rho]=1.
```

Combining this with `H_N=0` gives

```text
(N-2)(N-1)(rho-1)^2=0.
```

Thus there are exactly `N-3` admissible ordered ratios.  Moreover

```text
rho^N H_N(1/rho)=-H_N(rho),
```

so swapping the two double zeros acts by `rho<->1/rho`.  The only possible
fixed ratios are `+/-1`; `1` is forbidden and `-1` occurs exactly for even
`N`.  Hence the number of affine classes with an unordered pair of double
zeros is

```text
floor((N-2)/2).                                           (E)
```

The promised Chebyshev interpretation is exact.  Put

```text
rho=exp(2 i theta),                    0<theta<pi.
```

After multiplication by `exp(-iN theta)`, equation (D) becomes

```text
N sin((N-2)theta)-(N-2)sin(N theta)=0.                    (F)
```

But

```text
U_(N-2)(cos theta)=sin((N-1)theta)/sin(theta),
```

and differentiating with respect to `theta` shows that (F) is precisely

```text
d/dtheta U_(N-2)(cos theta)=0.
```

The `N-3` interior critical points of the real-rooted Chebyshev polynomial
`U_(N-2)` therefore give all admissible ratios.  In particular every ratio
lies on the unit circle and is simple; the endpoint `theta=0` is the
discarded triple root `rho=1`.  Reflection `theta<->pi-theta` is exactly
`rho<->rho^(-1)`, and the midpoint `theta=pi/2` gives `rho=-1` precisely
when `N` is even.  This supplies a second all-degree proof of the root count
and inversion-orbit count in (E).

The `N=4` class is split.  Every `N>=5` class is genuinely nonsplit by
THM-2796 criterion (4).

For the first open degree-26 chamber, `deg V=2N-2=26`, hence `N=14`.  Here

```text
H_14(rho)
=-2(rho-1)^3(rho+1)
  (6rho^10+5rho^9+10rho^8+8rho^7+12rho^6+9rho^5
   +12rho^4+8rho^3+10rho^2+5rho+6).
```

There are eleven ordered ratios and six unordered affine classes.
The THM-2796 dessin census independently gives class counts `1,1,2,2,3`
for `N=4,...,8`, agreeing with (E).

Remaining before canonization: standalone assertion-free companion,
normal/optimized/stored replay, exact equivalence scope, and independent
audit.  A useful optional refinement is to identify the reciprocal quotient
as a classical paraorthogonal/Jacobi-type polynomial and prove its unit-circle
root locus; that locus is not needed for (A)--(E).
