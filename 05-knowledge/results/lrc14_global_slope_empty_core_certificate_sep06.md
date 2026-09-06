# A universal slice slope and an explicit height-601 network reduction

**Status: PROVED ANALYTICALLY; global slope conditional only on the named
finite coefficient-box audit.** The arithmetic statements below do not use
a speed-height extrapolation. The final finite head is a separate explicit
proof obligation, pending completion by the root session.

This note extends the general affine-slice mechanism of the independently
audited [full-support note](lrc14_empty_core_certificate_sep06.md) and
[actual-zero-coordinate addendum](lrc14_pair_relation_empty_core_certificate_sep06.md).
The relevant consumer remains the exact complete THM-4414 network sum;
entry, synchronization, and LRC(14) are separate open obligations.

Let `w=(a,b,c)` be primitive, distinct, positive, odd, ternary-unit, and
sorted. Let `v.w=0` be a primitive nonzero integer relation. It has at least
two nonzero coordinates and even norm `S=||v||_1`. The impossible support-two
pattern `(1,1,0)` would force equal speeds. The norm-four pattern `(1,1,2)`
is already closed for the selected network target by THM-4422.

The new board is: projected error zonotope; even section width; lattice-rule
discrepancy; maximum coefficient; shortest integer relation; finite head.
The new connection preserves the *whole slice-width integral* and its peak.
These bound the arithmetic sampling loss without deleting any defect layer.
The inherited false alternative is to infer a universal count bound from
the one-ray average: the norm-four exception remains load-bearing throughout.

## 1. Slice widths are sections of one planar zonotope

Put `r=3/14` and choose an index `i` with `v_i!=0`. Apply the linear map

```text
e -> (v.e, (w cross e)_i/v_i)
```

to the cube `[-r,r]^3`, and call its planar image `Z`. This is a centrally
symmetric convex zonotope. Its section at first coordinate `delta` is the
exact scalar carrier interval of the earlier notes. Write its length as

```text
f(delta)=T_delta(v,w).
```

Central symmetry makes `f` even. Convexity makes section lengths concave on
their support: the convex combination of two section intervals lies inside
the intermediate section, so their lengths satisfy the concavity inequality.
Consequently `f` is nonincreasing on `[0,infinity)`, after extension by zero
outside its support. It is continuous in the interior of
`[-rS,rS]`, with all physically live defects satisfying the strict bounds.

At the support endpoints we set `f=0`, matching strict physical eligibility;
this changes no area integral. For a full-support relation the actual endpoint
width is already zero. For a relation with an actual zero coordinate the
closed zonotope can have vertical endpoint edges, so the chosen convention
can make a downward jump. Monotonicity and the lattice-rule estimate below
remain valid; endpoint continuity is neither asserted nor needed.

The three generator vectors of `Z` are the images of the three coordinate
unit vectors. Their pair determinants are, up to sign, exactly `a,b,c`.
For example, if `i=1`, they are

```text
g_1=(v_1,0), g_2=(v_2,-c/v_1), g_3=(v_3,b/v_1),
det(g_1,g_2)=-c, det(g_1,g_3)=b, det(g_2,g_3)=-a,
```

where the last identity uses `v.w=0`. The area of a planar zonotope generated
by three segments of coefficient range `[-r,r]` is the sum of the three
pair-parallelogram areas. This can be checked by ordering generator angles
and triangulating its hexagon; collinear limiting cases follow by continuity.
Thus Fubini gives the exact integral

```text
I=integral_R f(t)dt=area(Z)=4r^2(a+b+c)
                               =9(a+b+c)/49 <=27c/49. (1)
```

The area is independent of the relation coefficients. The coefficients
reappear in the arithmetic sampling and peak height, which must be retained.

## 2. A universal discrepancy bound from the central section

For any even nonnegative function decreasing on the positive half-line and
any `h>0`, rectangle comparison on consecutive length-`h` intervals gives

```text
|sum_(k in Z) f(hk)-I/h| <= f(0).                     (2)
```

Indeed, `h sum_(k>=1)f(hk) <= integral_0^infty f
<=h sum_(k>=0)f(hk)`; double and rearrange.

The earlier owner-residue decomposition defines the exact *slope load*

```text
F_v(w) = (2/3) sum_k f(3k)                       (unit relation),
F_v(w) = (1/3)(sum_k f(k)-sum_k f(3k))       (one-zero-mod3 relation).
```

The second formula also covers an actual zero integer coefficient. All
primitive relations fall into these two residue types: two zero residues
would force the third to vanish, contradicting primitivity. In both cases,
(2) yields

```text
F_v(w) <= (2/9)I+(2/3)f(0).                           (3)
```

Let `M=max_i |v_i|`, and choose a coordinate attaining it. At defect zero,
the image carrier is `t v`. Since it is `w cross e` for an error in the
cube, its chosen coordinate obeys

```text
|t| M <=r(sum(w)-w_i).
```

Thus

```text
f(0)<=2r(sum(w)-w_i)/M <=6c/(7M).
```

Equations (1)-(3) give, at every coefficient height,

```text
F_v(w)/c <=6/49+4/(7M).                               (4)
```

For `M>=19`, the right side is at most

```text
6/49+4/133=142/931 <15/98.                            (5)
```

Therefore the universal slope assertion

```text
F_v(w)<=15c/98                                        (6)
```

for every relation except the norm-four pattern reduces to the finite box
`max_i |v_i|<=18`. The box must include *all* primitive magnitude patterns,
of support two or three and even norm, with at most one zero residue modulo
three. It is not a norm-twenty or height census. The full-support norms can
reach fifty-four. The earlier convex-speed vertex compiler applies unchanged
after choosing a nonzero pivot coordinate. The support-two rectangle formula
provides a separate path for that boundary.

**Finite box obligation.** Check (6) for that complete box after excluding
the impossible `(0,1,1)` pattern and the inherited `(1,1,2)` exception.
The independent `three_ray_geometry` referee is compiling this exact box in
a separate artifact; this note does not treat a running computation as a
proved dependency.

## 3. Uniform intercept in terms of relation norm

Let `D` be the allowed integer defect list. Each unit-relation line has
strict count bound `(2/3)T_delta+4/3`, and each one-zero-mod-three line has
strict count bound `T_delta/3+1`. Summing gives

```text
N=|Lambda(w)| <F_v(w)+B_v.
```

In the unit case,

```text
|D|=|{delta in 3Z: |delta|<3S/14}| <S/7+1,
B_v=4|D|/3 <4S/21+4/3 <=2S/7+4/3.
```

In the one-zero case, the residue-deleted integer count in the same open
interval gives `B_v=|D|<2S/7+4/3`. Consequently (6) implies

```text
N <15c/98+2S/7+4/3.                                  (7)
```

The sufficient count gate `N<=2c/11`, and hence all three strict network
certificates, follows whenever

```text
c >= (308/31)S+4312/93.                               (8)
```

Strictness is inherited from the open carrier intervals. Neither a rounding
convention nor an endpoint equality can create an extra carrier here.

## 4. A short relation from an explicit planar area calculation

Project the relation lattice onto its first two coordinates:

```text
L={ (x,y) in Z^2 : ax+by=0 mod c }.
```

The map `(x,y)->ax+by mod c` is onto because `gcd(a,b,c)=1`; therefore this
lattice has determinant `c`. The third coordinate is
`z=-(ax+by)/c`. For `L0>0`, the projected `l1` ball

```text
K_(L0)={ (x,y): |x|+|y|+|(ax+by)/c|<=L0 }
```

has vertices

```text
+/- (L0*c/(a+c),0),
+/- (0,L0*c/(b+c)),
+/- (L0*b/(a+b),-L0*a/(a+b)).
```

A shoelace sum gives

```text
area(K_(L0)) = [2c(ab+ac+bc)/((a+b)(a+c)(b+c))] L0^2
                                                    >(3/4)L0^2. (9)
```

To check the inequality, scale `c=1`, so `0<a,b<1`. Its numerator difference
is the explicit nonnegative sum

```text
8(ab+a+b)-3(a+b)(a+1)(b+1)
 =3a(1-a)(b+1)+3b(1-b)(a+1)+2a(1-b)+2b(1-a)>0.
```

The strictness uses distinct sorted positive speeds, hence `a,b<c`.

Here is the needed planar lattice argument. If a centrally symmetric open
convex body has area greater than `4c`, its half-body has area greater than
one fundamental domain of `L`. Reducing translates of the half-body into
that domain forces positive overlap between two distinct lattice translates
(otherwise their total area could not exceed `c`). Their two points have a
nonzero lattice difference lying in the original open body, by symmetry and
convexity. This is the elementary fundamental-domain proof of the planar
Minkowski bound, with all hypotheses stated.

Apply it at `L0=4sqrt(c/3)`. By (9) there is a nonzero integer relation with

```text
S=||v||_1 <4sqrt(c/3).                                (10)
```

Dividing by its content makes it primitive without increasing its norm.
It still has even norm. If its pattern is norm four, invoke THM-4422;
otherwise (7)-(8) apply once the finite coefficient box is certified.

## 5. Exact high-height cutoff

Assume the finite coefficient box has passed. For `c>=603`, choose the
relation in (10). If `S<=56`, then

```text
(308/31)S+4312/93 <=56056/93 <603 <=c.
```

If `S>=58`, define

```text
g(S)=3S^2/16-(308/31)S-4312/93.
```

It has `g(58)=3023/372>0` and is increasing for `S>=58`. Equation (10) gives

```text
c>3S^2/16>(308/31)S+4312/93.
```

The two cases exhaust even norms. Thus every eligible triple with `c>=603`
has the selected network certificate, and every triple outside the inherited
norm-four class has all three projections strictly below the target.

The only remaining speed universe is therefore

```text
1<=a<b<c<=601,
a,b,c odd and nonzero modulo3,
gcd(a,b,c)=1.                                        (11)
```

The endpoint is `601` because `602` is even. An exact complete verification
of (11), with raw or independent literal projections and honest norm-four
handling, will close the universal degree-zero local projection target.
That finite verification is a separate obligation, not a consequence of
the coefficient-box audit. It still will not supply chart entry,
synchronization, or a proof of LRC(14).
