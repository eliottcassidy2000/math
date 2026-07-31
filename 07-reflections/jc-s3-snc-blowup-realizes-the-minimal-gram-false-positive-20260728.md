# An `S3`-equivariant SNC blow-up realizes the minimal Gram false positive

> **STATUS: PROVED GEOMETRIC CONTROL, NOT CANON.**  This constructs a smooth
> projective rational `S3`-surface with an `S3`-stable SNC boundary whose
> Hermitian Gram block has one standard zero mode while the exact boundary
> presentation has no mod-two kernel.  It is a geometric hostile for the
> implication “singular Gram block implies a Kummer plane.”  It is not a
> quartic-resolvent completion, a Keller map, or an `A4/S4/JC(2)` conclusion.

The abstract `I_(1,6)` hostile in
[`jc-boundary-relations-are-a-mod-two-presentation-not-a-gram-kernel-20260728.md`](jc-boundary-relations-are-a-mod-two-presentation-not-a-gram-kernel-20260728.md)
is actually realized by irreducible boundary curves, with the full reflection
already present.

## 1. Equivariant blow-up

Let `S3` act on `P^2` by permuting homogeneous coordinates and let

```text
p_0=[1:0:0],  p_1=[0:1:0],  p_2=[0:0:1].                 (1)
```

Blow up this orbit.  Write `E_i` for the total-transform exceptional class
over `p_i`.  At `p_i`, the order-two stabilizer swaps the other two local
coordinates, so it fixes the tangent direction in which those coordinates
are equal.  Let `q_i` be the corresponding point of the exceptional curve
over `p_i`.  The set `{q_0,q_1,q_2}` is `S3`-stable.  Blow up this second
orbit and write `F_i` for the new exceptional classes.

The action lifts at both stages.  On the resulting smooth rational surface
`X`, the strict transforms

```text
D_i = E_i-F_i                                             (2)
```

are three mutually disjoint smooth rational curves of self-intersection
`-2`.  Thus

```text
D=D_0 union D_1 union D_2                                 (3)
```

is a reduced `S3`-stable SNC divisor.  The three-cycle permutes the `D_i`,
and every transposition supplies the reflection completing that cycle to
`S3`.

## 2. The presentation is injective, even modulo two

In the standard total-transform basis,

```text
Pic(X)=Z H direct_sum (direct_sum_i Z E_i)
                   direct_sum (direct_sum_i Z F_i),
Q=diag(1,-1,-1,-1,-1,-1,-1).                            (4)
```

The boundary-class map is

```text
delta:Z^3 -> Pic(X),       e_i |-> E_i-F_i.              (5)
```

Its image is primitive: in each direct summand
`Z E_i direct_sum Z F_i`, the vector `E_i-F_i` is primitive, and the three
summands are disjoint.  Hence

```text
ker(delta)=0,             L_sat/L=0.                     (6)
```

Modulo two the three images are the independent vectors `E_i+F_i`, so

```text
ker(delta mod 2)=0.                                      (7)
```

For `U=X\D`, the exact boundary-presentation theorem therefore gives

```text
H^1_et(U,mu_2)=0.                                        (8)
```

In particular there is no unit, Picard, `C3`-standard, or `S3`-standard
Kummer class.

## 3. The Gram shadow nevertheless has a standard zero mode

The boundary Gram matrix is

```text
M=(D_i.D_j)=-2 I_3.                                      (9)
```

It is zero modulo two.  The permutation module `F_2^3` is the direct sum of
one trivial line and the standard plane

```text
W={(x_0,x_1,x_2):x_0+x_1+x_2=0}.                        (10)
```

Consequently the `F4` Hermitian standard block is the `1 by 1` zero matrix:

```text
B=[0],             nullity_F4(B)=1,                     (11)
```

although the exact standard presentation matrix is injective.  Equivalently,
if `A=[1;1]` and `G=diag(1,1)` are the standard coordinates of `(5)` and
`(4)`, then

```text
A^dagger G A=[0],   ker(A)=0,   rad(im A)=im A.          (12)
```

The mod-four pairing detects the false positive.  On the basis
`u=(1,1,0)`, `v=(0,1,1)` of `W`,

```text
lambda_M(x,y)=x^T M y/2 mod 2

lambda_M|_W = [0 1]
             [1 0].                                      (13)
```

Thus the Gram zero mode is symplectic, whereas every exact presentation
kernel is totally isotropic.

## 4. Exact boundary of the control

This construction removes two possible excuses for the Gram false positive:
the classes are effective irreducible SNC components, and the `C2`
reflection exists geometrically.  What it does **not** provide is the much
more restrictive realization relevant to the quartic-resolvent program:
there is no graph quartic, finite cover, Jelonek divisor, or Keller equation
attached to `(1)`--`(3)`.

The remaining geometric question should therefore be stated as

```text
which boundary presentations occur in the specified quartic-resolvent
completion?
```

not merely whether an `S3`-equivariant rational SNC realization exists.

QED.
