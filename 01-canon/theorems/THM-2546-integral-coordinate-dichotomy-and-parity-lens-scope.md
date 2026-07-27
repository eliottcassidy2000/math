---
id: THM-2546
title: "Integral-coordinate dichotomy, compactified double shadow, and the exact Jelonek-slice fibre law"
status: >
  PROVED for (1) the global constant-leading r- and z-coordinate
  equations, hence integrality of y and z and x as the unique escape
  coordinate; (2) the exact factorization and compactified two-sheet
  shadow on L=0; (3) the complete affine fibre law on L=0 (one explicit
  survivor when p != 0, empty precisely when p=0); (4) the a=0, q=0,
  and z-triple exceptional-locus audit; and (5) the pointwise parity
  rigidity lemma.  All polynomial and rational identities have an exact
  symbolic certificate.  This is a theorem about the fixed sporadic map
  only and does not exclude G1 or any planar Keller stratum.
source: opus-2026-07-27; codex-2026-07-27 exact L-slice audit
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - HYP-9030-keller-degree-semigroup (test (ii))
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
script: 04-computation/keller_integral_coordinate_parity_opus_20260727.py
output: 05-knowledge/results/keller_integral_coordinate_parity_opus_20260727.out
---

# THM-2546 -- one coordinate carries all escape

Throughout, `F` is the sporadic Keller map of THM-2473, with source
coordinates `(x,y,z)`, target coordinates `(a,b,c)`, and

```text
L = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
p = 4 - 3bc,                 q = b^2 - 12a.
```

The Jelonek set is `{L=0}`, and the x-core is

```text
E_x(X) = L X^3 + p X - 2c.                         (1)
```

## (1) Two global integral-coordinate equations [PROVED]

Put `r=b-y` and

```text
P_r(R) = 27a^2c - 18aR + 3bR^2 - 2R^3,
D_r(R) = R^2 - bR + 3a.                            (2)
```

On every fibre, direct substitution gives

```text
P_r(r)=0,                   x D_r(r)=r.             (3)
```

The sign convention is load-bearing:

```text
P_r'(R) = -6D_r(R),                                 (4)
```

not `+6D_r`.  (The positive sign belongs to the polynomial in `y`, before
the reversal `r=b-y`.)  Since `P_r` has constant leading coefficient
`-2`, `r`, hence `y=b-r`, is integral over `Q[a,b,c]`.  Its discriminant is

```text
disc_R(P_r) = -2916 a^2 L.                          (5)
```

There is also a global constant-leading cubic for the source coordinate
`z`.  With a new indeterminate `W`, set

```text
Q_z(W) = 8W^3 + q_2 W^2 + 6L S_z W + L T_z,         (6)

q_2 = 324a^2c^2 - 216abc + 408a - 15b^3c + 6b^2,

S_z = 27a^2c^2 - 18abc + 52a + b^3c + 14b^2,

T_z = 729a^4c^4 - 972a^3bc^3 + 2322a^3c^2
      + 54a^2b^3c^3 + 270a^2b^2c^2 - 3735a^2bc - 338a^2
      - 36ab^4c^2 + 122ab^3c + 1372ab^2
      + b^6c^2 - 2b^5c - 80b^4.
```

The companion verifies the polynomial identity

```text
Q_z(F_1(x,y,z),F_2(x,y,z),F_3(x,y,z); z) = 0        (7)
```

in `Q[x,y,z]`.  Dividing (6) by `8` gives a monic equation over
`Q[a,b,c]`, so `z` is integral.  This upgrades the former five-target
sampled claim to **PROVED**.

Consequently every sequence with bounded target has bounded `y` and `z`.
All Jelonek escape is carried by the one non-integral coordinate `x`.

The earlier gauge hostile remains sharp.  Viete applied to (2) gives
`Tr(r)=3b/2`, hence `Tr(y)=3b-Tr(r)=3b/2`, whereas (1) gives
`Tr(x)=0`.  The depressed trisection is a property of the selected escape
coordinate, not an invariant of every primitive coordinate of the map.

## (2) The generic double root is an infinity shadow [PROVED]

Assume first `L=0` and `q != 0`.  Define

```text
r_0 = 3a(b-9ac)/q,
r_1 = 3(b^3-16ab+36a^2c)/(2q).                      (8)
```

The following are global identities before imposing `L=0`:

```text
q^3 [P_r(R)+2(R-r_0)^2(R-r_1)]
  = 27a^2(b^3-108a^2c-6qR)L,                       (9)

q^2 D_r(r_0) = 27a^2 L.                            (10)
```

Therefore, in the localized coordinate ring of `{L=0,q!=0}`,

```text
P_r(R) = -2(R-r_0)^2(R-r_1),       D_r(r_0)=0.     (11)
```

If also `a!=0`, then `r_0!=0`: otherwise (10) would give
`D_r(0)=3a=0`.  Thus (3) forbids an affine point over `r_0`, because it
would say `0=xD_r(r_0)=r_0`.  The double root is not two finite fibre
points.

This algebraic obstruction is exactly the analytic escape picture.  In a
transverse approach to a point with `L=0`, `p!=0`, `q!=0`, and `a!=0`,
(1) has two
roots

```text
x_+- ~ +-sqrt(-p/L).
```

Their bounded `r`-coordinates coalesce at `r_0`.  From the third component
of `F`,

```text
z = 2/x^2 - 3y/x - c/x^3,                           (12)
```

so both escaping sheets satisfy `z -> 0`.  They converge only in the
partial `(r,z)` compactification, to the double shadow `(r_0,0)`.

The z-cubic sees precisely the same phenomenon.  The exact congruence

```text
q_2 - 12L = 9(b^2p-2q)                              (13)
```

gives

```text
Q_z(W)|_{L=0} = W^2(8W+q_2).                        (14)
```

The double root `W=0` is the projection of the two escaping x-sheets, not
two affine fibre points.

## (3) Complete finite-survivor law on `{L=0}` [PROVED]

For `p!=0`, define

```text
x_s = 2c/p,
r_s = 3cq/(2p),
y_s = b-r_s,
z_s = -(9/8)(b^2p-2q).                              (15)
```

This is not an interpolation from examples.  Direct substitution gives the
global identities

```text
(3bc-4)^3 [F_1(x_s,y_s,z_s)-a]
   = 4(54ac^2-27bc+28)L,

(3bc-4)^2 [F_2(x_s,y_s,z_s)-b] = -36cL,

F_3(x_s,y_s,z_s)-c = 0.                             (16)
```

Thus (15) is an affine fibre point at every target with `L=0,p!=0`.
Moreover,

```text
pq(r_s-r_1)=6bL,
p^2(D_r(r_s)-3q/4)=12L.                             (17)
```

On the generic locus it is the simple root `r_1`, with
`D_r(r_1)=3q/4` and `x_s=r_1/D_r(r_1)`.  Equations (13)-(15) also give
`z_s=-q_2/8` on `L=0`.

There is exactly one affine fibre point, including at the boundary pieces:

- If `a!=0` and `q!=0`, (11), (3), and (17) exclude the double root and retain only
  (15).  Once `x,y` are fixed, `z` is unique: if `x!=0`, use the third
  component of `F`; if `x=0`, then `u=1` and use the first component.
- If `a=0,q!=0`, then `b!=0` and `L=0` forces `bc=1`.  The explicit
  substitution leaves only `(2/b,-b/2,9b^2/8)`, which is (15).  The
  tempting double-root lift `(2/b,b,-9b^2/8)` maps instead to
  `(-3b^2/8,b/4,1/b)` and is not in the fibre.
- If `q=0,p!=0`, the identity in (20) below forces `a=b=0`; direct solution
  gives the sole point `(c/2,0,0)`, again (15).

If `p=0`, then `b,c!=0`, `c=4/(3b)`, and

```text
L = q^2/(3b^2).                                     (18)
```

Hence `L=p=0` is exactly the omitted rational curve

```text
(a,b,c) = (4/(27t^2), 4/(3t), t),      t!=0.        (19)
```

Here (1) is the nonzero constant `-2c`, so the fibre is empty.  Combining
the two cases:

> **Exact Jelonek-slice law.** Every target on `{L=0,p!=0}` has the one
> finite preimage (15).  Every target on `{L=p=0}` has empty fibre.  The
> remaining two generic sheets are represented by the double `(r,z)`
> shadow but live at `x=infinity`.

## (4) Exceptional projections, separated [PROVED]

### The `a=0` coordinate collision is not the Jelonek surface

Equation (5) has the extra square factor `a^2`.  Indeed,

```text
P_r|_{a=0}=r^2(3b-2r),          L|_{a=0}=b^2(bc-1).
```

For `b!=0,bc!=1`, the target lies off `L`: the full three-point fibre has
two distinct finite residual points with the same `y=b` (`r=0`) and one
contracted-divisor point with `r=3b/2`.  This is only loss of information
under the y-projection.  On `bc=1`, the two residual points escape and the
contracted point above is the sole survivor.  Thus `a=0` must not be folded
into the generic pole argument.

### The `q=0` triple-shadow locus

Substitution `a=b^2/12` gives

```text
L|_{q=0}=b^2p^2/48.                                 (20)
```

Thus `{L=q=0}` has exactly two pieces:

- `b=0`, hence `a=0`: the target z-axis.  Here
  `P_r=-2r^3` and `Q_z=8W^3`, but the affine fibre has only
  `(c/2,0,0)`; one root is finite and two are escape shadows.
- `p=0,b!=0`: the omitted curve (19).  Here
  `P_r=-2(r-b/2)^3` and `Q_z=8W^3`, but the fibre is empty.

So a triple coordinate root can encode one or zero affine points; it is not
a three-point fibre count.

### Where the survivor also has `z_s=0`

On `L=0`, the z-cubic is triple exactly when `q_2=0`, equivalently
`b^2p-2q=0`.  Writing

```text
a=b^2(3bc-2)/24
```

gives the exact restriction

```text
L = b^2p^2(9b^2c^2+12bc-28)/192.                    (21)
```

Besides the z-axis and the omitted curve, this produces two curves
`9(bc)^2+12bc-28=0` with `bp!=0`.  Their fibre still consists of the one
point (15), now with `z_s=0`; the other two copies of the triple z-root are
escape shadows.

Two hostile boundary controls make the distinction concrete:

```text
(a,b,c)=(2/27,1,1):
  survivor=(2,5/6,-7/8),
  P_r=-2(r-2/3)^2(r-1/6),       Q_z=W^2(8W+7);

c=0, L=0, b!=0:
  a=b^2/16,                     survivor=(0,b,-63b^2/16).
```

The first has two shadows at `(r,z)=(2/3,0)` and exactly one affine point;
the second shows that `x_s=0` is a genuine survivor, not an omitted-fibre
signal.

## (5) Pointwise parity rigidity and HYP-9030 test (ii) [PROVED]

Let `G:R^n->R^n` be a real polynomial map with everywhere nonzero constant
Jacobian and complex degree `d`.  At every value whose complex fibre is full
(`d` points), complex conjugation pairs all nonreal points and fixes exactly
the real points.  Therefore

```text
# real fibre = d (mod 2).                            (22)
```

An even-degree wild map has even real count at every full-fibre value.  A
pointwise odd-real-count witness cannot exist there; odd counts can arise
only after escape, where the fibre is no longer full.  HYP-9030 test (ii) is
therefore released as stated.  The surviving parity routes are
escape-coordinate subleading laws, sigma-fixed fibres, and sphere degree at
values over which the map is proper.

## Loss ledger

- This theorem concerns the fixed three-variable sporadic map `F`.  It says
  nothing directly about `JC(2)`, `GMC(2)`, G1, or the point-cap stratum.
- Roots of `P_r` and `Q_z` count coordinate-projection multiplicity, not
  affine fibre points.  Equations (3), (18), and the exceptional audit are
  required before interpreting them.
- The displayed square-root asymptotic is for a transverse generic approach
  with `p!=0`, `a!=0`, and `q!=0`.  The polynomial identities and the finite-survivor law are
  global on their stated localized slices.

## Independent referee (opus, 2026-07-27, second proof route)

The z-integrality and the exact core `8Z^3 + A_z Z^2 + 6 L B_z Z + L C_z`
are independently re-proved by degree-bounded exact interpolation
(Leibniz bound on the 8x8 Sylvester determinant: deg_a <= 16,
deg_b <= 24, deg_c <= 13; 5,950-node tensor Lagrange grid, 23,800
Bareiss determinants; a determinant is a polynomial identity in its
entries, so the interpolation is itself a proof), agreeing
coefficient-by-coefficient with the direct certificate. The sampled
lead statistic `{1,2,4,8}` is fully explained (symbolic lead `8`;
per-node integer content absorbs 2-powers). New uniform structure from
the referee run: ALL three coordinate cubics obey one discriminant law
`disc = -4 (square)^2 L` (`disc_x = -4(27ac^2-9bc+8)^2 L`,
`disc_r = -4(27a)^2 L`, `disc_z = -4(27M)^2 L`, `M` explicit weight-5),
`F` is quasi-homogeneous with weights `(a,b,c,z) = (2,1,-1,2)`, and on
`L = 0` the z-core degenerates to `Z^2(8Z + A_z)`.
Referee: `04-computation/keller_z_eliminant_exact_opus_20260727.py` ->
`05-knowledge/results/keller_z_eliminant_exact_opus_20260727.out`.
