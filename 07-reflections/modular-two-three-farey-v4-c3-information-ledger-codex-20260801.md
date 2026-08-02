# The same `2/3` grammar appears three times, but the intertwiner is still missing

**Synthesis, 2026-08-01.** This note separates three exact appearances of
the modular presentation

```text
PSL_2(Z) = <S,R | S^2=R^3=1> = Z_2 * Z_3
```

from the tempting but presently unproved claim that they are one physical
LRC/JC mechanism. The common object is an **information ledger**: an
order-two move changes a flank or an origin, an order-three move cycles a
triangle or a nonzero direction, and forgetting the order-two coordinate
leaves exactly the quotient on which the current quartic and LRC obstructions
stall.

## 1. The literal modular occurrence: THM-2056's Farey fan

`01-canon/theorems/THM-2056-kelvin-polar-farey-defect-certificate.md`
works with primitive lattice rays and acute unimodular cones
`|det(u,v)|=1`. These are the edges of the Farey tessellation. Its two
standard local views are:

- the binary Stern--Brocot subdivision of an edge by the mediant `u+v`;
- the trivalent dual tree whose vertex is the Farey triangle
  `(u,v,u+v)`.

The order-two generator reverses a Farey edge, while the order-three
generator cycles a Farey triangle. Thus the binary fraction tree and the
ternary triangle tree really are the two Bass--Serre faces of one modular
object here.

The load-bearing gate, however, is owner-dependent:

```text
F_p(m u+n v)
 =mF_p(u)+nF_p(v)
  +(m^2-m)||u||^2+(n^2-n)||v||^2+2mn u.v.             (1)
```

The determinant/Farey structure is modular, but the polar owner `p`, the
Euclidean norm, and the acuteness inequality are extra data. A triangle
rotation does not preserve `(1)` unless those data are transported too.
This is the first precise reason that "PSL2 symmetry" alone does not prove an
LRC gate.

## 2. The finite quartic shadow: affine points and three directions

On four labelled points identify the point torsor with `V_4=F_2^2`. Its
three nonzero directions are cycled by a `C_3` inside
`GL(2,2)=S_3`, and

```text
AGL(2,2)=V_4 semidirect S_3 is isomorphic to S_4,
S_4/V_4 is isomorphic to S_3.                            (2)
```

The exact referee
`07-reflections/modular-v4-c3-two-prime-point-direction-decoder-codex-20260801.md`
computes the two integral Fourier defects:

```text
point Walsh:      SNF(1,2,2,4), determinant 16;
direction C3:     determinant (1-zeta_3)^3, norm 27;
tensor index:     2^24 3^12.                            (3)
```

It reconstructs the twelve directed point--direction flags after choosing a
point origin and a `C3` orientation. Forgetting the origin is the `V4`
quotient in `(2)`; forgetting the orientation makes reflections semilinear.
The six undirected quartic edges are a two-to-one quotient of these twelve
flags. This is a finite quotient of the modular `2/3` grammar, not a proof
that the full free product acts on the physical carrier.

`01-canon/theorems/THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss.md`
makes the finite quotient exact.  The twelve flags are a regular
`A4=V4 semidirect C3` bitorsor.  Right edge reversal `S` and direction
rotation `R` satisfy `S^2=R^3=(SR)^3=1`; forgetting the point has kernel
`V4`, while forgetting edge orientation gives six edges on which `R` no
longer descends.  The odd involution that reverses the `C3` orientation is a
second sidecar, not the already-used translation `S`.

The identity

```text
disc(quartic)=disc(resolvent cubic)                      (4)
```

is now typed correctly: it retains the common sign double cover, but it does
not restore the lost affine point torsor or a branchwise cofactor.

## 3. The exact JC loss: product character versus pointed contrast

`01-canon/theorems/THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor.md`
inserts branch cofactors `c_i` by the natural vertex gauge
`B_ij=c_i c_j(z_i-z_j)`. Every perfect matching then receives the same
factor

```text
C=c_0c_1c_2c_3.                                        (5)
```

Hence the entire matching triple retains only the trivial multiplicative
character. Its invisible fibre is the three-dimensional product-one torus.
Keller constancy lives on that fibre because it asks for four separate
values `J_i=c_i f'(z_i)` to agree.

On a `C3 + fixed sheet` splitting, the unique `C3`-invariant missing line is

```text
(1,1,1,-3),       equivalently R_J=J_C3/J_fixed.         (6)
```

The exact hostile `(2,2,2,1/8)` moves along `(6)`, preserves `(5)` and every
matching jet, and breaks Keller equality.
`01-canon/theorems/THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md`
gives the pointed repair once the fixed sheet and true cofactor are supplied:

```text
Delta_J = Norm_(L/K)(R_J-1),
Delta_J=0 iff R_J=1                                    (7)
```

for a transitive cubic field `L/K`. Tame cubic ramification simultaneously
forces inverse-different cofactor valuation `-2`. Thus the ternary orbit is
not the missing datum by itself; the missing datum is its contrast with a
marked point.

## 4. The reflected LRC gain matroid

The cap-`5/2` addendum in
`01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md`
supplies an independent arithmetic shadow.
Represent a positive rational gain by its prime-exponent vector. For the
phase-zero alphabet

```text
4/3 -> ( 2,-1, 0),      3/2 -> (-1,1,0),
2   -> ( 1, 0, 0),      5/2 -> (-1,0,1),               (8)
```

the first three vectors have rank/nullity `(2,1)` because

```text
(3/2)(4/3)=2.                                           (9)
```

Adjoining `5/2` changes rank/nullity to `(3,1)`: a new prime coordinate
enters but no new circuit does. Consequently a unique `5/2` edge in a
realized zero-gain component must be a bridge. At the next cap, `3` lies in
the old span through `(3/2)2=3`, so nullity can increase instead. This is the
general bridge/circuit rule:

```text
new exponent direction outside the old span -> rank / bridge payload;
new direction inside the old span            -> nullity / circuit payload. (10)
```

The exact finite branch census reports the predicted bridge in all `65`
forced-full components.  MISTAKE-347 separates that structural fact from the
former cone claim: the split analytic tails covered one physical level
ordering twice, so `3m>=2D` is `AUDIT-REQUIRED`.  The last unaffected proved
cone is `3m>=4D`, leaving the `561`-body certificate-failure wedge
`D>=6,1<=m<4D/3`, not a physical-survivor census.

The audit-required cap-`3` scout makes the comparison with THM-3067 sharper.
Its `55` forced components contain the scale `K4`
`{1,3/2,2,3}={1,x,y,xy}`.  The vertices may be indexed by `F_2^2`, but the
gain decoration does not descend through point-forgetting: the two
`x`-matching edges both have gain `3/2`, the two `y`-matching edges both have
gain `2`, while the diagonal matching has gains `3` and `4/3`.  Their ratio
`9/4=x^2` is an exact origin defect.  Thus the tetrahedral incidence object
is visible in the finite CSP, but the `A4` action is not a gain symmetry and
the unfinished analytic cover cannot be supplied by quotienting the point.

## 5. The missing commuting square

The three appearances become one proof mechanism only after constructing a
map with the following ledger:

```text
source:       THM-2056 oriented Farey cone plus its polar owner;
target:       a physical V4 point torsor with a marked C3 direction/orbit;
map:          still OPEN;
must retain:  determinant sign, owner inequality, carrier positivity,
              pointed fixed sheet, and normalized cofactor contrast;
may forget:   only a proved kernel annihilated by the target predicate.
```

The cheapest decisive test is local: on every candidate Farey triangle,
transport the owner and evaluate whether its three flank defects determine
the pointed contrast `(6)`. If two triangles have the same `S3` direction
data and discriminant but different point origins or `R_J`, then the modular
bridge fails exactly as THM-3066 predicts. If the contrast is preserved,
`(7)` supplies the first genuine ternary decoder.

So the bold reframe survives in a precise form:

> the binary and ternary trees are two views of one modular incidence object,
> but current proofs retain only its `S3` direction quotient. The live
> sidecar is the affine point/origin contrast, not another scalar invariant.

This explains the repeated special role of `2` and `3` without claiming an
unconstructed PSL2 action, a Keller exclusion, or an LRC decrement.
