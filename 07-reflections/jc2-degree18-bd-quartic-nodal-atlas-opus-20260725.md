---
source: codex-2026-07-25-degree18-bd-quartic-nodes
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Within the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, the four roots of THM-2311's quartic B--D
  ratio polynomial are exactly the parameters for which the even spectral
  quotient is a geometrically irreducible nodal plane cubic.  The node is
  rational over the quartic ratio field and the quotient normalization is
  P1 over that field.  Restoring the discarded sign by y^2=x gives a
  double cover branched at three simple zeros and three simple poles of x,
  hence a geometrically irreducible genus-two original spectral curve.
  The genus/deck argument closes all four quartic ratios and, together
  with the two rational-ratio packets, empties the complete six-point
  B--D bank.  This remains a scoped degree-eighteen branch theorem, not
  JC(2) or DC(2).
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-bc-algebraic-finite-place-genus-floor-opus-20260725
  - jc2-degree18-remaining-two-sparse-finite-place-closure-opus-20260725
  - jc2-degree18-bd-4075-smooth-cubic-closure-opus-20260725
  - jc2-degree18-central-twall-closure-opus-20260725
script: 04-computation/jc2_degree18_bd_quartic_nodal_atlas_probe.py
output: 05-knowledge/results/jc2_degree18_bd_quartic_nodal_atlas_probe.out
script_sha256: c6a1d1a5b74f688568fdb63a7a12ce9b5fb24f69c165eccd3571927f58b945ed
output_sha256: c2dfd8dc957b21686aa65661d6545803aa3e78e1b28b4c9580c0d88d01a7d2ab
hash_basis: LF-normalized working-tree bytes
---

# The nodal B--D quotient lifts to genus two

## 1. Why this is the next decisive split

On the `B=1,C=W=0` chart, THM-2311 leaves six weighted ratios:

```text
D=25/126;

D=4075/85176;

and the four roots of

p(D)
 =22656250
  -772734375D
  +7600635000D^2
  -30805790400D^3
  +46376717184D^4.                                 (1)
```

The first ratio now lies on the canonically closed common-root wall of
THM-2345.  The work below is nonoverlapping: it closes the four
off-common-root quartic ratios by reconstructing the sign discarded by
the even quotient.  Together with THM-2345 and the separately audited
`4075/85176` packet, it empties this six-point bank.

The first ratio is removed by THM-2345's saturation of the full
common-root wall.  The second ratio has a smooth genus-one quotient and
is removed by the genus/deck argument.  The four roots in (1) behave in
the opposite way: they are precisely the parameters at which the
quotient cubic acquires one ordinary node.  Thus they are not four more
genus-one cases in disguise.

## 2. A coordinate in which the parameter is linear

Start with the even quotient coordinate

```text
x=y^2
```

and make the affine source change

```text
z=x-27u+120/7.                                     (2)
```

The spectral cubic becomes

```text
H(u,z;D)
 =27648000/7
  -907200z-37324800u+1959552Dz
  +20160z^2+2993760uz+88179840u^2
  +1127z^3-47628uz^2-3429216u^2z-61725888u^3.
                                                            (3)
```

The parameter occurs only in `1959552Dz`.  This is the useful hidden
normal form.

At a singular point `H_z=0`.  Put

```text
F=H-zH_z.                                           (4)
```

Then, on `H_z=0`, the equation `H=0` is equivalent to `F=0`.  It is
essential to retain the actual derivative `H_u`; differentiating `F`
instead would introduce a false elimination problem.

## 3. Complete singular-point elimination

Exact elimination in `u` gives

```text
primitive Res_u(F,H_u)
 =z^2 S(z),                                        (5)

S(z)
 =51840000+21168000z+123480z^3+55223z^4.           (6)
```

Both `S` and `p` from (1) are irreducible modulo `11`, hence irreducible
over `Q`; both are squarefree.

A linear subresultant of `F` and `H_u` has the form

```text
V(z)u-U(z).                                        (7)
```

The exact gcd control is

```text
gcd(S,V)=1.                                        (8)
```

Therefore on every root of `S`,

```text
u=U(z)/V(z).                                       (9)
```

Substitution in `H_z=0` gives

```text
D=N(z)/M(z),                  gcd(S,M)=1.           (10)
```

The parameter resultant is exactly

```text
primitive Res_z(S(z),N(z)-DM(z))=p(D).             (11)
```

Since the two irreducible polynomials in (6) and (1) both have degree
four, (10)--(11) identify their fields:

```text
Q[z]/(S)  isomorphic to  Q[D]/(p).                  (12)
```

Thus the node coordinates lie in the ratio field itself; no auxiliary
extension is needed.

The repeated factor `z^2` in (5) is also exact.  Its unique common
`u`-root is

```text
z=0,                    u=40/63,
```

and (10) specializes to

```text
D=25/126.                                           (13)
```

So the central point and the quartic bank exhaust the singular quotient
cubics in this one-parameter plane.  The smooth total-ramification point
`4075/85176` is correctly absent from (5).

## 4. Every quartic singularity is one ordinary node

At the point (9), consider the tangent Hessian determinant

```text
Theta=H_uu H_zz-H_uz^2.                            (14)
```

After clearing the square of `V` and reducing modulo `S`, its primitive
remainder is

```text
-59428800-20055000z-224273z^2+94325z^3.            (15)
```

The gcd of (15) with `S` is one.  Hence

```text
Theta!=0                                           (16)
```

at all four conjugate points.  The quadratic tangent cone has two
distinct lines over the algebraic closure, so each singularity is an
ordinary node.

There is exactly one such node on the cubic belonging to each root of
`p`.  Indeed, every affine singular point must appear in (5); `z=0`
forces the distinct central value (13); and the squarefree degree-four
resultant (11) makes the four roots of `S` map bijectively to the four
roots of `p`.

The homogeneous cubic at infinity is unchanged from the smooth-ratio
audit:

```text
1127-138915v+1607445v^2-26040609v^3,                (17)
```

with nonzero discriminant

```text
-153384762202971019112448.                          (18)
```

Thus the three points at infinity remain distinct and smooth.

## 5. Geometric irreducibility and rational normalization

A reducible reduced plane cubic is a line plus a conic or three lines.
Its component intersections are singular.  If a line and conic met only
at one point with total intersection multiplicity two, their tangents
there would coincide; that point would not have the nondegenerate tangent
cone (16).  A transverse ordinary node contributes only one to the
line--conic intersection number, so Bezout would force another singular
intersection.  Three lines force at least as much singular intersection
data.  A nonreduced cubic would have a positive-dimensional or
degenerate singular locus.

Because the projective cubic has exactly one ordinary node, all these
alternatives are impossible.  Therefore it is geometrically irreducible.

Its arithmetic genus is one and its node has delta-invariant one, so the
normalization has genus zero.  More concretely, the node is rational over
the quartic field (12), and the pencil of lines through it gives a
birational parameter

```text
t=(u-u_0)/(z-z_0).                                 (19)
```

Substitution into (3) removes the double intersection at the node and
solves the third intersection rationally in `t`.  Hence

```text
Normalization(H)=P^1 over Q[D]/(p).                 (20)
```

## 6. Restoring `y` raises the genus to two

The quotient normalization is rational, but the original trajectory did
not live merely on that quotient.  It retained

```text
y^2=x,

x=z+27u-120/7.                                     (21)
```

The divisor parity of `x` on the quotient normalization is therefore the
missing sidecar.

At `x=0`, equation (3) is the cubic

```text
-26040609u^3
+49601160u^2
+(-20995200-52907904D)u
+33592320D.                                        (22)
```

Its discriminant, as a polynomial in `D`, has primitive ascending
coefficients

```text
[-15625, 236250, -1190700, 2000376].                (23)
```

The gcd of (23) with `p(D)` is one.  Consequently, at every quartic
ratio, the fibre `x=0` consists of three distinct points.  In the
`(u,z)` coordinates, restricting to `x=0` means
`z=-27u+120/7`; hence the derivative of the restricted cubic is

```text
d/du H(u,-27u+120/7)=H_u-27H_z,                    (23a)
```

not the changed-coordinate partial `H_u`.  Nonvanishing of the fibre
discriminant says precisely that (23a) is nonzero at all three points.
Thus the curve meets `x=0` transversely, `x` is a local parameter with a
simple zero there, and in particular the affine node is not on `x=0`.

At infinity, the three roots in (17) are distinct.  The line at infinity
meets the cubic transversely at each, and `x` has a simple pole there.
Thus on the quotient normalization

```text
divisor(x):
  three simple zeros and three simple poles.        (24)
```

In particular `x` is not a square in its rational function field.  The
quadratic extension

```text
K(P^1)(y),                    y^2=x                  (25)
```

is geometrically irreducible and ramifies at exactly those six places.
Riemann--Hurwitz gives

```text
2g-2=2(-2)+6=2,

g=2.                                                (26)
```

This is the normalization of the original even-in-`y` spectral equation
`G(u,y)=H(u,y^2)=0`.

A hypothetical Keller trajectory therefore gives a rational map from
`P^1` to a genus-two curve and must be constant.  Hence `u` and `y` are
constant.  If `y!=0`, THM-2262 equation (11) makes `Z=T^2` constant,
then `T` and the nonzero deck coordinate `q` are constant, contradicting
the genuine deck involution.  If `y=0`, THM-2262 Section 4 already
eliminates the entire center by the third-flux/Keller/Faber argument.

Therefore all four roots of `p(D)` are empty in the retained branch.

## 7. What this proves and what it destroys

The six-point B--D bank is now classified by normalization type:

```text
D=25/126:
  common-root wall, canonically closed by THM-2345;

D=4075/85176:
  smooth genus-one cubic, closed by genus/deck;

p(D)=0:
  four conjugate irreducible nodal quotient cubics;
  restoring y gives genus two;
  closed by genus/deck.                             (27)
```

Thus the entire six-point `B`--`D` bank is empty.  The creative step is
that genus disappears under the quotient normalization and reappears
only after restoring the sign coordinate that the quotient discarded.

```text
source:
  THM-2311's quartic ratio factor;

map:
  expose the linear-D coordinate, eliminate the full singular equations,
  pass to the node-line normalization, then restore y^2=x;

preserved:
  the ratio field, the unique node, the rational spectral trajectory,
  the divisor of x, and the genuine deck;

destroyed by the first quotient:
  the sign of y and the six odd divisor places carrying the genus;

restoring sidecar:
  the Kummer divisor parity of x on the nodal normalization;

cheapest decisive next test:
  the finite-place packets now close every other two-sparse bank; move
  to the first genuinely three-sparse translation face.                 (28)
```

No claim about three- or four-coordinate strata, JC(2), or DC(2) is
made.

## 8. Exact reproduction

Run

```text
C:\Users\Eliott\.cache\codex-runtimes\codex-primary-runtime\
dependencies\python\python.exe
  04-computation/jc2_degree18_bd_quartic_nodal_atlas_probe.py
```

and repeat with `-O`.  The script verifies the coordinate change (2)--(3),
the full resultant (5), hostile interpolation values, the linear
subresultant, both denominator gcds, the exact field resultant (11),
mod-11 irreducibility, the central specialization (13), and the tangent
Hessian gcd (15)--(16).  It also verifies the exact zero-fibre
discriminant (23), its coprimality with `p`, the three distinct infinity
points, and the six-branch genus-two count.

Stored transcript:

```text
05-knowledge/results/jc2_degree18_bd_quartic_nodal_atlas_probe.out
```

Normal, optimized, and stored transcripts match.  LF-byte SHA-256:

```text
script  c6a1d1a5b74f688568fdb63a7a12ce9b5fb24f69c165eccd3571927f58b945ed
output  c2dfd8dc957b21686aa65661d6545803aa3e78e1b28b4c9580c0d88d01a7d2ab
```
