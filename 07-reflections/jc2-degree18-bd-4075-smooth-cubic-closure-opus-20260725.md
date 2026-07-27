---
source: codex-2026-07-25-degree18-bd-4075-closure
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Within the genuine nonsplit polynomial exact-square-prefix degree-eighteen
  branch of THM-2262/2297, the remaining rational B--D point
  D/B^2=4075/85176 has a smooth irreducible genus-one spectral quotient and
  is empty. This argument alone does not close the quartic B--D ratio
  bank; the later nodal/Kummer packet restores y^2=x and closes that bank
  by an independently audited genus-two argument. Other weighted planes,
  split/even descent, JC(2), and DC(2) remain open.
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2313-degree-eighteen-bd-linear-ratio-closure
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-bc-algebraic-finite-place-genus-floor-opus-20260725
  - jc2-degree18-remaining-two-sparse-finite-place-closure-opus-20260725
  - jc2-degree18-central-twall-closure-opus-20260725
  - jc2-degree18-bd-quartic-nodal-atlas-opus-20260725
script: 04-computation/jc2_degree18_bd_4075_closure_probe.py
output: 05-knowledge/results/jc2_degree18_bd_4075_closure_probe.out
script_sha256: 0992e7f383a74ecf8ed452cf16640919b31f2097c8a97c5f71aa62015daf2b62
output_sha256: 9d4f82a8170315960002a1bdaa328ed6463f893e726e8e2eb643bdd70f17f63d
hash_basis: LF-normalized working-tree bytes
---

# The remaining rational B--D point is a smooth plane cubic

> **FOLLOW-UP B--D BANK CLOSURE.** THM-2345 canonically closes the first
> rational point on the common-root wall.  This note closes the second.
> The later exact packet
> `07-reflections/jc2-degree18-bd-quartic-nodal-atlas-opus-20260725.md`
> identifies the quartic points as nodal quotients and restores the
> discarded sign to obtain genus-two original curves.  Its hostile audit
> has passed, so the quartic points are no longer a third-flux obligation.

## 1. Scope and inherited obstruction

THM-2311 reduces every two-sparse `B`--`D` survivor to six weighted
ratios. Its two rational linear factors are

```text
D/B^2=25/126,                  D/B^2=4075/85176.      (1)
```

The first ratio is already empty on its full `C,W` extension by
THM-2345. This note closes the second
ratio on the actual two-sparse chart

```text
B=1,             C=W=0,             D=4075/85176.    (2)
```

Here `B=1` is a constant weighted normalization: rescaling by a fixed
nonzero square root of `B` gives an isomorphic spectral curve and
preserves constancy of a Keller trajectory.  After the geometric
argument, the first-flux equation is read back in the original
coordinates, so this normalization does not identify a nonconstant
trajectory with a constant one.

The load-bearing correction is to quotient the even `y`-symmetry before
interpreting the repeated-root resultant. Multiplicity in
`Disc_y(delta(y^2))` mixes genuine branch collision with ramification of
the base change `y -> x=y^2`; it is not by itself a singularity of the
spectral normalization.

## 2. The exact quotient cubic

Put

```text
x=y^2.
```

The THM-2262/2297 spectral equation becomes the total-degree-three plane
curve

```text
H(u,x)
 =-26040609u^3
  +(49601160+1607445x)u^2
  +(-20995200-52907904D-2857680x-138915x^2)u
  +33592320D+(777600+1959552D)x+78120x^2+1127x^3,

D=4075/85176.                                      (3)
```

Let

```text
delta(x)=Disc_u H(u,x).                              (4)
```

Exact rational Euclidean division gives

```text
deg delta=6,

gcd(delta,delta')=x+1215/91,                         (5)

delta=K(x+1215/91)^2 Q_4(x),             K!=0,
```

where the primitive ascending coefficients of `Q_4` are

```text
[-3037628182500,
  1389005982000,
  -147517112925,
  203464170,
  1577224103].                                      (6)
```

Moreover,

```text
gcd(Q_4,Q_4')=1,

Q_4(-1215/91)!=0.                                   (7)
```

Thus there is one double discriminant value and four further simple
values. There are no hidden additional collisions.

## 3. The double discriminant value is smooth total ramification

At

```text
x_0=-1215/91,                  u_0=295/819,           (8)
```

the whole fibre is exactly

```text
H(u,x_0)=-26040609(u-u_0)^3.                         (9)
```

Nevertheless the curve is smooth there, because

```text
H_x(u_0,x_0)=-16329600/169!=0.                       (10)
```

Locally `x-x_0` is therefore a unit times `(u-u_0)^3`; the projection to
the `x`-line has ramification index three and contributes two to the
ramification divisor. This explains the square in (5) without creating a
singularity.

The leading coefficient of `H` as a cubic in `u` is the nonzero
constant `-26040609`.  Consequently the standard local
discriminant argument applies without a disappearing-leading-term
exception: at an affine singular point, the fibre has a multiple `u`
root and the linear base term of the corresponding double-root local
form vanishes (a triple root lies on the singular locus of the universal
discriminant), so its `x`-value is a multiple zero of `delta`.  The only
candidate is (8), and (10) excludes it. The four roots of `Q_4` are
simple branch values and contribute one each.

## 4. Infinity is smooth and unramified

The degree-three homogeneous part of (3), in the slope `v=u/x`, is

```text
L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (11)
```

Its discriminant is

```text
-153384762202971019112448!=0.                        (12)
```

Hence the projective completion has three distinct smooth points at
infinity. At each point `x` has a simple pole, so all three branches are
unramified for the degree-three map to the `x`-line.

The projective plane cubic is therefore smooth. A reducible projective
plane cubic over the algebraic closure would have intersecting components
and hence a singular point, so the curve is geometrically irreducible.
Equivalently, Riemann--Hurwitz sees

```text
Ram=2+4=6,

2g-2=-6+6=0,                    g=1.                 (13)
```

## 5. Keller contradiction

The rational functions `(u,x=y^2)` of a hypothetical Keller trajectory
define a rational map from `P^1` to the smooth projective curve (3).
Positive genus forces this map to be constant. Hence `x` and `u` are
constant.

If `x!=0`, then `y` is constant in the algebraically closed constant
field. THM-2262's first-flux square makes `T^2` constant, hence `T` and
then the nonzero deck coordinate `q` constant; the genuine deck fixes the
base constants and sends `q` to `-q`, a contradiction.

If `x=0`, then `y=0`, and THM-2262 Section 4 excludes the entire
exceptional center by the retained third flux, Keller one-form,
rational-primitive lemma, and the uncancellable whole-polynomial Faber
sidecar.

Therefore

```text
Within the genuine nonsplit polynomial exact-square-prefix branch,
there is no degree-eighteen survivor with

           C=W=0,        B!=0,        D/B^2=4075/85176.   (14)
```

Together with THM-2345, this closes both rational linear-factor points in
THM-2311's `B`--`D` bank.  The four roots of its remaining quartic ratio
polynomial are left open by the present genus-one argument; the
independently audited nodal/Kummer packet above gives their separate
genus-two closure.

## 6. Information ledger and exact reproduction

```text
source:
  the even-in-y two-sparse spectral cubic at the second rational B--D
  ratio;

map:
  quotient by x=y^2, then project the plane cubic to the x-line;

preserved:
  the rational Keller trajectory, branch multiplicities, the exceptional
  y=0 alternative, and the genuine nonsplit deck;

destroyed:
  the sign of y;

sidecar:
  x!=0 recovers y as a constant after the genus argument, while x=0 is
  discharged by THM-2262's third-flux/Faber center theorem;

hostile boundary:
  the double zero of delta is smooth e=3 ramification, showing why
  resultant multiplicity must not be called a curve singularity;

next target:
  the finite-place packets now close every other two-sparse bank; move
  to the first genuinely three-sparse translation face.                 (15)
```

Run

```text
python3 04-computation/jc2_degree18_bd_4075_closure_probe.py
python3 -O 04-computation/jc2_degree18_bd_4075_closure_probe.py
```

Both executions must match

```text
05-knowledge/results/jc2_degree18_bd_4075_closure_probe.out
```

after LF normalization.
