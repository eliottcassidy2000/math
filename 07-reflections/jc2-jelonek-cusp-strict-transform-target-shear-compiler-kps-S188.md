# The Jelonek cusp turns nonlinear target shears into a strict-transform curve compiler

**Status: DERIVED-EXACT identities + FINITE-EXACT diagnostic; OPEN construction
program.**  This note does not claim a planar counterexample and does not
reserve a theorem identifier.  It sharpens the nonlinear coordinate-surface
cell left by THM-3558--3560.

Companion:
[`jacobian_target_graph_cusp_strict_transform_kps_s188.py`](../04-computation/jacobian_target_graph_cusp_strict_transform_kps_s188.py),
with matching output
[`jacobian_target_graph_cusp_strict_transform_kps_s188.out`](../05-knowledge/results/jacobian_target_graph_cusp_strict_transform_kps_s188.out).

## 1. Inheritance pass

- **Closest proved mechanism:** THM-3560 integrates the exact `3/1/0` fibre
  stratification and gives the necessary coordinate-pullback equation
  `2 chi(D)+chi(e)=2`.
- **Canonical hostile:** `c=0` passes that Euler equation but pulls back to
  `F3=xC`, one affine plane plus the punctured Kummer surface of THM-3554.
  Passing Euler is necessary, never sufficient.
- **Corrected near miss:** THM-3559 closes affine target coordinates, and
  THM-3560 closes all monomial shears `c+lambda b^m`; neither statement
  classifies a mixed polynomial `phi(a,b)` or a coordinate factor of a
  reducible pullback.
- **Least-used sidecar:** THM-1335's cube-plus-square identity is not merely a
  description of the ambient Jelonek cusp.  It normalizes every target-graph
  section uniformly.

The live board is:

1. the target graph `T_phi: c+phi(a,b)=0`;
2. its Jelonek section `D_phi=T_phi intersect V(L)`;
3. the omitted-curve support `e_phi=T_phi intersect E`;
4. the strict-transform curve `Gamma_phi` below;
5. source coordinate status of the complete pullback and of each factor.

## 2. Universal strict-transform equation

Use the exact THM-1335 identity

```text
108a^2 L=P^3+E^2,
P=12a-b^2,
E=54a^2c-18ab+b^3.                                   (1)
```

Normalize the cusp `P^3+E^2=0` by

```text
P=-r^2,                      E=r^3.                    (2)
```

Then

```text
a=(b^2-r^2)/12.                                        (3)
```

On the graph `c=-phi(a,b)`, exact substitution and factorization give

```text
E-r^3
=-(b-r)^2/2 [ b+2r+(3/4)(b+r)^2
                    phi((b^2-r^2)/12,b) ].             (4)
```

The square `(b-r)^2` is the base component `a=0` created when `(1)` is
multiplied by `a^2`; it is not automatically a component of `D_phi`.
Dividing it produces the strict-transform equation

```text
G_phi(b,r)
=b+2r+(3/4)(b+r)^2 phi((b^2-r^2)/12,b)=0.             (5)
```

With `s=b+r`, this is the compact compiler

```text
a=s(2b-s)/12,

G_phi(b,s)
=2s-b+(3/4)s^2 phi(s(2b-s)/12,b)=0.                   (6)
```

Thus the mixed-target search is no longer an unstructured surface search.
For each `phi`, it is a plane-curve normalization/intersection calculation.
Conversely, one can try to design a curve in `(b,s)` and ask whether the
required value

```text
phi|Gamma=4(b-2s)/(3s^2)                               (7)
```

descends through the quadratic subring

```text
C[a,b]=C[s(2b-s)/12,b].                                (8)
```

Equation `(8)` is the new integrality gate.  It is the target-shear analogue
of the cusp-packet question whether a source-polynomial Bezout coefficient
descends to `C[L,T,U,S]`: rational recovery is cheap; polynomial descent is
load-bearing.

## 3. The omitted curve becomes the cusp origin

The omitted curve is

```text
E(t)=(4/(27t^2),4/(3t),t),               t in C*.     (9)
```

Its projection satisfies `12a-b^2=0`, so it is the `r=0` locus in `(2)`.
The strict-transform equation gives exactly

```text
G_phi(b,0)
=b+(3/4)b^2 phi(b^2/12,b)
=(3/4)b^2 [4/(3b)+phi(b^2/12,b)].                     (10)
```

Therefore the roots of `(10)` away from `b=0` are precisely `e_phi`.
The parity gate in THM-3560 becomes a visible root-support parity condition
on the `r=0` slice of one curve.

This is a useful design rule:

- avoiding `E` means `G_phi(b,0)` has no nonzero root;
- over `C`, that forces `1+(3/4)b phi(b^2/12,b)` to be a monomial unit on
  `C*`;
- the simplest way is `phi` divisible by `12a-b^2`, but its resulting
  strict-transform curve usually acquires hyperelliptic topology.

## 4. Why the ideal section is rational and why completion is hard

The formally ideal choice sets `E=0` on the target graph:

```text
phi_ideal=(b^3-18ab)/(54a^2).                          (11)
```

Then `(1)` collapses to

```text
L|T=P^3/(108a^2),                                      (12)
```

so the Jelonek section is the rational parabola `P=0` with triple contact.
This is exactly the topology one wants: `chi(D)=1`, and the ideal graph can
avoid the omitted curve.  But `(11)` has an unavoidable `a^(-2)` pole.  The
denominator is not cosmetic; it is the polynomial-completion problem.

This explains the incoming Danielewski architecture reserved in THM-3561.
The source, target, and map are:

```text
source: rational graph with an a^(-2) pole,
target: a polynomial Danielewski completion,
map:    retain the rational Keller collision while adding the missing divisor.
```

The preserved predicate is the constant tangential Jacobian and collision.
The destroyed information is affine-plane topology: the completion tends to
acquire nontrivial Picard/Euler data.  The needed sidecar is therefore not
another local Jacobian calculation but a divisor-class or unit cancellation
that makes the completed surface genuinely `A2`.

## 5. Quadratic hostile and what remains

As a cheap finite diagnostic, the companion exhausts every quadratic target
graph over `F_7`:

```text
g=c+l a+m b+A a^2+B ab+C b^2.
```

Exactly `19` parameter rows have uniformly `7^2` points in every fibre, and
all `19` have `(A,B,C)=(0,0,0)`.  No genuinely quadratic row survives.
This is **FINITE-EXACT over `F_7` only**.  It is not a characteristic-zero
no-go: a hypothetical complex coordinate and its inverse may have bad
reduction at seven, and coefficients need not be rational.

The sharp next targets are:

1. classify degree-two `phi` over characteristic zero through the curve
   `(6)`, not through finite-field balance alone;
2. search degree-three mixed shears satisfying
   `2chi(D_phi)+#e_phi=2`; pure `b^3` is already excluded by THM-3560;
3. test whether a reducible pullback has a coordinate factor even when its
   complete Euler characteristic fails;
4. formulate descent in `(8)` as an integral-closure/conductor problem and
   compare it directly with the cusp-square packet's proper birational
   subring and THM-3561's Danielewski completion.

The creative change of viewpoint is precise: do not guess a planar Keller
pair.  Design a target curve `Gamma_phi` with the right Euler/omitted support,
then solve the quadratic-subring descent `(8)`, and only then test whether the
source pullback is a coordinate plane.
