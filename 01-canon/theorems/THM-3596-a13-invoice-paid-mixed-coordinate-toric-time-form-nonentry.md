---
id: THM-3596
title: "A13 invoice-paid mixed-coordinate toric time-form nonentry"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / PENDING INDEPENDENT AUDIT.
  On the A13 three-arm surface c^3 e=b(b^2+1), an explicit one-parameter
  family of genuinely mixed coordinates (five weights when lambda is
  nonzero) has no rational constant-bracket mate.  The member lambda=3/2 is
  the fully arm-repaired
  coordinate from THM-3591.  Its time form is the interior Newton-polygon
  differential (1,-2); the sole degenerate face resolves into the three
  smooth affine arms, where the form is a unit.  No Darboux pair and no
  counterexample to JC(2) is constructed.
source: root / planar-Jacobian mixed-coordinate construction-hostile session, 2026-08-21
depends_on:
  - THM-3591-danielewski-arm-blind-separated-darboux-nonentry
related:
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
  - THM-3589-danielewski-central-arm-every-line-and-kummer-trace-darboux-gates
  - THM-3595-danielewski-separated-transverse-time-form-rational-nonentry
script: 04-computation/jc2_a13_mixed_toric_timeform_nonentry_thm3596.py
output: 05-knowledge/results/jc2_a13_mixed_toric_timeform_nonentry_thm3596.out
script_sha256: 0c8a9ae72dea21daf7d04a23527fd53717279f7055d322fd740e0b01dc05d4e1
output_sha256: 0bdc43c1e7816a770c7a664ec5737550abfc219963c2464e23dd213ba926dd92
hash_basis: raw LF bytes
---

# THM-3596 -- A13 invoice-paid mixed-coordinate toric time-form nonentry

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / PENDING INDEPENDENT
AUDIT.**  This theorem attacks the most developed explicit mixed coordinate,
not merely another bounded support box.  It proves that even an arbitrary
rational mate cannot repair that coordinate.  It does not exclude other
mixed coordinates or prove any case of `JC(2)`.

Work over `C`.  Put

```text
Sigma=b(b^2+1),
q=(225b^4+237b^2+8)/8,              v=105b/16,          (1)

Y: c^3 e=Sigma(b),
{b,c}=c^3,       {c,e}=-Sigma'(b),  {b,e}=-3c^2e.      (2)
```

For every `lambda in C`, define the genuinely mixed coordinate

```text
Q_lambda=e^3+q(b)e+c^2+v(b)c^4+lambda ce.              (3)
```

Then there are no `P in C(Y)` and `kappa in C*` such that

```text
{P,Q_lambda}=kappa.                                    (4)
```

The special value `lambda=3/2` is exactly `Q_ddagger` from THM-3591:
paired there with `P=e^2+c+c^3`, it pays both central-arm curves, the
labelled inverse-derivative arm interpolation, and the global scalar bracket
layer.  The present theorem shows that its remaining debt is not a bad choice
of polynomial mate: no rational mate exists at all.

## 1. The forced Laurent time form

On `c!=0`, eliminate `e=Sigma/c^3`.  The coordinate `(3)` becomes

```text
F_lambda(b,c)
 =Sigma^3 c^-9+q Sigma c^-3+c^2+v c^4+lambda Sigma c^-2.  (5)
```

The Poisson bracket is the standard Laurent Jacobian with density `c^3`:

```text
{R,S}=c^3(R_b S_c-R_c S_b).                            (6)
```

Let `w` be transcendental.  The generic fibre is integral: it is the
localization of the domain `C[Y]` from the inclusion `C[Q_lambda] subset
C[Y]`.  Let `C_w` be its smooth projective normalization.  On every smooth
affine chart the unique Hamiltonian time
form is

```text
eta=db/(c^3 F_c)=-dc/(c^3 F_b),             eta(D_F)=1. (7)
```

Thus `(4)` would imply

```text
dP=kappa eta.                                            (8)
```

on every generic component.  Exact computation gives

```text
gcd(16c^9 F_b,8c^10 F_c)=1 in C(lambda)[b,c],           (9)
```

so the Laurent presentation also exposes no fixed critical divisor.

## 2. The Newton polygon puts the time form strictly inside

The exponent support of `F_lambda-w` is

```text
(3,5,7,9)x{-9},   (1,3,5,7)x{-3},
{(1,-2),(3,-2),(0,2),(1,4),(0,0)},                     (10)
```

where the two `y=-2` points disappear when `lambda=0` without changing the
convex hull.  In counterclockwise order the hull vertices are

```text
(3,-9),(9,-9),(7,-3),(1,4),(0,2),(0,0).                (11)
```

For a Laurent plane curve, the residue differential attached to a lattice
point `(i,j)` is

```text
b^(i-1)c^(j-1) db/F_c.                                 (12)
```

The form `(7)` is therefore the point

```text
p=(1,-2).                                               (13)
```

Its primitive lattice distances from the six edges of `(11)` are

```text
7,17,36,6,1,1.                                         (14)
```

In particular `p` is strictly interior.  Toric adjunction now makes `eta`
holomorphic at every nondegenerate toric end; its boundary orders there are
the corresponding distances minus one.

The exact face polynomials, with monomial units suppressed, are

```text
bottom:  b^3(b^2+1)^3,
BC:      b^2+(225/8)c^6,
CD:      (225/8)b^6+(105/16)c^7,
DE:      1+(105/16)bc^2,
EF:      c^2-w,
FA:      u^3+u-w,                         u=bc^-3.      (15)
```

For generic `w`, every face except the bottom one is nondegenerate.  The
bottom face has only the two torus roots `b=+/-i`, each triple.  These are not
new infinity branches: they are the collapsed affine arms of `Y`.

## 3. Restoring the arm coordinate removes the apparent degeneracy

Let `beta` be any root of `Sigma`.  Because `Sigma'(beta)!=0`, `(c,e)` are
smooth local coordinates on `Y` near

```text
D_beta={b=beta,c=0}.                                    (16)
```

On that arm,

```text
Q_lambda|D_0=e^3+e,
Q_lambda|D_(+/-i)=e^3-e/2.                             (17)
```

The term `lambda ce` vanishes, so this calculation is uniform in `lambda`.
For transcendental `w`, all roots of the two cubics in `(17)` are simple.
On the branch through such a root `e_0`, `c` is a local parameter and the
time form is

```text
eta=-dc/[Sigma'(b) partial_e Q_lambda],                 (18)
```

which is a unit because

```text
Sigma'(0)=1,          Sigma'(+/-i)=-2,
partial_e Q|D_0=3e_0^2+1,
partial_e Q|D_(+/-i)=3e_0^2-1/2.                       (19)
```

Thus the bottom triple roots of `(15)` merely record that the Laurent chart
forgot `e`; after restoring that sidecar, every resulting point is smooth and
`eta` has order zero.  There is no hidden divergent-`e` branch: if
`b-beta` has order `r` in `c`, then `e` has order `r-3`; for `r<3`, `e^3`
is uniquely dominant in `Q_lambda-w`, while for `r>3`, the constant `-w` is
uniquely dominant.  Hence necessarily `r=3`, exactly the finite arm chart
above.  The root `b=0` lies at the `FA` corner and is already resolved by
`u=e+O(c^3)`; the same local calculation covers it.

## 4. Nonexactness

The generic fibre is smooth on the affine surface.  Sections 2--3 show that
`eta` extends to a nonzero holomorphic differential on its smooth projective
normalization.  In
characteristic zero, the differential of a rational function has a pole at
every pole of that function.  A rational function without poles on a
projective curve is constant.  Hence a nonzero holomorphic differential is
not exact, contradicting `(8)`.

This proves `(4)` for every `lambda`, including the fully invoice-paid value
`lambda=3/2`.

## 5. What this challenges

The coefficient `3/2` was introduced to repair labelled arm interpolation.
It contributes only the two **interior** monomials `(1,-2),(3,-2)` and changes
no Newton face.  Therefore it can repair the affine arm ledger without
changing the global first-kind obstruction.  Boundary interpolation and
global exactness are independent bills.

```text
source       the five-weight mixed coordinate Q_lambda on Y
target       the projective normalization of F_lambda=w
map          c!=0 Laurent chart followed by toric compactification
preserved    {P,Q}=kappa becomes dP=kappa eta
destroyed    the arm coordinate e at c=0
sidecar      the three smooth (c,e) arm charts
obstruction  p=(1,-2) is interior; eta is nonzero holomorphic
cheap test   six edge distances plus the two arm cubics in (17)
```

This theorem does **not** give a target-wide support floor.  It closes one
particularly well-calibrated coordinate, even against rational mates.  A
surviving construction must change a Newton face or force a controlled
meromorphic time form; adding only interior correction monomials cannot do
either.

## 6. Exact verification contract

The companion checks:

- the Laurent identity `(5)`, Poisson contraction `(6)--(8)`, and gcd `(9)`;
- the exact support, hull, primitive inequalities, and distances `(10)--(14)`;
- all six face identities and their generic nondegeneracy controls `(15)`;
- the arm residues, inverse-function expansions, and unit gates `(17)--(19)`;
- symbolic independence of every boundary invoice from `lambda`;
- `lambda=0` and `lambda=3/2` as active hostile/positive controls.

The proof is exact and universal in `lambda`; finite symbolic rows are
controls, not extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_a13_mixed_toric_timeform_nonentry_thm3596.py
python3 -O 04-computation/jc2_a13_mixed_toric_timeform_nonentry_thm3596.py
```
