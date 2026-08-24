---
id: THM-3918
title: "Rational moving-triple-root one-place decic and anchored color transition"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At rho^2=-6, a sparse
  tangent deformation of THM-3913 gives a primitive normal globally
  nonmonogenic S3 cubic order whose absolutely irreducible discriminant is a
  rational decic with one smooth infinity place.  Its complete finite packet
  is an eight-address origin of delta 34 and two external A2 cusps.  The
  quadratic resolvent still has constant units and nonzero class-group
  three-torsion.  The rationalization is exactly a simple tie between one
  conjugate radial color and the affine-infinity anchor; the unanchored
  two-color gain is constant and misses the transition.  This removes the
  elliptic/Jelonek obstruction of THM-3913 but does not construct a plane
  atlas or Keller map; JC(2) remains open.
source: root + jc_rational_deformation + jc_tournament_response, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  One route eliminated the repeated-root incidence to
  the genus-one family 18y^2=(w+rho)(w^2-rho^2-6) and normalized its nodal
  rho^2=-6 fibre.  A second route worked directly from the polynomial
  parametrization, proving the exact critical and collision factors and the
  complete 34+1+1 genus ledger.  A separate response audit factored the
  radial squareclass over k(sqrt(-6)) and proved the anchored simple-tie and
  unanchored no-go statements.  A separate referee reconstructed the
  Delone--Faddeev normality bridge, birationality, first-blowup contact-seven
  calculation, complete genus ledger, unit and Kummer arguments, and exact
  scope.  The canonical companion freezes 45 active gates; normal and
  optimized runs byte-match the frozen LF transcript.  No repair was needed.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3912-even-one-place-split-boundary-a2-three-torsion-design-sieve
related:
  - THM-3913-moving-triple-root-one-place-decic-normal-nonmonogenic-cubic
  - THM-3914-decic-boundary-three-class-degree-one-isotropic-divisor
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
  - THM-3916-positive-genus-collapsed-valuation-keller-obstruction
  - THM-3917-quintic-parameter-rational-collapsed-cubic
script: 04-computation/jc2_rational_moving_triple_root_decic_thm3918.py
output: 05-knowledge/results/jc2_rational_moving_triple_root_decic_thm3918.out
script_sha256: c8fa01f5c50f7dff6ffd31c844c418cc94fbdec7ed5c3c009e5fbb992493ae73
output_sha256: 426819660702d48c2a11b65268604c3a7198b6456193b996f5d0317c867537a9
semantic_sha256: 89d595a4a896120d32ad58b0db0a0357e673e8df466d59e85503b35183ef6d1b
hash_basis: raw LF bytes
---

# THM-3918 -- the elliptic moving-root decic has a rational boundary fibre

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an algebraically
closed field `k` of characteristic zero.  Fix `rho in k` with

```text
rho^2=-6.                                                   (1)
```

Put `R=k[A,C]`, `ell=AU+CV`, and

```text
Phi_rho=ell^3+C ell U^2+rho C ell V^2+A^2V^3.              (2)
```

Let `T_rho/R` be the Delone--Faddeev cubic algebra of `(2)`.  Then:

1. `T_rho` is a normal finite-flat domain of rank three, is globally
   nonmonogenic, and has generic Galois group `S3`;
2. its discriminant curve `Gamma_rho` is an absolutely irreducible decic
   with one smooth projective infinity point;
3. the normalization map is polynomial of coordinate degrees `(9,10)` and
   identifies the affine normalization with `A1`;
4. the complete finite singular packet is one eight-address origin of
   `delta=34` and two ordinary `A2` cusps;
5. for the affine quadratic resolvent

   ```text
   Q_rho=Spec R[W]/(W^2-Delta_rho),                         (3)
   ```

   one has `Q_rho^*=k^*` and `Cl(Q_rho)[3]!=0`.

This is a rational degeneration of the elliptic THM-3913 construction, not
a Jacobian counterexample.  No plane atlas or Keller map is constructed.

## 1. The order and its one-place discriminant

Expansion of `(2)` gives

```text
a=A^3+AC,
b=3A^2C+C^2,
c=3AC^2+rho AC,
d=C^3+rho C^2+A^2.                                        (4)
```

These coefficients are primitive.  Indeed, a common prime would divide
`c=AC(3C+rho)`, but `A`, `C`, and `3C+rho` each fail to divide the complete
row `(4)`.  All four coefficients vanish at `(A,C)=(0,0)`.  Hence the
intrinsic binary index form cannot represent a unit, so `T_rho` is globally
nonmonogenic.  At `C=0`,

```text
Phi_rho=A^2(AU^3+V^3).                                    (5)
```

Choose hypothetical factors primitive over the `C`-adic DVR.  The `U^3`
coefficient in `(4)` is not divisible by `C`, so both primitive factors have
nonzero positive-degree specializations.  The `A`-valuation of `-1/A` is
`-1`, not a multiple of three; thus `AU^3+V^3` is irreducible over `k(A)`,
and the generic cubic is irreducible.

Its exact discriminant is

```text
Delta_rho=
 -27A^10-54A^8C-12rho A^6C^3-27A^6C^2
 +36A^4C^5-12rho A^4C^4
 -12rho A^2C^7+44A^2C^6-4C^9-4rho C^8.                   (6)
```

The degree-ten row is `-27A^10`.  Thus the projective closure has only
`[0:1:0]` at infinity, and the `-4C^9Z` term in the homogenization gives

```text
partial_Z Delta_h(0,1,0)=-4.                              (7)
```

The infinity point is smooth.  Every positive-degree projective component
would meet the infinity line, so a geometric factorization would force two
components through that sole point and make it singular.  Hence `(6)` is
absolutely irreducible.  Its height-one valuation is one.  The discriminant-
index identity makes the finite-free order maximal at that height-one prime;
it is etale at all other height-one primes.  Finite freeness gives `S2`, so
Serre's criterion makes `T_rho` normal.  Irreducibility and the nonsquare
discriminant give generic Galois group `S3`.

## 2. Nodal incidence and polynomial normalization

For arbitrary parameter `r`, eliminating the repeated root of `Phi_r` gives
the birational incidence model

```text
18y^2=(w+r)(w^2-r^2-6).                                   (8)
```

At `(1)`, this becomes the nodal cubic

```text
18y^2=w^2(w+rho),                                         (9)
```

whose normalization is `w=18t^2-rho`, `y=t w`.  Pulling this normalization
back through the elimination gives

```text
L(t)=162t^6-18rho t^4+6t^2-rho,
G(t)=(18t^2-rho)L(t),
A(t)=tG(t),
C(t)=-(3t^2-rho/6)G(t).                                  (10)
```

Direct substitution puts `(10)` on `(6)`.  The degrees are `9` and `10`,
and `A/C` is a uniformizer at the sole infinity point.  Moreover

```text
C(t)/A(t)=-3t+rho/(6t).                                   (11)
```

Two parameters with equal ratio are either equal or related by

```text
iota(t)=-rho/(18t).                                       (12)
```

The exact collision factor below is not identically zero, so `(10)` is
generically injective.  It is therefore the normalization, and the affine
normalization is `A1`.

## 3. Complete finite singular packet

The polynomial `G` has eight simple roots; exactly,

```text
Res(G,G')=-58535441882591744284714139375370240000.        (13)
```

They are the eight normalization addresses over the origin.  The only
nonimmersive addresses satisfy

```text
gcd(A',C')=t^2+rho/18.                                    (14)
```

At those two addresses the second/third derivative determinant reduces to
`-276480rho`, so they are ordinary `A2` cusps.  Their distinct targets are

```text
A=-10t,                         C=-10rho/3.                (15)
```

Among the eight origin branches, `(12)` pairs only the two roots with
`t^2=rho/18`.  The exact response factors are

```text
num(A(t)-A(iota(t)))
 =153055008(t^2-rho/18)^6(t^2+rho/18)^3,

num(C(t)-C(iota(t)))
 =-459165024(t^2-rho/18)^7(t^2+rho/18)^3.                 (16)
```

In the first blowup chart `q=C/A`, one has `dq/dt=-6` at either address, so
`q` is a local parameter on both strict transforms.  Their radial `A`
coordinates differ to order six by `(16)`.  Hence their strict transforms
intersect with multiplicity six and the original smooth branches have
intersection multiplicity `1+6=7`.  Every other origin pair has distinct
tangent and intersection multiplicity one.  Since all eight individual
branches are smooth,

```text
delta_origin=binom(8,2)+(7-1)=34.                         (17)
```

The two external cusps contribute one each.  A plane decic has arithmetic
genus `36`, while the normalization has genus zero, and

```text
34+1+1=36.                                                 (18)
```

The genus budget is exhausted, proving both the stated anatomy and the
absence of hidden finite singularities.

## 4. The three-class survives rationalization

The surface `Q_rho` is normal: it is a hypersurface and hence `S2`, while
the reduced plane branch `(6)` has only the finite singular packet in
Section 3, so the double plane is regular in codimension one.

The split infinity curves on a resolution of `(3)` have Gram matrix

```text
K_10=[-4  5],                    Smith(K_10)=(1,9).        (19)
     [ 5 -4]
```

Its nonzero determinant forces every unit on `Q_rho` to be constant by the
boundary-divisor argument.  Over the quadratic resolvent, the connected
`S3` splitting field becomes a connected cyclic cubic cover.  Generic
transposition inertia is killed by the quadratic base change, and purity
extends the cover over `Q_rho,reg`.  Since `k^*` is three-divisible, Kummer
theory sends this nontrivial torsor to

```text
0 != xi in Pic(Q_rho,reg)[3]=Cl(Q_rho)[3].                (20)
```

This proves existence, not purity: THM-3919 shows on the original elliptic
fibre why the visible boundary block alone does not determine the support.
The full removed lattice of the rational fibre remains a separate invoice.

## 5. The faithful discrete carrier is an anchored tie

The radial quadratic discriminant for the family is the cube of

```text
h_r(z)=(r^2+6)z^4+2rz^2+1.                               (21)
```

Put `q=z^2`, choose `eta^2=-6`, and homogenize:

```text
H_r(q,s)=(s+(r+eta)q)(s+(r-eta)q).                       (22)
```

The two colors never collide: their binary-`q` discriminant is the constant
`-24`.  Introduce the affine-infinity anchor `I=(0,1)` and color vectors

```text
v_+=(r+eta,1),              v_-=(r-eta,1),
omega(p,q)=det(p,q).                                        (23)
```

Then

```text
omega(I,v_+)=-(r+eta),
omega(I,v_-)=-(r-eta),
omega(v_+,v_-)=2eta,
omega(I,v_+)omega(I,v_-)=r^2+6.                           (24)
```

At `r^2=-6` exactly one anchored edge ties.  Its derivative in `r` is `-1`,
and the next `q` coefficient `2r` is nonzero.  The intrinsic natural depth
multiset is therefore `{1,0,0}`.  Encoding the event by the natural number
`1` retains its mathematical content; ranking the two colors alone loses it,
because their full unanchored gain `2eta` is constant.  At the critical
fibre the faithful object is a valued graph with a first tie response, not a
tournament.

Geometrically, one factor in `(22)` becomes `s`: a color crosses the affine
boundary instead of colliding with the other color.  This is exactly why the
elliptic quartic becomes a rational conic.

## 6. Boundary of the result

THM-3913 excluded a plane atlas because its affine normalization was an
elliptic curve minus one point.  Here the affine normalization is `A1`, so
that obstruction disappears.  Its disappearance is not an atlas theorem:
one must still construct or exclude a plane open in the cubic normalization,
compute the actual rational-fibre removed lattice and boundary tariff, and
meet the Keller condition.  Those tasks and `JC(2)` remain **OPEN**.

## 7. Exact replay

```powershell
python -B 04-computation/jc2_rational_moving_triple_root_decic_thm3918.py
python -B -O 04-computation/jc2_rational_moving_triple_root_decic_thm3918.py
```

Both streams byte-match the frozen LF transcript.  The companion verifies
the cubic expansion, exact family discriminant, unique smooth infinity,
radial cube, nodal incidence normalization, polynomial parametrization,
complete collision and genus ledger, anchored response, and boundary Smith
packet in 45 active gates.  **QED.**
