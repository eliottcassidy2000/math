# Planar Jacobian: the Delta-zero endpoint tower and the weight-ten hidden tail

## Status firewall

- **PROVED RELATIVE + VERIFIED-EXACT + independently source-audited:**
  [THM-4217](../01-canon/theorems/THM-4217-complete-mixed-off-antidiagonal-delta-zero-planar-jacobian-exclusion.md)
  closes the mixed off-antidiagonal exact-`M=9` face `Delta=0` inside the
  inherited `b=d=0` reduced `(2,3)` seam.
- **PROVED RELATIVE + VERIFIED-EXACT + independently audited:**
  [THM-4218](../01-canon/theorems/THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion.md)
  closes the dense-top exact-`M=10` chamber for arbitrary lower coefficients.
- **PROVED RELATIVE + VERIFIED-EXACT + independently audited:**
  [THM-4220](../01-canon/theorems/THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion.md)
  closes the entire `zeta=0` exact-`M=10` top-unit wall, including its
  split-conic strict transform.
- **PROVED RELATIVE + VERIFIED-EXACT + independently audited:**
  [THM-4222](../01-canon/theorems/THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion.md)
  closes the dense-side exact-`M=11` chamber by primitive CM and degree zero.
- **FINITE-EXACT / TRACKED SCRATCH:** one independent critical-open control
  in [the scratch report](../.scratch/jc_m10_dense_20260826/REPORT.md).
- **OPEN:** the three remaining exact-`M=10` top walls, five dense-`M=11`
  walls, nonzero-`Delta`/nonzero-`K` exact-`M=9` critical walls, seam entry,
  other cells, exact `M>=12`, `JC(2)`, and `DC(2)`.

This reflection preserves the session's operation choices and next questions.
The theorem files, not this reflection, are the truth sources for
THM-4217/4218/4220/4222.

## Portfolio and inheritance

| lane | closest proved mechanism | hostile or corrected near miss | least-used sidecar |
|---|---|---|---|
| Anchor: mixed `Delta=0` | THM-4147 packet/carrier response | a vanished `Delta` face cannot keep its generic label | replacement `W=sp` face plus source-first endpoint tower |
| Niche: exact `M=10` | THM-4147 critical length and prime carrier | THM-4045's hidden positive-genus tail defeats highest-face-only reasoning | component Jacobian plus attachment difference |
| Wildcard: scalable wall compiler | THM-4159/4180 strict transforms | a projected resultant discriminant is not critical-scheme loss | low-`p` Hasse rows before the full resultant |

The repeated-top route was deliberately retired at inheritance: THM-4176
already closes all of it, and THM-4205 closes the final exact-`M=9`, `K=0`
row. MISTAKE-519 records why theorem-title archaeology and ambient quantifiers
must precede a new derivation.

## 1. Why the `Delta=0` face closes completely

The generic mixed source resultant has two endpoint factors with different
roles:

```text
bottom:  D=4K0^2 Theta-27zeta^2,
top:     eta^3 zeta^2(eta+zeta)^4.                    (1)
```

The top factor is exactly the structural scope: P-only, Y-only, and
repeated-top are separate faces. The bottom factor is an affine critical
wall, not a carrier degeneration. Once `D=0` is parameterized without a sign
quotient, the successive bottom rows are

```text
D -> J -> S -> T0 -> unit.                            (2)
```

The mixed coefficient `eta` first appears in `S`, one row later than the
Y-only analogue. This delayed visibility is useful structure: low-order
critical jets see `Theta,zeta,Phi` before they see the second mixed top owner.
The five residual degrees `21,20,19,18,17`, after four universal points are
restored, give exact lengths

```text
25,24,23,22,21.                                       (3)
```

At the boundary the old `Delta` face disappears, but the diagonal face

```text
Qp^4(-1376/135+eta sp)                                (4)
```

replaces it with one simple rational index-four place. The length-two top
face remains split because `eta*zeta*(eta+zeta)!=0`; its two rational places
both have index eight. Together with the prime cubic carrier, the packet is

```text
(8,8,4,2,2,2,1),       genus=11,       defect=20.     (5)
```

The response inequalities have large slack. Full-response commutator caps
are `4,6,8,10,12<20`. In the finite response, deleting the cubic carrier
leaves origin index `17`; the carrier-orbit caps are `3,5,7,9<17`, while the
length-25 row fails even the required handle-support union. Thus every
critical strict transform remains on the safe side of the same boundary
obstruction.

## 2. The weight-ten face that almost gives a different proof

Let the two new top coefficients be typed as

```text
upsilon=[p^5]H,       xi=[p^2y^2]H,
upsilon*xi*zeta*(upsilon+xi)!=0.                      (6)
```

Exact enumeration over all `256` lower-support subsets retains only two lower
planes and gives

```text
polygon=(0,1),(2,0),(5,3),(4,4),(0,6),
(2Area,boundary,g)=(36,12,13),
packet=(9,9,6,2,2,2,1).                               (7)
```

The hull itself needs only `upsilon*xi*zeta!=0`: its three outer vertices
are uniquely owned by those coefficients. The extra `upsilon+xi!=0` is
exactly the gate separating the two rational index-nine top places, not a
hidden hull hypothesis.

The main nonrational face normalizes to

```text
xi Y^2=1-upsilon P^5,                                 (8)
```

a genus-two curve with an order-five automorphism. Its Jacobian has a
`Q(zeta_5)` action and is simple, so it has no elliptic quotient. This first
suggests a face-Hom obstruction against the elliptic target.

The complete lower model supplies the hostile. A side face normalizes to

```text
xi W^2=1-zeta T^3,                                    (9)
```

a `j=0` elliptic curve. Its attachments `(0,1)` and `(0,-1)` differ by exact
three-torsion. Hence the highest face forgets the one component capable of
mapping to the target, and the naive attachment repair is itself killed by
the Eisenstein endomorphism `1-zeta_3`. What survives is a divisibility
sidecar: any compatible specialized map must have degree divisible by three.
That is information, not a contradiction.

## 3. A second route independently finds a weight-ten obstruction

At the exact rational control

```text
Delta=1, K=5591/90, Phi=2, Theta=5, eta=7, zeta=11,
upsilon=13, xi=17,                                     (10)
```

two source projections give the same squarefree degree-25 residual, and a
disjoint normalized projection gives another squarefree degree-25 residual.
All endpoint and chart-loss checks are nonzero. Restoring the four universal
Morse points gives `L=29`. With packet `(7)`, the inherited response arithmetic
is

```text
full:    2(31-29)=4<24,
finite:  2*25-29-1+3=23<24.                           (11)
```

Thus the exact critical-open part of this chamber is nonempty before the
Keller obstruction is imposed. This calculation is an independent positive
control, not the reason the full coefficient chamber closes. The exact scout
and frozen output remain routed through the scratch report.

## 4. Degree-three closure from the hidden tail

The hostile tail may carry a stronger obstruction than the critical length.
On a compatible regular specialization, rational components are constant by
Riemann--Hurwitz and the genus-two component is constant by simplicity. The
`j=0` tail is then the only possible nonconstant component. Because its two
attachments differ by nonzero three-torsion, every map identifying them has
degree divisible by three. Degree conservation would force

```text
3 | degree(generic fibre).                             (12)
```

But the two carrier responses have degrees `31` and `25`, both congruent to
one modulo three. THM-4218 proves this closes the whole packet chamber
independently of the critical length. Its regular-total-space certificate
after `Q=sigma^30` checks primitive face multiplicity one, every compactified
face/edge root, rationality of the resolution chains, and specialization
degree conservation with the labelled attachments. An independent
reconstruction recovered the same model and degree obstruction.

## 5. Rational side faces and primitive CM

On `zeta=0`, the old tail disappears but the `y^2` row creates a rational
side face. With `D_V=Theta^2-4Kxi`, the three polygons have genera `12,12,11`
according as `K!=0`, `K=0,Theta!=0`, or `K=Theta=0`. Off `D_V=0` the side
normalizes to the conic

```text
Y^2=K+Theta P+xi P^2.                                  (13)
```

The genus-two main component is then the sole positive-genus component, so
degree conservation is zero. On `D_V=0,K!=0`, the first strict-transform row
has order `sigma^15`. If `A_0=Phi+eta alpha` is nonzero it gives an `A_29`
rational path; if `A_0=0` normalization gives two rational sheets. Thus no
hidden tail survives, and THM-4220 closes the entire top-unit wall.

At exact weight eleven, the new coefficients

```text
A=[p^4y]H,       B=[py^3]H
```

together with `U=[p^5]H,Z=[y^3]H` force three lower faces. Under
`A*B*U*Z*(A+B)!=0`, the global polygon has Pick genus `15` and packet
`(10,10,5,4,2,2,2,1)`. Both side faces are rational. The main component is

```text
P^11=(B/A^3)x^3/(1-x),                                (14)
```

a genus-five cyclic cover with CM type `{4,5,8,9,10}` and trivial
stabilizer. Milne's primitive-CM classification makes its Jacobian simple,
so it has no elliptic quotient. The `Q=sigma^330` model has only rational
side, toric, and `A_329` components; degree conservation again gives zero.

This also records a genuine method boundary: the full and finite carrier
degrees are `36,30`, both divisible by three, so THM-4218's congruence does
not extend. The order-eleven attachment difference is exact torsion but is
not a genus-one degree invoice; primitive CM is the stronger coordinate.

## 6. New operation suggested by the endpoint tower

Full symbolic resultants become expensive at `M>=10`, but the obstruction
uses only four pieces of their information:

```text
p-valuation, residual degree, first nonzero p-row, last nonzero p-row. (15)
```

This suggests a scalable **endpoint Hasse compiler**:

1. derive the source pair over the complete coefficient ring;
2. compute successive subresultant or Sylvester Hasse rows only at `p=0`;
3. stratify by the first nonzero row, retaining an explicit parameter chart;
4. compute the reciprocal last row independently at `p=infinity`;
5. call the full resultant only on the final small quotient algebras.

The source is the unsaturated critical ideal; the target is exact affine
length. The map forgets individual projected roots but preserves intersection
length once the Hessian bridge and endpoint units are retained. Required
sidecars are the second projection, coordinate fibres, source infinity, and
terminal quotient gcds. THM-4217 is a positive control; projected-root
collisions are the hostile showing why squarefreeness must not enter the
compiler.

## 7. Ranked continuation

1. Resolve the exact-`M=10` contractions `upsilon=0`, `xi=0`, and
   `upsilon+xi=0`; each changes a vertex, component, or contact divisor.
2. Build the finite exact-`M=11` wall atlas `A=0,B=0,U=0,Z=0,A+B=0`; the
   replacement side faces are rational but the main CM type may change.
3. Enumerate exact `M=12`, where new face moduli may replace cyclotomic
   rigidity; retain all positive-genus components before testing `Hom`.
4. Build the endpoint Hasse compiler and reproduce the exact `M=9`
   `D,J,S,T0` tower before trusting it on higher critical walls.
5. On remaining `M=9` critical walls, use carrier orbits before elementary
   support caps; projected-root squarefreeness is not a source-length gate.
6. Keep other cells and seam entry separate. No seam theorem is evidence for
   global entry without an explicit target-preserving reduction.

This is the session's main conceptual shift: coefficient walls are best
organized by what the source endpoint tower can still see, while boundary
components are organized by the maps and attachment classes they can carry.
Neither organization alone is complete; their response intersection is the
actual obstruction.
