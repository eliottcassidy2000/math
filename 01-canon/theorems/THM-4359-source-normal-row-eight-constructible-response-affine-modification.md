---
id: THM-4359
title: "Source-normal row-eight constructible response affine modification"
status: >
  PROVED FINITE-ROW COROLLARY RELATIVE TO THM-4308 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The effective three-parameter row-eight response has
  a smooth hypersurface closure but a strictly smaller constructible image.
  It is an affine modification that is an isomorphism over zeta_3!=0 and has
  an affine-line fibre over each point of one exceptional line at zeta_3=0.
  Specializing the global elimination kernel misses one additional equation.
  The global and special-fibre algebraic matroids and the separate complete
  zero-wall nonface clutters are determined exactly. This is finite response
  geometry only; no all-row lift, seam entry, Keller pair, JC(2), or DC(2)
  conclusion is asserted.
source: root + row-eight algebraic-matroid scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
related:
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-527
primary_script: 04-computation/jc2_source_normal_row8_constructible_response_affine_modification_thm4359.py
primary_output: 05-knowledge/results/jc2_source_normal_row8_constructible_response_affine_modification_thm4359.out
primary_script_sha256: 3311b86ea74c83fcd9cb5b3d5ffc8bc6d79b487ebbabb7b9952c8343017284bb
primary_output_sha256: 11adca06a2e89c7218bd5db55ad6c59f216f3c2790e3d8c60a40b37cc2fa776d
referee_script: 04-computation/jc2_source_normal_row8_constructible_response_affine_modification_independent_referee_thm4359.py
referee_output: 05-knowledge/results/jc2_source_normal_row8_constructible_response_affine_modification_independent_referee_thm4359.out
referee_script_sha256: 39e64790b87568accb441aa979123456ad4b82c454d7c36f8941239d29978eba
referee_output_sha256: 5bf2daae017bddbb539c7e5caf9210c21690b557277115bbd185091e7095ff63
hash_basis: raw LF bytes
audit: >
  PASS. The 84-check primary and 93-check import-free clean-room referee independently
  derive both elimination kernels, inverse maps, constructible strata,
  affine-modification identity, fibre dimensions, matroid circuits, all
  zero-wall unit certificates, and an exhaustive witness bank. Normal and
  optimized runs byte-match both frozen outputs.
---

# THM-4359 -- Source-normal row-eight constructible response affine modification

**PROVED FINITE-ROW COROLLARY RELATIVE TO THM-4308 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS DISTINGUISHES AN ACTUAL CONSTRUCTIBLE RESPONSE
IMAGE FROM ITS ZARISKI CLOSURE. IT IS NOT AN ALL-ROW LIFT, SEAM-ENTRY,
KELLER-PAIR, `JC(2)`, OR `DC(2)` THEOREM.**

## 1. Effective response and silent fibres

Work over an algebraically closed field of characteristic zero. Put

```text
P=Phi,       e=eta,       x=xi_10,       z=zeta_3.       (1)
```

THM-4308's effective row-eight response is the polynomial map

```text
rho: A^3_(P,e,x) -> A^4_(z,U,W,Z),

z=-3P/2,
U=475515904/200475-6x/11,
W=-13P^2/12+424x/99-35956576256/1002375,
Z=3126529518592/27064125-130P^2/81
  -65Pe/36-22886x/2673.                                 (2)
```

The full finite THM-4308 solution family also has silent source parameters
`alpha_11,beta_11` and an affine seven-dimensional terminal tangent fibre.
Every isomorphism statement below concerns the effective map `(2)`, or the
full map after quotienting those nine silent directions.

The inheritance pass is:

- closest formula: THM-4308's explicit source response; closest geometric
  mechanism: THM-3406's `B[I/g]` affine modification;
- canonical hostile: a point satisfying the global closure equation but
  absent from the actual image;
- corrected near miss: elimination and then specialization need not equal
  specialization and then elimination;
- least-used sidecar: the constructible stratum and fibre address, retained
  separately from the generic algebraic matroid.

The live board was

```text
closure kernel | actual image | special fibre | affine modification
algebraic circuit | zero-wall nerve | silent fibre | next consumer.        (3)
```

## 2. Smooth closure and exact constructible image

Define

```text
F=1184625z^2+19318500U+2460375W+42434609152.             (4)
```

The full kernel of

```text
Q[z,U,W,Z] -> Q[P,e,x]
```

induced by `(2)` is the prime principal ideal `(F)`. Indeed, `F` vanishes
under `(2)`, is linear in `W`, and the minor

```text
det d(z,U,Z)/d(P,x,e)=-65P/44                           (5)
```

is nonzero. Thus the image has dimension three and `(F)` is the whole
height-one kernel. Since `dF/dW=2460375`, its closure `V(F)` is smooth.

Over `z!=0`, the effective response is an isomorphism onto `V(F) cap D(z)`.
An exact inverse is

```text
P=-2z/3,
x=237757952/54675-11U/6,
e=54/(65z)*(Z-5200877686784/66430125
             -125873U/8019+520z^2/729).                 (6)
```

At `z=0`, hence `P=0`, the coordinate `e` disappears. The actual image is
the line in `V(F)` cut out by

```text
G0=-1042743375U+66430125Z-5200877686784=0.              (7)
```

Conversely, `(4)` and `(7)` recover exactly the `W,Z` values parametrized by
the unique `x` determined by `U`. Therefore

```text
Image(rho)=V(F) cap (D(z) union V(G0)),                  (8)

closure(Image(rho)) minus Image(rho)=V(F,z) cap D(G0).  (9)
```

For example,

```text
(z,U,W,Z)=(0,0,-42434609152/2460375,0)                 (10)
```

lies in `V(F)` but has `G0=-5200877686784`, so it is not a response.

## 3. Affine-modification structure and fibres

Use `(4)` to identify the closure coordinate ring with

```text
A=Q[z,U,Z].                                               (11)
```

In the source ring, exact substitution gives

```text
2G0=z(159924375e-94770000z).                             (12)
```

Thus

```text
Q[P,e,x]=A[G0/z].                                        (13)
```

The response is the affine modification of smooth `A^3=Spec(A)` along the
center `(z,G0)` with divisor `z=0`. Its exceptional `A^2` maps to the center
line `V(z,G0)` with affine-one `e`-fibres. Hence the effective fibres are

```text
one point over D(z),       A^1 over the attained special line,
empty over the closure defect.                           (14)
```

Restoring `alpha_11,beta_11` and the seven terminal directions makes the
full finite-family fibres `A^9` generically and `A^10` on the exceptional
line. In particular, the full THM-4308 response map itself is not an
isomorphism.

This is exactly THM-3406's ring pattern `B[I/g]` at its first denominator
level, with `B=A`, `g=z`, and the second center generator `G0`. THM-3406's
principal-part tower suggests a higher transverse-jet question along `z=0`;
no such higher-level statement is asserted here.

Equation `(12)` also gives the exact specialization failure

```text
(z,F) strictly contained in (z,F,G0):                   (15)
```

the left side is the specialized global kernel, while the right side is the
kernel of the specialized map. The hostile `(10)` proves strictness. This is
the same algebraic warning used in THM-4255—specialization can create a new
kernel direction—but it implies nothing about the external density claim
audited there.

## 4. Retaining `xi_10`

If `x=xi_10` is kept as an observed target coordinate, put

```text
R1=109350x+200475U-475515904,
R2=482625z^2-4293000x+1002375W+35956576256,
R7=231720750x+27064125Z-3126529518592.                  (16)
```

The closure kernel in `Q[z,x,U,W,Z]` is the prime smooth codimension-two
complete intersection `(R1,R2)`. It equals `(R1,F)` because

```text
11F=27R2+1060R1.                                        (17)
```

Its exact image is

```text
V(R1,R2) cap (D(z) union V(R7)).                        (18)
```

The inverse over `z!=0` keeps the observed `x`, takes `P=-2z/3`, and uses

```text
e=54/(65z)*(Z-3126529518592/27064125
             +520z^2/729+22886x/2673).                 (19)
```

At `z=0`, equation `R7=0` is exactly the additional special-line condition.

## 5. Algebraic matroids: generic and special

For the closure function field on `{z,U,W,Z}`, the algebraic-matroid rank is
three and its complete circuit list is

```text
{z,U,W},       circuit polynomial F;                    (20)
```

`Z` is a coloop. With `x` retained, the rank remains three and the complete
circuits are

```text
{x,U}: R1,       {z,x,W}: R2,       {z,U,W}: F;         (21)
```

again `Z` is a coloop. Passing to the dense open `z!=0` does not change the
function field, so it does not change either matroid.

On the actual special image line, `{z}` is a loop. After deleting it, the
rank is one. All three pairs in `{U,W,Z}` are circuits, with primitive
equations

```text
19318500U+2460375W+42434609152=0,
-1042743375U+66430125Z-5200877686784=0,
115860375W+57955500Z-2539122786304=0.                  (22)
```

With `x` retained, all six pairs in `{x,U,W,Z}` are circuits; the three new
equations are

```text
109350x+200475U-475515904=0,
-4293000x+1002375W+35956576256=0,
231720750x+27064125Z-3126529518592=0.                  (23)
```

These are algebraic matroids of the closure or the actual special image.
They are not the zero-wall feasibility system considered next.

## 6. Complete zero-wall nonface clutters

For a ground set of observed affine coordinates, call a subset a zero-wall
nonface when their simultaneous vanishing has empty inverse image under
`rho`, and minimal when every proper subset occurs. On `{z,U,W,Z}`, the
complete list is

```text
{z,U,W},       {z,U,Z},       {z,W,Z}.                  (24)
```

On `{z,x,U,W,Z}`, the complete list is

```text
{x,U},
{z,x,W}, {z,x,Z}, {z,U,W}, {z,U,Z}, {z,W,Z}.           (25)
```

These lists form a pointed affine feasibility clutter, equivalently the
minimal nonfaces of the zero-wall nerve. They are not algebraic-matroid
circuits. For instance, `{z,x,Z}` is globally algebraically independent even
though its three zero walls have no common response.

Here are exact source-ring unit certificates, in the order displayed in
`(25)`:

```text
1=(109350x+200475U)/475515904;

1=-(482625z^2-4293000x+1002375W)/35956576256;

1=((38610000z-65154375e)z+463441500x+54128250Z)
  /6253059037184;

1=-(1184625z^2+19318500U+2460375W)/42434609152;

1=(-2085486750U+132860250Z+(94770000z-159924375e)z)
  /10401755373568;

1=(115860375W+57955500Z+(97124625z-69761250e)z)
  /2539122786304.                                       (26)
```

In each identity, the displayed response coordinates mean their source
polynomials from `(2)`. Setting the named subset to zero gives `1=0`, so the
subset is infeasible.

## 7. Minimality and exhaustiveness

On `P=z=0`, the four zero addresses in the source coordinate `x` are

```text
x: 0,
U: 237757952/54675,
W: 4494572032/536625,
Z: 1563264759296/115860375.                             (27)
```

They are pairwise distinct. Hence every pair `{z,q}` occurs for
`q in {x,U,W,Z}`, while `z` and two different such coordinates do not.
The singletons `{x}` and `{U}` also occur.

The only triple controls still needed away from `z=0` are

```text
x=W=Z=0:
  P^2=-143826305024/4343625,
  Pe=6086390091776/65154375;

U=W=Z=0:
  x=237757952/54675,
  P^2=-13056802816/820125,
  Pe=707514056704/12301875.                              (28)
```

Both have `P!=0` and satisfy `(2)` exactly over the algebraic closure. Any
zero set avoiding `z` and `{x,U}` is a subset of one of these two triples.
Any zero set containing `z` and avoiding `(25)` has at most one of
`x,U,W,Z`, hence occurs at `(27)`. This proves minimality and exhaustiveness
without relying on enumeration alone. Omitting `x` gives `(24)`.

The first, second, and fourth nonfaces in `(25)` already follow from global
closure circuits. The other three involve `Z` and appear only after
specializing to `z=0`. Thus the global algebraic matroid forgets three real
owner-wall incompatibilities.

## 8. Audit and scope

The 84-check primary and 93-check independent referee each derive the formulas
rather than sampling points. They separately check lexicographic elimination, primeness,
smoothness, both inverse maps, `(8)--(19)`, generic and specialized Jacobian-
rank circuit enumeration, every unit identity `(26)`, every address and
hostile witness `(27)--(28)`, and exhaustive zero-wall enumeration. The
referee imports no primary code. Normal and optimized runs byte-match both
frozen outputs in the frontmatter.

This theorem classifies one finite row-eight response map. It does not show
that any response extends to row nine or all rows, that a hypothetical
counterexample enters the source-normal or exact seam, or that a finite jet
terminates polynomially. It produces no Keller pair and does not prove
`JC(2)` or `DC(2)`, which remain open.
