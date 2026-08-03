---
id: THM-3280
title: "Coprime selector-geodesic potentials and integral augmentation completion"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT.  THM-3278's selector cut potential and
  THM-3277's even-geodesic-length potential have augmentation indices 35 and
  81 in their canonically root-normalized abstract C12 coordinates.  Their
  orbit lattices therefore generate the full integral augmentation lattice;
  exact Smith computation gives the stronger statement directly, including
  explicit sparse unimodular bases.  Promotion awaits independent hostile
  audit of the saturated-basis calculation and the cross-coordinate scope.
source: root/2026-08-03
audit: >
  The assertion-independent exact companion pins all four dependency
  surfaces, reconstructs both words from their promoted transcripts, and
  checks the uncentered and augmentation Smith forms in the saturated basis
  x^j-1.  It verifies the coprime completion, every one-vector extension,
  all 792 sparse 10+1 trials and all 220 reverse 2+9 trials.  Normal,
  optimized and stored output agree.  Independent audit is pending.
depends_on:
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
  - THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization
  - THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas
  - THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary
related:
  - THM-3274-norm-fibre-constrained-phase-transfer-and-refinement-invoice
script: 04-computation/gmc_coprime_selector_geodesic_potentials_thm3280.py
output: 05-knowledge/results/gmc_coprime_selector_geodesic_potentials_thm3280.out
script_sha256: 85d86e7707a852425a33108182559437a73f92979d9491f1221b19296336c017
output_sha256: 97065ac5eeaaaa12136cf1f4c9d2e7aed665e6d10f9f37aaed27f7a41fc6a10c
hash_basis: LF-normalized bytes
---

# THM-3280 -- coprime selector-geodesic potentials and integral augmentation completion

**PROVED CANDIDATE + VERIFIED-EXACT.  INDEPENDENT AUDIT PENDING.**

Two canonical binary potentials have appeared on the order-twelve response
coordinate, but each one separately misses part of the integral lattice.
Their missing primes are complementary.  This theorem turns that numerical
signal into an exact integral completion and transports it through the
THM-3260 holotopy atlas.

## 1. Coordinate contract and the two words

Put

```text
R=Z[C_12]=Z[x]/(x^12-1),
N=1+x+...+x^11,
Aug(R)=ker(epsilon:R-->Z).                              (1)
```

THM-3278 supplies the nonlinear rank-label coordinate on the twelve response
rows.  In its root gauge

```text
ell(17)=0,       ell(16)=1,                            (2)
```

the selector small-class word is

```text
f_A=x+x^2+x^4+x^5+x^10.                               (3)
```

Independently, THM-3277 uses the genuine critical quotient gauge

```text
j_12(17)=0,     j_12(16)=1.                            (4)
```

Its twelve minimum path lengths, indexed by target residue, are

```text
(2,2,1,1,2,1,1,1,2,1,1,2).                           (5)
```

The even-length marker is therefore

```text
f_B=1+x+x^4+x^8+x^11.                                 (6)
```

Equations `(2)` and `(4)` fix a unique identity between the two abstract
rooted copies of `C_12`: residue `r` is sent to residue `r`.  This is the
coordinate identification used below.  It is canonical relative to the
promoted gauges, but it is not a vertex-equivariant identification.  In
particular, `ell` is a bijection on vertices whereas `j_12` has nontrivial
vertex fibres, and `(6)` originally marks target-path lengths rather than
vertices.

## 2. The uncentered obstruction

Both words have mass five.  Exact circulant Smith calculations give

```text
det Circ(f_A)=175,       SNF Circ(f_A)=(1^10,5,35),
det Circ(f_B)=405,       SNF Circ(f_B)=(1^10,9,45).     (7)
```

The combined uncentered lattice has

```text
SNF [Circ(f_A) Circ(f_B)]=(1^11,5).                    (8)
```

Indeed every element of `f_A R+f_B R` has augmentation divisible by five,
so it lies in `epsilon^{-1}(5Z)`, a sublattice of index five.  Equation `(8)`
shows equality:

```text
f_A R+f_B R=epsilon^{-1}(5Z).                          (9)
```

Thus the two words do **not** generate `R`.  Centering is load-bearing; the
common mass-five obstruction is the sharp failure boundary.

## 3. Coprime completion of the saturated augmentation lattice

Use the integral basis

```text
b_j=x^j-1,       1<=j<=11,                             (10)
```

of `Aug(R)`.  This basis is saturated.  Deleting one coordinate from a
twelve-vector orbit instead would introduce the familiar spurious factor
twelve.

Define

```text
L_A=f_A Aug(R)=(x-1)f_A R,
L_B=f_B Aug(R)=(x-1)f_B R.                             (11)
```

In basis `(10)`, the multiplication matrices have Smith forms

```text
SNF(L_A -> Aug(R))=(1^10,35),
SNF(L_B -> Aug(R))=(1^9,9,9).                          (12)
```

Consequently

```text
[Aug(R):L_A]=35,       [Aug(R):L_B]=81.                (13)
```

Since `gcd(35,81)=1`, the finite quotient
`Aug(R)/(L_A+L_B)` has order dividing both indices and is trivial.  Hence

```text
L_A+L_B=Aug(R).                                        (14)
```

The direct stacked Smith calculation is

```text
SNF [L_A L_B]=(1^11),                                  (15)
```

so `(14)` is an exact integral equality, not merely rational spanning.

## 4. Sparse unimodular certificates and asymmetry

Let

```text
g_A=(x-1)f_A=-x+x^3-x^4+x^6-x^10+x^11,
g_B=(x-1)f_B=x^2-x^4+x^5-x^8+x^9-x^11.                (16)
```

The full orbit of `g_A` has index 35 and the full orbit of `g_B` has index
81.  Nevertheless, fix the unshifted `g_B`.  Each of the following eleven
element sets is an integral basis of `Aug(R)`:

```text
{x^i g_A : i in {0,1,2,3,5,6,7,8,10,11}} union {g_B}, det=+1,

{x^i g_A : i in {0,2,4,5,6,7,8,9,10,11}} union {g_B}, det=-1.             (17)
```

For every fixed shift `x^t g_B`, exactly two of the `binom(12,10)=66`
ten-shift subsets of the `g_A` orbit are unimodular.  Thus `(17)` belongs to
an exact bank of 24 sparse bases.

There is an even cheaper extension if the whole selector lattice is already
available.  For

```text
h_{B,j}=f_B(x^j-1),       1<=j<=11,                    (18)
```

the indices of `L_A+Z h_{B,j}` are, in order,

```text
(1,1,7,5,1,7,1,5,7,1,1).                              (19)
```

Thus one geodesic relative potential kills the entire selector defect for
`j in {1,2,5,7,10,11}`.

The reverse direction is genuinely less economical.  Because

```text
Aug(R)/L_B = Z/9Z direct-sum Z/9Z,                     (20)
```

one added vector cannot kill its two-generator 3-primary quotient.  Exact
calculation gives the one-selector-vector residual indices

```text
(9,9,9,27,9,9,9,27,9,9,9).                            (21)
```

Two selector relative vectors are sufficient and minimal.  For example,

```text
f_A(x-1), f_A(x^2-1),
{x^i g_B : i in {0,1,2,4,5,7,8,10,11}}                (22)
```

form an integral basis of determinant `-1`.  Exactly seven of the 220
nine-shift choices in `(22)` are unimodular.

## 5. Holotopy completion of the response cycle lattice

THM-3260 constructs, from a cyclic vertex label and an ordered bispanning
chart, an integral isomorphism

```text
Psi:Aug(Z[C_12]) --> H_1(G_0;Z).                       (23)
```

The two formerly external choices are now available: THM-3278 supplies the
canonical rank-label bijection `(2)`, while THM-3269 supplies the uniquely
weighted ordered chart.  Therefore `(23)` is fixed relative to the promoted
labeled response data.  Applying it to `(14)` gives

```text
Psi(L_A)+Psi(L_B)=H_1(G_0;Z),                          (24)
```

and either basis in `(17)` gives an explicit eleven-cycle integral basis.

Every symmetric exchange in THM-3260 changes the cycle frame by an element
of `GL_11(Z)`.  Hence `(24)` and unimodularity of `(17)` transport through
the connected chart atlas.  Around a closed exchange loop the transition
matrices still telescope: the two-potential completion creates no new
holonomy.  Its content is instead a **holotopy completion**--two locally
nonsaturated canonical marker systems have coprime torsion defects and
together furnish a globally transportable integral frame.

## 6. Map and loss contract

The connection proved here is

```text
source A: selector face cut in the nonlinear rank C12 coordinate;
source B: parity of canonical minimum-path lengths in genuine target J12;
map:      root- and primitive-normalized identity of the abstract C12 copies,
          then cyclic difference and the selected THM-3260 chart;
target:   the integral response cycle lattice H_1(G_0;Z);
preserved: cyclic translation, augmentation, lattice index and unimodularity;
destroyed: vertex identity, path endpoints, path realization, face legality,
           response signs, owner labels and every positive cone.             (25)
```

The coordinate identity in `(25)` is not induced by a graph automorphism or
a physical response transition.  It does not make the selector cut a
critical-group character, identify a geodesic target with a response row, or
show that the two marker events can occur simultaneously.  The sparse bases
are coefficient/cycle bases, not physical walks.

Accordingly `(24)` proves no row exclusion, no positive moment functional,
no `FC(3)` or `SFC(3)`, no `LRC(14)` decrement, and no new Gaussian Moment
Conjecture case.

## 7. Exact verification

Run

```text
python 04-computation/gmc_coprime_selector_geodesic_potentials_thm3280.py
python -O 04-computation/gmc_coprime_selector_geodesic_potentials_thm3280.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_coprime_selector_geodesic_potentials_thm3280.out.
```

The companion uses exact integer, Smith-normal-form and finite enumeration
only.  It has no assertion node, floating literal, randomness or fitted
recurrence.

QED.
