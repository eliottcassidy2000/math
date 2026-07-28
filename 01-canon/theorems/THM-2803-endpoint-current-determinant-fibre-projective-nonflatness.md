---
id: THM-2803
title: "Endpoint-current determinant-fibre projective arcs and marked-power boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonzero central direction of the canonical
  THM-2625 endpoint current, the thirteen determinant-fibre profiles are
  pairwise maximally separated after every relative cyclic origin, scalar,
  and orientation choice: all 78 projective minors survive in both certified
  fields.  Every pair is therefore a 13-point projective arc and an
  [13,2,12] MDS coefficient code.  A separate exact pointwise-function lemma
  identifies the first-forgotten digit as a persistent marked-power
  separator.  This is coefficient-side anti-quotient evidence, not physical
  endpoint allocation, a Special Image theorem, a runner-row exclusion, or
  LRC(14).
source: current-frontier-synthesis-2026-07-28
depends_on:
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
related:
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
  - THM-2801-sharp-special-image-boundary-and-beta-shift-witness
  - THM-2802-normalized-unit-bispectrum-and-projective-cyclic-orbit-reconstruction
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
script: 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.py
output: 05-knowledge/results/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.out
independent_script: 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.py
independent_output: 05-knowledge/results/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.out
script_sha256: d3c011f3583a6f58e07be3a4f170f10f641dc35ba986e9d502de7f6e11e6958d
output_sha256: 9e216edd09aa460e5a683556e1a5a47bbcd502519813d44f403648a424b5a8ed
independent_script_sha256: d5d0ddeb6e1367f120562364fc8265b728345427258124e3f40619bd90f554c5
independent_output_sha256: 5d5a59dededc066b0bda9ee7ce5026bf8effc8083771a212e753321c5f477b09
hash_basis: LF-normalized bytes
---

# THM-2803 -- endpoint-current determinant-fibre projective arcs and marked-power boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2790 proves that every spatial coefficient and every Fourier coordinate
of every canonical endpoint cycle is nonzero.  THM-2792 consequently gives
an abstract full-unit circulant intertwiner from the THM-2791 semantic chain
to each endpoint cycle.  The unresolved question is whether those separate
intertwiners hide a much smaller transport law between determinant fibres.

They do not.  Every two-fibre coefficient matrix is as far from projectively
flat as a `2 x 13` matrix can be: every maximal minor survives, after every
origin choice and even after orientation reversal.  Thus the endpoint bank
has a flat transport groupoid only in the **full convolution-unit group**;
its structure group cannot be reduced to scalar monomials.

The same coefficient vector also carries a second multiplication,
coefficientwise rather than convolutional.  Keeping those two algebra
structures typed reveals the exact role of the first forgotten address
digit: a fibre-difference functional kills every pointwise power of a coarse
pullback, while multiplication by that digit makes every marked power
survive.  This is an elementary finite-function Mathieu-subspace analogue,
not an instance or consequence of the Special Image Conjecture.

## 1. Two algebra structures on one endpoint cycle

Let

```text
K=Q(zeta_N),                 N=50,334,435,734,703,120,
V=F_13^2.
```

For every `s in V\{0}`, THM-2625 gives cyclotomic-integral endpoint factors
`P_R,Q_R`.  Put

```text
J_s(R)=P_(R+s) Q_R.                                           (1)
```

Choose `t in V` with `det(s,t)=1`.  The determinant-`delta` cycle profile is

```text
j_(s,delta)(w)=J_s(w s+delta t),       w,delta in F_13.        (2)
```

There are

```text
168*13=2,184                                                (3)
```

such profiles.  Each lives simultaneously in two commutative algebras:

| algebra | multiplication | unit test |
|---|---|---|
| `A_pt=K^(C_13)` | pointwise | all thirteen coefficients are nonzero |
| `A_conv=K[C_13]` | cyclic convolution | all thirteen Fourier coordinates are nonzero |

THM-2790 gives the twelve nontrivial Fourier coordinates and all spatial
coefficients.  THM-2625 gives the nonzero zero-mode.  Hence every endpoint
profile is a unit in **both** algebras.  We call this only the local
`biunit` locus; it does not identify the two multiplications.

THM-2791's normalized semantic chain `z^6(1+z)` is a convolution unit but
has eleven zero coefficients, so it is not a pointwise unit.  This already
explains why the THM-2792 convolution intertwiner need not be a pointwise
physical identification.

## 2. Coefficient logarithmic derivative

The pointwise-unit torus has a simple complete projective invariant which is
dual to THM-2802's proved Fourier bispectrum.  For any field `L` and
`f in (L^x)^(C_n)`, define

```text
(partial_x f)(w)=f(w+1)/f(w).                               (4)
```

Then

```text
product_w (partial_x f)(w)=1.                               (5)
```

Conversely, every nonzero word `rho` with product one is `(4)` for a word
`f`, unique after choosing the single scalar `f(0)`.  Thus `(4)` gives the
exact sequence

```text
1 -> L^x -> (L^x)^(C_n)
  -> {rho in (L^x)^(C_n): product rho=1} -> 1,              (6)
```

and is equivariant for cyclic translation.  Consequently the cyclic orbit
of `partial_x f` is a complete invariant of `f` modulo nonzero scalar and
cyclic origin.

For reversal `f^vee(w)=f(-w)`,

```text
partial_x(f^vee)(w)=1/partial_x f(-w-1).                    (7)
```

Thus the dihedral quotient additionally identifies a ratio word with its
reversed inverse.  Equations `(4)--(7)` require coefficientwise nonvanishing,
not convolution invertibility.  THM-2802's bispectrum requires Fourier
nonvanishing instead.  The endpoint cycles admit both descriptions because
they are biunits.

## 3. Maximal projective separation

Fix `s`, two distinct determinant fibres `delta!=epsilon`, a relative
origin shift `c`, and an orientation `eta in {+1,-1}`.  Write

```text
f(w)=j_(s,delta)(w),
g(w)=j_(s,epsilon)(eta*w+c).                                (8)
```

For every pair `0<=u<v<=12`, form the projective minor

```text
D_(u,v)=f(u)g(v)-f(v)g(u).                                  (9)
```

The exact theorem is

```text
D_(u,v) != 0                                                (10)
```

for **every** choice in `(8)--(9)`.

The primary computation checks `(10)` in both certified THM-2625 fields.
Per field its exact universe is

```text
168 directions;
78 unordered determinant-fibre pairs per direction;
13 relative cyclic origins;
2 orientations;
78 minors per comparison.                                  (11)
```

Thus each field has

```text
170,352 oriented comparisons,
170,352 reversed comparisons,
13,287,456 nonzero oriented minors,
13,287,456 nonzero reversed minors.                         (12)
```

There are zero equality, scalar-projective, cyclic-projective, or
dihedral-projective matches.  More sharply, every one of the

```text
26,574,912
```

minors in `(12)` is nonzero in each field.

Each finite-field map is a certified specialization of the cyclotomic
integer construction.  If a characteristic-zero minor in `(9)` vanished,
its image would vanish in every good specialization.  Either certified
nonzero image therefore proves `(10)` over `K`; the second field is an
independent arithmetic control.

## 4. Projective-arc, MDS, and decoder consequences

Because every coefficient of `f` is nonzero, `(10)` is equivalent to the
thirteen pointwise ratios

```text
g(w)/f(w),                  w in C_13,                       (13)
```

being pairwise distinct.  Equivalently, the thirteen columns

```text
[f(w):g(w)] in P^1(K)                                       (14)
```

form a projective `13`-arc.

The `2 x 13` matrix with rows `f,g` therefore generates an exact

```text
[13,2,12] MDS code.                                         (15)
```

Indeed, a nonzero combination `alpha f+beta g` cannot vanish at two
coordinates, since those two zeros would force the corresponding minor to
vanish.  Conversely one can cancel any chosen coordinate, so the minimum
support is exactly twelve.  Thus every two-fibre linear blend retains at
least `12/13` coefficient support.  This is an algebraic support floor, not
a positive real magnitude floor.

Equation `(13)` also recovers the central cycle coordinate from any chosen
pair of determinant profiles and relative gauge.  It is a coefficient
address decoder, but it is not yet a physical endpoint allocation: the
input already assumes two coefficient profiles on the THM-2625 endpoint
plane.

## 5. Frame and origin audit

Every other normalized transverse vector is

```text
t'=t+a s.
```

Substitution in `(2)` gives

```text
j'_(s,delta)(w)=j_(s,delta)(w+a*delta).                    (16)
```

Thus a frame shear only changes each fibre origin.  In a two-fibre
comparison it changes the relative shift by `a(epsilon-delta)`, already
quantified over in `(8)`.  The primary companion explicitly checks all

```text
168*13*13*13=369,096                                      (17)
```

profile entries under all frame shears in each field.

An independent hostile audit avoids a transverse vector altogether.  It
constructs the thirteen affine lines

```text
{R in V:det(s,R)=delta}
```

directly, chooses their least points, and tries every equivariant or
antiequivariant image of that origin.  Every distinct-fibre comparison has

```text
12/12 nonzero anchored cross-products,
78/78 nonzero pairwise cross-products.                      (18)
```

This independently proves that `(10)` is not an origin convention, frame
formula, or orientation artifact.

No covariance under replacing the increment `s` by a different nonzero
direction is asserted.  The current `P_(R+s)Q_R` itself changes with `s`;
a linear reindexing of the endpoint plane is not automatically a physical
current symmetry.

## 6. What “nonflat” means

Because every `j_(s,delta)` is a convolution unit, define the unique
coefficient transition

```text
T_(epsilon<-delta)
 =j_(s,epsilon) * j_(s,delta)^(-1) in K[C_13]^x.             (19)
```

These full-unit transitions satisfy the exact pair-groupoid law

```text
T_(zeta<-epsilon) * T_(epsilon<-delta)
 =T_(zeta<-delta).                                          (20)
```

So the full convolution-unit transport is flat.  The content of
`(10)` is instead the sharp nonreduction

```text
T_(epsilon<-delta) notin K^x <z>                            (21)
```

after every origin choice; the reversed comparison also fails.  In other
words, the transition structure group cannot be reduced from the full
circulant unit group `K[C_13]^x` to scalar times translation.

This corrects a potentially misleading use of “nonflatness.”  There is no
curvature obstruction in the full coefficient groupoid.  The obstruction
is to a smaller pointwise/projective transport law.  THM-2792's full
circulant intertwiners are therefore not dispensable bookkeeping; they are
the least coefficient-level structure currently known to transport the
bank.

## 7. First-forgotten-digit marked-power separator

The SIC prompt suggests a second, differently typed operation.  In general,
let `pi:Omega->Y` be a quotient of finite sets, choose two points
`x_0,x_1` in one fibre over `y`, and put

```text
A=K^Omega,                         Lambda=ev_(x_0)-ev_(x_1)
```

with **pointwise** multiplication.  For a marker `H in A`, a low-sheet
function `a in K^Y`, and `f=a o pi`, the exact discrete Leibniz
transgression is

```text
Lambda(f^m)=0,
Lambda(H f^m)=(H(x_0)-H(x_1))a(y)^m             for every m>=1. (22)
```

Thus the functional form of the `H`-drift is completely rigid: the
unmarked coarse power disappears and only the marker difference remains.

For the chosen two-digit chart, specialize to

```text
Omega=F_13 x F_13,             pi(r,h)=r,
A=K^Omega
```

with **pointwise** multiplication.  Put

```text
H(r,h)=h,
Lambda_r(F)=F(r,0)-F(r,1).                                  (23)
```

For any low-sheet function `a in K^(F_13)`, set `f=a o pi`.  Then for every
integer `m>=1`,

```text
Lambda_r(f^m)=0,
Lambda_r(H f^m)=-a(r)^m.                                    (24)
```

This is just the discrete Leibniz transgression: `f` is constant on the
chosen fibre while `H` changes by one.  If `a(r)!=0`, the second value in
`(24)` survives for every `m`.

For THM-2791's scalar-normalized two-point chord

```text
a=1_{6}+1_{7},
```

either `r=6` or `r=7` gives the persistent value `-1`.  The primary script
checks exponents `1,...,26`; the independent script exhausts all `8,192`
Boolean low-sheet functions at exponents `1,2,13`, for

```text
319,488 exact formula checks per field.                      (25)
```

The identity `(22)`, rather than the finite range, proves every exponent.

The proper unital subalgebra `ker Lambda_r` is therefore not a Mathieu
subspace of the pointwise function algebra; indeed the presence of `1`
already makes that non-Mathieu conclusion elementary.  What is useful here
is the **typed persistent witness** with the proved two-point semantic
support and the explicitly chosen first-forgotten-digit model.  No physical
attachment between those two data is inferred.

This is structurally parallel to THM-2801's marked-power failure of the
Special Image Conjecture, but it is not SIC.  It has:

- no Weyl contraction `E_n`;
- no identification with `Im(xi_i-partial_(z_i))`;
- no `xi/z` bihomogeneous grading;
- no polynomial-coordinate multiplier; and
- no demonstrated nilpotent-Jacobian sector.

Thus `(22)` proves no statement about the Jacobian conjecture.  It says
exactly that the high digit destroyed by the quotient can remain
load-bearing under every marked pointwise power.

## 8. Source-to-target ledger and LRC boundary

The projective-arc statement has the following type:

| item | exact content |
|---|---|
| source | canonical THM-2625 endpoint-current profiles on thirteen determinant fibres |
| target | projective cyclic/dihedral coefficient orbits, relative-ratio arcs, and `[13,2,12]` codes |
| map | profile pair to all minors, equivalently to its pointwise ratio word |
| preserved | coefficient support, relative origin gauge, orientation when requested, and central-cycle coordinate |
| destroyed | high address digits, rail ancestry, interval atom, positivity, root/Cech data, and endpoint allocation |
| needed sidecar | THM-2806's reserved natural map from THM-2791's fixed literal rail sheet to one allocated THM-2625 endpoint atom |
| cheapest decisive test | audit THM-2806's fixed-sheet scalar law against `(19)` and the high-digit/carry commutator before any further mode census |

THM-2791 now supplies a genuine same-ancestry partial rail-sheet germ with
the full graded chord.  This theorem shows that the remaining map to the
endpoint-origin current cannot be represented merely by a scalar and a
cycle-origin choice.  It must retain at least the full circulant transition
or an equivalent high-digit marker.  THM-2806 reserves exactly the
fixed-sheet allocation test, while THM-2807 reserves the complementary
positive graded-address simplex; this theorem neither duplicates nor
settles either reservation.

That is an obstruction and a target correction, not the missing map.  It
does not exclude any of the `165` LRC relation rows, prove a positive
magnitude floor, attach a semantic atom to an endpoint atom, or prove
LRC(14).

## 9. Exact replay and validity controls

Run

```bash
python 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.py
python -O 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_thm2803.py
python 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.py
python -O 04-computation/lrc14_endpoint_current_determinant_fibre_nonflatness_hostile_audit_thm2803.py
```

Normal and optimized outputs byte-match their stored transcripts.  Both
companions use explicit exception gates and no truth-bearing Python
assertions.  The primary path uses normalized determinant coordinates, all
projective minors, all frame shears, Fourier controls, and a fixed semantic
chord.  The independent path instead uses literal affine-line sets,
equivariant origin images, anchored cross-products, and an exhaustive
Boolean marked-power hostile.

The finite-field implication is one-way in the lawful direction: a nonzero
specialized minor proves the characteristic-zero cyclotomic minor nonzero.
No converse lifting from a zero specialization is used.  Every quantified
profile entry is nonzero in the same fields, so the scalar/projective and
ratio tests never divide by a specialized zero.

QED.
