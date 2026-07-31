---
id: THM-2874
title: "Endpoint Kummer--Galois clutch and Bockstein seam transgression"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The THM-2868 projective frequency orbit is carried to the THM-2857
  centered endpoint Galois orbit by an explicit Q(zeta_91)-rational
  scalar, and recombination recovers exactly the THM-2861 Hermitian
  edge.  Independently, the thirteen frequency sections realize the
  character-three coordinate on the sharp C169 carry fibre over q3.
  Collapsing that fibre along the tempting q3/q11/q7 seam has exact
  Cech defect omega^3-1, the THM-2851 Bockstein generator.  The full
  endpoint-pattern complement supplies the coefficient-support half of
  a q7 chart at the zero address, but its formal frequency triangle is
  flat and its literal factor projection is a three-corner horn.  The
  q11 signed row still cancels, q7 remains E3-zero, and no physical
  common-support extension or LRC decrement is proved.
source: root/lrc-holotopy-allocation-2026-07-28
depends_on:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
  - THM-2861-endpoint-hermitian-edge-holonomy-and-semilinear-clutch-test
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2870-prime-power-convolution-versus-physical-diagonal-intertwiner-obstruction
related:
  - THM-2863-endpoint-prony-splitter-and-carry-character-three-intertwiner
script: 04-computation/lrc14_endpoint_kummer_galois_bockstein_thm2874.py
output: 05-knowledge/results/lrc14_endpoint_kummer_galois_bockstein_thm2874.out
script_sha256: 3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0
output_sha256: 90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455
addendum_script: 04-computation/lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.py
addendum_output: 05-knowledge/results/lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.out
addendum_script_sha256: 44dfdefbf5392e7840f74e63d190a96a484af71ba9bd31df3ce62a22b827d67e
addendum_output_sha256: f914d934c40ef58ea5df0f0df0c61c357c9ab9073db3ba7cbb044d8564886cab
hash_basis: LF-normalized bytes
---

# THM-2874 -- Endpoint Kummer--Galois clutch and Bockstein seam transgression

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact base-field clutch

Put

```text
xi=zeta_2366,       omega=xi^182,       F=Q(xi^13)=Q(zeta_91).
```

THM-2857 writes the centered relative-Galois endpoint orbit as

```text
c_r-mean(c)=-B omega^(3r),             B=xi^1020.       (1)
```

THM-2868 gives the projective physical-frequency orbit

```text
t_r=xi^(955+546r)=xi^955 omega^(3r).                    (2)
```

The two character-three lines have an explicit normalization over the
base field:

```text
t_r=kappa [c_r-mean(c)],

kappa=-xi^(-65)=xi^1118=(xi^13)^86 in F.               (3)
```

Indeed, retaining both displayed signs gives
`(-xi^(-65))(-xi^1020)=xi^955`; equivalently the exponent calculation is
`1118+1020+1183=955 mod2366`.  Thus (3) is not merely an equality of
character labels: relative to the pinned normalizations, it is an explicit
`F`-linear clutch between the projective multiplier-frequency coordinate
and the centered coefficient-field Galois coordinate.  It cancels the
non-rational common source scale of THM-2868.

This remains an algebraic coefficient reference.  The coordinate `t_r`
comes from nonlinear projectivization of two signed split branches, and
THM-2857 proves that the canonical physical ancestry action is
coefficient-linear rather than the required Galois-semilinear action.

### 1a. Recombination is exactly the Hermitian edge

The two THM-2868 Prony branches retain more than the projective ratio.
Writing their common nonzero source scalar as `P_src`, their exact
normalization is

```text
V_0=P_src A,                U_r=-P_src B omega^(3r),
S_r=U_r+V_r=P_src c_r,      1+t_r=c_r/A.               (H1)
```

Moreover `1+t_r^(-1)=A bar(c_r)`.  Therefore the normalized adjacent
projective product is not just of the same character type as THM-2861:

```text
(1+t_(r+1))(1+t_r^(-1))=c_(r+1) bar(c_r),              (H2)

S_(r+1) bar(S_r)
 =P_src bar(P_src)c_(r+1)bar(c_r).                     (H3)
```

Thus the signed frequency-charted coefficient edge differs from the
canonical Hermitian edge only by the fixed nonzero norm scalar in (H3).
The edge has Fourier support exactly `{0,3,10}` and thirteen distinct
nonzero values.  It is not a physical full current.

The symmetric trace does not create a new unoriented measurement.  Once
the forward frequency orientation is independently retained, the law
`E_r=omega^3 bar(E_r)` recovers

```text
E_r=omega^3 (1+omega^3)^(-1)(E_r+bar(E_r)).             (H4)
```

Reversing the edge has the same trace, however.  No physical polarization
measurement or adjacent common support has been constructed.

## 2. The frequency atlas realizes one sharp carry fibre

Do not identify the frequency index with a target residue.  Instead use

```text
S=C13_r x C13_q,

iota(r,q)=13r+q mod169,                                 (4)
```

with natural lifted residue action

```text
L_h(r,q)
 =(r+floor((q+h)/13), q+h mod13).                       (5)
```

Then

```text
iota(L_h(r,q))=iota(r,q)+h mod169,                       (6)

L_k L_h=T^epsilon(h,k) L_(h+k mod13),

T(r,q)=(r+1,q),       epsilon(h,k)=floor((h+k)/13),      (7)
```

and `L_1^13=T`.  Hence (4)--(7) are exactly the sharp first-carry
`C169 -> C13` torsor of THM-2851, with `r` as the carry coordinate.

THM-2868's split-left frequency coefficients satisfy

```text
U_(r+1)=omega^3 U_r.                                    (8)
```

Consequently

```text
U_(r+floor((q+h)/13))
 =omega^(3 floor((q+h)/13)) U_r.                        (9)
```

All thirteen values in (8) are nonzero and are reconstructed from 26
actual 91-unit multiplier probes.  In this precise coefficient-level
sense, the frequency measurement atlas supplies a faithful
character-three coordinate on the entire carry fibre over the populated
target row `q=3`.  It realizes the correct thirteen-state carrier
externally.

This does not put thirteen internal states into one positive packet.
Only the `q3` fibre is populated in the signed table used here:
THM-2868's origin-resolved `q11` rows are nonzero but cancel under that
selector, while its `q7` row is zero by `E3`.  Installing an
origin-resolved `q11` fibre and then extending (8)--(9) to `q7` with one
common physical support remains a separate problem.

## 3. The tempting seam is the forbidden section

The two repaired THM-2868 frequency charts give the numerical seam

```text
q3 <-> r=0,       q11 <-> r=8,       q7 <-> r=4.        (10)
```

This is the affine assignment `r(q)=q-3`.  It collapses the independent
carry coordinate in (4) onto the target residue and therefore cannot be
a global section of the nonsplit lift (5).

The THM-2851 arrows are

```text
q3 --8--> q11 --9--> q7,
q3 --------4-------> q7.                                (11)
```

On the frequency line, section displacement `h` acts by `omega^(3h)`.
The three edge scalars induced by (10) are

```text
E_3,11=omega^11,       E_11,7=omega,       E_3,7=omega^12,

E_11,7 E_3,11=E_3,7.                                  (12)
```

Thus the collapsed frequency seam is flat.  The ancestry carry line
instead uses (5):

```text
A_3,11=1,       A_11,7=omega^3,       A_3,7=1,

A_11,7 A_3,11=omega^3 A_3,7.                          (13)
```

Suppose nonzero vertex scalars extended the normalized one-line
comparison of THM-2863 over all three edges:

```text
A_i,j phi_i=phi_j E_i,j.                               (14)
```

Propagating from `q3` to `q7` directly and through `q11` gives

```text
phi_7^via=omega^3 phi_7^direct,

phi_7^via-phi_7^direct
 =(omega^3-1)phi_7^direct !=0.                         (15)
```

Therefore no nonzero global edgewise gauge exists.  The defect in (15)
is the exact Cech/mapping-cone transgression of the collapsed seam.  It
lies on the same character-three line as THM-2851's oriented derivative

```text
449(omega^(-3)-1),                                     (16)
```

and the companion computes their nonzero normalization in both certified
endpoint fields.  There is no further representation-theoretic mismatch:
the failure is precisely the Bockstein generator.

There is also an exact cotangent shadow on the semantic carry edge.  From
the exponents in THM-2868,

```text
t_4/t_8=xi^(773-591)=xi^182=omega.                     (17)
```

Evaluating THM-2851's first-residue derivative
`449(Z-1)` at `Z=omega` therefore gives

```text
449(t_4-t_8)/t_8=449(omega-1).                         (18)
```

In the cyclotomic cotangent quotient
`Z[omega]/(13,(omega-1)^2)`, (18) is

```text
7(omega-1),                                            (19)
```

exactly matching `partial_T h_L=7e mod e^2` under
`e <-> omega-1`.  Equations (17)--(19) are a normalized cotangent
shadow, not an equality of physical coefficient vectors: the assignment
(10) and the character evaluation are still external choices.

## 4. Sharp physical boundary

At the physical target boundary, the THM-2868 signed current has mask

```text
h=delta_q3:       q3 !=0,       q11=0,       q7=0.      (20)
```

The two zeros respectively come from origin cancellation and `E3`
absence.  Invertible edge transport would be nonzero at all three seam
vertices, so vertexwise gauge cannot turn (20) into the local system of
Section 3.

The literal deep-bit-only endpoint repair also fails.  The THM-2868
endpoint word has the deep factor dangerous and every other factor safe;
flipping only that bit makes all nine factors safe.  Exact restriction at
both origins for each of `q3,q7,q11` is empty.  The 20-cell **macro**
`E3` complement of THM-2847 is therefore not the same object as an
endpoint deep-bit complement.

### 4a. The full-pattern complement closes the cheapest test negatively

At the zero endpoint chart, the underlying half-open interval locus of
the entire nine-factor `PAT_E3` pattern—not a single-bit
approximation—is exactly the underlying THM-2847 macro `E3` locus.  This
is an extensional equality, not a typed endpoint-to-macro functor.
Consequently its full set-theoretic complement gives a Boolean
coefficient-support chart.  The three choices

```text
(r,q,truth,offset)
 =(0,3,E3,0), (4,7,not-E3,1), (8,11,E3,2)              (21)
```

are full translated weight-one atoms on the exact 20-cell horn.  Their six
Prony multipliers

```text
1, 2, 170, 171, 339, 340
```

and frequencies

```text
38, 64, 4432, 4458, 8826, 8852
```

are all `91`-units.  The translated charts have the same two endpoint
nodes, and their transported split coefficients agree with the full
THM-2868 atlas.  Character `b=10` is the unique target character for which
`3r+bq` is constant on (21), so `U_r omega^(10q)` is constant.

This executes the coefficient-support half of the proposed cheap
complement test at the same zero address, and its formal coefficient seam
is still flat.  No typed transition or Cech morphism has been built.  The
`q11 -> q7` projective ratio is `omega`, but

```text
rho(9)rho(8)rho(4)^(-1)=1,       rho(h)=omega^(3h),     (22)
```

whereas the aligned ancestry carry has holonomy `omega^3`.  Thus the
complement supplies the q7 coefficient vertex but not the Bockstein
attachment.

There is also a parallel literal Boolean boundary explaining what the
binary truth projection forgets.  The q7 pattern is `DSSSSDSSD`:
complement occupancy is paid for by losing both guard and q5 safety while
retaining the `C3` danger.  Along the complete target-residue orbit, the
guard/q5 factor projection assumes exactly

```text
(safe,safe), (dangerous,safe), (dangerous,dangerous),  (23)
```

and never `(safe,dangerous)`.  This two-factor projection—not the entire
nine-bit word—is uniform over all `169` endpoint addresses because every
canonical representative fixes both corresponding endpoint coordinates to
zero.  The q7 chart is the
`(dangerous,dangerous)` corner of this punctured square.  Passing only to
the extensional `E3/not-E3` truth value forgets the two atomic safeties and
cannot manufacture the missing fourth corner.  A typed macro
`E3/complement` contraction compatible with `QA/QAB` ancestry is still
absent.

THM-2870 makes the same boundary linear and sharp.  On `C13`,
convolution by `delta_3` is an invertible cyclic permutation, while
physical multiplication by the same mask has rank one.  Its convolution
`1`-eigenspace also has dimension one.  Every same-mask intertwiner has
rank at most one; Fourier inversion cannot manufacture the missing
physical `q11/q7` rows or turn the projective clutch (3) into a linear
physical coefficient/character reference.

## 5. Exact evidence

The primary companion pins the promoted THM-2851, THM-2857, THM-2868, and
THM-2870 exact scripts.  In both certified endpoint fields it verifies:

1. all thirteen instances of the `F`-rational clutch (3);
2. all `2,197` actions and `28,561` compositions of (5)--(7), together
   with the coefficient response (9);
3. the cotangent seam (17)--(19), every frequency and ancestry edge
   scalar in (12)--(13), the nonzero transgression (15), and its nonzero
   normalization to (16);
4. the rank `13` versus rank `1` same-mask hostile and the
   one-dimensional convolution fixed space; and
5. all six literal all-nine-safe endpoint restrictions.

Run from the repository root:

```text
python 04-computation/lrc14_endpoint_kummer_galois_bockstein_thm2874.py
python -O 04-computation/lrc14_endpoint_kummer_galois_bockstein_thm2874.py
```

Both modes LF-normalized byte-match the stored 13-line transcript.

```text
script 3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0
output 90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455
```

The addendum independently pins the primary companion and THM-2847,
THM-2851, THM-2861, and THM-2868.  It verifies (H1)--(H4), the complete
thirteen-residue literal factor orbit, the uniform three-corner projection
(23), the exact Boolean partition, all six `91`-unit probes, common Prony
nodes, the unique `(3,10)` compensated character, and the hostile equality
(22).  Normal and optimized modes byte-match its stored transcript:

```text
python 04-computation/lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.py
python -O 04-computation/lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.py

script 44dfdefbf5392e7840f74e63d190a96a484af71ba9bd31df3ce62a22b827d67e
output f914d934c40ef58ea5df0f0df0c61c357c9ab9073db3ba7cbb044d8564886cab
```

The independent immutable audit rederived the base-field clutch, the
transported `C169` group law and carry cocycle, the signed-table population
boundary, the flat-versus-Bockstein seam holonomy, the absence of a global
vertex gauge, the cotangent-shadow typing, the literal all-nine-safe
hostile, and the THM-2870 rank-one obstruction.  Normal, optimized, and
stored output agree exactly, and both declared LF-normalized hashes match.
The audit also confirmed that the origin-resolved `q11` values are nonzero
and equal; only their signed origin difference vanishes.  No remaining
mathematical, typing, or evidence defect was found.

The addendum was separately replayed in normal and optimized modes and
independently rederived at the interval, factor, Prony, Hermitian, character,
and carry levels.  The audit confirmed that only the guard/q5 projection
is uniform over the `169` addresses, that the endpoint/macro equality is
only an underlying zero-chart locus equality, and that (22) is a formal
coefficient seam rather than a physical Cech transporter.

## 6. Connection contract

```text
source:
  THM-2868 projective frequency and split-left chi_3 atlases,
  THM-2857 centered endpoint Galois torsor, and
  THM-2851 first ancestry-carry quotient;

target:
  an F-rational frequency/coefficient clutch, the sharp C169 carrier
  over q3, and the q3/q11/q7 Bockstein seam;

map:
  apply kappa in (3), retain frequency r as the independent carry
  coordinate in (5), then compare the collapsed section r=q-3 with the
  three natural ancestry lifts;

preserved:
  projective gauge cancellation, F-structure, character 3, all thirteen
  carry-fibre states over q3, target labels, arrow orientation, and the
  exact Bockstein derivative line;

destroyed / missing:
  a single positive packet, internal rather than external-measurement
  states, q11 signed mass, q7 E3 support, the guard/q5 (safe,dangerous)
  corner, a typed macro-to-ancestry contraction, and a linear physical
  coefficient/character map;

needed sidecar:
  within the canonical address-orbit/full-pattern route, a lawful
  guard/q5 fourth corner or equally faithful sidecar, followed by a typed
  macro-E3/complement contraction aligned with QA/QAB ancestry and
  carrying the Bockstein edge phase; a different transporter may bypass
  this literal chamber;

cheapest decisive test:
  CLOSED NEGATIVELY by (21)--(23): the full-pattern complement gives a
  nonzero q7 coefficient at the common chart but has holonomy one.
  Within this route, the next test must retain the guard and q5 truth
  coordinates separately, realize the missing (safe,dangerous) corner or
  an equally faithful sidecar, and only then compare its QA/QAB-typed
  q11-via-q7 phase with the direct phase.
```

No row is excluded and the LRC(14) ledger remains `165`.
