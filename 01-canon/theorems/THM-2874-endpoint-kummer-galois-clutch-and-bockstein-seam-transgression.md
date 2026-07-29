---
id: THM-2874
title: "Endpoint Kummer--Galois clutch and Bockstein seam transgression"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  The THM-2868 projective frequency orbit is carried to the THM-2857
  centered endpoint Galois orbit by an explicit Q(zeta_91)-rational
  scalar.  Independently, the thirteen frequency sections realize the
  character-three coordinate on the sharp C169 carry fibre over q3.
  Collapsing that fibre along the tempting q3/q11/q7 seam has exact
  Cech defect omega^3-1, the THM-2851 Bockstein generator.  The q11
  signed row still cancels, q7 remains E3-zero, and no physical
  common-support extension or LRC decrement is proved.
source: root/lrc-holotopy-allocation-2026-07-28
depends_on:
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2870-prime-power-convolution-versus-physical-diagonal-intertwiner-obstruction
related:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2863-endpoint-prony-splitter-and-carry-character-three-intertwiner
script: 04-computation/lrc14_endpoint_kummer_galois_bockstein_thm2874.py
output: 05-knowledge/results/lrc14_endpoint_kummer_galois_bockstein_thm2874.out
script_sha256: 3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0
output_sha256: 90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455
hash_basis: LF-normalized bytes
---

# THM-2874 -- Endpoint Kummer--Galois clutch and Bockstein seam transgression

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

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

Indeed, substituting (1) into the right side of (3) gives
`xi^(1118+1020) omega^(3r)=xi^(955+546r)`, because the
exponents differ by `2366`.  Thus (3) is not merely an equality of
character labels: relative to the pinned normalizations, it is an explicit
`F`-linear clutch between the projective multiplier-frequency coordinate
and the centered coefficient-field Galois coordinate.  It cancels the
non-rational common source scale of THM-2868.

This remains an algebraic coefficient reference.  The coordinate `t_r`
comes from nonlinear projectivization of two signed split branches, and
THM-2857 proves that the canonical physical ancestry action is
coefficient-linear rather than the required Galois-semilinear action.

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
Only the `q3` fibre is populated: THM-2868's signed `q11` row cancels and
its `q7` row is zero by `E3`.  Extending (8)--(9) over those target fibres
with common physical support is still open.

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

THM-2870 makes the same boundary linear and sharp.  On `C13`,
convolution by `delta_3` is an invertible cyclic permutation, while
physical multiplication by the same mask has rank one.  Its convolution
`1`-eigenspace also has dimension one.  Every same-mask intertwiner has
rank at most one; Fourier inversion cannot manufacture the missing
physical `q11/q7` rows or turn the projective clutch (3) into a linear
physical coefficient/character reference.

## 5. Exact evidence

The companion pins the promoted THM-2851, THM-2857, THM-2868, and
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
  states, q11 signed mass, q7 E3 support, a common endpoint chart over
  all target fibres, and a linear physical coefficient/character map;

needed sidecar:
  a lawful common-support extension of the populated q3 frequency fibre
  over q11 and q7, with the Bockstein edge phase and a genuine macro-E3
  complement transporter;

cheapest decisive test:
  construct one nonzero q7 complement coefficient at the same source,
  target, and endpoint chart as the q3 split, retain one uncancelled q11
  origin coefficient, and compare the q11-via-q7 phase with the direct
  q3-to-q7 phase.  Holonomy omega^3 realizes the ancestry attachment;
  holonomy one reproduces the flat hostile.
```

No row is excluded and the LRC(14) ledger remains `165`.
