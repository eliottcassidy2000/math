---
id: THM-2446
title: "The 2-jet z-graded Jacobian decomposition and its cone system"
status: >
  PROVED (the six-bracket z-graded decomposition of det J for
  z-quadratic maps F = A z^2 + B z + C and its Keller reformulation;
  Euler slice reductions E5 vacuous, E4 = -2[a_x,a_y,b],
  E3 = (4/3)[b_x,b_y,a] - 4[a_x,a_y,g]; the cone law
  [Phi_x,Phi_y,Phi] = (z/4) det J(Phi) on leading forms) +
  VERIFIED-EXACT (generic 27- and 54-symbol residuals zero; canon
  THM-1310 map passes with exact D1 cancellation and D0 = -2;
  hostile xy^3 perturbation localized to D1'/D0'; tame witness
  (x+z^2, y, z)). This opens the 2-jet d=4 door of THM-1335 SS5 /
  backlog HYP-8130(a): the identities exist and the design staircase
  is one quadratic layer deeper than THM-1310 SS7. The theorem does
  not itself exhibit or exclude a wild 2-jet Keller map; THM-3438 later
  supplies an explicit wild degree-four member.  It does not prove the
  order-{1,3} conjecture or any box emptiness, and system-consistency
  in the 2-jet box is explicitly NOT a wildness certificate.
source: kind-pasteur-2026-07-26-S131
depends_on:
  - THM-1310-conic-pair-fibers-and-design-equations
related:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-1335-trisection-modulus-master-identities-trace-polynomiality
  - THM-1340-engine-trichotomy-zaffine-keller
  - THM-2451-twojet-plane-map-cone-prediction-hostile
script: 04-computation/jacobian_twojet_zgraded_identity_kps_S131.py
output: 05-knowledge/results/jacobian_twojet_zgraded_identity_kps_S131.out
script_sha256: 38122318f4086f4c14b570c65189bd1aca0e04179cad7fe6cabb698728cc45ae
output_sha256: 4a439f183b9eb2b717080f64057528eba7fb505457d02585f7cb20b8596944b3
hash_basis: working-tree bytes (LF)
---

# THM-2446 -- the 2-jet Jacobian is six brackets; its cone is Euler's

**PROVED + VERIFIED-EXACT** as itemized in the status.

> **EXISTENCE UPDATE (THM-3438 / MISTAKE-396).**  The graded identities here
> were never an existence obstruction.  The weighted quartic now gives a
> wild `S_4` 2-jet Keller map in the high-degree line-cap tail.  All identities,
> tame controls, and scoped negative statements below remain unchanged.

The dormant "2-jet architecture" door (THM-1335 SS5, backlog item
HYP-8130(a): "derive the z-graded det J identities; one hard
session") asks what replaces THM-1310 SS7's design equations when the
seed is z-quadratic. This theorem supplies the complete graded
system, its gauge bookkeeping, and its leading-form cone law.

## 1. The decomposition (PROVED)

Let `F(x,y,z) = A(x,y) z^2 + B(x,y) z + C(x,y)` with
`A, B, C : C^2 -> C^3` (C^1 components suffice; polynomial in every
application here). Take `J = [F_x | F_y | F_z]` columnwise and write
`[u,v,w] = det(u|v|w)`. Then, identically,

```text
det J = D0 + D1 z + D2 z^2 + D3 z^3 + D4 z^4 + D5 z^5,          (1)

D5 = 2[A_x,A_y,A]
D4 = [A_x,A_y,B] + 2[A_x,B_y,A] + 2[B_x,A_y,A]
D3 = [A_x,B_y,B] + [B_x,A_y,B]
     + 2([A_x,C_y,A] + [C_x,A_y,A] + [B_x,B_y,A])
D2 = [B_x,B_y,B] + [A_x,C_y,B] + [C_x,A_y,B]
     + 2([B_x,C_y,A] + [C_x,B_y,A])
D1 = [B_x,C_y,B] + [C_x,B_y,B] + 2[C_x,C_y,A]
D0 = [C_x,C_y,B].                                                (2)
```

*Proof.* The three columns split by z-degree:
`F_x = A_x z^2 + B_x z + C_x`, `F_y` likewise, and `F_z = 2A z + B`
(every factor `2` in (2) traces to this one derivative).
Trilinearity of the determinant expands `det J` into exactly
`3 x 3 x 2 = 18` brackets `[P_x, Q_y, R]` with `P, Q in {A,B,C}`,
`R in {2A, B}`; collecting by total z-degree gives (2), each bracket
occurring exactly once, with coefficient 2 precisely when the third
slot is `A`. Since `1, z, ..., z^5` are `C[x,y]`-independent:

```text
F Keller  <=>  D5 = D4 = D3 = D2 = D1 = 0 in C[x,y]
              and D0 = [C_x,C_y,B] = nonzero constant.           (3)
```

QED. The companion additionally verifies (1)-(2) with fully generic
degree-1 components (27 independent symbols -- itself proof-grade,
since both sides depend only on first jets and jet evaluation is
surjective onto `C^27`, so any wrong integer multiplicity would
leave a nonzero universal trilinear form) and with generic degree-2
components (54 symbols).

At `A = 0` the system degenerates to the z-affine triple
`[B_x,B_y,B] = 0`, `[B_x,C_y,B] + [C_x,B_y,B] = 0`,
`[C_x,C_y,B] = const`, the THM-1310 setting.

## 2. Canon grounding and controls (VERIFIED-EXACT)

Positive control: the THM-1310 counterexample-fiber map (`u = 1+xy`,
`det J = -2`) embeds as `A = 0`, `B = (u^3, 3xu^2, -x^3)`, and
passes (3) through a nontrivial exact cancellation:

```text
[B_x,C_y,B] = -6x^2(xy+1),   [C_x,B_y,B] = +6x^2(xy+1),
D1 = 0,   D2 = 0,   D0 = -2.                                     (4)
```

Hostile control: perturbing `C_1 -> C_1 + xy^3` breaks Keller and
the harness localizes the failure to the correct graded slots:

```text
D2' = 0   (C never enters D2 when A = 0),
D1' = -6x^4 y^2 (xy+1)(xy+3),
D0' = -2(9x^5y^5 + 30x^4y^4 + 5x^3y^3 - 9x^2y^2 + 1).            (5)
```

The d=3 ruling identity `B_2^3 + 27 B_1^2 B_3 = 0` (cone over the
cuspidal cubic) is re-verified on the same `B`.

## 3. The cone system on leading forms (PROVED)

Let `alpha = A^(2)`, `beta = B^(3)`, `gamma = C^(4)` be the leading
forms in the minimal quartic box `deg(A,B,C) <= (2,3,4)` and
`Phi = alpha z^2 + beta z + gamma` the homogeneous leading map
(weighted so that `Phi` is quartic for `z` of weight one... the
grading is by total (x,y)-degree per z-slice). Euler's relation
`4 Phi = x Phi_x + y Phi_y + z Phi_z` gives the exact cone law

```text
4 [Phi_x, Phi_y, Phi] = z det J(Phi),                            (6)
```

so Keller forces `[Phi_x,Phi_y,Phi] = 0`: the leading map is
non-dominant -- its image lies in a cone hypersurface. This is the
d=4 analogue of the d=3 cuspidal-cubic ruling cone. The top slices
`E_k := D_k^(9-k)` of the graded system reduce (verified exactly on
generic homogeneous forms):

```text
E5 = 0 identically             (the z^5 top slice is Euler-vacuous;
                                the constraining D5 slices are lower),
E4 = -2 [alpha_x, alpha_y, beta],
E3 = (4/3)[beta_x, beta_y, alpha] - 4 [alpha_x, alpha_y, gamma]. (7)
```

`E4 = 0` says `beta` lies in the tangent span of the alpha-cone;
`E3 = 0` is the first genuinely coupled balance law
`[beta_x,beta_y,alpha] = 3 [alpha_x,alpha_y,gamma]`.

Full-strength `D5 = 0` says the direction map `[A] : C^2 -> P^2` is
non-dominant: image a point, a line, or a conic. The cusp-carrying
capacity of the direction map is exhausted at d=3; at d=4 any
swallowtail/A_3 geometry must be carried by the next flag layer, the
gauge-covariant plane map `w = A x B` (degree <= 5).

## 4. Gauge, counts, staircase (DERIVED-EXACT)

Gauge group of the 2-jet shape has dimension 22: z-translation
`z -> z + phi`, `phi` affine (3; acts
`(A,B,C) -> (A, B + 2A phi, C + B phi + A phi^2)`), z-scaling (1),
source `Aff(C^2)` (6), target `Aff(C^3)` (12). Parameters in the
minimal quartic box: `18 + 30 + 45 = 93`, effective `71`.
Coefficient equations: at most `194` genuine closed conditions --
an upper bound only: measured per-box ranks are strictly smaller
(in the `deg C <= 2` sub-box the degree-4 part of `D2` cancels
identically, leaving 10 conditions). Any hunt must measure ranks
per box, never quote the naive bounds. Given a conical `A`
(`D5 = 0`): `D4` is linear in `B`; `D3, D2` linear in `C` with
B-quadratic sources; `D1` linear-plus-quadratic in `C` (the term
`2[C_x,C_y,A]` has no d=3 analogue); `D0` quadratic in `C` -- one
quadratic layer deeper than THM-1310 SS7's staircase.

**Hunt-design warning (established by witness).** The 2-jet box
contains tame Keller maps: `F = (x + z^2, y, z)` satisfies the full
graded system with `D0 = 1` and field degree 1 (companion leg 6).
Unlike the z-affine d=3 box, "system consistent + det nonzero" is
NOT a wildness certificate; any Groebner/Rabinowitsch hunt must
stratify away the tame locus (alpha-conic nondegeneracy plus
non-planar `w`-curve, or a field-degree-4 certificate) and report
detection floors accordingly.

## 5. Open layer (labelled; P2 now REFUTED)

Recorded as predictions only: (P1) a swallowtail master identity
`mu(P) L_4 = disc(z^4 + P~ z^2 + Q~ z + R~)` generalizing
`108 a^2 L = (12a-b^2)^3 + E^2`, with Jelonek set `{L_4 = 0}` the
swallowtail hypersurface ("A_3 at both ends"); (P3) the `S_4 > V_4`
tower matches resolvent-cubic / fold-double-cover layers; (P4) a
Veronese-projection generative ansatz for the staircase.

**(P1) update (S132):** the classical layer is now PROVED as
THM-2455 (`27 Delta_4 = 4I^3 - J^2`; exact resolvent-disc equality,
unit 1; ideal-exact swallowtail strata), which also RETIRES the
"A_3 at both ends" slogan (the quadrisection pencil provably avoids
the cuspidal edge; endpoints A1 and Z/2-symmetric A1A1) and the
naive `unit x square x Jelonek` disc shape (`mu` is a non-unit on
explicit degree-4 2-jet maps). The corrected hunt target is
HYP-9027's odd-exponent Jelonek law.  THM-3438 removes the former blocking
hypothesis G1 by supplying a field-degree-four 2-jet Keller map; its full
finite Jelonek/infinity-inertia audit is now the cheapest nonempty test of the
Keller-restricted successor.

The former prediction (P2), that `D5=D4=0` plus the alpha-part of
`D3` forces `[w_x,w_y,w]=0`, is **REFUTED twice over, concurrently**
(codex THM-2451 and kind-pasteur S131b, same day, independent code
paths). THM-2451's separated-ruling family

```text
A=a+xn,            B=b+yn,            C=0
```

has `D5=...=D0=0` but `[w_x,w_y,w]=[a,b,n]^2`. The S131b companions
`04-computation/jacobian_twojet_p2_refutation_core_kps_S131b.py` and
`..._conic_kps_S131b.py` (outputs in 05-knowledge/results, sha256 LF
`2c18ec3a...`/`415b5833...` and `615d7b62...`/`85a86847...`) extend
the refutation to every stratum and settle THM-2451's stated
residual ("a repaired conic-direction version remains open"):

- universal bracket lemma, identically in first jets:
  `[w_x,w_y,w] = [A,B,A_x][A,B,B_y] - [A,B,A_y][A,B,B_x]`
  (THM-2451's `[a,b,n]^2` is its evaluation on the ruling family);
- no identity `E = P D5 + Q D4 + R S3` exists with ANY polynomial
  multipliers (explicit first-jet point with all constraints zero
  and `E = 1` -- non-membership even in the radical);
- line-cap counterexample `A = (x,1,0), B = (y,0,1)` extends to a
  full tame KELLER map `F = (xz^2+yz+xy-x^3+x, z^2+3x^2-2y, x+z)`
  with `det J = -2` and explicit polynomial inverse: even the whole
  graded system does not force the w-cone;
- **the conic-direction version is also refuted**: the nondegenerate
  conic cap `A = (x^2, xy, y^2)`,
  `B = (-2x^3, x^2-2x^2y, x+2xy-2xy^2)` has `D5 = D4 = S3 = 0` and
  `E = -x^8`;
- P2 survives exactly on the point-cap stratum (`E = 0` for rank-1
  `A`, arbitrary `B`) and for `A, B` valued in a common plane.

Reading: the plane map `w = A x B` genuinely can carry dominant
(A_3-capable) geometry on Keller solutions -- the flag layer is
richer than the d=3 analogy suggested, which is favorable for (P1)
and removes a false constraint from any seed hunt.

## 6. Scope

Nothing proved in this theorem itself exhibits or excludes a wild 2-jet
Keller map; THM-3438 supplies one by a different construction.  Nothing here proves
the order-{1,3} conjecture, or decides any box emptiness. Groebner
box probes run during derivation were inconclusive within their
watchdogs and are deliberately not cited as negatives.

## 7. Reproduction

```bash
python 04-computation/jacobian_twojet_zgraded_identity_kps_S131.py
```

Fully symbolic and deterministic (sympy); every check is a hard
failure; final line `ALL CHECKS PASSED`.
