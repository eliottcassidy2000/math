---
id: THM-3312
title: "Exceptional-quadratic trace, norm, and cofactor anti-descent"
status: >
  PROVED (universal quadratic algebra) + INDEPENDENTLY VERIFIED-EXACT
  (THM-3309's two fixed slices).  For B=A[t]/(P t^2-Q t+R), trace, norm,
  and conjugate-difference square give the complete unordered passport of
  every a+bt.  In the
  exceptional critical fibre the critical y-root, first-normal velocity, and
  elimination cofactor ratio are genuinely C2-conjugate.  Although the
  elimination pair U,W is unimodular, its projective ratio does not descend:
  its alternating cross-square is 256 P^4 delta.  It generates B/A, while
  its unordered trace/norm passport and conjugate-difference square descend.
  The gradient pair itself vanishes, so no
  Keller Bezout row exists and the mate-integrability class is not entered.
  This proves no inverse, Keller mate, JC(2), or DC(2).
source: root/creative-passports/2026-08-03
depends_on:
  - THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus
  - THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction
related:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go
script: 04-computation/jc_exceptional_quadratic_cofactor_passport_thm3312.py
output: 05-knowledge/results/jc_exceptional_quadratic_cofactor_passport_thm3312.out
script_sha256: f2c87512b3d67c005822e934bb8ba8f8011f6ea4a0bda42f6202c10723cd0106
output_sha256: 2ce88681395c40081f5bb68c9de81a5bac2e7d7f56626834bbfac6c304f547d0
imported_response_script: 04-computation/jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py
imported_response_script_sha256: 6f050a583004172f812c3f7729427079d5df45c3a985c2e470b2a0d34ad8f337
hash_basis: LF-normalized bytes
---

# THM-3312 -- exceptional-quadratic trace, norm, and cofactor anti-descent

**PROVED + INDEPENDENTLY VERIFIED-EXACT.**

THM-3309 proves the stronger fixed-slice pointed-deck and true-gradient
theorem over THM-3306's transverse base ideal.  The present theorem extracts
its universal quadratic trace/norm passport and supplies an independent exact
replay of the branch anti-descent identities.  It does not replace THM-3309's
direct gradient computation.

## 1. Universal quadratic passport

Let `A` be a field of characteristic zero and

```text
F(t)=P t^2-Q t+R,       P != 0,
delta=Q^2-4PR != 0,     B=A[t]/(F),                       (1)
```

with `F` irreducible.  Then `B/A` is quadratic and

```text
sigma(t)=Q/P-t.                                           (2)
```

For every `z=a+bt`, direct reduction modulo `(1)` gives

```text
Tr(z)=2a+bQ/P,
Nm(z)=a^2+abQ/P+b^2R/P,
(z-sigma(z))^2=b^2 delta/P^2.                             (3)
```

Thus trace, norm, and conjugate-difference square descend to `A`.  The element
`z` itself descends exactly when `b=0`; otherwise it generates `B/A`.

On the exceptional divisor the strict linear row is `y+t=0`.  Therefore

```text
Tr(y)=-Q/P,       Nm(y)=R/P,
(y-sigma(y))^2=delta/P^2.                                 (4)
```

If the first-normal velocity is `v_0+v_1t` with `v_1!=0`, formula `(3)` gives
its descended unordered passport and proves that the two velocities remain
distinct.

## 2. The elimination ratio carries the lost branch label

Put

```text
D=6P-4Q,       h:=D/4,
J=D-4Pt,       W=-4J,       U=P^2-hJ.                     (5)
```

These are the exceptional specializations of the inverse-graph elimination
pair.  They satisfy

```text
U-(h/4)W=P^2.                                             (6)
```

Hence `(U,W)` is unimodular in `B`.  Irreducibility of `F` also makes `W`
nonzero: `W=0` would put `t=D/(4P)` in `A`.

In the basis `(1,t)`, the coefficient determinant of `U,W` is

```text
(P^2-hD)(16P)-(4hP)(-4D)=16P^3.                          (7)
```

Consequently

```text
(U sigma(W)-sigma(U)W)^2=256 P^4 delta != 0.              (8)
```

For `rho=U/W`, division by `Nm(W)^2` gives

```text
(rho-sigma(rho))^2=256 P^4 delta/Nm(W)^2 != 0.            (9)
```

Thus the projective elimination-cofactor ratio does not descend to `A`; it is
a primitive generator of `B/A`.  Its trace and norm do descend and determine
the unordered pair.  Passing to them loses exactly the quadratic branch label.

## 3. The Keller gate is a separate, earlier condition

The passport algebra above concerns the elimination pair.  On the fixed slice,
THM-3309 independently evaluates the true gradients and proves

```text
P_x=P_z=0 in B_i.                                         (10)
```

Its localized cubics obey the invertible triangular change

```text
R_1=V P_z,
R_2=V^3 P_x-(V'y/2)R_1,       V in B_i^*.                 (11)
```

Thus the new companion's check `R_1=R_2=0` agrees with, but is not substituted
for, THM-3309's direct gradient check.  A hypothetical gradient Bezout row

```text
S_1P_x+S_2P_z=1                                           (12)
```

would specialize to `0=1`; relative trace would give `0=2`.  No primitive
Keller cofactor exists on this fibre.  The mate-integrability/divergence class,
which is defined only after gradient unimodularity, is therefore not a new
obstruction here: it is not entered.

Identity `(6)` does not contradict this.  It makes the **elimination** pair
unimodular; it does not express `1` as a combination of the true gradients.

## 4. Independent exact replay of the THM-3309 specialization

THM-3309 pins the independently frozen nonsquare certificates and proves that
the exceptional quadratic is irreducible in each of the two THM-3212 accessory
fields on the fixed `C=c+x,d=k=1` slice.  The new companion then reconstructs

```text
A_i=K_i[x]/(linear_a),       deg(linear_a)=36,
B_i=A_i[t]/(P t^2-Q t+R).                                (13)
```

It keeps `linear_a` distinct from the degree-32 pivot `P`, verifies the
quadratic relation without choosing a root, and checks `(3)--(9)` for the
critical `y`-root, first-normal velocity, `U`, `W`, and `U/W`, together with
the localized pair in `(11)`.  Every displayed
trace, norm, and difference square is a nonzero unit in `A_i`; none of these
packets collapses further to the accessory field `K_i`.  In both cases the
cross-square is exactly `256P^4 delta` and nonzero.

The ordinary and optimized exact replays both complete with
`ALL EXACT PASSPORT CHECKS PASSED`.  The stored output is the optimized replay;
the ordinary replay has the same deterministic invariant lines.

## 5. Scope and reproduction

The universal passport theorem is `(1)--(9)`; the gradient gate `(10)--(12)`
is inherited from and independently cross-checked against THM-3309.  The
finite specialization is sharply limited to the two named fields and the
fixed slice `C=c+x,d=k=1`.  It releases no
deformation parameter, supplies no `A_i`-rational root, and constructs no
Keller mate or inverse.  `JC(2)` and `DC(2)` remain open.

Run

```text
python 04-computation/jc_exceptional_quadratic_cofactor_passport_thm3312.py
python -O 04-computation/jc_exceptional_quadratic_cofactor_passport_thm3312.py
```

from the repository root.  The computation uses exact polynomial and finite
extension arithmetic and imports only the frozen exceptional-quadratic
response constructor pinned in the frontmatter.

**QED.**
