---
id: THM-3505
title: "The r=5 projective slope pencil and anchored mixed-Haar basis"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED for the one canonical
  THM-2594 r=5 common-base table.  All 936 primitive C7 x C13 x C13
  coefficients with nonzero word and root-label frequencies survive at the
  certified split embedding in F_547, 72/72 on each affine slope.  The 144
  anchored mixed-Haar coordinates form another exact basis of the nonaxial
  quotient, but only 134/144 survive for each word character.  This is a
  retained-table theorem, not a U_full ancestry/address map, physical current,
  K4 bridge, row exclusion, or LRC(14) result.
source: candidate commit 5662617ac plus independent audit session 2026-08-16
depends_on:
  - THM-2594 (canonical r=5 common-base table and fixed-root hostile)
  - THM-2512 (mixed block and toothpick factorization conventions)
related:
  - MISTAKE-293 (one ancestry base is not one circle point)
  - MISTAKE-295 (constant-column erasure is not the fixed-root hostile)
candidate_script: 04-computation/lrc14_stage2_root_shear_contrast_probe_20260816.py
candidate_output: 05-knowledge/results/lrc14_stage2_root_shear_contrast_probe_20260816.out
audit_script: 04-computation/lrc14_r5_slope_pencil_independent_audit_20260816.py
audit_output: 05-knowledge/results/lrc14_r5_slope_pencil_independent_audit_20260816.out
audit_script_sha256: 74f53e6d633d2862f64013387d54bf0d0bc8f66dcb694841b2b3381665e71814
audit_output_sha256: 1e09ac04818325183d32f57bff7560fac4f172051b353580112822ca7c7da7a7
hash_basis: working-tree bytes after CRLF-to-LF normalization
---

# THM-3505 -- the r=5 projective slope pencil

## 1. Exact object and independent reconstruction

Let `N(u,q,ell,theta)` be the nonnegative integer numerator table of
[THM-2594](THM-2594-realized-theta-slaved-contraction-at-the-r5-window.md),
with `u,q,theta in F_13` and `ell in F_7`, and retain the owner and derived-root
coordinates:

```text
T_ell(u,theta) = sum_q N(u,q,ell,theta).                 (1)
```

The independent audit imports neither the THM-2594 constructor nor the
root-shear/Fourier helper.  Starting from the defining circle inequalities, it
rebuilds `E_1`, `Q_{1,{b}}`, and the source-safe-deleted word `T_b` by generic
interval intersection and difference.  A separately written transfer-profile
engine then reconstructs all `13*13*7*13 = 15,379` entries of `N`, including
the ten zero `theta` columns.  Its exact digest is

```text
18463c7af393c6090c7419c76d48b7ed40f9ac8b17d8c1be247062493e334ce8. (2)
```

This agrees with a direct digest of the canonical constructor.  Summing `q`
gives (1), digest

```text
1ba72585187c014e482f89959ff4f19e2185e6e4fcd982d497923a304e0d37d8. (3)
```

As a second exact route, the audit folds `sum_q F(y,q)` before multiplying by
the owner profile and obtains every entry of (1) again.  Thus (2)--(3) compare
the mathematical tables before any Fourier transform, not merely their final
support counts.

## 2. Complete affine slope pencil

In `F_547`, the element `w=64` has exact order `91`; hence

```text
eta = w^13 = 81  has order 7,
xi  = w^7  = 475 has order 13.                              (4)
```

For `beta in F_7`, `r,s in F_13`, define

```text
T_hat(beta,r,s)
  = sum_(ell,u,theta) T_ell(u,theta)
      eta^(-beta ell) xi^(-r u-s theta).                    (5)
```

Then, for every `beta != 0`, every `s != 0`, and every `r`,

```text
T_hat(beta,r,s) != 0 in F_547.                              (6)
```

There are `6*13*12=936` coefficients in (6), all nonzero.  Writing
`lambda=r/s`, each of the thirteen affine slopes has exactly `72/72`
nonzero coefficients.  The ordered ledger digest is

```text
e92e3f1b072db16ada1daa28925803ebd9e11658deb3532680911ed637dee85d. (7)
```

This is exact finite-field certification of nonvanishing over the integral
table: a nonzero reduction cannot come from an exact zero.  It does not claim
that reductions at other primes have the same coordinate support.

### Why slopes zero and two are the inherited controls

At `r=0`, equation (5) first sums out `u`, so it is the `C_7 x C_13` transform
of THM-2594's aggregate response.  For `beta,s != 0`, double centring changes
the table only by one-coordinate terms annihilated by these characters.
Consequently slope zero is exactly the already-proved THM-2594/THM-2512 mixed
block, not a new bridge theorem.  All `72/72` entries survive; their digest is

```text
31475920623317779b3c6de6334f258309256fb71734cc6ecd7ea6a6476b3e68. (8)
```

For the absolute root `t=theta+2u`, direct change of variables gives

```text
A_hat(beta,r,s) = T_hat(beta,r+2s,s).                       (9)
```

Therefore marginalizing the owner in the absolute-root chart (`r=0` on the
left) is slope `r=2s` on the right.  Slope two is precisely THM-2594's genuine
fixed-absolute-root hostile, again `72/72` nonzero, with digest

```text
808df0ceb8616773cff1b5c12de2333d1495c07d119c57b9b49e68aa9bb2627f. (10)
```

The other eleven slopes localize owner/root interaction that is pure in
neither of these two affine charts.

## 3. Anchored mixed-Haar coordinates are only partially supported

For `u,t in {1,...,12}`, define the anchored rectangle

```text
Delta_ell(u,t)
 = T_ell(u,t)-T_ell(u,0)-T_ell(0,t)+T_ell(0,0).             (11)
```

This is the checkerboard or mixed-Haar pairing on the `2 x 2` subtable with
rows `{0,u}` and columns `{0,t}`.  Its raw nonzero counts for
`ell=1,...,6` are

```text
(132,132,133,4,4,4),                                       (12)
```

and the first nonzero entry is

```text
Delta_1(1,2)=258805184656089173356800.                      (13)
```

After the `C_7` transform, each `beta=1,...,6` has `134/144` nonzero anchored
rectangles, hence `804/864` in total.  The first is `(beta,u,t)=(1,1,2)` with
value `180` in `F_547`; the ordered digest is

```text
5f284f6693c1f3441f51d9a86fddf6bd5d5fcaa1f76b0a9bd8add437b5e3752f. (14)
```

For every `beta`, the same ten coordinates vanish:

```text
(u,t)=(1,1),(2,1),(3,1),(4,1),(5,1),
      (7,1),(8,1),(9,1),(10,1),(11,1).                     (15)
```

Equations (6) and (15) are compatible.  The matrix

```text
M_(r,u)=xi^(-ru),       r,u in F_13^x,                      (16)
```

has rank `12` over `F_547`; its tensor square has rank `144`.  Directly,

```text
T_hat(beta,r,s)
 = sum_(u,t != 0) Delta_hat_beta(u,t) xi^(-ru-st),
                  r,s != 0.                                (17)
```

Thus the `144` anchored rectangles and the `144` nonaxial Fourier characters
are two exact bases of the same mixed quotient.  All `864/864` Fourier-basis
coordinates fire, while only `804/864` anchored coordinates fire.  Coordinate
support is not invariant under an invertible change of basis.

## 4. Exact comparison with the U_full minimum joint-address gate

The minimum U_full gate proves that a `2 x 2` marginal kernel is generated by

```text
C = [[ 1,-1],
     [-1, 1]],                                              (18)
```

and that one joint mixed-Haar scalar detects it.  Equation (11) shows that the
retained THM-2594 table has many algebraic target coordinates of exactly this
shape.  For any fixed `u,t`, the formal label map

```text
(i,j) -> (u_i,t_j),   u_0=t_0=0, u_1=u, t_1=t              (19)
```

is independent of `ell` and carries the abstract checkerboard pairing to
`Delta_ell(u,t)`.

Map (19) is not a lawful U_full transport.  The current source and target
contract is:

| field | audited content |
|---|---|
| source | actual pre-merge U_full factor-labelled endpoint cells |
| target | THM-2594 linked-node atoms underlying `T_ell(u,theta)` |
| formal map | binary labels `(i,j)` to one chosen rectangle `(0,u)x(0,t)` |
| preserved | abstract row/column margins and checkerboard pairing |
| lost | circle-cell measure, E/Q lineage, source root `q`, common-base identity, owner/word/source sheets, horizons, and `F_13^3` address |
| needed sidecar | an actual atom map preserving those fields and one fixed address before character selection |
| hostile | (18), invisible to every one-sided marginal; the calibrated positive and flat U_full couplings have identical margins but different bridge |
| verdict | current API supplies no lawful map; mathematical existence remains open |

The visible U_full record has only circle interval, E/Q lineage, component,
wrap, and boundary provenance.  It has none of `base`, `root`, `owner_sheet`,
`word_sheet`, `source_sheet`, `left_horizon`, `right_horizon`, or `address`.
Renaming `(i,j)` as `(u,theta)` would therefore discard exactly the typing
data the U_full gate says are missing.  The ten zeros in (15) additionally
show that even formal rectangle choices are not uniformly active.

Hence (11)--(17) are the right finite **mixed-Haar capacity and basis**
statement for the canonical r=5 table.  They are not a mixed-Haar bridge from
U_full.

## 5. Scope and reproduction

This theorem proves complete primitive slope support and the exact
coordinate-basis boundary for one realized THM-2594 candidate table.  It does
not identify a THM-2449 one-point response, construct the U_full atom/address
map, produce a physical phase or current, evaluate the K4 bridge, exclude a
scalar row, or change the LRC(14) ledger of `165`.

Run:

```text
python -B 04-computation/lrc14_r5_slope_pencil_independent_audit_20260816.py
python -B -O 04-computation/lrc14_r5_slope_pencil_independent_audit_20260816.py
```

Normal and optimized outputs are byte-identical and match the stored output.
The candidate and audit agree on (7)--(10); the two independent table
constructors agree on (2)--(3).  **QED.**
