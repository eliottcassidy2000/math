---
id: HYP-3110
title: LRC14 De Moivre, Jacobi theta, and crystallographic proof-carrier frontier
status: EVIDENCE / synthesis and Lean sidecar ledger; not a proof
source: codex-2026-06-27-S263
tangent: T1186
technique: LTI-247
tournament_technique: LTT-145
script: 04-computation/lrc14_de_moivre_jacobi_crystallographic_frontier_codex_s263.py
result: 05-knowledge/results/lrc14_de_moivre_jacobi_crystallographic_frontier_codex_s263.out
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCCrystallographicThetaFrontier.lean
reflection: 07-reflections/lrc14-de-moivre-jacobi-crystallographic-frontier-codex-s263.md
related:
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-3073
  - HYP-3063
  - HYP-2614
  - HYP-2613
  - HYP-2309
  - OPEN-Q-108
---

# HYP-3110: LRC14 De Moivre, Jacobi Theta, And Crystallographic Proof-Carrier Frontier

## Claim

This lane tests the user's De Moivre quintic, Jacobi theta,
17-wallpaper-group, and 230-space-group prompt against the live LRC14 proof
frontier.

It is not a proof.  The intended contribution is a controlled-forgetting
interface:

```text
de_moivre_quintic_fold
  -> exact finite-depth quintic cancellation detector;
jacobi_theta_tail
  -> signed residue-cusp / support-six theta tail carrier;
wallpaper_17_orbifold
  -> finite 2D crystallographic quotient audit;
space_group_230_orbifold
  -> finite 3D crystallographic quotient audit;
lee_yang_root_curve
  -> HYP-3109 zero-locus / zero-real ear sidecar;
finite_address_packet
  -> HYP-3107 terminal proof interface.
```

The preserved LRC predicate is: a residual row must retain enough coverage,
observer, endpoint-owner, residue-cusp, and finite-address data to imply
`LRC14Statement` through the HYP-3107 frontier.  The destroyed coordinates are
raw runner labels, raw time, and any scalarized root/moment/lattice count.

## Exact Scout Results

Artifacts:

```text
04-computation/lrc14_de_moivre_jacobi_crystallographic_frontier_codex_s263.py
05-knowledge/results/lrc14_de_moivre_jacobi_crystallographic_frontier_codex_s263.out
```

The scout verifies De Moivre's solvable quintic fold as exact
Laurent-polynomial cancellation.  For `x = u - a/u`,

```text
x^5 + 5*a*x^3 + 5*a^2*x = u^5 - a^5/u^5.
```

Equivalently, `x^5+5a*x^3+5a^2*x+b=0` reduces to
`y^2+b*y-a^5=0` with `y=u^5`.  This matters for LRC only as a
finite-depth cancellation detector: if a degree-5 interaction does not expose
the `1,5,5` De Moivre fold and the substitution sidecar, no quintic radical
shortcut has been found.

The crystallographic finite catalogs are recorded as sidecar budgets:

```text
wallpaper_group_count = 17
space_group_3d_count = 230
space_group_counts_by_crystal_system =
  {triclinic:2, monoclinic:13, orthorhombic:59, tetragonal:68,
   trigonal:25, hexagonal:27, cubic:36}
bravais_3d_count = 14
jacobi_theta_channels = theta1_odd, theta2_shifted,
  theta3_lattice, theta4_alternating
```

Readout: the 17 and 230 are finite quotient audits.  They do not prove LRC14
by analogy.  A row quotient that uses them must name translation lattice,
stabilizer word, glide/screw/torsion sidecar, preserved LRC predicate, and
destroyed coordinate.

Tournament Analysis uses proof-carrier sidecars as vertices:

```text
finite_address_exit
observer_gluing_certificate
jacobi_theta_tail
lee_yang_root_curve
de_moivre_quintic_fold
space_group_230_orbifold
wallpaper_17_orbifold
raw_scalar_shadow
```

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
edge_flips_vs_novelty_first_gauge = 10
tie_path =
  finite_address_exit > observer_gluing_certificate > jacobi_theta_tail
  > lee_yang_root_curve > de_moivre_quintic_fold
  > space_group_230_orbifold > wallpaper_17_orbifold > raw_scalar_shadow
```

The `10` edge flips under a novelty-first gauge are useful: exact finite
catalogs are attractive, but proof readiness still puts observer-gluing and
finite-address production first.

## Lean Sidecar Ledger

The Lean module

```text
TournamentH7.LRCCrystallographicThetaFrontier
```

records exact finite counts:

```text
wallpaperGroup_count = 17
spaceGroup3D_count = 230
bravais3D_count = 14
jacobiThetaChannel_count = 4
crystallographicThetaCarrier_count = 8
```

It also proves the De Moivre fold over `Rat`:

```text
deMoivreQuintic_fold_rat
```

and gives a conservative conditional theorem:

```text
lrc14_from_theta_crystallographic_residuals
```

This theorem says HYP-3110 closes against `LRC14Statement` only after the
theta/crystallographic residual producer emits an `ObserverGluingCertificate`
or `FiniteAddressBranchPacket`.  In other words, the new sidecars may feed
HYP-3107; they do not bypass it.

## Proof-Frontier Consequences

1. De Moivre is the solvable-normal-form lane, not the earlier A5
   Abel-Ruffini wall.  It should be used to detect exact degree-5 folds in
   finite expansions, then demoted if the coefficient pattern is absent.
2. Jacobi theta belongs directly to HYP-2614's signed residue-cusp/support-six
   tail after HYP-2616-style low-height wall deletion.  The likely theorem is a
   theta tail bound with deleted anti-cosets and boundary faces named.
3. The 17 wallpaper groups and 230 space groups are finite orbifold quotient
   ledgers.  They are useful precisely when a residual packet has periodic
   local structure and the quotient records the sidecar that prevents
   illegal forgetting.
4. HYP-3109 root curves remain upstream.  A crystallographic quotient that
   forgets the zero-collision/root-locus sidecar is not proof-legal.
5. The next proof-facing data structure is a HYP-2963 residual-row schema with
   theta/orbifold columns that can emit concrete `ObserverGluingCertificate`
   rows and then `FiniteAddressBranchPacket` rows when compression is
   available.

## Assumption Challenge

Do not assume the useful vertices are runners, arcs, roots, or group names.
The candidate vertices are proof-carrier sidecars: quintic-fold normal forms,
theta tails, wallpaper orbifold quotients, space-group orbifold quotients,
root-curve packets, observer-gluing certificates, and finite-address exits.
