# Session 144: Jacobian classification, ETH and MRB

**2026-11-03. Checkpoint 1: independently audited Jacobian composition and correction.** The ETH and MRB continuations are in progress and are not claimed as results of this checkpoint.

## What is established now

[The composition theorem](jacobian/equivariant_composition_tower.md) constructs an equivariant Keller tower `G^m` with determinant `2^m`, generic degree `3^m`, and orbit-scale degree `(3*5^(m-1)-1)/2`. The first nonlinear scale already occurs at `m=2`, with degree seven. Its reduced scale divisor has at least m disjoint smooth irreducible components. This is an infinite structural family, not a finite-box inference.

The tower corrects the literal higher-scale existence question in [THM-4727](../../../01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md). Affine-scale tangent-sweep seeds do not exhaust their composition closure. The divisor component count obstructs equivalence to an affine-scale map under equivariant source/target automorphisms; scale degree is not asserted invariant under arbitrary coordinate changes.

[The independent audit](jacobian/independent/INDEPENDENT_AUDIT.md) also gives a complete complex-fiber classification for the fixed degree-three map, with explicit boundary equations, a missing-curve certificate, and hostile examples for the discarded-wall and global-deformation claims. This classifies one map's fibers, not all Keller maps.

[Correction lineage](jacobian/CORRECTION_LINEAGE.md) and [source pins](jacobian/SOURCE_PINS.json) preserve the inherited evidence. No complete classification of dimension-three or higher counterexamples, primitive solutions, or arbitrary Aut-by-Aut equivalence classes is claimed. The planar Jacobian conjecture remains open.

## Verification

The producer's 37 explicit checks and the auditor's 46 checks use exact algebra and remain active under Python optimization. Seven corrupted certificates are rejected; independent source enumeration agrees with reconstructed fibers at all 1,799 targets over `F5`, `F7`, and `F11`. Those finite checks do not replace the universal proofs.

The [frozen Jacobian manifest](jacobian/MANIFEST.json) excludes this mutable overview and navigation. Run from the `jacobian` directory:

```
python3 -B verify_frozen_packet.py
python3 -B -O verify_frozen_packet.py
```

Both modes verify every frozen file hash and replay the two exact programs against pinned output bytes. SymPy is required. No new Isabelle or Lean build was run. Full repository documentation checks were not run because they read the restricted central ledger; scoped owned-file/link/hash checks and unchanged startup-surface identities are recorded with the publication evidence.
